// FAST5 reader. Compiled only when ENABLE_HDF5=ON.
// Uses the HDF5 C API directly (no hdf5_tools / fast5 C++ wrapper), so this
// file stands on its own against plain libhdf5.a.
//
// Supports both FAST5 variants:
//   - Single-read: top-level "Raw/Signal" + "UniqueGlobalKey/channel_id"
//   - Multi-read:  top-level "read_<UUID>" groups, each with its own
//                  "Raw/Signal" and "channel_id"

#include "signal_reader.hpp"

#include <hdf5.h>

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace segmentation {

namespace {

struct H5FileGuard {
    hid_t id = -1;
    ~H5FileGuard() { if (id >= 0) H5Fclose(id); }
};
struct H5GroupGuard {
    hid_t id = -1;
    ~H5GroupGuard() { if (id >= 0) H5Gclose(id); }
};
struct H5AttrGuard {
    hid_t id = -1;
    ~H5AttrGuard() { if (id >= 0) H5Aclose(id); }
};
struct H5DsetGuard {
    hid_t id = -1;
    ~H5DsetGuard() { if (id >= 0) H5Dclose(id); }
};
struct H5SpaceGuard {
    hid_t id = -1;
    ~H5SpaceGuard() { if (id >= 0) H5Sclose(id); }
};
struct H5TypeGuard {
    hid_t id = -1;
    ~H5TypeGuard() { if (id >= 0) H5Tclose(id); }
};

// Read a double-valued attribute, whatever its on-disk type.
bool read_attr_double(hid_t parent, const char *name, double *out) {
    if (H5Aexists(parent, name) <= 0) return false;
    H5AttrGuard a; a.id = H5Aopen(parent, name, H5P_DEFAULT);
    if (a.id < 0) return false;
    return H5Aread(a.id, H5T_NATIVE_DOUBLE, out) >= 0;
}

// Read a string-valued attribute (variable- or fixed-length).
bool read_attr_string(hid_t parent, const char *name, std::string *out) {
    if (H5Aexists(parent, name) <= 0) return false;
    H5AttrGuard a; a.id = H5Aopen(parent, name, H5P_DEFAULT);
    if (a.id < 0) return false;
    H5TypeGuard t; t.id = H5Aget_type(a.id);
    if (t.id < 0) return false;

    if (H5Tis_variable_str(t.id) > 0) {
        char *buf = nullptr;
        H5TypeGuard mem; mem.id = H5Tcopy(H5T_C_S1); H5Tset_size(mem.id, H5T_VARIABLE);
        if (H5Aread(a.id, mem.id, &buf) < 0) return false;
        if (buf) { out->assign(buf); free(buf); }
        return true;
    } else {
        size_t sz = H5Tget_size(t.id);
        std::vector<char> buf(sz + 1, 0);
        if (H5Aread(a.id, t.id, buf.data()) < 0) return false;
        out->assign(buf.data());
        return true;
    }
}

// Read the int16 "Signal" dataset at <raw_path>/Signal into `out`.
bool read_signal_dataset(hid_t file, const std::string &raw_path,
                         std::vector<int16_t> *out) {
    const std::string sig_path = raw_path + "/Signal";
    H5DsetGuard d; d.id = H5Dopen2(file, sig_path.c_str(), H5P_DEFAULT);
    if (d.id < 0) return false;
    H5SpaceGuard s; s.id = H5Dget_space(d.id);
    if (s.id < 0) return false;
    hssize_t npts = H5Sget_simple_extent_npoints(s.id);
    if (npts < 0) return false;
    out->resize(static_cast<size_t>(npts));
    return H5Dread(d.id, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   out->data()) >= 0;
}

struct CalibrationAttrs { double digitisation = 0, offset = 0, range = 0; };

bool read_calibration(hid_t parent_group_id, CalibrationAttrs *c) {
    return read_attr_double(parent_group_id, "digitisation", &c->digitisation) &&
           read_attr_double(parent_group_id, "offset",       &c->offset) &&
           read_attr_double(parent_group_id, "range",        &c->range);
}

std::vector<float> calibrate_and_filter(const std::vector<int16_t> &raw,
                                        const CalibrationAttrs &c) {
    const double scale = c.range / c.digitisation;
    std::vector<float> pA;
    pA.reserve(raw.size());
    for (int16_t v : raw) {
        float x = static_cast<float>((v + c.offset) * scale);
        if (x > 30.0f && x < 200.0f) pA.push_back(x);
    }
    return pA;
}

// Detect single-read FAST5: has top-level "Raw" group.
bool is_single_read_fast5(hid_t file) {
    return H5Lexists(file, "/Raw", H5P_DEFAULT) > 0;
}

// Collect names of all "read_*" groups at "/".
struct ReadGroupCollector { std::vector<std::string> names; };

herr_t collect_read_groups_cb(hid_t /*loc*/, const char *name,
                              const H5L_info_t * /*info*/, void *op_data) {
    auto *c = static_cast<ReadGroupCollector*>(op_data);
    if (std::strncmp(name, "read_", 5) == 0)
        c->names.emplace_back(name);
    return 0;
}

void read_one_multi(hid_t file, const std::string &read_group,
                    std::vector<Read> *out) {
    const std::string raw_path = "/" + read_group + "/Raw";
    const std::string ch_path  = "/" + read_group + "/channel_id";

    H5GroupGuard ch; ch.id = H5Gopen2(file, ch_path.c_str(), H5P_DEFAULT);
    if (ch.id < 0) return;
    CalibrationAttrs cal;
    if (!read_calibration(ch.id, &cal)) return;

    // read_id attribute lives on the <read_group>/Raw group.
    H5GroupGuard raw; raw.id = H5Gopen2(file, raw_path.c_str(), H5P_DEFAULT);
    if (raw.id < 0) return;
    std::string read_id;
    if (!read_attr_string(raw.id, "read_id", &read_id))
        read_id = read_group; // fallback

    std::vector<int16_t> sig;
    if (!read_signal_dataset(file, raw_path, &sig)) return;

    Read r;
    r.read_id   = read_id;
    r.signal_pA = calibrate_and_filter(sig, cal);
    out->push_back(std::move(r));
}

void read_single(hid_t file, std::vector<Read> *out) {
    // Calibration lives at /UniqueGlobalKey/channel_id
    H5GroupGuard ch; ch.id = H5Gopen2(file, "/UniqueGlobalKey/channel_id", H5P_DEFAULT);
    if (ch.id < 0) return;
    CalibrationAttrs cal;
    if (!read_calibration(ch.id, &cal)) return;

    // Each read is at /Raw/Reads/Read_<N>
    H5GroupGuard reads; reads.id = H5Gopen2(file, "/Raw/Reads", H5P_DEFAULT);
    if (reads.id < 0) return;

    struct Collector { std::vector<std::string> names; } coll;
    H5Literate(reads.id, H5_INDEX_NAME, H5_ITER_INC, nullptr,
               [](hid_t /*loc*/, const char *name, const H5L_info_t * /*info*/,
                  void *op_data) -> herr_t {
                   static_cast<Collector*>(op_data)->names.emplace_back(name);
                   return 0;
               }, &coll);

    for (const auto &read_name : coll.names) {
        const std::string raw_path = "/Raw/Reads/" + read_name;
        H5GroupGuard rg; rg.id = H5Gopen2(file, raw_path.c_str(), H5P_DEFAULT);
        if (rg.id < 0) continue;

        std::string read_id;
        if (!read_attr_string(rg.id, "read_id", &read_id))
            read_id = read_name;

        std::vector<int16_t> sig;
        if (!read_signal_dataset(file, raw_path, &sig)) continue;

        Read r;
        r.read_id   = read_id;
        r.signal_pA = calibrate_and_filter(sig, cal);
        out->push_back(std::move(r));
    }
}

} // anonymous namespace

std::vector<Read> read_fast5(const std::string &path) {
    H5FileGuard f; f.id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f.id < 0)
        throw std::runtime_error("FAST5: failed to open " + path);

    std::vector<Read> reads;
    if (is_single_read_fast5(f.id)) {
        read_single(f.id, &reads);
    } else {
        ReadGroupCollector coll;
        H5Literate(f.id, H5_INDEX_NAME, H5_ITER_INC, nullptr,
                   collect_read_groups_cb, &coll);
        for (const auto &rg : coll.names)
            read_one_multi(f.id, rg, &reads);
    }
    return reads;
}

} // namespace segmentation
