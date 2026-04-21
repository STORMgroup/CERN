// POD5 reader. Compiled only when ENABLE_POD5=ON.
// Uses the official pod5_format C API (same pattern as
// RawHash2/src/rsig.c::ri_read_sig_pod5).

#include "signal_reader.hpp"

#include "pod5_format/c_api.h"

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace segmentation {

namespace {

// RAII guard for pod5_init / pod5_terminate.
struct Pod5Init {
    Pod5Init()  { pod5_init(); }
    ~Pod5Init() { pod5_terminate(); }
};

struct Pod5FileGuard {
    Pod5FileReader_t *p = nullptr;
    ~Pod5FileGuard() { if (p) pod5_close_and_free_reader(p); }
};

struct Pod5BatchGuard {
    Pod5ReadRecordBatch_t *p = nullptr;
    ~Pod5BatchGuard() { if (p) pod5_free_read_batch(p); }
    void reset(Pod5ReadRecordBatch_t *next = nullptr) {
        if (p) pod5_free_read_batch(p);
        p = next;
    }
};

} // anonymous namespace

std::vector<Read> read_pod5(const std::string &path) {
    static Pod5Init pod5_lifetime;

    Pod5FileGuard file;
    file.p = pod5_open_file(path.c_str());
    if (!file.p) {
        throw std::runtime_error(std::string("POD5: failed to open ") + path +
                                 ": " + pod5_get_error_string());
    }

    size_t batch_count = 0;
    if (pod5_get_read_batch_count(&batch_count, file.p) != POD5_OK) {
        throw std::runtime_error(std::string("POD5: failed to query batch count: ") +
                                 pod5_get_error_string());
    }

    std::vector<Read> reads;

    for (size_t batch_idx = 0; batch_idx < batch_count; ++batch_idx) {
        Pod5BatchGuard batch;
        Pod5ReadRecordBatch_t *raw_batch = nullptr;
        if (pod5_get_read_batch(&raw_batch, file.p, batch_idx) != POD5_OK) {
            throw std::runtime_error(std::string("POD5: failed to get batch: ") +
                                     pod5_get_error_string());
        }
        batch.p = raw_batch;

        size_t row_count = 0;
        if (pod5_get_read_batch_row_count(&row_count, batch.p) != POD5_OK) {
            throw std::runtime_error(std::string("POD5: failed to get row count: ") +
                                     pod5_get_error_string());
        }

        for (size_t row = 0; row < row_count; ++row) {
            uint16_t read_table_version = 0;
            ReadBatchRowInfo_t read_data;
            if (pod5_get_read_batch_row_info_data(batch.p, row,
                                                  READ_BATCH_ROW_INFO_VERSION,
                                                  &read_data,
                                                  &read_table_version) != POD5_OK) {
                throw std::runtime_error("POD5: failed to read row info");
            }

            // Calibration formula: pA = (raw + offset) * scale — same as
            // signal_reader.py / extern/RawHash2/src/rsig.c.
            const size_t n = read_data.num_samples;
            std::vector<int16_t> raw(n);
            if (pod5_get_read_complete_signal(file.p, batch.p, row, n, raw.data()) != POD5_OK) {
                throw std::runtime_error(std::string("POD5: failed to read signal: ") +
                                         pod5_get_error_string());
            }

            std::vector<float> pA;
            pA.reserve(n);
            const float scale  = read_data.calibration_scale;
            const float offset = read_data.calibration_offset;
            for (size_t i = 0; i < n; ++i) {
                float v = (static_cast<float>(raw[i]) + offset) * scale;
                if (v > 30.0f && v < 200.0f)
                    pA.push_back(v);
            }

            char read_id_buf[37] = {0};
            pod5_format_read_id(read_data.read_id, read_id_buf);

            Read r;
            r.read_id   = std::string(read_id_buf);
            r.signal_pA = std::move(pA);
            reads.push_back(std::move(r));
        }
    }

    return reads;
}

} // namespace segmentation
