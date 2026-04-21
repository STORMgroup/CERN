#include "signal_reader.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace segmentation {

std::vector<float> filter_pA(const float *samples, size_t n) {
    std::vector<float> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        float v = samples[i];
        if (v > 30.0f && v < 200.0f)
            out.push_back(v);
    }
    return out;
}

static std::string lower_ext(const std::filesystem::path &p) {
    std::string ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return ext;
}

// --- Stubs for disabled backends -------------------------------------------
// When a format toggle is OFF the corresponding .cpp is excluded from the
// build, so the linker would fail without these. They throw at runtime so
// the error message points the user at the right CMake option.

#ifndef SEG_HAVE_HDF5
std::vector<Read> read_fast5(const std::string &path) {
    throw std::runtime_error(
        "FAST5 support not compiled in (file: " + path + "). "
        "Rebuild with -DENABLE_HDF5=ON");
}
#endif

#ifndef SEG_HAVE_POD5
std::vector<Read> read_pod5(const std::string &path) {
    throw std::runtime_error(
        "POD5 support not compiled in (file: " + path + "). "
        "Rebuild with -DENABLE_POD5=ON");
}
#endif

#ifndef SEG_HAVE_SLOW5
std::vector<Read> read_slow5(const std::string &path) {
    throw std::runtime_error(
        "SLOW5/BLOW5 support not compiled in (file: " + path + "). "
        "Rebuild with -DENABLE_SLOW5=ON");
}
#endif

// --- Dispatch ---------------------------------------------------------------

std::vector<Read> read_signal_file(const std::string &path) {
    std::filesystem::path p(path);
    const std::string ext = lower_ext(p);

    if (ext == ".fast5")                    return read_fast5(path);
    if (ext == ".pod5")                     return read_pod5(path);
    if (ext == ".slow5" || ext == ".blow5") return read_slow5(path);

    throw std::runtime_error("Unsupported file format: " + ext +
                             " (supported: .fast5, .pod5, .slow5, .blow5)");
}

std::vector<Read> read_signal_path(const std::string &path) {
    std::filesystem::path p(path);
    if (!std::filesystem::exists(p))
        throw std::runtime_error("Path does not exist: " + path);

    if (std::filesystem::is_regular_file(p))
        return read_signal_file(path);

    if (!std::filesystem::is_directory(p))
        throw std::runtime_error("Not a file or directory: " + path);

    std::vector<std::filesystem::path> files;
    for (const auto &entry : std::filesystem::directory_iterator(p)) {
        if (!entry.is_regular_file()) continue;
        const std::string ext = lower_ext(entry.path());
        if (ext == ".fast5" || ext == ".pod5" ||
            ext == ".slow5" || ext == ".blow5")
            files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());

    std::vector<Read> all_reads;
    for (const auto &f : files) {
        try {
            auto reads = read_signal_file(f.string());
            all_reads.insert(all_reads.end(),
                             std::make_move_iterator(reads.begin()),
                             std::make_move_iterator(reads.end()));
        } catch (const std::exception &e) {
            std::cerr << "Warning: failed to read " << f.string()
                      << ": " << e.what() << "\n";
        }
    }
    return all_reads;
}

} // namespace segmentation
