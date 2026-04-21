#ifndef SEGMENTATION_SIGNAL_READER_HPP
#define SEGMENTATION_SIGNAL_READER_HPP

#include <string>
#include <vector>

namespace segmentation {

struct Read {
    std::string read_id;
    std::vector<float> signal_pA; // calibrated and 30-200 pA filtered
};

// Apply the 30-200 pA hardware filter that rsig.c uses in RawHash/RawHash2.
// Samples outside that band are discarded. In place is not safe; returns a new
// vector so callers don't need to worry about iterator invalidation.
std::vector<float> filter_pA(const float *samples, size_t n);

// Read every read from one signal file. Extension (lowercased) picks the
// backend. Throws std::runtime_error on I/O failure or unsupported format.
std::vector<Read> read_signal_file(const std::string &path);

// Read every read from a single file, or every supported file in a directory.
// Warnings for individual file failures are printed to stderr.
std::vector<Read> read_signal_path(const std::string &path);

// --- Backend entry points, implemented in fast5_reader.cpp / pod5_reader.cpp /
// slow5_reader.cpp. Each is compiled only when the corresponding feature is
// enabled; signal_reader.cpp stubs them out otherwise. ---

std::vector<Read> read_fast5(const std::string &path);
std::vector<Read> read_pod5(const std::string &path);
std::vector<Read> read_slow5(const std::string &path);

} // namespace segmentation

#endif
