// SLOW5/BLOW5 reader. Compiled only when ENABLE_SLOW5=ON.
// Uses the slow5lib C API (same pattern as the ri_read_sig_slow5 path in RawHash2),
// but streaming single reads via slow5_get_next
// since this tool processes every read sequentially.

#include "signal_reader.hpp"

extern "C" {
#include "slow5/slow5.h"
}

#include <stdexcept>
#include <string>
#include <vector>

namespace segmentation {

namespace {

struct Slow5FileGuard {
    slow5_file_t *p = nullptr;
    ~Slow5FileGuard() { if (p) slow5_close(p); }
};

struct Slow5RecGuard {
    slow5_rec_t *p = nullptr;
    ~Slow5RecGuard() { if (p) slow5_rec_free(p); }
};

} // anonymous namespace

std::vector<Read> read_slow5(const std::string &path) {
    Slow5FileGuard f;
    f.p = slow5_open(path.c_str(), "r");
    if (!f.p) {
        throw std::runtime_error("SLOW5: failed to open " + path);
    }

    std::vector<Read> reads;

    while (true) {
        slow5_rec_t *rec = nullptr;
        int ret = slow5_get_next(&rec, f.p);
        if (ret < 0) {
            // Normal EOF path (SLOW5_ERR_EOF). Any other negative value is
            // treated as EOF too; slow5lib doesn't differentiate in a way
            // that matters for this read-all-or-fail tool.
            if (rec) slow5_rec_free(rec);
            break;
        }
        Slow5RecGuard guard; guard.p = rec;

        // Calibration: pA = (raw + offset) * (range / digitisation)
        const double scale  = rec->range / rec->digitisation;
        const double offset = rec->offset;

        std::vector<float> pA;
        pA.reserve(rec->len_raw_signal);
        for (uint64_t i = 0; i < rec->len_raw_signal; ++i) {
            float v = static_cast<float>((rec->raw_signal[i] + offset) * scale);
            if (v > 30.0f && v < 200.0f)
                pA.push_back(v);
        }

        Read r;
        r.read_id   = std::string(rec->read_id);
        r.signal_pA = std::move(pA);
        reads.push_back(std::move(r));
    }

    return reads;
}

} // namespace segmentation
