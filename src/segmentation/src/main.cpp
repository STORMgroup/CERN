// generate_events — CLI port of the old generate_events.py.
// Reads FAST5/POD5/SLOW5/BLOW5 signal files, runs RawHash or RawHash2 event
// detection, and writes a TSV of `read_id\tevent1\tevent2\t...\n` per read.

#include "rawhash_events.h"
#include "rawhash2_events.h"
#include "signal_reader.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace {

// --- Default parameters (match python-legacy/event_detector.py) -------

struct RH1Params {
    uint32_t window_length1 = 3;
    uint32_t window_length2 = 6;
    float    threshold1     = 4.30265f;
    float    threshold2     = 2.57058f;
    float    peak_height    = 1.0f;
};

struct RH2Params {
    uint32_t window_length1     = 3;
    uint32_t window_length2     = 9;
    float    threshold1         = 4.0f;
    float    threshold2         = 3.5f;
    float    peak_height        = 0.4f;
    uint32_t min_segment_length = 0;
    uint32_t max_segment_length = 500;
};

// R10.4.1 preset (used when --r10 is passed). Matches
// python-legacy/event_detector.py::RH2_R10_DEFAULTS
RH2Params rh2_r10_defaults() {
    RH2Params p;
    p.window_length1 = 3;
    p.window_length2 = 6;
    p.threshold1     = 6.5f;
    p.threshold2     = 4.0f;
    p.peak_height    = 0.2f;
    p.min_segment_length = 0;
    p.max_segment_length = 500;
    return p;
}

// --- CLI parsing -----------------------------------------------------------

struct Args {
    std::string input;
    std::string output;
    std::string method = "rawhash2";
    bool        r10    = false;

    // Per-method optional overrides. nullopt = use default.
    std::optional<uint32_t> rh1_w1, rh1_w2, rh2_w1, rh2_w2,
                             rh2_min_seg, rh2_max_seg;
    std::optional<float>    rh1_t1, rh1_t2, rh1_ph,
                             rh2_t1, rh2_t2, rh2_ph;
};

void print_usage(const char *prog) {
    std::cout <<
        "Usage: " << prog << " -i <input> -o <output> [options]\n"
        "\n"
        "Generate precomputed event values from nanopore signal files using\n"
        "RawHash or RawHash2 segmentation algorithms.\n"
        "\n"
        "Required:\n"
        "  -i, --input  <path>    Signal file (FAST5/POD5/SLOW5/BLOW5) or\n"
        "                         a directory of signal files.\n"
        "  -o, --output <path>    Output TSV file path.\n"
        "\n"
        "Common:\n"
        "  -m, --method {rawhash,rawhash2}   Algorithm (default: rawhash2)\n"
        "      --r10                         Use RawHash2 R10.4.1 preset\n"
        "  -h, --help                        Show this message\n"
        "\n"
        "RawHash v1.0 parameters (override defaults 3/6/4.30265/2.57058/1.0):\n"
        "      --rh1-window1 <int>\n"
        "      --rh1-window2 <int>\n"
        "      --rh1-threshold1 <float>\n"
        "      --rh1-threshold2 <float>\n"
        "      --rh1-peak-height <float>\n"
        "\n"
        "RawHash2 parameters (defaults 3/9/4.0/3.5/0.4/0/500; --r10 changes them):\n"
        "      --rh2-window1 <int>\n"
        "      --rh2-window2 <int>\n"
        "      --rh2-threshold1 <float>\n"
        "      --rh2-threshold2 <float>\n"
        "      --rh2-peak-height <float>\n"
        "      --rh2-min-seg-len <int>\n"
        "      --rh2-max-seg-len <int>\n";
}

// Pull the next argv entry as a string, or die.
const char *next_arg(int argc, char **argv, int &i, const char *flag) {
    if (++i >= argc) {
        std::cerr << "Error: missing value for " << flag << "\n";
        std::exit(2);
    }
    return argv[i];
}

Args parse_args(int argc, char **argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string f = argv[i];
        if (f == "-h" || f == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (f == "-i" || f == "--input") {
            a.input = next_arg(argc, argv, i, "--input");
        } else if (f == "-o" || f == "--output") {
            a.output = next_arg(argc, argv, i, "--output");
        } else if (f == "-m" || f == "--method") {
            a.method = next_arg(argc, argv, i, "--method");
        } else if (f == "--r10") {
            a.r10 = true;
        } else if (f == "--rh1-window1") {
            a.rh1_w1 = static_cast<uint32_t>(std::atoi(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh1-window2") {
            a.rh1_w2 = static_cast<uint32_t>(std::atoi(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh1-threshold1") {
            a.rh1_t1 = static_cast<float>(std::atof(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh1-threshold2") {
            a.rh1_t2 = static_cast<float>(std::atof(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh1-peak-height") {
            a.rh1_ph = static_cast<float>(std::atof(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh2-window1") {
            a.rh2_w1 = static_cast<uint32_t>(std::atoi(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh2-window2") {
            a.rh2_w2 = static_cast<uint32_t>(std::atoi(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh2-threshold1") {
            a.rh2_t1 = static_cast<float>(std::atof(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh2-threshold2") {
            a.rh2_t2 = static_cast<float>(std::atof(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh2-peak-height") {
            a.rh2_ph = static_cast<float>(std::atof(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh2-min-seg-len") {
            a.rh2_min_seg = static_cast<uint32_t>(std::atoi(next_arg(argc, argv, i, f.c_str())));
        } else if (f == "--rh2-max-seg-len") {
            a.rh2_max_seg = static_cast<uint32_t>(std::atoi(next_arg(argc, argv, i, f.c_str())));
        } else {
            std::cerr << "Error: unknown argument: " << f << "\n\n";
            print_usage(argv[0]);
            std::exit(2);
        }
    }

    if (a.input.empty() || a.output.empty()) {
        std::cerr << "Error: --input and --output are required\n\n";
        print_usage(argv[0]);
        std::exit(2);
    }
    if (a.method != "rawhash" && a.method != "rawhash2") {
        std::cerr << "Error: --method must be 'rawhash' or 'rawhash2'\n";
        std::exit(2);
    }
    return a;
}

// Write "read_id\tevent1\tevent2\t...\n" using the same %.6f format as
// python-legacy/generate_events.py so outputs stay byte-for-byte comparable.
void write_tsv_line(std::ostream &out, const std::string &read_id,
                    const float *events, size_t n) {
    out << read_id;
    char buf[32];
    for (size_t i = 0; i < n; ++i) {
        std::snprintf(buf, sizeof(buf), "%.6f", events[i]);
        out << '\t' << buf;
    }
    out << '\n';
}

} // anonymous namespace

int main(int argc, char **argv) {
    Args args = parse_args(argc, argv);

    std::cout << "Reading signals from: " << args.input << std::endl;
    std::vector<segmentation::Read> reads;
    try {
        reads = segmentation::read_signal_path(args.input);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    std::cout << "  Found " << reads.size() << " read(s)" << std::endl;
    if (reads.empty()) {
        std::cerr << "No reads found. Exiting.\n";
        return 1;
    }

    // Create the parent directory of the output path, matching the Python
    // wrapper's `output_path.parent.mkdir(parents=True, exist_ok=True)`.
    std::filesystem::path out_path(args.output);
    if (out_path.has_parent_path())
        std::filesystem::create_directories(out_path.parent_path());
    std::ofstream fout(out_path);
    if (!fout) {
        std::cerr << "Error: cannot open output file: " << args.output << "\n";
        return 1;
    }

    size_t n_written = 0, n_skipped = 0;

    if (args.method == "rawhash") {
        RH1Params p;
        if (args.rh1_w1) p.window_length1 = *args.rh1_w1;
        if (args.rh1_w2) p.window_length2 = *args.rh1_w2;
        if (args.rh1_t1) p.threshold1     = *args.rh1_t1;
        if (args.rh1_t2) p.threshold2     = *args.rh1_t2;
        if (args.rh1_ph) p.peak_height    = *args.rh1_ph;

        for (auto &r : reads) {
            uint32_t n_events = 0;
            float *ev = rh1_detect_events(
                static_cast<uint32_t>(r.signal_pA.size()),
                r.signal_pA.data(), &n_events,
                p.window_length1, p.window_length2,
                p.threshold1, p.threshold2, p.peak_height);
            if (!ev || n_events == 0) { ++n_skipped; if (ev) rh1_free(ev); continue; }
            write_tsv_line(fout, r.read_id, ev, n_events);
            rh1_free(ev);
            ++n_written;
        }
    } else {
        RH2Params p = args.r10 ? rh2_r10_defaults() : RH2Params{};
        if (args.rh2_w1)      p.window_length1     = *args.rh2_w1;
        if (args.rh2_w2)      p.window_length2     = *args.rh2_w2;
        if (args.rh2_t1)      p.threshold1         = *args.rh2_t1;
        if (args.rh2_t2)      p.threshold2         = *args.rh2_t2;
        if (args.rh2_ph)      p.peak_height        = *args.rh2_ph;
        if (args.rh2_min_seg) p.min_segment_length = *args.rh2_min_seg;
        if (args.rh2_max_seg) p.max_segment_length = *args.rh2_max_seg;

        // Cumulative normalization state — each call to rh2_detect_events
        // updates these, and subsequent reads see the running statistics.
        // Matches detect_events_rawhash2() in python-legacy/event_detector.py
        // where the state is re-initialized per-call (so, per-read here to
        // preserve parity).
        for (auto &r : reads) {
            double   mean_sum     = 0.0;
            double   std_dev_sum  = 0.0;
            uint32_t n_events_sum = 0;
            uint32_t n_events     = 0;
            float *ev = rh2_detect_events(
                static_cast<uint32_t>(r.signal_pA.size()),
                r.signal_pA.data(), &n_events,
                p.window_length1, p.window_length2,
                p.threshold1, p.threshold2, p.peak_height,
                p.min_segment_length, p.max_segment_length,
                &mean_sum, &std_dev_sum, &n_events_sum);
            if (!ev || n_events == 0) { ++n_skipped; if (ev) rh2_free(ev); continue; }
            write_tsv_line(fout, r.read_id, ev, n_events);
            rh2_free(ev);
            ++n_written;
        }
    }

    std::cout << "Done. Wrote " << n_written << " read(s) to "
              << args.output << " (" << n_skipped << " skipped)\n";
    return 0;
}
