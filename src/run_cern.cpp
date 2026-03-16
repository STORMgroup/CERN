#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <chrono>
#include <cstdint>

using namespace std;

struct StayRemoved {
    std::vector<float> events;
    std::vector<float> states;
};

struct HMM {
    vector<float> means;
    vector<float> variances;

    vector<float> log_norm;   // -0.5*(log(2pi)+log(var))
    vector<float> inv2_var;   // 1/(2*var)

    vector<vector<float>> transitions;
    vector<vector<float>> log_transitions;
    vector<vector<int>> to_adj_list;
};

/*
 * Utility: split space-separated line into floats
 */
static vector<float> parse_floats(const string& line) {
    vector<float> values;
    stringstream ss(line);
    float x;
    while (ss >> x) {
        values.push_back(x);
    }
    return values;
}

/*
 * Utility: split space-separated line into ints
 */
static vector<int> parse_ints(const string& line) {
    vector<int> values;
    stringstream ss(line);
    int x;
    while (ss >> x) {
        values.push_back(x);
    }
    return values;
}

/*
 * Load HMM from file
 */
HMM load_hmm(const string& filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        throw runtime_error("Could not open HMM file: " + filename);
    }

    HMM hmm;

    string line;

    // ---- means ----
    if (!getline(fin, line)) {
        throw runtime_error("Failed to read means line");
    }
    hmm.means = parse_floats(line);
    int N = hmm.means.size();

    // ---- variances ----
    if (!getline(fin, line)) {
        throw runtime_error("Failed to read variances line");
    }
    hmm.variances = parse_floats(line);
    if ((int)hmm.variances.size() != N) {
        throw runtime_error("Variances size mismatch");
    }

    // ---- emission precompute ----
    const float LOG_2PI = log(2.0f * M_PI);

    hmm.log_norm.resize(N);
    hmm.inv2_var.resize(N);

    for (int s = 0; s < N; ++s) {
        float var = hmm.variances[s];

        if (var <= 0.0f) {
            throw runtime_error("Variance must be > 0");
        }

        hmm.log_norm[s] = -0.5f * (LOG_2PI + log(var));
        hmm.inv2_var[s] = 1.0f / (2.0f * var);
    }


    // ---- transitions ----
    if (!getline(fin, line)) {
        throw runtime_error("Failed to read transitions line");
    }
    vector<float> flat_trans = parse_floats(line);
    if ((int)flat_trans.size() != N * N) {
        throw runtime_error("Transition matrix size mismatch");
    }

    hmm.transitions.resize(N, vector<float>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            hmm.transitions[i][j] = flat_trans[i * N + j];
        }
    }

    // ---- log transitions ----
    hmm.log_transitions.resize(N, vector<float>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float p = hmm.transitions[i][j];
            if (p > 0.0f)
                hmm.log_transitions[i][j] = log(p);
            else
                hmm.log_transitions[i][j] = -INFINITY;
        }
    }

    // ---- error-free transitions (adjacency) ----
    if (!getline(fin, line)) {
        throw runtime_error("Failed to read error-free transitions line");
    }
    vector<int> valid = parse_ints(line);
    if ((int)valid.size() != N * N) {
        throw runtime_error("Error-free transition matrix size mismatch");
    }

    hmm.to_adj_list.resize(N);

    // If valid transition goes from i -> j, append i to jth list
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int v = valid[i * N + j];
            if (v == 1) {
                hmm.to_adj_list[j].push_back(i);
            }
        }
    }

    fin.close();
    return hmm;
}


vector<int> viterbi_fast(
    const HMM& hmm,
    const vector<float>& events
) {
    int T = events.size();          // number of observations
    int N = hmm.means.size();       // number of states

    vector<int> path;

    // ================================
    // Log-space emission probabilities
    // Gaussian model
    // log N(x | μ, σ²)
    // ================================
    vector<vector<float>> log_emit(T, vector<float>(N));

    for (int t = 0; t < T; ++t) {
        float x = events[t];
        for (int s = 0; s < N; ++s) {
            float diff = x - hmm.means[s];
            log_emit[t][s] =
                hmm.log_norm[s]           // precomputed -0.5*(log(2π)+log(var))
                - diff * diff * hmm.inv2_var[s];  // precomputed 1/(2*var)
        }
    }

    // ================================
    // DP + backtracking matrices
    // ================================
    vector<vector<float>> dp(T, vector<float>(N, -INFINITY));
    vector<vector<int>> backptr(T, vector<int>(N, N));

    // ================================
    // Initialization (t = 0)
    // No start probabilities — only emissions
    // ================================
    int last_best_state = N;
    float last_best_score = -INFINITY;

    for (int s = 0; s < N; ++s) {
        dp[0][s] = log_emit[0][s];
        backptr[0][s] = N;   // no predecessor at t=0

        if (dp[0][s] > last_best_score) {
            last_best_score = dp[0][s];
            last_best_state = s;
        }
    }

    // ================================
    // Dynamic Programming
    // ================================
    for (int t = 1; t < T; ++t) {

        int new_last_best_state = N;
        float new_last_best_score = -INFINITY;

        for (int s = 0; s < N; ++s) {

            float best_score = -INFINITY;
            int best_prev = N;

            // ---- adjacency predecessors ----
            const vector<int>& preds = hmm.to_adj_list[s];
            for (int k = 0; k < (int)preds.size(); ++k) {
                int p = preds[k];

                float score = dp[t-1][p] + hmm.log_transitions[p][s];
                if (score > best_score) {
                    best_score = score;
                    best_prev = p;
                }
            }
            // Also consider self
            float score = dp[t-1][s] + hmm.log_transitions[s][s];
            if (score > best_score) {
                best_score = score;
                best_prev = s; 
            }

            // ---- also consider last_best_state ----
            if (last_best_state < N) {
                int p = last_best_state;
                float score = dp[t-1][p] + hmm.log_transitions[p][s];

                if (score > best_score) {
                    best_score = score;
                    best_prev = p; 
                }
            }

            // ---- emission + transition ----
            dp[t][s] = best_score + log_emit[t][s];
            backptr[t][s] = best_prev;

            // ---- track best state for next timestep ----
            if (dp[t][s] > new_last_best_score) {
                new_last_best_score = dp[t][s];
                new_last_best_state = s;
            }
        }

        // update for next iteration
        last_best_state = new_last_best_state;
    }
    // ================================
    // Backtracking
    // ================================
    path.resize(T);

    // start from best state at final time
    int cur_state = last_best_state;
    path[T-1] = cur_state;

    for (int t = T - 1; t > 0; --t) {

        int prev = backptr[t][cur_state];

        path[t-1] = prev;
        cur_state = prev;
        
    }
    
    return path;
}

std::pair<std::vector<float>, std::vector<int>>
remove_stays(const std::vector<int>& path,
             const std::vector<float>& events,
             const HMM& hmm)
{
    std::vector<float> out_events;
    std::vector<int>   out_states;
    out_events.reserve(path.size());
    out_states.reserve(path.size());

    if (path.empty()) return {out_events, out_states};

    int prev_state = INT32_MIN;
    double accum = 0.0;
    size_t count = 0;

    // initial state push (same as original)
    out_states.push_back(path[0]);

    for (size_t i = 0; i < path.size(); ++i) {
        int state = path[i];
        int abs_state = (state < 0) ? -state : state;

        if (abs_state == prev_state && hmm.transitions[abs_state][abs_state] > 0.001) {
            accum += events[i];
            count++;
        } else {
            if (count > 0) {
                out_events.push_back(accum / count);
                out_states.push_back(state);
            }
            prev_state = abs_state;
            accum = events[i];
            count = 1;
        }
    }

    if (count > 0) {
        out_events.push_back(accum / count);
    }

    return {out_events, out_states};
}

std::vector<float> remove_noise(
    const std::vector<float>& events,
    const std::vector<int>& states,
    const HMM& hmm,
    int window_size)
{
    int N = events.size();
    std::vector<float> corrected_events = events;  // start as copy

    if (N == 0 || window_size == 0) return corrected_events;
    if (window_size > N) return corrected_events;

    // enforce odd window size
    if (window_size % 2 == 0) {
        window_size++;
    }
    int halfW = (window_size - 1) / 2;

    float sum;
    float avg;
    // slide window
    sum = 0;

    // build windowed noise
    for (int j = 0; j < window_size; ++j)
    {
        sum += events[j] - hmm.means[std::abs(states[j])];
    }

    // Slide window
    for (int start = 0; start < N - window_size - 1; ++start)
    {
        avg = sum / window_size;
        int center = start + halfW;
        corrected_events[center] = events[center] - avg;
        sum -= events[start] - hmm.means[std::abs(states[start])];
        sum += events[start + window_size] - hmm.means[std::abs(states[start + window_size])];
    }
    avg = sum / window_size;
    int center = N - window_size + halfW;
    corrected_events[center] = events[center] - avg;

    return corrected_events;
}

std::vector<float> remove_noise_with_stays(
    const std::vector<float>& events,
    const std::vector<int>& states,
    const HMM& hmm,
    int window_size)
{
    int N = events.size();
    auto result = remove_stays(states, events, hmm);
    std::vector<float> stay_removed_events = std::move(result.first);
    std::vector<int> stay_removed_path = std::move(result.second);
    int N_stay_removed = stay_removed_events.size();
    int cur_events_idx = 0;
    int cur_stay_removed_events_idx = 0;


    std::vector<float> corrected_events = events;  // start as copy

    if (N_stay_removed == 0 || window_size == 0) return corrected_events;
    if (window_size > N_stay_removed) return corrected_events;

    // enforce odd window size
    if (window_size % 2 == 0) {
        window_size++;
    }
    int halfW = (window_size - 1) / 2;

    float sum;
    float avg;

    // prepare the indexes
    while (cur_stay_removed_events_idx < halfW) {
        while (cur_events_idx < N && std::abs(states[cur_events_idx]) == std::abs(stay_removed_path[cur_stay_removed_events_idx])) {
            cur_events_idx++;
        }
        cur_stay_removed_events_idx++;
    }

    // slide window
    for (int start = 0; start < N_stay_removed - window_size; ++start)
    {
        sum = 0;

        // build windowed noise
        for (int j = 0; j < window_size; ++j)
        {
            sum += stay_removed_events[start + j] - hmm.means[std::abs(stay_removed_path[start + j])];
        }

        avg = sum / window_size;

        int center = start + halfW;


        while (cur_events_idx < N && std::abs(stay_removed_path[center]) == std::abs(states[cur_events_idx])) {
            corrected_events[cur_events_idx] = events[cur_events_idx] - avg;
            cur_events_idx++;
        }

    }

    return corrected_events;
}


int main(int argc, char* argv[]) {

    auto start = std::chrono::steady_clock::now();

    // Command line args

    // ---- defaults for optional args ----
    int window_size = 20;
    bool removing_noise = true;
    bool removing_stays = true;

    // ---- required positional args ----
    if (argc < 3) {
        cerr << "Usage: " << argv[0]
             << " <hmm_file> <event_file> "
             << "[--window-size N] "
             << "[--no-noise-removal] "
             << "[--no-stay-removal]\n";
        return 1;
    }

    string hmm_filename = argv[1];
    string event_filename = argv[2];

    // ---- parse optional flags ----
    for (int i = 3; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "--window-size") {
            if (i + 1 >= argc) {
                cerr << "Error: --window-size requires a value\n";
                return 1;
            }
            window_size = stoi(argv[++i]);
        }
        else if (arg == "--no-noise-removal") {
            removing_noise = false;
        }
        else if (arg == "--no-stay-removal") {
            removing_stays = false;
        }
        else {
            cerr << "Warning: unknown argument ignored: " << arg << "\n";
        }
    }

    // ---- load HMM ----
    HMM hmm = load_hmm(hmm_filename);

    // ---- open event file ----
    ifstream fin(event_filename);
    if (!fin.is_open()) {
        cerr << "Could not open event file: " << event_filename << "\n";
        return 1;
    }

    string line;

    int lines = 0;
    // ---- process each read ----
    while (getline(fin, line)) {
        if (line.empty()) continue;
        lines++;

        stringstream ss(line);

        // first block = read ID
        string read_id;
        ss >> read_id;

        // remaining blocks = events
        vector<float> events;
        float val;
        while (ss >> val) {
            events.push_back(val);
        }


        // ---- run viterbi ----
        vector<int> path = viterbi_fast(hmm, events);

        
        // ---- remove stays ----
        if (removing_stays) {
            auto result = remove_stays(path, events, hmm);
            events = std::move(result.first);
            path = std::move(result.second);
        } 


        // ---- remove noise ----
        if (removing_noise) {
            if(removing_stays) {
                events =
                remove_noise(events,
                                path,
                                hmm,
                                window_size);
            } else {
                events =
                remove_noise_with_stays(events,
                                path,
                                hmm,
                                window_size);
            }
        } 

        // ---- output ----
        cout << read_id;
        for (float v : events) {
            cout << "\t" << v;
        }
        cout << "\n";
    }

    fin.close();

    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> duration = end - start;

    float time_per_read = duration.count() / lines;
    float total_seconds = duration.count() / 1000.0f;

    std::cerr << "Time per read: " << time_per_read << " ms\n";
    std::cerr << "Total time: " << total_seconds << " s\n";


    return 0;
}
