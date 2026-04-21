/**
 * rawhash2_events.c
 *
 * Adapted from RawHash2's event detection logic.
 *
 * Default parameters (from roptions.c):
 *   window_length1 = 3,  threshold1 = 4.0
 *   window_length2 = 9,  threshold2 = 3.5
 *   peak_height    = 0.4
 *   min_segment_length = 0,  max_segment_length = 500
 *
 * Compiled as a shared library for Python ctypes FFI.
 */

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ---------- prefix sums -------------------------------------------------- */

static void comp_prefix_prefixsq(const float *sig, uint32_t s_len,
                                  float *prefix_sum, float *prefix_sum_square)
{
    assert(s_len > 0);
    prefix_sum[0] = 0.0f;
    prefix_sum_square[0] = 0.0f;
    for (uint32_t i = 0; i < s_len; ++i) {
        prefix_sum[i + 1] = prefix_sum[i] + sig[i];
        prefix_sum_square[i + 1] = prefix_sum_square[i] + sig[i] * sig[i];
    }
}

/* ---------- t-statistic (RawHash2 formula) ------------------------------- */

static float *comp_tstat(const float *prefix_sum,
                         const float *prefix_sum_square,
                         uint32_t s_len, uint32_t w_len)
{
    const float eta = FLT_MIN;
    float *tstat = (float *)calloc(s_len + 1, sizeof(float));
    if (s_len < 2 * w_len || w_len < 2) return tstat;
    memset(tstat, 0, w_len * sizeof(float));

    for (uint32_t i = w_len; i <= s_len - w_len; ++i) {
        float sum1 = prefix_sum[i];
        float sumsq1 = prefix_sum_square[i];
        if (i > w_len) {
            sum1 -= prefix_sum[i - w_len];
            sumsq1 -= prefix_sum_square[i - w_len];
        }
        float sum2 = prefix_sum[i + w_len] - prefix_sum[i];
        float sumsq2 = prefix_sum_square[i + w_len] - prefix_sum_square[i];
        float mean1 = sum1 / w_len;
        float mean2 = sum2 / w_len;
        /* RawHash2: combined_var is divided by w_len here */
        float combined_var = (sumsq1 / w_len - mean1 * mean1
                            + sumsq2 / w_len - mean2 * mean2) / w_len;
        combined_var = fmaxf(combined_var, eta);
        const float delta_mean = mean2 - mean1;
        tstat[i] = fabs(delta_mean) / sqrt(combined_var);
    }
    memset(tstat + s_len - w_len + 1, 0, w_len * sizeof(float));
    return tstat;
}

/* ---------- signal normalization (z-score + 3-sigma filter) -------------- */

/**
 * Normalize signal: z-score then keep only values within +/-3 sigma.
 * Supports cumulative statistics across reads via mean_sum/std_dev_sum/n_events_sum.
 */
static float *normalize_signal(const float *sig, uint32_t s_len,
                               double *mean_sum, double *std_dev_sum,
                               uint32_t *n_events_sum, uint32_t *n_sig)
{
    double sum = *mean_sum, sum2 = *std_dev_sum;
    for (uint32_t i = 0; i < s_len; ++i) {
        sum += sig[i];
        sum2 += sig[i] * sig[i];
    }

    *n_events_sum += s_len;
    *mean_sum = sum;
    *std_dev_sum = sum2;

    double mean = sum / (*n_events_sum);
    double std_dev = sqrt(sum2 / (*n_events_sum) - mean * mean);

    float *events = (float *)calloc(s_len, sizeof(float));
    int k = 0;
    for (uint32_t i = 0; i < s_len; ++i) {
        float norm_val = (float)((sig[i] - mean) / std_dev);
        if (norm_val < 3.0f && norm_val > -3.0f)
            events[k++] = norm_val;
    }

    *n_sig = (uint32_t)k;
    return events;
}

/* ---------- detector struct ----------------------------------------------- */

typedef struct {
    int   DEF_PEAK_POS;
    float DEF_PEAK_VAL;
    float *sig;
    uint32_t s_len;
    float threshold;
    uint32_t window_length;
    uint32_t masked_to;
    int   peak_pos;
    float peak_value;
    int   valid_peak;
} rh2_detect_t;

/* ---------- peak detection (RawHash2 style: n_detectors loop) ------------ */

static uint32_t gen_peaks(rh2_detect_t **detectors, uint32_t n_detectors,
                          float peak_height, uint32_t *peaks)
{
    uint32_t cur = 0;
    for (uint32_t i = 0; i < detectors[0]->s_len; i++) {
        for (uint32_t k = 0; k < n_detectors; k++) {
            rh2_detect_t *d = detectors[k];
            if (d->masked_to >= i) continue;
            float cv = d->sig[i];
            if (d->peak_pos == d->DEF_PEAK_POS) {
                if (cv < d->peak_value) {
                    d->peak_value = cv;
                } else if (cv - d->peak_value > peak_height) {
                    d->peak_value = cv;
                    d->peak_pos = i;
                }
            } else {
                if (cv > d->peak_value) {
                    d->peak_value = cv;
                    d->peak_pos = i;
                }
                if (d->peak_value > d->threshold) {
                    for (uint32_t n_d = k + 1; n_d < n_detectors; n_d++) {
                        detectors[n_d]->masked_to = d->peak_pos + detectors[0]->window_length;
                        detectors[n_d]->peak_pos = detectors[n_d]->DEF_PEAK_POS;
                        detectors[n_d]->peak_value = detectors[n_d]->DEF_PEAK_VAL;
                        detectors[n_d]->valid_peak = 0;
                    }
                }
                if (d->peak_value - cv > peak_height &&
                    d->peak_value > d->threshold) {
                    d->valid_peak = 1;
                }
                if (d->valid_peak &&
                    (i - (uint32_t)d->peak_pos) > d->window_length / 2) {
                    peaks[cur++] = d->peak_pos;
                    d->peak_pos = d->DEF_PEAK_POS;
                    d->peak_value = cv;
                    d->valid_peak = 0;
                }
            }
        }
    }
    return cur;
}

/* ---------- IQR-filtered mean -------------------------------------------- */

static int compare_floats(const void *a, const void *b)
{
    const float *da = (const float *)a;
    const float *db = (const float *)b;
    return (*da > *db) - (*da < *db);
}

static float calculate_mean_of_filtered_segment(float *segment,
                                                 uint32_t segment_length)
{
    qsort(segment, segment_length, sizeof(float), compare_floats);
    float q1 = segment[segment_length / 4];
    float q3 = segment[3 * segment_length / 4];
    float iqr = q3 - q1;
    float lower_bound = q1 - iqr;
    float upper_bound = q3 + iqr;

    float sum = 0.0f;
    uint32_t count = 0;
    for (uint32_t i = 0; i < segment_length; i++) {
        if (segment[i] >= lower_bound && segment[i] <= upper_bound) {
            sum += segment[i];
            ++count;
        }
    }
    return count > 0 ? sum / count : 0.0f;
}

/* ---------- event generation (RawHash2: IQR-filtered + segment filter) --- */

static float *gen_events(float *sig, const uint32_t *peaks, uint32_t peak_size,
                         uint32_t s_len, uint32_t min_seg_len,
                         uint32_t max_seg_len, uint32_t *n_events)
{
    uint32_t n_ev = 0;
    for (uint32_t pi = 0; pi < peak_size; ++pi)
        if (peaks[pi] > 0 && peaks[pi] < s_len) n_ev++;

    if (!n_ev) { *n_events = 0; return NULL; }

    float *events = (float *)malloc(n_ev * sizeof(float));
    uint32_t start_idx = 0, n_written = 0;

    for (uint32_t pi = 0; pi < peak_size && n_written < n_ev; pi++) {
        if (!(peaks[pi] > 0 && peaks[pi] < s_len)) continue;
        uint32_t segment_length = peaks[pi] - start_idx;
        if (segment_length >= min_seg_len && segment_length < max_seg_len)
            events[n_written++] = calculate_mean_of_filtered_segment(
                sig + start_idx, segment_length);
        start_idx = peaks[pi];
    }

    if (n_written == 0) {
        free(events);
        *n_events = 0;
        return NULL;
    } else if (n_written < n_ev) {
        events = (float *)realloc(events, n_written * sizeof(float));
    }

    *n_events = n_written;
    return events;
}

/* ========================================================================= */
/* Public API                                                                */
/* ========================================================================= */

/**
 * Detect events using the RawHash2 algorithm.
 *
 * @param s_len            Length of the input signal (pA, already 30-200 filtered).
 * @param sig              Input signal array (float32, pA values).
 * @param n_events         [out] Number of events produced.
 * @param window_length1   Short window length (default 3).
 * @param window_length2   Long window length (default 9).
 * @param threshold1       Short detector threshold (default 4.0).
 * @param threshold2       Long detector threshold (default 3.5).
 * @param peak_height      Minimum peak rise/fall (default 0.4).
 * @param min_seg_len      Minimum segment length (default 0).
 * @param max_seg_len      Maximum segment length (default 500).
 * @param mean_sum         [in/out] Cumulative sum for normalization (init 0).
 * @param std_dev_sum      [in/out] Cumulative sum-of-squares (init 0).
 * @param n_events_sum     [in/out] Cumulative sample count (init 0).
 * @return Pointer to float array of events (caller must free via rh2_free).
 */
float *rh2_detect_events(uint32_t s_len, const float *sig,
                         uint32_t *n_events,
                         uint32_t window_length1, uint32_t window_length2,
                         float threshold1, float threshold2,
                         float peak_height,
                         uint32_t min_seg_len, uint32_t max_seg_len,
                         double *mean_sum, double *std_dev_sum,
                         uint32_t *n_events_sum)
{
    *n_events = 0;
    if (s_len == 0) return NULL;

    /* Step 1: Normalize signal (z-score + 3-sigma filter) */
    uint32_t n_signals = 0;
    float *norm_signals = normalize_signal(sig, s_len, mean_sum, std_dev_sum,
                                           n_events_sum, &n_signals);
    if (n_signals == 0) { free(norm_signals); return NULL; }

    /* Step 2: Prefix sums on normalized signal */
    float *prefix_sum = (float *)calloc(n_signals + 1, sizeof(float));
    float *prefix_sum_sq = (float *)calloc(n_signals + 1, sizeof(float));
    comp_prefix_prefixsq(norm_signals, n_signals, prefix_sum, prefix_sum_sq);

    /* Step 3: Compute t-statistics */
    float *tstat1 = comp_tstat(prefix_sum, prefix_sum_sq, n_signals, window_length1);
    float *tstat2 = comp_tstat(prefix_sum, prefix_sum_sq, n_signals, window_length2);

    rh2_detect_t short_det = {
        .DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = tstat1, .s_len = n_signals,
        .threshold = threshold1, .window_length = window_length1,
        .masked_to = 0, .peak_pos = -1, .peak_value = FLT_MAX, .valid_peak = 0
    };
    rh2_detect_t long_det = {
        .DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = tstat2, .s_len = n_signals,
        .threshold = threshold2, .window_length = window_length2,
        .masked_to = 0, .peak_pos = -1, .peak_value = FLT_MAX, .valid_peak = 0
    };

    /* Step 4: Peak detection */
    uint32_t *peaks = (uint32_t *)malloc(n_signals * sizeof(uint32_t));
    rh2_detect_t *detectors[2] = {&short_det, &long_det};
    uint32_t n_peaks = gen_peaks(detectors, 2, peak_height, peaks);

    /* Step 5: Generate events from peaks */
    float *events = NULL;
    if (n_peaks > 0)
        events = gen_events(norm_signals, peaks, n_peaks, n_signals,
                            min_seg_len, max_seg_len, n_events);

    free(tstat1);
    free(tstat2);
    free(prefix_sum);
    free(prefix_sum_sq);
    free(norm_signals);
    free(peaks);

    return events;
}

/** Free memory returned by rh2_detect_events. */
void rh2_free(void *ptr) { free(ptr); }
