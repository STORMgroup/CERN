/**
 * rawhash_events.c
 *
 * Adapted from RawHash v1.0's event detection logic.
 *
 * Default parameters (from roptions.c):
 *   window_length1 = 3,  threshold1 = 4.30265
 *   window_length2 = 6,  threshold2 = 2.57058
 *   peak_height    = 1.0
 */

#include "rawhash_events.h"

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

/* ---------- t-statistic (RawHash v1.0 formula) --------------------------- */

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
        /* RawHash v1.0: combined_var is NOT divided by w_len here */
        float combined_var = sumsq1 / w_len - mean1 * mean1
                           + sumsq2 / w_len - mean2 * mean2;
        combined_var = fmaxf(combined_var, eta);
        const float delta_mean = mean2 - mean1;
        /* Division by w_len happens inside the sqrt */
        tstat[i] = fabs(delta_mean) / sqrt(combined_var / w_len);
    }
    memset(tstat + s_len - w_len + 1, 0, w_len * sizeof(float));
    return tstat;
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
} rh_detect_t;

/* ---------- peak detection (identical to RawHash v1.0) -------------------- */

static uint32_t gen_peaks(rh_detect_t *short_det, rh_detect_t *long_det,
                          float peak_height, uint32_t *peaks)
{
    assert(short_det->s_len == long_det->s_len);
    uint32_t cur = 0;
    uint32_t ndetector = 2;
    rh_detect_t *detectors[2] = {short_det, long_det};

    for (uint32_t i = 0; i < short_det->s_len; i++) {
        for (uint32_t k = 0; k < ndetector; k++) {
            rh_detect_t *d = detectors[k];
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
                /* short detector masks long detector */
                if (d == short_det && d->peak_value > d->threshold) {
                    long_det->masked_to = d->peak_pos + d->window_length;
                    long_det->peak_pos = long_det->DEF_PEAK_POS;
                    long_det->peak_value = long_det->DEF_PEAK_VAL;
                    long_det->valid_peak = 0;
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

/* ---------- event generation (RawHash v1.0: simple mean + z-norm) --------- */

static float *gen_events(const uint32_t *peaks, uint32_t peak_size,
                         const float *prefix_sum,
                         uint32_t s_len, uint32_t *n_events)
{
    /* count valid peaks */
    uint32_t n_ev = 1; /* +1 for the last segment */
    for (uint32_t i = 1; i < peak_size; ++i)
        if (peaks[i] > 0 && peaks[i] < s_len) n_ev++;

    if (!n_ev) { *n_events = 0; return NULL; }

    float *events = (float *)malloc(n_ev * sizeof(float));
    double sum = 0, sum2 = 0;
    float l_prefixsum = 0;
    float l_peak = 0;

    for (uint32_t pi = 0; pi < n_ev - 1; pi++) {
        events[pi] = (prefix_sum[peaks[pi]] - l_prefixsum) / (peaks[pi] - l_peak);
        sum += events[pi];
        sum2 += events[pi] * events[pi];
        l_prefixsum = prefix_sum[peaks[pi]];
        l_peak = peaks[pi];
    }
    /* last event: ends at s_len */
    events[n_ev - 1] = (prefix_sum[s_len] - l_prefixsum) / (s_len - l_peak);
    sum += events[n_ev - 1];
    sum2 += events[n_ev - 1] * events[n_ev - 1];

    /* post-hoc z-score normalization of events */
    double mean = sum / n_ev;
    double std_dev = sqrt(sum2 / n_ev - mean * mean);
    for (uint32_t i = 0; i < n_ev; ++i)
        events[i] = (events[i] - mean) / std_dev;

    *n_events = n_ev;
    return events;
}

/* ========================================================================= */
/* Public API                                                                */
/* ========================================================================= */

/**
 * Detect events using the RawHash v1.0 algorithm.
 *
 * @param s_len          Length of the input signal (pA, already 30-200 filtered).
 * @param sig            Input signal array (float32, pA values).
 * @param n_events       [out] Number of events produced.
 * @param window_length1 Short window length (default 3).
 * @param window_length2 Long window length (default 6).
 * @param threshold1     Short detector threshold (default 4.30265).
 * @param threshold2     Long detector threshold (default 2.57058).
 * @param peak_height    Minimum peak rise/fall (default 1.0).
 * @return Pointer to float array of events (caller must free via rh1_free).
 */
float *rh1_detect_events(uint32_t s_len, const float *sig,
                         uint32_t *n_events,
                         uint32_t window_length1, uint32_t window_length2,
                         float threshold1, float threshold2,
                         float peak_height)
{
    *n_events = 0;
    if (s_len == 0) return NULL;

    float *prefix_sum = (float *)calloc(s_len + 1, sizeof(float));
    float *prefix_sum_sq = (float *)calloc(s_len + 1, sizeof(float));
    comp_prefix_prefixsq(sig, s_len, prefix_sum, prefix_sum_sq);

    float *tstat1 = comp_tstat(prefix_sum, prefix_sum_sq, s_len, window_length1);
    float *tstat2 = comp_tstat(prefix_sum, prefix_sum_sq, s_len, window_length2);

    rh_detect_t short_det = {
        .DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = tstat1, .s_len = s_len,
        .threshold = threshold1, .window_length = window_length1,
        .masked_to = 0, .peak_pos = -1, .peak_value = FLT_MAX, .valid_peak = 0
    };
    rh_detect_t long_det = {
        .DEF_PEAK_POS = -1, .DEF_PEAK_VAL = FLT_MAX,
        .sig = tstat2, .s_len = s_len,
        .threshold = threshold2, .window_length = window_length2,
        .masked_to = 0, .peak_pos = -1, .peak_value = FLT_MAX, .valid_peak = 0
    };

    uint32_t *peaks = (uint32_t *)malloc(s_len * sizeof(uint32_t));
    uint32_t n_peaks = gen_peaks(&short_det, &long_det, peak_height, peaks);

    float *events = NULL;
    if (n_peaks > 0)
        events = gen_events(peaks, n_peaks, prefix_sum, s_len, n_events);

    free(tstat1);
    free(tstat2);
    free(prefix_sum);
    free(prefix_sum_sq);
    free(peaks);

    return events;
}

/** Free memory returned by rh1_detect_events. */
void rh1_free(void *ptr) { free(ptr); }
