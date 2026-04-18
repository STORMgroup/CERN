#ifndef RAWHASH2_EVENTS_H
#define RAWHASH2_EVENTS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Detect events using the RawHash2 algorithm (ScrappieR10 compatible).
 *
 * mean_sum / std_dev_sum / n_events_sum are in/out accumulators for
 * cumulative z-score normalization across reads. Initialize each to 0
 * before the first call and pass the same three pointers across all
 * reads belonging to the same normalization cohort.
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
                         uint32_t *n_events_sum);

/** Free memory returned by rh2_detect_events. */
void rh2_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif
