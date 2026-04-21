#ifndef RAWHASH_EVENTS_H
#define RAWHASH_EVENTS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Detect events using the RawHash v1.0 algorithm (ScrappieR9 compatible).
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
                         float peak_height);

/** Free memory returned by rh1_detect_events. */
void rh1_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif
