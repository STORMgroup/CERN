#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define M_PI 3.14159265358979323846

/* ------------------------------------------------------------------ */
/* De Bruijn sequence via Eulerian path                               */
/* ------------------------------------------------------------------ */

static const char BASES[4] = {'A', 'C', 'G', 'T'};

static void decode_kmer(int val, int k, char *buf)
{
    buf[k] = '\0';
    for (int i = k - 1; i >= 0; i--) {
        buf[i] = BASES[val & 3];
        val >>= 2;
    }
}

static char *debruijn(int k)
{
    if (k == 1) {
        char *seq = malloc(5);
        if (!seq) { fprintf(stderr, "OOM\n"); exit(1); }
        seq[0]='A'; seq[1]='C'; seq[2]='G'; seq[3]='T'; seq[4]='\0';
        return seq;
    }
    int node_count = 1 << (2 * (k - 1));
    int edge_count = node_count * 4;
    int *ptr     = calloc(node_count, sizeof(int));
    int *stack   = malloc((edge_count + 1) * sizeof(int));
    int *circuit = malloc((edge_count + 1) * sizeof(int));
    if (!ptr || !stack || !circuit) { fprintf(stderr, "OOM\n"); exit(1); }
    int stack_top = 0, circuit_top = 0;
    stack[stack_top++] = 0;
    while (stack_top > 0) {
        int v = stack[stack_top - 1];
        if (ptr[v] < 4) {
            int b = ptr[v]++;
            int u = ((v << 2) | b) & (node_count - 1);
            stack[stack_top++] = u;
        } else {
            circuit[circuit_top++] = v;
            stack_top--;
        }
    }
    int seq_len = edge_count + k - 1;
    char *seq = malloc((seq_len + 1) * sizeof(char));
    if (!seq) { fprintf(stderr, "OOM\n"); exit(1); }
    char kbuf[64];
    decode_kmer(circuit[circuit_top - 1], k - 1, kbuf);
    memcpy(seq, kbuf, k - 1);
    for (int i = k - 1; i < seq_len; i++) {
        int node = circuit[circuit_top - 2 - (i - (k - 1))];
        seq[i] = BASES[node & 3];
    }
    seq[seq_len] = '\0';
    free(ptr); free(stack); free(circuit);
    return seq;
}

/* ------------------------------------------------------------------ */
/* K-mer table (open-addressing, FNV-1a)                             */
/* ------------------------------------------------------------------ */

typedef struct { char *kmer; double value; } KmerEntry;
typedef struct { KmerEntry *entries; int size, capacity, k; } KmerTable;

static inline unsigned int hash_str(const char *s, int capacity)
{
    unsigned int h = 2166136261u;
    while (*s) { h ^= (unsigned char)(*s++); h *= 16777619u; }
    return h & (unsigned int)(capacity - 1);
}

static KmerTable *kmer_table_new(int cap)
{
    int p = 1; while (p < cap) p <<= 1;
    KmerTable *t = malloc(sizeof(KmerTable));
    t->entries = calloc(p, sizeof(KmerEntry));
    t->size = 0; t->capacity = p; t->k = 0;
    return t;
}

static void kmer_table_insert(KmerTable *t, const char *kmer, double val)
{
    if (t->size * 10 >= t->capacity * 7) {
        int new_cap = t->capacity * 2;
        KmerEntry *ne = calloc(new_cap, sizeof(KmerEntry));
        for (int i = 0; i < t->capacity; i++) {
            if (t->entries[i].kmer) {
                unsigned int h = hash_str(t->entries[i].kmer, new_cap);
                while (ne[h].kmer) h = (h + 1) & (new_cap - 1);
                ne[h] = t->entries[i];
            }
        }
        free(t->entries); t->entries = ne; t->capacity = new_cap;
    }
    unsigned int h = hash_str(kmer, t->capacity);
    while (t->entries[h].kmer && strcmp(t->entries[h].kmer, kmer) != 0)
        h = (h + 1) & (t->capacity - 1);
    if (!t->entries[h].kmer) {
        t->entries[h].kmer = strdup(kmer);
        t->entries[h].value = val;
        t->size++;
    } else {
        t->entries[h].value = val;
    }
}

static inline double kmer_table_lookup(const KmerTable *t, const char *kmer)
{
    unsigned int h = hash_str(kmer, t->capacity);
    while (t->entries[h].kmer) {
        if (strcmp(t->entries[h].kmer, kmer) == 0) return t->entries[h].value;
        h = (h + 1) & (t->capacity - 1);
    }
    return NAN;
}

static KmerTable *load_kmer_table(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (!f) { fprintf(stderr, "Error: cannot open '%s'\n", filename); exit(1); }
    KmerTable *t = kmer_table_new(1024);
    char line[256]; int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        if (line[0] == '#' || line[0] == '\n') continue;
        char kmer[64]; double val;
        if (sscanf(line, "%63s %lf", kmer, &val) != 2) {
            fprintf(stderr, "Warning: bad line %d\n", lineno); continue;
        }
        if (t->k == 0) t->k = (int)strlen(kmer);
        kmer_table_insert(t, kmer, val);
    }
    fclose(f);
    fprintf(stderr, "Loaded %d k-mers (k=%d)\n", t->size, t->k);
    return t;
}

/* ------------------------------------------------------------------ */
/* GMM: k-means++ seeding + EM with two-pass variance                */
/* ------------------------------------------------------------------ */

#define GMM_MAX_ITER    500
#define GMM_CONV_THRESH 1e-6

static inline unsigned long lcg_next(unsigned long *s)
{
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return *s;
}
static inline double lcg_double(unsigned long *s)
{
    return (double)(lcg_next(s) >> 11) / (double)(1ULL << 53);
}

static double *extract_kmer_values(const KmerTable *t, int *out_m)
{
    double *vals = malloc(t->size * sizeof(double));
    if (!vals) { fprintf(stderr, "OOM\n"); exit(1); }
    int m = 0;
    for (int i = 0; i < t->capacity; i++)
        if (t->entries[i].kmer) vals[m++] = t->entries[i].value;
    *out_m = m;
    return vals;
}

static void kmeans_pp_seed(const double *vals, int M, int N,
                           double *centres, unsigned long *rng)
{
    centres[0] = vals[(int)(lcg_double(rng) * M) % M];
    double *dist2 = malloc(M * sizeof(double));
    if (!dist2) { fprintf(stderr, "OOM\n"); exit(1); }
    for (int c = 1; c < N; c++) {
        double total = 0.0;
        for (int i = 0; i < M; i++) {
            double best = INFINITY;
            for (int k = 0; k < c; k++) {
                double d = vals[i] - centres[k];
                double d2 = d * d;
                if (d2 < best) best = d2;
            }
            dist2[i] = best;
            total   += best;
        }
        double r = lcg_double(rng) * total, cum = 0.0;
        int chosen = M - 1;
        for (int i = 0; i < M; i++) {
            cum += dist2[i];
            if (cum >= r) { chosen = i; break; }
        }
        centres[c] = vals[chosen];
    }
    free(dist2);
}

static void gmm_fit(const double *vals, int M, int N,
                    double *means, double *vars, double *weights,
                    unsigned long *rng)
{
    double *resp = malloc((size_t)M * N * sizeof(double));
    if (!resp) { fprintf(stderr, "OOM\n"); exit(1); }

    kmeans_pp_seed(vals, M, N, means, rng);

    double global_mean = 0.0;
    for (int i = 0; i < M; i++) global_mean += vals[i];
    global_mean /= M;
    double global_var = 0.0;
    for (int i = 0; i < M; i++) {
        double d = vals[i] - global_mean; global_var += d * d;
    }
    global_var /= M;
    double init_var = (global_var > 0.0 ? global_var : 1.0) / N;
    if (init_var < 1e-4) init_var = 1e-4;
    for (int k = 0; k < N; k++) { vars[k] = init_var; weights[k] = 1.0 / N; }

    double prev_ll = -INFINITY;
    for (int iter = 0; iter < GMM_MAX_ITER; iter++) {
        /* E-step */
        #pragma omp parallel for schedule(static)
        for (int m = 0; m < M; m++) {
            double x = vals[m], row_sum = 0.0;
            for (int k = 0; k < N; k++) {
                double var = vars[k];
                double r = weights[k] / sqrt(2.0*M_PI*var)
                           * exp(-0.5*(x-means[k])*(x-means[k])/var);
                resp[m*N+k] = r;
                row_sum += r;
            }
            double inv = (row_sum > 0.0) ? 1.0/row_sum : 1.0/N;
            for (int k = 0; k < N; k++)
                resp[m*N+k] = (row_sum > 0.0) ? resp[m*N+k]*inv : inv;
        }
        /* M-step: two-pass variance (avoids catastrophic cancellation) */
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < N; k++) {
            double sk = 0.0, sx = 0.0;
            for (int m = 0; m < M; m++) { double r = resp[m*N+k]; sk += r; sx += r*vals[m]; }
            weights[k] = sk;
            double mu = (sk > 0.0) ? sx/sk : means[k];
            means[k] = mu;
            double svar = 0.0;
            for (int m = 0; m < M; m++) {
                double d = vals[m] - mu; svar += resp[m*N+k] * d * d;
            }
            double var = (sk > 0.0) ? svar/sk : vars[k];
            vars[k] = (var < 1e-6) ? 1e-6 : var;
        }
        double wsum = 0.0;
        for (int k = 0; k < N; k++) wsum += weights[k];
        if (wsum > 0.0) {
            double inv = 1.0/wsum;
            for (int k = 0; k < N; k++) weights[k] *= inv;
        }
        /* log-likelihood */
        double ll = 0.0;
        #pragma omp parallel for reduction(+:ll) schedule(static)
        for (int m = 0; m < M; m++) {
            double x = vals[m], mix = 0.0;
            for (int k = 0; k < N; k++) {
                double var = vars[k];
                mix += weights[k]/sqrt(2.0*M_PI*var)
                       * exp(-0.5*(x-means[k])*(x-means[k])/var);
            }
            ll += (mix > 0.0) ? log(mix) : -700.0;
        }
        double delta = ll - prev_ll;
        if (iter > 0 && fabs(delta) < GMM_CONV_THRESH * fabs(ll)) break;
        prev_ll = ll;
    }
    free(resp);
}

static void gmm_init(const KmerTable *table, int N, double *means, double *vars)
{
    int M;
    double *vals    = extract_kmer_values(table, &M);
    double *weights = malloc(N * sizeof(double));
    if (!weights) { fprintf(stderr, "OOM\n"); exit(1); }
    if (M < N)
        fprintf(stderr, "Warning: fewer k-mers (%d) than states (%d)\n", M, N);
    fprintf(stderr, "Fitting GMM with %d components to %d k-mer values...\n", N, M);
    unsigned long rng = (unsigned long)time(NULL) ^ 0xdeadbeefcafeULL;
    gmm_fit(vals, M, N, means, vars, weights, &rng);
    /* Sort by mean for natural state ordering */
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (means[j] < means[i]) {
                double tmp;
                tmp = means[i]; means[i] = means[j]; means[j] = tmp;
                tmp = vars[i];  vars[i]  = vars[j];  vars[j]  = tmp;
            }
    fprintf(stderr, "GMM init complete. Mean range: [%.4f, %.4f]\n",
            means[0], means[N-1]);
    free(vals); free(weights);
}

/* ------------------------------------------------------------------ */
/* Single-chunk forward-backward                                      */
/*                                                                    */
/* Runs one complete Baum-Welch E-step on obs[0..T-1] and            */
/* accumulates sufficient statistics into the provided arrays.        */
/* All arrays are private to the calling thread — no locking needed. */
/*                                                                    */
/* Memory passed in (all caller-allocated, size noted):              */
/*   alpha [T*N]   forward variables (scratch)                        */
/*   beta  [2*N]   two ping-pong rows for backward (scratch)          */
/*   emit  [T*N]   emission probabilities (scratch)                   */
/*   scale [T]     forward scaling factors (scratch)                  */
/*   xi    [N*N]   transition sufficient statistics  (accumulated)    */
/*   gsum  [N]     gamma sum                         (accumulated)    */
/*   gx    [N]     gamma-weighted obs sum            (accumulated)    */
/*   gvar  [N]     gamma-weighted variance sum       (accumulated)    */
/*                 (gvar is filled on pass 2 only, needs means[])     */
/*                                                                    */
/* Returns the log-likelihood contribution of this chunk.            */
/* ------------------------------------------------------------------ */
static double chunk_forward_backward(
    int T, int N,
    const double *obs,
    const double *means,
    const double *vars,
    const double *trans,
    double *alpha, double *beta, double *emit, double *scale,
    double *xi, double *gsum, double *gx, double *gvar,
    int do_gvar)
{
    /* ---- Emission probabilities ---- */
    for (int i = 0; i < N; i++) {
        double mu    = means[i];
        double var   = vars[i];
        double coeff = 1.0 / sqrt(2.0 * M_PI * var);
        double inv2v = -0.5 / var;
        for (int t = 0; t < T; t++) {
            double d = obs[t] - mu;
            emit[t*N + i] = coeff * exp(d * d * inv2v);
        }
    }

    /* ---- Forward pass with scaling ---- */
    /* t = 0: uniform start */
    double s = 0.0;
    double inv_N = 1.0 / N;
    for (int i = 0; i < N; i++) {
        alpha[i] = inv_N * emit[i];
        s += alpha[i];
    }
    scale[0] = (s == 0.0) ? 1e-300 : s;
    double inv = 1.0 / scale[0];
    for (int i = 0; i < N; i++) alpha[i] *= inv;

    for (int t = 1; t < T; t++) {
        const double *a_prev = alpha + (t-1)*N;
        double       *a_cur  = alpha + t*N;
        const double *e_cur  = emit  + t*N;
        s = 0.0;
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int i = 0; i < N; i++)
                sum += a_prev[i] * trans[i*N + j];
            a_cur[j] = sum * e_cur[j];
            s += a_cur[j];
        }
        scale[t] = (s == 0.0) ? 1e-300 : s;
        inv = 1.0 / scale[t];
        for (int j = 0; j < N; j++) a_cur[j] *= inv;
    }

    double loglik = 0.0;
    for (int t = 0; t < T; t++) loglik += log(scale[t]);

    /* ---- Backward pass + accumulate statistics ---- */
    /*
     * We use two ping-pong rows for beta (beta[cur*N] and beta[next*N])
     * rather than storing the full T*N backward matrix.
     * At each step we have alpha[t] available (full matrix kept) and
     * beta[t] freshly computed, so we can immediately accumulate gamma
     * and xi without storing them.
     */
    int cur = 0, nxt = 1;

    /* beta[T-1] = 1/scale[T-1] */
    inv = 1.0 / scale[T-1];
    for (int i = 0; i < N; i++) beta[cur*N + i] = inv;

    /* gamma at T-1 */
    {
        const double *a   = alpha + (T-1)*N;
        const double *bc  = beta  + cur*N;
        double norm = 0.0;
        for (int i = 0; i < N; i++) norm += a[i] * bc[i];
        if (norm == 0.0) norm = 1e-300;
        double inv_norm = 1.0 / norm;
        double ot = obs[T-1];
        if (!do_gvar) {
            for (int i = 0; i < N; i++) {
                double g = a[i] * bc[i] * inv_norm;
                gsum[i] += g;
                gx[i]   += g * ot;
            }
        } else {
            for (int i = 0; i < N; i++) {
                double g = a[i] * bc[i] * inv_norm;
                double d = ot - means[i];
                gvar[i] += g * d * d;
            }
        }
    }

    for (int t = T-2; t >= 0; t--) {
        const double *e_next  = emit  + (t+1)*N;
        const double *a_cur_p = alpha + t*N;
        const double *bc      = beta  + cur*N;   /* beta[t+1] */
        double       *bn      = beta  + nxt*N;   /* will hold beta[t] */
        double inv_c = 1.0 / scale[t];

        /* Compute beta[t] */
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            const double *tr_row = trans + i*N;
            for (int j = 0; j < N; j++)
                sum += tr_row[j] * e_next[j] * bc[j];
            bn[i] = sum * inv_c;
        }

        /* xi: uses alpha[t] and beta[t+1]=bc */
        if (!do_gvar) {
            for (int i = 0; i < N; i++) {
                const double *tr_row = trans + i*N;
                double        ai     = a_cur_p[i];
                double       *xi_row = xi + i*N;
                for (int j = 0; j < N; j++)
                    xi_row[j] += ai * tr_row[j] * e_next[j] * bc[j];
            }
        }

        /* gamma at t: uses alpha[t] and beta[t]=bn */
        double norm = 0.0;
        for (int i = 0; i < N; i++) norm += a_cur_p[i] * bn[i];
        if (norm == 0.0) norm = 1e-300;
        double inv_norm = 1.0 / norm;
        double ot = obs[t];

        if (!do_gvar) {
            for (int i = 0; i < N; i++) {
                double g = a_cur_p[i] * bn[i] * inv_norm;
                gsum[i] += g;
                gx[i]   += g * ot;
            }
        } else {
            for (int i = 0; i < N; i++) {
                double g = a_cur_p[i] * bn[i] * inv_norm;
                double d = ot - means[i];
                gvar[i] += g * d * d;
            }
        }

        /* ping-pong */
        { int tmp = cur; cur = nxt; nxt = tmp; }
    }

    return loglik;
}

/* ------------------------------------------------------------------ */
/* Baum-Welch: chunk-parallel                                         */
/*                                                                    */
/* The observation sequence is split into n_threads contiguous        */
/* chunks.  Each thread runs a complete serial forward-backward on    */
/* its chunk, with uniform start at the chunk boundary.              */
/*                                                                    */
/* This eliminates all inter-thread communication during the E-step.  */
/* The only synchronisation is a single implicit barrier at the end   */
/* of the omp parallel for, after which the main thread reduces the   */
/* per-thread statistics and does the M-step serially.               */
/*                                                                    */
/* Memory layout: each thread gets its own private alpha/beta/emit/   */
/* scale/xi/gsum/gx/gvar arrays sized to its chunk.  All are         */
/* allocated up front in one block per thread.                        */
/* ------------------------------------------------------------------ */
static double baum_welch_full(
    int N, double *means, double *vars, double *trans,
    const double *obs, int T, double thresh, int max_iter,
    int n_threads)
{
    /* Split T observations into n_threads contiguous chunks */
    int *chunk_start = malloc(n_threads * sizeof(int));
    int *chunk_len   = malloc(n_threads * sizeof(int));
    if (!chunk_start || !chunk_len) { fprintf(stderr, "OOM\n"); exit(1); }
    {
        int base = T / n_threads, rem = T % n_threads, off = 0;
        for (int p = 0; p < n_threads; p++) {
            chunk_start[p] = off;
            chunk_len[p]   = base + (p < rem ? 1 : 0);
            off += chunk_len[p];
        }
    }

    /*
     * Per-thread private arrays.  Each is allocated to its thread's
     * actual chunk length — no max_chunk needed, no offset arithmetic.
     */
    double **t_alpha = malloc(n_threads * sizeof(double *));
    double **t_beta  = malloc(n_threads * sizeof(double *));
    double **t_emit  = malloc(n_threads * sizeof(double *));
    double **t_scale = malloc(n_threads * sizeof(double *));
    double **t_xi    = malloc(n_threads * sizeof(double *));
    double **t_gsum  = malloc(n_threads * sizeof(double *));
    double **t_gx    = malloc(n_threads * sizeof(double *));
    double **t_gv    = malloc(n_threads * sizeof(double *));
    if (!t_alpha||!t_beta||!t_emit||!t_scale||!t_xi||!t_gsum||!t_gx||!t_gv)
        { fprintf(stderr, "OOM\n"); exit(1); }
    for (int p = 0; p < n_threads; p++) {
        int Tp = chunk_len[p];
        t_alpha[p] = malloc((size_t)Tp * N * sizeof(double));
        t_beta[p]  = malloc(2 * N * sizeof(double));
        t_emit[p]  = malloc((size_t)Tp * N * sizeof(double));
        t_scale[p] = malloc(Tp * sizeof(double));
        t_xi[p]    = malloc((size_t)N * N * sizeof(double));
        t_gsum[p]  = malloc(N * sizeof(double));
        t_gx[p]    = malloc(N * sizeof(double));
        t_gv[p]    = malloc(N * sizeof(double));
        if (!t_alpha[p]||!t_beta[p]||!t_emit[p]||!t_scale[p]||
            !t_xi[p]||!t_gsum[p]||!t_gx[p]||!t_gv[p])
            { fprintf(stderr, "OOM\n"); exit(1); }
    }

    /* Shared reduction targets */
    double *xi   = malloc((size_t)N * N * sizeof(double));
    double *gsum = malloc(N * sizeof(double));
    double *gx   = malloc(N * sizeof(double));
    double *gvar = malloc(N * sizeof(double));
    if (!xi || !gsum || !gx || !gvar) { fprintf(stderr, "OOM\n"); exit(1); }

    double prev = -INFINITY, final_ll = -INFINITY;

    for (int iter = 0; iter < max_iter; iter++) {

        /* ---- Pass 1: each thread runs forward-backward on its chunk ---- */
        double total_ll = 0.0;

        #pragma omp parallel for schedule(static,1) num_threads(n_threads) \
                                 reduction(+:total_ll)
        for (int p = 0; p < n_threads; p++) {
            int Tp  = chunk_len[p];
            int off = chunk_start[p];
            memset(t_xi[p],   0, (size_t)N * N * sizeof(double));
            memset(t_gsum[p], 0, N * sizeof(double));
            memset(t_gx[p],   0, N * sizeof(double));
            total_ll += chunk_forward_backward(
                Tp, N, obs + off, means, vars, trans,
                t_alpha[p], t_beta[p], t_emit[p], t_scale[p],
                t_xi[p], t_gsum[p], t_gx[p], NULL,
                /*do_gvar=*/0);
        }

        /* ---- Reduce pass-1 statistics ---- */
        memset(xi,   0, (size_t)N * N * sizeof(double));
        memset(gsum, 0, N * sizeof(double));
        memset(gx,   0, N * sizeof(double));
        for (int p = 0; p < n_threads; p++)
            for (int i = 0; i < N; i++) {
                gsum[i] += t_gsum[p][i];
                gx[i]   += t_gx[p][i];
                double *xd = xi + i*N, *xs = t_xi[p] + i*N;
                for (int j = 0; j < N; j++) xd[j] += xs[j];
            }

        /* ---- Update means ---- */
        for (int i = 0; i < N; i++)
            if (gsum[i] != 0.0) means[i] = gx[i] / gsum[i];

        /* ---- Pass 2: gvar with updated means ---- */
        #pragma omp parallel for schedule(static,1) num_threads(n_threads)
        for (int p = 0; p < n_threads; p++) {
            int Tp  = chunk_len[p];
            int off = chunk_start[p];
            memset(t_gv[p], 0, N * sizeof(double));
            chunk_forward_backward(
                Tp, N, obs + off, means, vars, trans,
                t_alpha[p], t_beta[p], t_emit[p], t_scale[p],
                NULL, NULL, NULL, t_gv[p],
                /*do_gvar=*/1);
        }

        /* ---- Reduce gvar ---- */
        memset(gvar, 0, N * sizeof(double));
        for (int p = 0; p < n_threads; p++)
            for (int i = 0; i < N; i++) gvar[i] += t_gv[p][i];

        /* ---- Update variances ---- */
        for (int i = 0; i < N; i++) {
            if (gsum[i] != 0.0) {
                double v = gvar[i] / gsum[i];
                vars[i] = (v < 1e-6) ? 1e-6 : v;
            }
        }

        /* ---- Update transition matrix ---- */
        for (int i = 0; i < N; i++) {
            double row = 0.0;
            double *xi_row = xi    + i*N;
            double *tr_row = trans + i*N;
            for (int j = 0; j < N; j++) row += xi_row[j];
            if (row == 0.0) continue;
            double inv_r = 1.0 / row;
            row = 0.0;
            for (int j = 0; j < N; j++) {
                double p = xi_row[j] * inv_r;
                if (iter >= 100) {
                    tr_row[j] = (p < 0.001) ? 0.0 : p;
                } else {
                    tr_row[j] = p;  // no pruning before iteration 100
                }
                row += tr_row[j];
            }
            if (row > 0.0) {
                inv_r = 1.0 / row;
                for (int j = 0; j < N; j++) tr_row[j] *= inv_r;
            }
        }

        final_ll = total_ll;
        double delta = total_ll - prev;
        fprintf(stderr, "Iter %d | loglik=%.4f | delta=%.4e\n",
                iter+1, total_ll, delta);
        if (iter > 1000 && fabs(delta) < thresh) break;
        prev = total_ll;
    }

    for (int p = 0; p < n_threads; p++) {
        free(t_alpha[p]); free(t_beta[p]); free(t_emit[p]); free(t_scale[p]);
        free(t_xi[p]);    free(t_gsum[p]); free(t_gx[p]);   free(t_gv[p]);
    }
    free(t_alpha); free(t_beta); free(t_emit); free(t_scale);
    free(t_xi);    free(t_gsum); free(t_gx);   free(t_gv);
    free(chunk_start); free(chunk_len);
    free(xi); free(gsum); free(gx); free(gvar);
    return final_ll;
}


/* ------------------------------------------------------------------ */
/* Main                                                               */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    /*
     * Keep OpenMP threads spinning rather than sleeping between iterations.
     * Must be set before the first parallel region.
     * OMP_WAIT_POLICY : standard OpenMP 3.0+
     * GOMP_SPINCOUNT  : GCC libgomp
     * KMP_BLOCKTIME   : Intel / LLVM OpenMP
     */
    setenv("OMP_WAIT_POLICY", "active",   1);
    setenv("GOMP_SPINCOUNT",  "INFINITE", 1);
    setenv("KMP_BLOCKTIME",   "infinite", 1);

    if (argc < 4 || argc > 5) {
        fprintf(stderr, "Usage: %s <n_states> <kmer_file> <k> <n_threads>\n", argv[0]);
        return 1;
    }
    int N         = atoi(argv[1]);
    int seq_k     = atoi(argv[3]);
    int n_threads = (argc == 5) ? atoi(argv[4]) : omp_get_max_threads();
    omp_set_num_threads(n_threads);
    fprintf(stderr, "Using %d OpenMP thread(s)\n", n_threads);

    KmerTable *table = load_kmer_table(argv[2]);
    int obs_k = table->k;
    if (seq_k < 1) { fprintf(stderr, "Error: de Bruijn k must be >= 1\n"); return 1; }
    if (obs_k < 1) { fprintf(stderr, "Error: invalid k-mer table\n"); return 1; }

    char *db  = debruijn(seq_k);
    int   len = (int)strlen(db);
    if (len < obs_k) {
        fprintf(stderr, "Error: de Bruijn sequence of order %d is too short "
                        "for %d-mer lookup\n", seq_k, obs_k);
        free(db); return 1;
    }
    int     n_obs = len - obs_k + 1;
    double *obs   = malloc(n_obs * sizeof(double));
    char    buf[64];
    for (int t = 0; t < n_obs; t++) {
        memcpy(buf, db + t, obs_k); buf[obs_k] = '\0';
        double v = kmer_table_lookup(table, buf);
        obs[t] = isnan(v) ? 0.0 : v;
        if (isnan(v))
            fprintf(stderr, "Warning: missing %d-mer '%s' at position %d\n",
                    obs_k, buf, t);
    }
    free(db);

    double *means = malloc(N * sizeof(double));
    double *vars  = malloc(N * sizeof(double));
    double *trans = malloc((size_t)N * N * sizeof(double));
    if (!means || !vars || !trans) { fprintf(stderr, "OOM\n"); exit(1); }

    gmm_init(table, N, means, vars);

    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < N; j++) {
            trans[i*N+j] = (i == j) ? 0.0 : 1.0;
            sum += trans[i*N+j];
        }
        for (int j = 0; j < N; j++) trans[i*N+j] /= sum;
    }

    double wall0 = omp_get_wtime();
    baum_welch_full(N, means, vars, trans, obs, n_obs, 0.0001, 10000, n_threads);
    double wall1 = omp_get_wtime();

    fprintf(stderr, "Wall time: %.3f seconds\n", wall1 - wall0);
    fprintf(stderr, "Number of states: %d\n", N);
    fprintf(stderr, "DeBruijn sequence generated with all possible %d-mers\n", seq_k);

    for (int i = 0; i < N; i++) { if (i) putchar('\t'); printf("%.6f", means[i]); }
    putchar('\n');
    for (int i = 0; i < N; i++) { if (i) putchar('\t'); printf("%.6f", vars[i]); }
    putchar('\n');
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            if (i || j) putchar('\t');
            printf("%.6f", trans[i*N+j]);
        }
    putchar('\n');

    free(obs); free(means); free(vars); free(trans);
    return 0;
}