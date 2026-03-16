import sys
import random
import numpy as np
from hmmlearn.hmm import GaussianHMM
import time

### Fine-tune an HMM on real segmented data

# This will train a starting HMM on real nanopore data, the data must already be segmented

# The format of nanopore data must be 1 read per line, starting with read ID, events seperated by whitespace:
# read1 event1 event2...
# read2 event1 event2...

### Parameters ###

STARTING_HMM = "latent_hmm_256_10_finished.txt"             # HMM to load in
EVENT_FILE = "../d9_3_campolina_events_nonan.sorted.txt"    # File to load real segmented nanopore data from
HMM_FILENAME = "hmm_L256_ft_t_db_camp_2000000.txt"          # File to save HMM to
EST_SKIP = 0.05                 # An estimation of how often a skip error occurs, helpful for training
EST_STAY = 0.3                  # An estimation of how often a stay error occurs, helpful for training
DATA_LEN = 2000000              # Amount of events to load
PARAMETERS = "t"                # Add t for transitions, c for covariances, m for means.
DEBIAS = True                   # Whether or not to force the original non-zero transitions to be proportionally restored after fine-tuning
ERRORS_IN_HMM = False           # True if there are already errors in the HMM we're loading

##################

def load_hmm(hmm_file):
    """
    Load an HMM saved as:
      line 1: means
      line 2: covariances
      line 3: transitions (row-major)
    """
    with open(hmm_file) as f:
        means = np.array([float(x) for x in f.readline().split()])
        covars = np.array([float(x) for x in f.readline().split()])
        trans_flat = np.array([float(x) for x in f.readline().split()])

    n_states = len(means)
    transmat = trans_flat.reshape((n_states, n_states))

    hmm = GaussianHMM(
        n_components=n_states,
        covariance_type="diag",
        init_params="",
        params=PARAMETERS,
        random_state=0
    )

    hmm.startprob_ = np.ones(n_states) / n_states
    hmm.means_ = means.reshape(-1, 1)
    hmm.covars_ = covars.reshape(-1, 1)
    hmm.transmat_ = transmat

    return hmm


def load_random_sequences(event_file, total_events):
    """
    Randomly sample event sequences from file.
    Each line is treated as an independent sequence.

    Returns:
        X: concatenated events, shape (N, 1)
        lengths: list of sequence lengths
    """
    with open(event_file) as f:
        lines = f.readlines()

    random.shuffle(lines)

    all_events = []
    lengths = []
    idx = 0

    while sum(lengths) < total_events:
        if idx >= len(lines):
            random.shuffle(lines)
            idx = 0

        parts = lines[idx].strip().split()
        idx += 1

        # discard ID
        seq = [float(x) for x in parts[1:]]
        if not seq:
            continue

        remaining = total_events - sum(lengths)
        seq = seq[:remaining]

        all_events.extend(seq)
        lengths.append(len(seq))

    X = np.array(all_events, dtype=float).reshape(-1, 1)
    return X, lengths

def add_error_to_hmm(hmm):
    """
    Return a new HMM identical to `hmm` except with modified transitions:
      - self-transition with probability EST_STAY
      - skip transitions (two-step reachable states) with total probability EST_SKIP,
        distributed proportionally to the two-step path probabilities
    """
    n = hmm.n_components
    old_T = hmm.transmat_

    remaining = 1.0 - EST_STAY - EST_SKIP

    # Uniform skip matrix: zero on diagonal, uniform elsewhere
    uniform_skip = np.ones((n, n))
    np.fill_diagonal(uniform_skip, 0.0)
    uniform_skip /= (n - 1)

    new_T = (
        EST_STAY * np.eye(n)
        + remaining * old_T
        + EST_SKIP * uniform_skip
    )

    # --- construct new HMM ---
    new_hmm = GaussianHMM(
        n_components=n,
        covariance_type="diag",
        init_params="",
        params=PARAMETERS,
        random_state=0
    )
    new_hmm.n_features = 1
    hmm.n_features = 1
    new_hmm.startprob_ = hmm.startprob_.copy()
    new_hmm.means_ = hmm.means_.copy()
    covars = hmm.covars_.copy()
    covars = covars.reshape(-1, 1)
    new_hmm.covars_ = covars
    new_hmm.transmat_ = new_T

    return new_hmm

def save_hmm(hmm, old_T):

    means = hmm.means_

    covars = hmm.covars_.reshape(hmm.n_components)

    trans = hmm.transmat_

    error_free_trans = []
    for t in old_T.flatten():
        if t > 0:
            error_free_trans.append("1")
        else:
            error_free_trans.append("0")

    with open(HMM_FILENAME, "w") as file:
        # Line 1: means
        file.write(" ".join(str(m) for m in means[:, 0]) + "\n")

        # Line 2: covariances
        file.write(" ".join(str(c) for c in covars) + "\n")

        # Line 3: transitions (row-major)
        file.write(" ".join(str(t) for t in trans.flatten()) + "\n")

        # Line 4: Error free transitions
        file.write(" ".join(error_free_trans) + "\n")
    
    return

def debias(err_hmm, old_T):
    """
    Reweights only originally-allowed transitions to match old proportions.
    All other transitions are left untouched.

    For each state i:
      X = {j | old_T[i,j] != 0}
      mass = sum_j∈X err_T[i,j]
      err_T[i,j] = mass * old_T[i,j] / sum_k∈X old_T[i,k]   for j∈X
      all other transitions unchanged

    Parameters
    ----------
    err_hmm : HMM object
        Fine-tuned HMM
    old_T : np.ndarray (N, N)
        Original transition matrix
    eps : float
        Numerical stability constant

    Returns
    -------
    err_hmm : HMM object
    """

    T_new = err_hmm.transmat_.copy()
    N = T_new.shape[0]

    for i in range(N):
        # Indices of transitions that originally existed
        X = np.nonzero(old_T[i])[0]

        if X.size == 0:
            continue

        # Mass currently assigned to these transitions
        mass = T_new[i, X].sum()


        old_probs = old_T[i, X]
        old_sum = old_probs.sum()


        # Old proportions
        old_props = old_probs / old_sum

        # Redistribute ONLY within X
        T_new[i, X] = mass * old_props

        # everything else is untouched

    err_hmm.transmat_ = T_new
    return err_hmm


if __name__ == "__main__":

    hmm_file = STARTING_HMM
    event_file = EVENT_FILE

    print("="*60)
    print("HMM Fine-tuning (Training Start)")
    print("="*60)
    print(f"HMM input file      : {hmm_file}")
    print(f"Event input file    : {event_file}")
    print(f"Output HMM filename : {HMM_FILENAME}")
    print("")
    print("Training configuration:")
    print(f"  DATA_LEN              = {DATA_LEN}")
    print(f"  PARAMETERS            = {PARAMETERS}")
    print(f"  EST_SKIP              = {EST_SKIP}")
    print(f"  EST_STAY              = {EST_STAY}")
    print(f"  ERRORS_IN_HMM         = {ERRORS_IN_HMM}")
    print(f"  DEBIAS= {DEBIAS}")
    print("-"*60)

    t0 = time.time()

    print("Loading HMM...")
    hmm = load_hmm(hmm_file)
    old_T = hmm.transmat_.copy()
    print(f"HMM loaded with {hmm.transmat_.shape[0]} states")

    if not ERRORS_IN_HMM:
        print("Injecting error model into HMM...")
        err_hmm = add_error_to_hmm(hmm)
    else:
        print("Using HMM as-is (no error injection)")
        err_hmm = hmm

    print(f"Loading training data (DATA_LEN={DATA_LEN})...")
    X, lengths = load_random_sequences(event_file, DATA_LEN)
    print(f"Loaded {len(lengths)} sequences")
    print(f"Total observations: {X.shape[0]}")

    print("Starting training (Baum-Welch / EM)...")
    t_train_start = time.time()
    err_hmm.fit(X, lengths)
    t_train_end = time.time()
    print("Training complete")

    print(f"Training time: {t_train_end - t_train_start:.2f} seconds")

    if DEBIAS:
        print("Restoring original non-zero transition structure (debiasing)...")
        err_hmm = debias(err_hmm, old_T)
        print("Debiasing complete")

    print("Saving fine-tuned HMM...")
    save_hmm(err_hmm, old_T)

    t1 = time.time()

    print("="*60)
    print("Training Summary")
    print("="*60)
    print(f"Total runtime        : {t1 - t0:.2f} seconds")
    print(f"Training runtime     : {t_train_end - t_train_start:.2f} seconds")
    print(f"States               : {hmm.transmat_.shape[0]}")
    print(f"Sequences trained on : {len(lengths)}")
    print(f"Observations         : {X.shape[0]}")
    print(f"Parameters trained   : {PARAMETERS}")
    print(f"Output file          : {HMM_FILENAME}")
    print("="*60)