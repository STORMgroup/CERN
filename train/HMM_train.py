import numpy as np
import pandas as pd
import random
from hmmlearn.hmm import GaussianHMM
import time

### Train an HMM using the Baum-Welch algorithm on nanopore sequencing data

# This runs multiple rounds of Baum-Welch on different batches of data, 
# which helps to get the training process un-stuck when in local minima 

# Transitions near 0 probability are set to 0 in the second half of training batches by default
# If training finishes before then, we start trimming earlier


### Parameters ###

S = 196                             # Number of states
DATA_LEN = 100000                   # Amount of events in each training batch
TRAINING_ITERS = 30                 # Training batches to run, will stop here or when model converges
HMM_FILENAME = f"hmm_{S}_{DATA_LEN}_{TRAINING_ITERS}_phase1.txt"       # HMM file to save to
NANOPORE_MODEL_FILE = "../extern/kmer_models/uncalled_r1041_model_only_means.txt"     # Nanopore file to get expected k-mer levels from

###################

def gen_data(length):
    """
    Generate a sequence of events of given length.

    Steps:
    1. Generate a random DNA sequence of length (length + 8)
    2. Load kmer -> level_mean mapping from model_file
    3. Convert DNA sequence into events using 9-mers

    Returns:
        dna_seq (str): DNA sequence of length (length + 8)
        events (np.ndarray): array of floats of length `length`
    """

    # -----------------------
    # Generate DNA sequence
    # -----------------------
    bases = ["A", "C", "G", "T"]
    dna_len = length + 8
    dna_seq = "".join(random.choices(bases, k=dna_len))

    # -----------------------
    # Load k-mer model
    # -----------------------
    kmer_df = pd.read_csv(
        NANOPORE_MODEL_FILE,
        sep="\t",
        header=None,
        names=["kmer", "level_mean"],
        dtype={"kmer": str, "level_mean": float}
    )
    kmer_dict = dict(zip(kmer_df["kmer"], kmer_df["level_mean"]))

    # -----------------------
    # Convert DNA to events
    # -----------------------
    events = np.zeros(length, dtype=float)

    for i in range(length):
        kmer = dna_seq[i:i+9]
        if kmer not in kmer_dict:
            raise KeyError(f"K-mer {kmer} not found in model file")
        events[i] = kmer_dict[kmer]

    events = events.reshape(-1, 1)

    return dna_seq, events

def build_latent_hmm():
    """
    Build a Gaussian HMM with S
    
    States are ordered as:
    AA (K states), AC (K states), AG (K states), AT (K states),
    CA (K states), ..., TT (K states)
    """

    # -----------------------
    # Parameters
    # -----------------------
    n_states = S

    # -----------------------
    # Transition matrix
    # -----------------------
    transmat = np.zeros((n_states, n_states))

    for i in range(n_states):
        for j in range(n_states):
            if i != j:
                transmat[i][j] = 1/(n_states-1)

    # -----------------------
    # Emission parameters
    # -----------------------
    # One-dimensional Gaussian emissions
    means = np.zeros((n_states, 1))
    covars = np.zeros((n_states, 1))

    # Base emission centers: AA = 0, AC = 1000, AG = 2000, ...
    for i in range(n_states):
        means[i][0] = (random.random()) - 0.5
        covars[i][0] = 2.0  # std = 2

    # -----------------------
    # Initial state distribution
    # -----------------------
    startprob = np.ones(n_states) / n_states

    # -----------------------
    # Build HMM
    # -----------------------
    model = GaussianHMM(
        n_components=n_states,
        covariance_type="diag",
        init_params="",      # we set everything manually
        params="tmc",
        random_state=0
    )

    model.startprob_ = startprob
    model.transmat_ = transmat
    model.means_ = means
    model.covars_ = covars

    return model

def eval_hmm(hmm, length):
    sequence, events = gen_data(length)

    log_likelihood = hmm.score(events)

    return log_likelihood

def save_hmm(hmm):

    means = hmm.means_

    covars = hmm.covars_.reshape(hmm.n_components)

    trans = hmm.transmat_

    with open(HMM_FILENAME, "w") as file:
        # Line 1: means
        file.write(" ".join(str(m) for m in means[:, 0]) + "\n")

        # Line 2: covariances
        file.write(" ".join(str(c) for c in covars) + "\n")

        # Line 3: transitions (row-major)
        file.write(" ".join(str(t) for t in trans.flatten()) + "\n")

        # Line 4: error-free transitions (1 where transition > 0, 0 otherwise)
        error_free_trans = ["1" if t > 0 else "0" for t in trans.flatten()]
        file.write(" ".join(error_free_trans) + "\n")

    return

def trim_transitions(hmm, threshold=0.001):
    """
    Zero out transitions with probability < threshold,
    then renormalize each row to sum to 1.
    """

    trans = hmm.transmat_.copy()
    n_states = trans.shape[0]

    for i in range(n_states):
        row = trans[i]

        # Zero out small probabilities
        row[row < threshold] = 0.0

        row_sum = row.sum()

        if row_sum == 0.0:
            raise ValueError("Transitions zeroed out.")
        else:
            # Renormalize
            row /= row_sum

    hmm.transmat_ = trans


if __name__ == "__main__":

    hmm = build_latent_hmm()
    prev_lik = -9999999
    num_skips = 0
    print("Beginning latent HMM training with these params:")
    print(f"States: {S}")
    print(f"Data batch size size: {DATA_LEN}")
    print(f"Training iterations: {TRAINING_ITERS}")
    print(f"Out file: {HMM_FILENAME}")
    start = time.perf_counter()

    for i in range(TRAINING_ITERS):
    
        sequence, events = gen_data(DATA_LEN)

        hmm.fit(events)
        if i > TRAINING_ITERS / 2: # Don't want to trim transitions too early, will break the HMM
            trim_transitions(hmm)

        lik = eval_hmm(hmm, 200000)
        print(f"Done training iteration {i}")
        print(f"Current log likelihood on 200000 events: {lik}")
        if lik < prev_lik:
            if (i <= TRAINING_ITERS / 2):
                trim_transitions(hmm)
                continue
            else:
                break
        save_hmm(hmm) # Save as we go
        prev_lik = lik
    print("Done training.")
    end_time = time.perf_counter()
    totaltime = end_time - start
    print(f"Took {totaltime:.4f} seconds.")
    
