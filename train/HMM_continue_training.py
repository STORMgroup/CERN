import numpy as np
import pandas as pd
import random
from hmmlearn.hmm import GaussianHMM
import time

### Train an already created HMM using the Baum-Welch algorithm on nanopore sequencing data

# This is very similar to HMM_train.py, but this begins training from an HMM already made

# A notable difference is that near-zero transitions are set to zero starting from the first training iteration


### Parameters ###

S = 256                             # Number of states
DATA_LEN = 150000                   # Amount of events in each training batch
TRAINING_ITERS = 20                 # Training batches
HMM_TO_LOAD = "latent_hmm_256_150000_40.txt"                                        # Starting HMM
HMM_FILENAME = f"latent_hmm_{S}_{DATA_LEN}_{TRAINING_ITERS}_continued.txt"          # HMM file to save to
NANOPORE_MODEL_FILE = "../extern/kmer_models/uncalled_r1041_model_only_means.txt"   # Nanopore file to get expected k-mer levels from

###################



def load_hmm():
    """
    Load an HMM saved as:
      line 1: means
      line 2: covariances
      line 3: transitions (row-major)
    """
    with open(HMM_TO_LOAD) as f:
        means = np.array([float(x) for x in f.readline().split()])
        covars = np.array([float(x) for x in f.readline().split()])
        trans_flat = np.array([float(x) for x in f.readline().split()])

    n_states = len(means)
    transmat = trans_flat.reshape((n_states, n_states))

    hmm = GaussianHMM(
        n_components=n_states,
        covariance_type="diag",
        init_params="",
        params="mtc",
        random_state=0
    )

    hmm.startprob_ = np.ones(n_states) / n_states
    hmm.means_ = means.reshape(-1, 1)
    hmm.covars_ = covars.reshape(-1, 1)
    hmm.transmat_ = transmat

    return hmm

def gen_data(length, model_file="uncalled_r1041_model_only_means.txt"):
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
        model_file,
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
        events[i] = kmer_dict[kmer] + random.random() / 10 - 0.05

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

    hmm = load_hmm()
    prev_lik = -9999999
    print("Continuing latent HMM training with these params:")
    print(f"States: {S}")
    print(f"Data batch size size: {DATA_LEN}")
    print(f"Training iterations: {TRAINING_ITERS}")
    print(f"Starting HMM file: {HMM_TO_LOAD}")
    print(f"Out file: {HMM_FILENAME}")
    start = time.perf_counter()

    for i in range(TRAINING_ITERS):
    
        sequence, events = gen_data(DATA_LEN)

        hmm.fit(events)
        trim_transitions(hmm)

        lik = eval_hmm(hmm, 200000)
        print(f"Done training iteration {i}")
        print(f"Current likelihood on 100000 events: {lik}")
        if lik < prev_lik:
            break
        save_hmm(hmm)
        prev_lik = lik
    print("Done training.")
    end_time = time.perf_counter()
    totaltime = end_time - start
    print(f"Took {totaltime:.4f} seconds.")
    
