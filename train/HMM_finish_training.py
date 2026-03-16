import numpy as np
import pandas as pd
import random
from hmmlearn.hmm import GaussianHMM

### Train an HMM once on all possible k-mers

# The purpose of this is to do one last training batch on a completely balanced dataset
# This should reduce any biases incurred by training on randomly generated data

### Parameters ###

KMER_LEN = 10 # The length of the kmers that we want to see all of once in the training data
HMM_TO_LOAD = "hmm_256_150000_20_continued.txt" # HMM that we are loading
HMM_FILENAME = f"hmm_256_{KMER_LEN}_finished.txt" # File to save the resulting HMM to

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


def de_bruijn_sequence(alphabet, k):
    """
    Construct a De Bruijn sequence for alphabet and subsequences of length k.
    Returns a linear string (not cyclic) of length |alphabet|^k + k - 1
    """
    n = len(alphabet)
    a = [0] * (n * k)
    sequence = []

    def db(t, p):
        if t > k:
            if k % p == 0:
                for i in range(1, p + 1):
                    sequence.append(alphabet[a[i]])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, n):
                a[t] = j
                db(t + 1, t)

    db(1, 1)
    return "".join(sequence) + "".join(sequence[:k-1])  # linearized


def gen_data(model_file="uncalled_r1041_model_only_means.txt"):
    """
    Generate events from a De Bruijn DNA sequence that contains
    every possible KMER_SIZE-mer exactly once.

    Returns:
        dna_seq (str): DNA De Bruijn sequence
        events (np.ndarray): array of floats, one per k-mer
    """

    bases = ["A", "C", "G", "T"]
    k = 9

    # -----------------------
    # Generate De Bruijn DNA
    # -----------------------
    dna_seq = de_bruijn_sequence(bases, KMER_LEN)

    length = len(dna_seq) - k + 1   # number of k-mers (should be 4^k)

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
        kmer = dna_seq[i:i+k]
        if kmer not in kmer_dict:
            raise KeyError(f"K-mer {kmer} not found in model file")
        events[i] = kmer_dict[kmer] + random.random() / 10 - 0.05

    events = events.reshape(-1, 1)

    return dna_seq, events

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

if __name__ == "__main__":

    hmm = load_hmm()
    
    sequence, events = gen_data()

    hmm.fit(events)

    save_hmm(hmm)

    print(f"Done training on all possible {KMER_LEN}-mers.")
    
