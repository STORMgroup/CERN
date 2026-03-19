# Error Correction Read Mapping Evaluation

In order to evaluate CERN, make sure you have followed the steps in [evaluation](../), and have ran the script [create_indexes.sh](../create_indexes.sh).

In this directory, run each script in this order:

```bash
conda activate cern-env
bash map_all.sh [THREAD_COUNT]
bash annotate_pafs.sh
bash summarize_annotations.sh
```

Substitute [THREAD_COUNT] with the number of wanted threads, and make sure that the cern-env conda environment is setup correctly and activated.

Ensure `map_all.sh` is ran with sufficient memory and CPU resources to avoid errors and long runtimes.