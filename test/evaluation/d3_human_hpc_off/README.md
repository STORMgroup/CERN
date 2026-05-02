# Error Correction Read Mapping Evaluation

In order to evaluate CRANE, make sure you have followed the steps in [evaluation](../).

In this directory, just run the following command:

```bash
conda activate crane-env
bash map_all.sh [THREAD_COUNT]
```

Substitute [THREAD_COUNT] with the number of wanted threads, and make sure that the crane-env conda environment is setup correctly and activated.

Ensure `map_all.sh` is ran with sufficient memory and CPU resources to avoid errors and long runtimes.