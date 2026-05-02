# Error Correction Read Mapping Evaluation

In order to evaluate CRANE, follow these steps:

1. Make the segmenter and CERN executable in [../../src/](../../src/)

2. Install the necessary tools by following the steps in [evaluation](../)
   
3. Download then segment and correct the datasets in [../data/cern_datasets/](../data/cern_datasets/)

4. Run RawHash2 with and without the corrected events, using the script map_all.sh in a dataset's evaluation directory
```bash
cd d1_ecoli_hpc_off # Navigate to the wanted directory
bash ./map_all.sh [THREAD_COUNT]
```

1. View the results:
```bash
cat summary.txt
```