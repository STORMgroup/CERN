# CRANE Evaluation Data

Get started by downloading the CRANE evaluation datasets:

```bash
bash download_cern_data.sh
```


Then prepare the datasets using the preparation scripts below.

 - Before running these, we recommend inspecting them to ensure the needed tools have been installed and/or compiled.
 - Ensure that you have compiled CRANE and segmentation code in [../../src](../../src)
 - Running these will require access to GPUs due to the the dorado basecalling, and campolina segmentation steps.

```bash
bash prepare_d1.sh
bash prepare_d2.sh
bash prepare_d3.sh
```

If you wish to train your own HMM models using the provided configuration, you will need to prepare the training dataset as well:

```bash
bash prepare_training.sh
```