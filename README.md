# TESS Injection + Recovery Pipeline

Welcome to the injection+recovery pipeline!

## Overview

The pipeline is organized into 3 main parts:
1. **Getting the TESS sample** to run the injection+recovery
2. **Injecting the Lightcurves** and sending them to Geryon for processing
3. **Turning the results into heatmaps** that account for false positives

The folders are numbered in chronological order of the steps to take, and inside each folder, the files/folders are also numbered in order of execution (e.g., start with `1_get_sample/`, and inside that folder start with `1_TOI_CTL.ipynb`).

---

## Part One: Sample Selection (`1_get_sample/`)

### Steps:

1. **Open** `1_TOI_CTL.ipynb`

2. **Add your input catalog**: Place a CSV file containing the TIC IDs you'd like to get a representative sample for with the relevant parameters in the `1_get_sample/` folder
   - See `TOI_Mar2025_1pt5to4_R.csv` as an example

3. **Ensure the CTL catalog is present**: Make sure `CTL_2025April29.npz` (or whatever date) is in the same folder, and modify cell 8 to reflect the correct date

4. **Run the notebook**: This will generate the sample and automatically add it to the relevant folders

### Customization:

- In `histogram_sample_2d.py`, in the `make_control` function, you can modify `n_targets` to get a different number of targets from the CTL
  - Currently set to get **10x** the number of inputted targets
- The default is set up to use **stellar radius** and **TESS magnitude** as the parameters, but this can be changed

---

## Part Two: Injection (`2_injection/`)

### Step 1: Generate Injected Lightcurves

1. **Open and run** `1_output_injected_flattened_lcs_randomized_periods_radii.ipynb`
   - This will put the injected lightcurves into the relevant folder
   - If you want to modify the period/radius space, you can do that towards the bottom of the cell
   - The notebook should run automatically
   - **This can take a few hours** for ~1000 lightcurves

2. **Progress checking**: 
   - Progress checking is now implemented! The notebook will track completed injections in `progress.csv` files
   - Set `resume=True` to skip already-completed work if you need to re-run
   - **Internet connection required** - lightcurves can't be accessed without it

3. **Verify output**: Use the last cell to ensure you generated the number of lightcurves that you intended before uploading to Geryon

### Step 2: Generate PBS Job Files

1. **Open** `2_generate_pbs_files.ipynb`
2. **Edit** to fit the way you want to run the jobs in Geryon
3. **Run** - This will generate the PBS files in the relevant folder
   - You can modify:
     - `n_jobs`: How many jobs to split the pipeline into
     - `chunk_size`: How many lightcurves are run by each job

### Step 3: Prepare for Geryon

Open the folder `3_import_to_geryon/` and ensure you have the following files:
- `in_geryon_run_tls.py`
- `stellar_params_CTL.csv`
- `run_tls_part*.pbs` (however many you generated)
- The folder `injected_flattened_lcs/` with all the TIC ID folders and lightcurve CSV files inside

---

## Part Geryon: Remote Processing

### Upload to Geryon

Send the entire folder `3_import_to_geryon/` to Geryon:

```bash
rsync -avz --progress /Users/danayaptangco/local_code/Mulders/inj_rec_yaptangco/2_injection/3_import_to_geryon \
  dyaptangco@geryon2.astro.puc.cl:/home/dyaptangco/
```

### Run Jobs on Geryon

1. Log into Geryon
2. Navigate into the `3_import_to_geryon` folder
3. Submit jobs:
   ```bash
   qsub run_tls_part*.pbs
   ```

### Retrieve Results

Once the jobs are done running, copy the results back to your local computer:

```bash
scp -r dyaptangco@geryon2.astro.puc.cl:'/home/dyaptangco/3_import_to_geryon/tls_results_per_tic' \
  /Users/danayaptangco/local_code/Mulders/inj_rec_yaptangco/3_recovery_heatmaps/
```

---

## Part Three: Recovery Analysis (`3_recovery_heatmaps/`)

### Step 1: False Positive Detection

**While Geryon is running**, you can run the code in `1_rec_pipeline_CTL_sample.ipynb`:
- This will catch false positive cells
- **Can take up to a full day** to run ~1000 stars
- **Has progress checking** - feel free to start and stop again
- **Internet connection required** - lightcurves can't be accessed without it

### Step 2: Generate Heatmaps

Once Geryon has finished running and you've copied the results back to your local machine, open `2_heatmaps_analysis.ipynb`:

- **Cell 1**: Creates one map per TIC
- **Cell 2**: Creates one map per specific temperature bin
- **Cell 3**: Adds contours
- **Cell 4**: Contour analysis

---

## Important Notes

**Recommendations:**
- I would recommend starting with a **fresh folder** every time you run this pipeline so that the file/folder generation doesn't get mixed up
- **Be very careful** that you are sending and receiving data to the **correct folder** every time from Geryon! This is an easy mistake to make

---

## Contact

If you have any questions, you can reach me at **dsy24@ic.ac.uk**.

---

## Requirements

- Python 3.x with the following packages:
  - `numpy`
  - `pandas`
  - `matplotlib`
  - `lightkurve`
  - `batman-package`
  - `transitleastsquares`
  - `seaborn`
  - `scipy`
  
- Internet connection for downloading TESS lightcurves
