# TESS Injection + Recovery Pipeline

Welcome to the injection+recovery pipeline!

## Overview

The way that it's organized is into 3 parts:
1. Getting the TESS sample to run the injection+recovery
2. Injecting the Lightcurves
Sending the Lightcurves to Geryon, running them, and retrieving the results
3. Turning the results into heatmaps that also take into account false positives

The folders are numbered in chronological order of what steps to take, and inside each folder, the files/folders you need to use are also numbered in order of what to do. (e.g. start with `1_get_sample/`, and inside that folder start with `1_TOI_CTL.ipynb`)

---

## Part One: `1_get_sample/`

1. **Open** `1_TOI_CTL.ipynb`

2. **Add to the same folder** (`1_get_sample/`) a csv that contains the TIC IDs that you'd like to get a representative sample for with the relevant parameters
   - See `TOI_Mar2025_1pt5to4_R.csv` as an example

3. **Make sure** that `CTL_2025April29.npz` is in the same folder as well, or from whatever date, and modify cell 8 to reflect the correct date

4. **Run the cell**, which will generate the sample and automatically add it to the relevant folders

**Note:**
- In `histogram_sample_2d.py`, in the `make_control` function, you can modify `n_targets` to get a different number of targets from the CTL. It's currently set to get 10x the number of inputted targets.
- The default is set up to use stellar radius and tess magnitude as the parameters, but this can be changed

---

## Part Two: `2_injection/`

1. **Open and run** `1_output_injected_flattened_lcs_randomized_periods_radii.ipynb`. This will put the injected lightcurves into the relevant folder.
   - If you want to modify the period/radius space, you can do that towards the bottom of the cell. Otherwise, the notebook should run automatically.
   - **This can take a few hours** for ~1000 lightcurves
   - **Be careful to only run this cell once**, because previous progress checking is not implemented. It's non-trivial since the injections are randomized in the cell. Feel free to take a stab at adding previous progress checking though!! If you have to re-run part way through, you'll have to go into `3_import_to_geryon/` and delete the existing TIC folders, or you can manually define which TICs to run again.
   - Also be aware that **you need to have an internet connection** while it's running or else the lightcurves can't be accessed.
   - **Use the last cell** to ensure that you generated the number of lightcurves that you intended before uploading to Geryon

2. **Open** `2_generate_pbs_files.ipynb` and edit to fit the way you want to run the jobs in Geryon, and run. This will generate the pbs files in the relevant folder.
   - You can modify `n_jobs` (how many jobs to split the pipeline into) and `chunk_size` (which is how many lightcurves are run by each job)

3. **Open the folder** `3_import_to_geryon/` and ensure that you have the following files:
   - `in_geryon_run_tls.py`
   - `stellar_params_CTL.csv`
   - `run_tls_part*.pbs` (how ever many you generated)
   - AND the folder `injected_flattened_lcs/` with all the TIC ID folders and lightcurve csv files inside

---

## Part in Geryon:

1. **Now we need to send** the entire folder `3_import_to_geryon/` to Geryon (surprise)
   - Example command (in local terminal):
     ```bash
     rsync -avz --progress /Users/danayaptangco/local_code/Mulders/inj_rec_yaptangco/2_injection/3_import_to_geryon \
       dyaptangco@geryon2.astro.puc.cl:/home/dyaptangco/
     ```

2. **Now, log into Geryon** and navigate into the `3_import_to_geryon` folder

3. **Run job** by running the relevant .pbs files
   ```bash
   qsub run_tls_part*.pbs
   ```

4. **Once the jobs are done running**, send the results back to local computer
   - Example command (in local terminal):
     ```bash
     scp -r dyaptangco@geryon2.astro.puc.cl:'/home/dyaptangco/3_import_to_geryon/tls_results_per_tic' \
       /Users/danayaptangco/local_code/Mulders/inj_rec_yaptangco/3_recovery_heatmaps/
     ```

---

## Part Three: `3_recovery_heatmaps/`

1. **While Geryon is running**, you can run the code in `1_rec_pipeline_CTL_sample.ipynb` which will catch false positive cells. **This can take up to a full day** to run ~1000 stars. This code does have previous progress checking, so feel free to start and stop again.
   - However, be aware that **you need to have an internet connection** while it's running or else the lightcurves can't be accessed.

2. **Once Geryon has finished running** and you've copied the results back to your local machine, you can open `2_heatmaps_analysis.ipynb` to see the results of the injections in heatmap form.
   - The first cell creates one map per TIC
   - The second cell creates one map per specific temperature bin
   - The third cell adds contours
   - The fourth cell is contour analysis

---

## Important Notes

**That's everything!** If you have any questions you can reach me at dsy24@ic.ac.uk.

**I would recommend starting with a fresh folder** every time running this pipeline so that the file/folder generation doesn't get mixed up. However, **be very careful** that you are sending and receiving data to the correct folder every time from Geryon!! I could foresee that being an easy mistake to make.
