import os
import pandas as pd
import numpy as np
import glob
from transitleastsquares import transitleastsquares
import sys

def generate_filename(tic, period, rp_earth, trial):
    """Generate filename for individual TIC result CSV"""
    tic_id = str(tic).replace(' ', '_')
    # Round to match the precision used in the data row creation
    period_rounded = round(float(period), 4)
    rp_earth_rounded = round(float(rp_earth), 2)
    trial_int = int(trial)
    return f"tls_results_per_tic/{tic_id}_P{period_rounded}_Rp{rp_earth_rounded}_trial{trial_int}.csv"

def check_if_already_processed(tic, period, rp_earth, trial):
    """Check if this specific combination has already been processed"""
    filename = generate_filename(tic, period, rp_earth, trial)
    return os.path.exists(filename)

def parse_slice(slice_str):
    if slice_str.startswith('[') and slice_str.endswith(']'):
        slice_str = slice_str[1:-1]
    parts = slice_str.split(':')
    start = int(parts[0]) if parts[0] else None
    stop = int(parts[1]) if len(parts) > 1 and parts[1] else None
    return slice(start, stop)

def run_tls(tic, time, flux, true_period):

    tls = transitleastsquares(time, flux)
    result = tls.power(period_max = 25)

    powers = result.power
    periods = result.periods

    periodogram = pd.DataFrame({
        'period': periods,
        'power': powers
    })

    power_7_periods = periodogram[periodogram['power'] > 7]['period'].values
    print(f"TLS Periods with power > 7: {power_7_periods}")
    power_7_power = periodogram[periodogram['power'] > 7]['power'].values


    print(f"{tic}, FAP: {result.FAP}, SDE: {result.SDE}, True Period: {true_period}")
    # Mark as found if period matches true period or a simple alias (1/2x, 2x, 1/3x, 3x, etc.)
    aliases = [true_period, true_period / 2, true_period * 2, true_period / 3, true_period * 3]
    #found = result.FAP < 1e-4 and any(0.99 * alias < result.period < 1.01 * alias for alias in aliases)
    found = any(0.99 * alias < hp < 1.01 * alias for hp in power_7_periods for alias in aliases)


    alias_found = any(0.99 * alias < result.period < 1.01 * alias for alias in aliases)
    print(f"Detection found: {found}")
    if alias_found == True and found == False:
        print(f"True transit or alias found, but high FAP: {alias_found}")
    #print("TLS Transit duration (hours):", result.duration * 24)
    

    return found, alias_found, result.period, result.SDE, power_7_periods, power_7_power, result.period_uncertainty, result.depth, result.rp_rs, result.FAP, result.snr


# Create output directory for individual TIC results
output_dir = 'tls_results_per_tic'
os.makedirs(output_dir, exist_ok=True)

if len(sys.argv) > 1:
    tic_slice_str = sys.argv[1]
    #save_file = sys.argv[2]
else:
    print("Usage: python3 script.py '[start:stop]' output.csv")
    sys.exit(1)

tic_slice = parse_slice(tic_slice_str)

lc_directories = sorted(glob.glob('injected_flattened_lcs/*/'))
all_TICs = [dir.rstrip('/').split('/')[-1] for dir in lc_directories]
TICs = all_TICs[tic_slice]
print(len(TICs))



for tic in TICs:
    files = glob.glob(f'injected_flattened_lcs/{tic}/*.csv')
    #files = files[2:]
    for file_name in files:
        
        
        tic = file_name.split('/')[1].replace('_', ' ')
        parameters = file_name.split('/')[2]
        period = float(parameters.split('_')[0].replace('P', ''))
        rp_earth = float(parameters.split('_')[1].replace('Rp', ''))
        trial = parameters.split('_')[2].replace('trial', '').replace('.csv', '')

        print("tic:", tic, type(tic))
        print("period:", period, type(period))
        print("rp_earth:", rp_earth, type(rp_earth))
        print("trial:", trial, type(trial))
        
        # Check if this combination has already been processed
        if check_if_already_processed(tic, period, rp_earth, trial):
            print(f"Skipping {tic} P={period} Rp={rp_earth} trial={trial} (already processed)")
            continue
	

        print(f'{tic}, period: {period}, rp_earth: {rp_earth}, trial: {trial}')
        true_period = float(period)

        data = pd.read_csv(file_name)
        time = data['time'].values
        flux = data['flux'].values
        true_t0 = data['true_t0'].values[0]
        true_inc = data['true_inc'].values[0]
        print(type(time), type(flux))

        found, alias_found, tls_period, SDE_strongest,  power_array_7, SDE_array_7, tls_period_uncertainty, tls_depth, tls_rp_rs, tls_FAP, tls_snr = run_tls(tic,time, flux, true_period)

        star_params = pd.read_csv("stellar_params_CTL.csv").query(f"id == {tic.split()[1]}")
        if star_params.empty:
            print(f"⚠️ No stellar params for {tic}")
            continue

        st_radius = star_params['RAD'].values[0]
        st_mass = star_params['MASS'].values[0] 
        st_teff = star_params['teff'].values[0]
        st_tmag = star_params['tmag'].values[0]

        row = {
        'TIC': tic,
        'True Radius (Earth Radii)': round(float(rp_earth), 2),
        'True Period (Days)': round(float(period), 4),
        'True t0': round(float(true_t0), 4),
        'True Inclination': round(float(true_inc), 4),        
        'Stellar Radius': round(float(st_radius), 2),
        'Stellar Temperature': st_teff,
        'Stellar Magnitude': st_tmag,
        'Detection': found,
        'TLS Period': tls_period,
        'TLS SDE strongest': SDE_strongest,
        'TLS SDE > 7 array': SDE_array_7,
        'TLS Periods array': power_array_7,
        'TLS Period Uncertainty': tls_period_uncertainty,
        'TLS Depth': tls_depth,
        'TLS Rp/Rs': tls_rp_rs,
        'TLS FAP': tls_FAP,
        'TLS SNR': tls_snr,
        #'Detection %': detections / trials,
        '# trials': trial
            }

        row_df = pd.DataFrame([row])

        # Save to individual CSV file
        output_filename = generate_filename(tic, period, rp_earth, trial)
        row_df.to_csv(output_filename, index=False)
        print(f"Saved results to {output_filename}")
