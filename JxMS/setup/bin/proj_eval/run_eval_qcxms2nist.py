#!/usr/bin/env python3
# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                        — by Jesi Lee —                         ║
# ║ Project     : Blood Exposome CID-MS/MS Prediction via QCxMS    ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : run_eval_qcxms2nist.py                           ║
# ║ Description : automated evaluation of QCxMS results            ║ 
# ║               creates evaluate/                                ║
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝


import os, sys
sys.path.append(os.path.abspath(".."))
from jxms_config import (
    jxms_projhome, jxms_infohome, 
    jxms_nisthome, jxms_calcshome, 
    jxms_codeBIN, jxms_evalBIN, 
    jxms_initBIN, jxms_paramBIN
)
sys.path.append(jxms_initBIN)
import Daphnis.methods

from datetime import datetime
import numpy as np
import pandas as pd
import glob

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

import scipy
from scipy.spatial.distance import cosine
from collections import defaultdict
import re

working_dir = os.getcwd()
sys.path.append(os.getcwd())




def parse_filename( theo_jdx_in):

    cid_id = theo_jdx_in.rsplit('.', 1)[0]
    base_name = theo_jdx_in.rsplit('.', 1)[0]
    parts = base_name.split('_')
    bexp_id = parts[0] + "_" + parts[1]  
    
    protomer_index = parts.index("protomer") 
    mol_name = "_".join(parts[2:protomer_index])  # everything between bexp_id and "protomer"
    
    small_parts = parts[protomer_index+1].split('-')
    prot_num = small_parts[0]
    prot_id = f"prot_{prot_num}"
    param_id = "-".join(small_parts[1:]) 
    short_id = f"{prot_id}-{param_id}"

    return cid_id, bexp_id, mol_name, prot_id, param_id, short_id




def read_spectrum_file(filename):
    
    """used in : def read_info_files()
    returns a dictionary"""
    
    data = {}
    with open(filename, 'r') as file:
        
        for line in file:
            match = re.match(r'([\w_]+)\s*=\s*([\d\.]+)', line)
            if match:
                ion_type = match.group(1).strip()
                value = float(match.group(2))
                data[ion_type] = value
    return data




def read_info_files(file1,file2):

    """used in : def lets_eval() for ### getting match info """ 
    
    spectrum1 = read_spectrum_file(file1)
    spectrum2 = read_spectrum_file(file2)

    return spectrum1, spectrum2



def find_matching_peaks(spectrum1, spectrum2):
    
    """used in : def lets_eval() for ### getting match info """ 
    
    mol_ion_1 = spectrum1.get('Molecular_ion', 0)
    mol_ion_2 = spectrum2.get('Molecular_ion', 0)
    
    mol_ion_diff = abs(mol_ion_1 - mol_ion_2)
    
    if mol_ion_diff <= 0.009:
        print(mol_ion_diff)
        print("Correct comparison: Molecular ion values are within the threshold.")
        pass
    else:
        print("!!! Incorrect comparison:") 
        print("!!! Molecular ion values differ beyond the threshold.")
        exit()

    spectrum1.pop('Molecular_ion', None)
    spectrum2.pop('Molecular_ion', None)
    matching_peaks = {}
    
    for key1, val1 in spectrum1.items():
        
        for key2, val2 in spectrum2.items():
            if abs(val1 - val2) <= 0.009:
                matching_peaks[key1, key2] = (val1, val2)
    
    return mol_ion_1, mol_ion_2, matching_peaks




def save_matching_peaks( mol_ion_1, mol_ion_2, matching_peaks, theo_jdx_in, nist_jdx_in):
    
    """used in : def lets_eval() for ### getting match info """ 
    
    nist_id = nist_jdx_in.split('_')[-1].split('.')[0]
    output_file="match_"+nist_id+".out"
    
    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in)
    
    with open(output_file, 'w') as f:
        f.write("="*65 + "\n")
        f.write(f"{cid_id}\n")
        # this is without the mol_ion_peak)
        f.write(f"Matching_peak_num : {len(matching_peaks)}\n\n") 
        f.write(f"{'Peak Type 1':<15}{param_id:>15}{nist_id:>15}{'Peak Type 2':>20}\n")
        f.write("-"*65 + "\n")
        f.write(f"{'mol_ion':<15}{mol_ion_1:>15.5f}{mol_ion_2:15.5f}{'mol_ion':>20}\n")
        
        for (key1, key2), (value1, value2) in matching_peaks.items():
            f.write(f"{key1:<15}{value1:>15.5f}{value2:>15.5f}{key2:>20}\n")
        f.write(f"\n")
        matching_peak_num = len(matching_peaks)
    
    return cid_id, nist_id, matching_peak_num




def parse_mol_info(file_path):
    
    """used in : def lets_eval() for ### gethering mol_info (step1)""" 

    mol_name = mol_id = mol_cid = mol_inchi = molion = mol_xm = mol_mf = nistnumber = None
    try:
        with open(file_path, 'r') as file:
            
            for line in file:
                line = line.strip()
                # Skip header or non-data lines
                
                if line.startswith("###__") or not line:
                    continue
                
                if '=' in line:
                    key, value = map(str.strip, line.split('=', 1))
                    # Assign values to the corresponding variables
                    
                    if key == 'mol_name':
                        mol_name = value
                    elif key == 'mol_id':
                        mol_id = value
                    elif key == 'mol_cid':
                        mol_cid = value
                    elif key == 'mol_inchi':
                        mol_inchi = value
                    elif key == 'molion':
                        try:
                            molion = float(value)
                        except ValueError:
                            molion = None
                    elif key == 'mol_xm':
                        try:
                            mol_xm = float(value)
                        except ValueError:
                            mol_xm = None
                    elif key == 'mol_mf':
                        mol_mf = value
                    elif key == 'nistnumber':
                        nistnumber = value
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return mol_name, mol_id, mol_cid, mol_inchi, molion, mol_xm, mol_mf, nistnumber




def process_molion(molion, dp):
    
    """used in : def lets_eval() for ### comparing qcxms2nist (step?)""" 

    if molion is not None:
        mol_ion = round(molion, dp)
    else:
        print("Error: molion is None")
        mol_ion = 0.0  # Default value or handle the error
    
    return mol_ion




def parse_QCxMSjdx_df(infilename):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist (step?) 

    Parses one .jdx file and 
    Returns a pandas DataFrame with m/z and intensity columns."""
    
    spectrum = []
    reading_data = False
    
    with open(infilename, 'r', encoding='utf-8') as file:
        
        for line in file:
            line = line.strip()
            
            if "##PEAK TABLE" in line:
                reading_data = True
                continue
            
            if reading_data:
                if line.startswith('##') or line == '':
                    break
                parts = line.split()
                if len(parts) == 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])*10
                        spectrum.append({'mz': mz, 'intensity': intensity})
                    except ValueError:
                        print(f"Skipping line due to value error: {line}")
                        continue
    
    df = pd.DataFrame(spectrum)
    return df




def sort_nistfiles_by_instrument(file_paths):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist (step?) 
    
    Sorts a list of file paths by instrument type then by collision energy
    (HCD first, ITFT second, QTOF last)(CE_*) """
    
    #by instrument type
    files_by_instrument = defaultdict(list)

    for file_path in file_paths:
        
        if '-HCD' in file_path:
            files_by_instrument['HCD'].append(file_path)
        
        elif '-ITFT' in file_path:
            files_by_instrument['ITFT'].append(file_path)
        
        elif '-QTOF' in file_path:
            files_by_instrument['QTOF'].append(file_path)


    def extract_collision_energy(file_name):

        """Sort files within each instrument group by collision energy """

        match = re.search(r'CE(\d+)', file_name)
        
        # Return infinity if no match is found
        return int(match.group(1)) if match else float('inf')  


    sorted_files = []
    
    for instrument_type in ['HCD', 'QTOF', 'ITFT']:
        # Sort files within each instrument type by collision energy value
        sorted_files.extend(sorted(files_by_instrument[instrument_type], key=extract_collision_energy))
    
    return sorted_files




def parse_NISTjdx_df(infilename):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist (step?) 
    
    Parses a .jdx file and returns a pandas DataFrame 
    with mz and intensity columns."""
    
    spectrum = []
    reading_data = False
    
    try:
        with open(infilename, 'r') as file:
            
            for line in file:
                line = line.strip()
                
                if "##XYDATA" in line:
                    reading_data = True
                    continue
                
                if reading_data:
                    if line.startswith('##') or line == '':
                        break
                    #spectrum data
                    parts = line.split()
                    
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            spectrum.append({'mz': mz, 'intensity': intensity})
                        except ValueError:
                            print(f"Skipping line due to value error: {line}")
                            continue
                    else:
                        print(f"Unexpected line format: {line}")
    
    except FileNotFoundError:
        print(f"Error: The file '{infilename}' was not found.")
    
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return pd.DataFrame(spectrum)




def clean_dfspectra( df1, df2, decimal_place_included):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist -#clean 
    
    # include til
    # 1 decimal places  -> round(1)  -> tolereance = 0.0999
    # 2 decimal places  -> round(2)  -> tolereance = 0.00999
    # dp decimal places -> round(dp) -> tolerance = 0.99999/(10**dp)"""

    dp = decimal_place_included
    
    df1['mz'] = df1['mz'].round(dp)
    df2['mz'] = df2['mz'].round(dp)
    
    tolerance = 0.99999/(10**dp)
    
    return df1, df2, tolerance




def merge_with_tolerance( df1, df2, tolerance):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist-#clean""" 

    merged_data = []
    df2_used = np.zeros(len(df2), dtype=bool)
    
    for i, row1 in df1.iterrows():
        mz1 = row1['mz']
        intensity1 = row1['intensity']
        matched = False
        
        for j, row2 in df2.iterrows():
            
            if not df2_used[j] and abs(mz1 - row2['mz']) <= tolerance:
                avg_mz = round((mz1 + row2['mz']) / 2, 2)
                merged_data.append([avg_mz, intensity1, row2['intensity']])
                df2_used[j] = True
                matched = True
                break
        
        if not matched:
            merged_data.append([mz1, intensity1, 0])
    
    for j, row2 in df2.iterrows():
        
        if not df2_used[j]:
            merged_data.append([row2['mz'], 0, row2['intensity']])
    

    merged_df = pd.DataFrame(merged_data, columns=['mz', 'intensity_1', 'intensity_2'])
    
    # Remove duplicates and average intensities if necessary
    merged_df = merged_df.groupby('mz').agg({'intensity_1': 'sum', 'intensity_2': 'sum'}).reset_index()
    return merged_df




def normalize_peaks(merged_df):

    """used in : def lets_eval() for ### comparing qcxms 2 nist-#normalize the peaks""" 

    # Extract intensity vectors in array
    mz = merged_df.iloc[:, 0].values
    intensities_1 = merged_df.iloc[:, 1].values
    intensities_2 = merged_df.iloc[:, 2].values

    # Normalize the intensity vectors
    norm_1 = np.linalg.norm(intensities_1)
    norm_2 = np.linalg.norm(intensities_2)
    #if norm_1 == 0 or norm_2 == 0:
    #    return 0.0
    
    normalized_intensities_1 = intensities_1 / norm_1
    normalized_intensities_2 = intensities_2 / norm_2
    
    merged_df['normal_1']= normalized_intensities_1.round(3)
    merged_df['normal_2']= normalized_intensities_2.round(3)
    
    return merged_df




def calculate_cosine_similarity(merged_df):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore 
    
    Calculates the cosine similarity score using a tolerance for m/z matching."""

    intensities_1 = merged_df['intensity_1'].values
    intensities_2 = merged_df['intensity_2'].values
    
    norm_1 = np.linalg.norm(intensities_1)
    norm_2 = np.linalg.norm(intensities_2)
    normalized_intensities_1 = intensities_1 / norm_1
    normalized_intensities_2 = intensities_2 / norm_2
    
    merged_df['normal_1']= normalized_intensities_1.round(3)#*1000
    merged_df['normal_2']= normalized_intensities_2.round(3)#*1000
    
    cosine_similarity = np.dot(normalized_intensities_1, normalized_intensities_2)
    
    # Check if the similarity score is close to 1
    if np.isclose(cosine_similarity, 1, atol=1e-8):
        print("Similarity score is close to 1")
    
    return cosine_similarity*1000




def calculate_mass_weighted_dot_product(merged_df):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore"""
    
    merged_df['weight_1'] = (merged_df['intensity_1'] ** 0.6) * (merged_df['mz'] ** 0.3)
    merged_df['weight_2'] = (merged_df['intensity_2'] ** 0.6) * (merged_df['mz'] ** 0.3)

    numerator = np.sum(merged_df['weight_1'] * merged_df['weight_2']) ** 2
    denominator = np.sum(merged_df['weight_1'] ** 2) * np.sum(merged_df['weight_2'] ** 2)
    
    merged_df['weight_1'] = merged_df['weight_1'].round(3)
    merged_df['weight_2'] = merged_df['weight_2'].round(3)
    
    WDot = np.sqrt(numerator / denominator)
    
    return WDot*1000




def read_ms_spec( infilename, normalize_to=0):
    
    """used in : def generate_predicted_peaks_tohave() """

    spec_pd = pd.read_csv(infilename, comment='#', header=None, names=("mz", "intensity"), dtype={"mz": np.float64, "intensity": np.float64}, skipinitialspace=True, sep=r'\s+')
    spec_mz_np = spec_pd.mz.to_numpy()
    spec_intensity_np = spec_pd.intensity.to_numpy()
    num_mzs = len(spec_mz_np)
    spec_mz = spec_mz_np.reshape(num_mzs,1)
    spec_intensity = spec_intensity_np.reshape(num_mzs,1)
    
    if normalize_to != 0:
        spec_intensity /= spec_intensity.max()
        spec_intensity *= normalize_to
    
    return spec_pd, spec_mz, spec_intensity




def calc_entropy_simscore(merged_df):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore"""
    
    theo_array = np.array(merged_df[['mz','intensity_1']])
    nist_array = np.array(merged_df[['mz','intensity_2']])
    
    entropy_simscore = 1 - Daphnis.methods.distance(theo_array, nist_array, method="weighted_entropy", ms2_da=0.05)
    
    entropy_simscore_scaled = entropy_simscore*1000
    
    return entropy_simscore_scaled




def write_comparison_df(merged_df, clean_merged_df, theo_jdx_in, nist_id):

    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore"""

    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in)
    
    raw_outfilename = "compare_raw."+cid_id+"."+nist_id+"."+"out"
    
    mz          = merged_df.iloc[:, 0].values
    raw_theo    = merged_df.iloc[:, 1].values
    raw_nist    = merged_df.iloc[:, 2].values
    normal_theo = merged_df.iloc[:, 3].values*1000
    normal_nist = merged_df.iloc[:, 4].values*1000
    weight_theo = merged_df.iloc[:, 5].values
    weight_nist = merged_df.iloc[:, 6].values
    
    raw_merged_df =pd.DataFrame({'mz': mz,'raw_theo': raw_theo, 'raw_nist': raw_nist, 'normal_theo':normal_theo, 'normal_nist':normal_nist, 'weight_theo':weight_theo, 'weight_nist':weight_nist})
    
    with open(raw_outfilename, 'w') as file:
        file.write(raw_merged_df.to_string())
        file.write('\n')
    
    clean_outfilename = "compare_clean."+cid_id+"."+nist_id+"."+"out"
    mz          = clean_merged_df.iloc[:, 0].values
    raw_theo    = clean_merged_df.iloc[:, 1].values
    raw_nist    = clean_merged_df.iloc[:, 2].values
    normal_theo = clean_merged_df.iloc[:, 3].values*1000
    normal_nist = clean_merged_df.iloc[:, 4].values*1000
    weight_theo = clean_merged_df.iloc[:, 5].values
    weight_nist = clean_merged_df.iloc[:, 6].values
    
    clean_merged_df =pd.DataFrame({'mz': mz,'raw_theo': raw_theo, 'raw_nist': raw_nist, 'normal_theo':normal_theo, 'normal_nist':normal_nist, 'weight_theo':weight_theo, 'weight_nist':weight_nist})

    with open(clean_outfilename, 'w') as file:
        file.write(clean_merged_df.to_string())
        file.write('\n')
        #file.write(df.to_string(index=False))




def get_peaks_tohave_minmax(peaks_tohave, peaks_tohave_thresh ):
    
    """used in : def generate_predicted_peaks_tohave() """

    peaks_tohave_min=[]
    peaks_tohave_max=[]
    
    for i in peaks_tohave:
        check_peak = float(i)
        check_peak_min = i - peaks_tohave_thresh
        check_peak_max = i + peaks_tohave_thresh
        
        peaks_tohave_min.append(check_peak_min)
        peaks_tohave_max.append(check_peak_max)
    
    return peaks_tohave_min, peaks_tohave_max



"""not used right now !!! 
def get_predicted_peaks(peaks_tohave, peaks_tohave_min, peaks_tohave_max, theo_jdx  ):
    predicted = []
    for j in range(len(peaks_tohave)):
        #print("checking ", peaks_tohave[j])
        for i in range(len(theo_jdx)):
            if peaks_tohave_min[j] < theo_jdx[i,0] < peaks_tohave_max[j]:
                predicted.append(theo_jdx[i,0])
    return predicted
"""



def check_peaks(theo_spec_pd, peaks_tohave, peak_min, peak_max ):
    
    """used in : def generate_predicted_peaks_tohave() """

    theo_spec_np = theo_spec_pd.to_numpy()
    theo_spec_np = np.round(theo_spec_np, 5)
    #theo_spec_intensity = theo_spec_pd.intensity
    #theo_spec_intensity *= 1000
    keep_theo_peaks=[]
    keep_nist_peaks=[]
    #print('Nist peaks , Predicted')
    
    for i in range(len(theo_spec_np)) :
        
        for j in range(len(peaks_tohave)):
            
            if peak_min[j] <  theo_spec_np[i,0] < peak_max[j]:
                #print(peaks_tohave[j] , theo_spec_np[i,0] )
                #print(peaks_tohave[j] , theo_spec_np[i,0], theo_spec_intensity[i])
                
                keep_theo_peaks.append(theo_spec_np[i])
                keep_nist_peaks.append(peaks_tohave[j])
    
    np.set_printoptions(suppress=True)
    theo_peaks = np.round(keep_theo_peaks, 5)
    nist_peaks = np.round(keep_nist_peaks, 5)
    missing_peaks = np.setdiff1d(peaks_tohave, keep_nist_peaks)
    
    return theo_peaks, nist_peaks, missing_peaks




def save_peaks( peaks_outfilename, peaks_tohave, missing_peaks, keep_nist_peaks, keep_theo_peaks):
    
    """used in : def generate_predicted_peaks_tohave() """
    
    f = open(peaks_outfilename,'w')
    f.write('# nist_peaks_tohave    : mz = %s \n' % peaks_tohave )
    f.write('# theo_peaks_missing  : mz = %s \n' % missing_peaks )
    f.write('#  nist_mz ,  theo_mz , theo_intensity\n')
    
    for nist_peak, (mz, intensity) in zip(keep_nist_peaks, keep_theo_peaks):
        f.write(f"{nist_peak:10.2f} {mz:10.2f} {intensity*1000:10.2f}\n")
    
    f.close()
    #print('\n* Filename:', peaks_outfilename, ' is created.\n')




def save_simscores( simscores_outfilename, theo_jdx_in, nist_jdx_in, cos_score, dot_score, ent_score ):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore"""
    
    f= open(simscores_outfilename, 'a')
    f.write(f"{theo_jdx_in:25s} {nist_jdx_in:12s} {cos_score:>5.4f} {dot_score:>5.4f} {ent_score:>4.5f}\n")
    
    f.close()




def prettyprint_simscores(raw_simscores_outfilename, prettyprint_outfilename, theo_jdx_in):
    
    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore"""

    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in) 
    
    with open(raw_simscores_outfilename, 'r') as infile:
        infilelines = infile.readlines()
    
    headers = ["QCxMS_spec", "NIST23_spec", "Cos", "Wdot", "Entropy"]
    
    with open(prettyprint_outfilename, 'w') as outfile:
        outfile.write("=" * 68 + "\n\n")
        outfile.write(cid_id+'\n\n')
        outfile.write(f"{headers[0]:22} {headers[1]:>11} {headers[2]:>9} {headers[3]:>9} {headers[4]:>9}\n")
        outfile.write("-" * 68 + "\n")
        
        for line in infilelines:
            columns = line.split()
            outfile.write(f"{short_id:<22} {columns[1]:<11} {columns[2]:>9} {columns[3]:>9} {columns[4]:>9}\n")




def prettyprint_result_out(raw_simscores_outfilename, prettyprint_result_out, theo_jdx_in, matching_peak_num_list):

    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore"""
    
    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in)
    
    with open(raw_simscores_outfilename, 'r') as infile:
        infilelines = infile.readlines()
    
    headers = ["QCxMS_spec", "NIST23_spec", "Cos", "Wdot", "Entropy", "Match"]
    with open(prettyprint_result_out, 'w') as outfile:
        outfile.write("=" * 68 + "\n\n")
        outfile.write(cid_id+'\n\n')
        outfile.write(' ')
        outfile.write(f"{headers[0]:^21} {headers[1]:<11} {headers[2]:^8} {headers[3]:^8} {headers[4]:>8} {headers[5]:<4}\n")
        outfile.write("-" * 68 + "\n")
        i = 0
        for line in infilelines:
            if not line.strip():
                continue  # Skip empty lines
            columns = line.split()
            outfile.write(f"{param_id:<21} {columns[1]:11} {columns[2]:>8} {columns[3]:>8} {columns[4]:>8} {matching_peak_num_list[i]:>4}\n")
            i += 1
        outfile.write(f"\n")




def generate_predicted_peaks_tohave( peaks_infilename, peak_thresh, peaks_outfilename ):

    """used in : def save_peaks() """

    peaks_tohave = np.genfromtxt(peaks_infilename, 
                                 delimiter=',', dtype=float)
    peak_min, peak_max = get_peaks_tohave_minmax(peaks_tohave, 
                                                 peak_thresh)
    
    keep_theo_peaks, keep_nist_peaks, keep_missing_peaks = check_peaks(
            theo_spec, peaks_tohave, peak_min, peak_max)
    
    save_peaks(peaks_outfilename, peaks_tohave, keep_missing_peaks, 
               keep_nist_peaks, keep_theo_peaks)





""" now plotting from HERE """

def plot_headtotail(theo_jdx_in, nist_jdx_in, saveto=None):
    qcjdx, theo_mz, theo_intensity = read_ms_spec(theo_jdx_in, normalize_to=100) #top aka target = theo 
    nijdx, nist_mz, nist_intensity = read_ms_spec(nist_jdx_in, normalize_to=-100)  #bottom aka ref.= nist

    fig = plt.figure(figsize=(20,14), dpi=350)
    matplotlib.rc("font", **{ "size": 28 })
    ax = fig.add_subplot( 111)

    qc_spec = ax.bar( qcjdx.mz, qcjdx.intensity, width=0.8, label='Theoretical spectrum', color='magenta')
    ni_spec = ax.bar( nijdx.mz, nijdx.intensity, width=0.8, label='Experimental spectrum', color='blue')

    #ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    #ax.spines['left'].set_smart_bounds(True)
    #ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_ylabel('Intensity')
    ax.set_xlabel('m/z')
    ax.xaxis.set_label_coords(0.96, 0.48)

    #ax.legend(fontsize=20,loc='upper left')
    ax.legend(fontsize=20,loc='best')

    cid_id = theo_jdx_in.split('.')[0]
    nist_id = nist_jdx_in.split('_')[-1].replace('.JDX','')
    ax.set_title("Molecule " + cid_id+".h2t."+nist_id, fontsize=14)
#    ax.set_title("Molecule " + 
#            theo_jdx.split('/')[-1].replace('.jdx',''), 
#            fontsize=14)

    ax.set_xticks([])
    ax.tick_params(axis='y', labelsize=20)
    label_peaks_headtotail(ax, qc_spec)
    label_peaks_headtotail(ax, ni_spec)
    ax.set_yticks(list(range(-100,101,25)))
    
    ticks = ax.get_yticks()
    #ax.set_yticklabels([int(abs(tick)) for tick in ticks])

    ax.yaxis.set_major_locator(FixedLocator(ticks))  # Use FixedLocator to set tick positions
    ax.set_yticklabels([int(abs(tick)) for tick in ticks])  # Set labels accordingly

    trim_xaxis = True
    if trim_xaxis:
        #min_mz = min( qcjdx.mz.min(), nijdx.mz.min() )
        min_mz = 0.
        # defined as the mz with the maximum intensity
        max_mz = max( qcjdx.mz.max(), nijdx.mz.max() )
        #print("Trimming plot to mz = [",min_mz, max_mz, "]")
        #ax.set_xlim(5+min_mz, 1.1* max_mz)
        #print("Trimming plot to mz = [ 0.,", max_mz, "]")
        ax.set_xlim(min_mz, 1.1* max_mz)
        #ax.set_xbound(lower=0.0, upper=1.1* max_mz)
    fig.tight_layout()
    if saveto:
        fig.savefig(saveto)
    matplotlib.pyplot.close()
    return fig




def peaks_headtotail( x, y, separation=1., ratio=1.01, threshold=0.0):
    """
    separation means the peak should be separated in mz by +- this number
    ratio means the peak should be this many times greater than the other points
    threshhold means the intensity should be at least this number to be a peak
    """
    assert separation > 0 and ratio > 0
    x = np.asarray(x)
    y = np.asarray(y)

    # sort them to make sure we always prefer the largest peaks first
    srt = np.abs(y).argsort()[::-1]

    xp = []
    yp = []

    for i in range( len(x)):
        i = srt[i]
        mask = np.abs(x - x[i])
        mask = mask < separation
        if not mask.any():
            #print("skipping due to empty mask")
            continue
        group = y[mask]

        # quick way to just get the 2 largest values
        if len(group) > 1:
            max_indices = np.argpartition(np.abs(group), len(group) - 2)[-2:]
            xval = x[mask][ max_indices[1]]
            yval = group[ max_indices[1]]
            second_max = group[ max_indices[0]]
        else:
            idx = group.argmax()
            xval = x[mask][ idx]
            yval = group[ idx]
            second_max = 0

        if second_max == 0:
            ratio_i = ratio * 2
        else:
            ratio_i = yval / second_max
        uniq = [ abs(xpi - xval) > separation for xpi in xp]
        #print("******ADD", xval, yval)
        xp.append( xval)
        yp.append( yval)

    return xp, yp




def label_peaks_headtotail(ax, rects, peaks_only=True):

    """
    Attach a text label above each bar displaying its intensity
    """
    dat = [[rect.get_x() + rect.get_width()/2., rect.get_height()] for rect in rects]
    if peaks_only:
        x = [v[0] for v in dat]
        y = [v[1] for v in dat]
        xp, yp = peaks_headtotail( x, y)
    else:
        xp = [v[0] for v in dat if v[1] > 0.0]
        yp = [v[1] for v in dat if v[1] > 0.0]

    for mz, intensity in zip( xp, yp):
        ax.annotate('%g' % mz, (mz, np.sign(intensity)*7 + intensity), ha='center', va='center', fontsize=18)




def plot_headtotail_slide(theo_jdx_in, nist_jdx_in, saveto=None):
    
    qcjdx, theo_mz, theo_intensity = read_ms_spec(theo_jdx_in, normalize_to=100) 
    nijdx, nist_mz, nist_intensity = read_ms_spec(nist_jdx_in, normalize_to=-100)  
    
    fig = plt.figure(figsize=(20,14), dpi=350)
    matplotlib.rc("font", **{ "size": 48 })
    ax = fig.add_subplot( 111)
    
    qc_spec = ax.bar( qcjdx.mz, qcjdx.intensity, width=0.8, label='Theoretical spectrum', color='magenta')
    ni_spec = ax.bar( nijdx.mz, nijdx.intensity, width=0.8, label='Experimental spectrum', color='blue')

    
    #ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    #ax.spines['left'].set_smart_bounds(True)
    #ax.spines['bottom'].set_smart_bounds(True)
    
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylabel('Intensity')
    ax.set_xlabel('m/z')
    ax.xaxis.set_label_coords(0.96, 0.48)
    ax.legend(fontsize=40,loc='best')

    cid_id = theo_jdx_in.split('.')[0]
    nist_id = nist_jdx_in.split('_')[-1].replace('.JDX','')
    ax.set_title("Molecule " + cid_id+".h2t."+nist_id, fontsize=14)


    ax.set_xticks([])
    ax.tick_params(axis='y', labelsize=30)
    label_peaks_headtotail_slide(ax, qc_spec)
    label_peaks_headtotail_slide(ax, ni_spec)
    #ax.set_yticks(list(range(-100,101,25))) 
     
    ticks = ax.get_yticks()
    #ax.set_yticklabels([int(abs(tick)) for tick in ticks])
    
    ax.yaxis.set_major_locator(FixedLocator(ticks)) 
    ax.set_yticklabels([int(abs(tick)) for tick in ticks])  
    
    trim_xaxis = True
    if trim_xaxis:
        #min_mz = min( qcjdx.mz.min(), nijdx.mz.min() )
        min_mz = 0.
        # defined as the mz with the maximum intensity
        max_mz = max( qcjdx.mz.max(), nijdx.mz.max() )
        #print("Trimming plot to mz = [",min_mz, max_mz, "]")
        #ax.set_xlim(5+min_mz, 1.1* max_mz)
        #print("Trimming plot to mz = [ 0.,", max_mz, "]")
        ax.set_xlim(min_mz, 1.1* max_mz)
        #ax.set_xbound(lower=0.0, upper=1.1* max_mz)
    fig.tight_layout()
    if saveto:
        fig.savefig(saveto)
    matplotlib.pyplot.close()
    return fig




def peaks_headtotail_slide( x, y, separation=3., ratio=1.01, threshold=0.0):
    """
    separation means the peak should be separated in mz by +- this number
    ratio means the peak should be this many times greater than the other points
    threshhold means the intensity should be at least this number to be a peak
    """
    
    assert separation > 0 and ratio > 0
    x = np.asarray(x)
    y = np.asarray(y)

    srt = np.abs(y).argsort()[::-1]
    xp = []
    yp = []
    
    for i in range( len(x)):
        i = srt[i]
        mask = np.abs(x - x[i])
        mask = mask < separation
   
        if not mask.any():
            #print("skipping due to empty mask")
            continue
        
        group = y[mask]
        # quick way to just get the 2 largest values
        
        if len(group) > 1:
            max_indices = np.argpartition(np.abs(group), len(group) - 2)[-2:]
            xval = x[mask][ max_indices[1]]
            yval = group[ max_indices[1]]
            second_max = group[ max_indices[0]]
        
        else:
            idx = group.argmax()
            xval = x[mask][ idx]
            yval = group[ idx]
            second_max = 0
        #print("Considering xval=", xval, "yval=", yval, "second_max=", second_max)
        
        if i > 0:
            if np.abs( y[ i]) < np.abs( y[ i-1]):
##                print("SKIPPING", xval, y[ i], y[ i-1])
                continue
        
        if i < len(y) - 1:
            if np.abs( y[ i]) < np.abs( y[ i+1]):
##                print("SKIPPING", xval, y[ i], y[ i+1])
                continue
        
        if second_max == 0:
            ratio_i = ratio * 2
        
        else:
            ratio_i = yval / second_max

##        print("RATIO=", ratio_i, "THRESH=", np.abs(yval) > threshold)
        
        if ratio_i > ratio and np.abs(yval) > threshold:
            uniq =  [ abs(xpi - xval) > separation for xpi in xp]
            
            if all(uniq):
##                print("******ADD", xval, yval)
                xp.append( xval)
                yp.append( yval)
        
        uniq = [ abs(xpi - xval) > separation for xpi in xp]
        #print("******ADD", xval, yval)
        xp.append( xval)
        yp.append( yval)
    
    #print(" DONE")
    #print("xp = ", xp)
    #print("yp = ", yp)
    return xp, yp




def label_peaks_headtotail_slide(ax, rects, peaks_only=True):
    """
    Attach a text label above each bar displaying its intensity
    """

    dat = [[rect.get_x() + rect.get_width()/2., rect.get_height()] for rect in rects]
    
    if peaks_only:
        x = [v[0] for v in dat]
        y = [v[1] for v in dat]
        xp, yp = peaks_headtotail_slide( x, y)
    else:
        xp = [v[0] for v in dat if v[1] > 0.0]
        yp = [v[1] for v in dat if v[1] > 0.0]
    
    for mz, intensity in zip( xp, yp):
        ax.annotate('%g' % mz, (mz, np.sign(intensity)*7 + intensity), ha='center', va='center', fontsize=30)



def generate_headtotail(theo_jdx_in, nist_jdx_in):
    
    cid_id = theo_jdx_in.split('.')[0]
    nist_id = nist_jdx_in.split('_')[-1].replace('.JDX','')
    h2t_outgraph = cid_id+".h2t."+nist_id+".png"
    
    fig = plot_headtotail(theo_jdx_in, nist_jdx_in, saveto= cid_id+".h2t."+nist_id+".png")   # must read it in here due to the range(100, -100)
    
    print('\n* Filename:', h2t_outgraph, '     (raw ver.  ) is created.\n')
    return fig




def generate_headtotail_slide(theo_jdx_in, nist_jdx_in):
    
    cid_id = theo_jdx_in.split('.')[0]
    nist_id = nist_jdx_in.split('_')[-1].replace('.JDX','')
    
    fig = plot_headtotail_slide(theo_jdx_in, nist_jdx_in, saveto= cid_id+".h2t_slide."+nist_id+".png")
    #print('\nHead-to-tail graphs (slide ver.) are generated.\n')
    
    h2t_outgraph = cid_id+".h2t_slide."+nist_id+".png"
    print('\n* Filename:', h2t_outgraph, '(slide ver.) is created.\n')
    
    return fig


""" done plotting to HERE """





def save_check_match_to_file( outfilename):
    
    match_files = glob.glob(f'match_CE*.out')
    match_files.sort()
    found_covered = False  # track if any file is covered
    
    with open(outfilename, 'w') as f:
        
        for match_file in match_files:
            
            if check_all_matched(match_file):
                #print("Found covered: ", match_file)
                f.write(f"FOUND COVERED: {match_file}\n")
                found_covered = True  # Set flag if we find any covered file
            else:
                #print("Not covered: ", match_file)
                f.write(f"Not covered: {match_file}\n")
    
    if found_covered:
        
        with open('COVERED.out', 'w') as covered_file:
            covered_file.write("Covered files were found.\n")
            
            for match_file in match_files:
                
                if check_all_matched(match_file):
                    covered_file.write(f"FOUND COVERED: {match_file}\n")
                    print("\nCOVERED!!\n")





def get_best_entropy_out(infilename, outfilename, theo_jdx_in):

    """used in : def lets_eval() for ### comparing qcxms 2 nist-#calc simscore"""
    
    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in)
    
    with open(infilename, 'r') as file:
        lines = file.readlines()
    
    start_line = lines[0]
    emp_1 = lines[1]
    cid_id = lines[2]
    emp_2 = lines[3]
    title = lines[4]
    chart = lines[5]
    
    data_pd = pd.read_csv(infilename, header=3, sep=r'\s+',names=['QCxMS_spec', 'NIST23_spec', 'Cos', 'Wdot', 'Entropy', 'Match'],dtype={'QCxMS_spec': str, 'NIST24_spec': str, 'Cos': np.float64, 'Wdot': np.float64, 'Entropy': np.float64, 'Match': np.int32})
    
    highest_entropy_row = data_pd.loc[data_pd['Entropy'].idxmax()]
    
    with open(outfilename, 'w') as f:
        f.write(f"{start_line}")
        f.write(f"{emp_1}")
        f.write(f"{cid_id}\n")
        f.write(f"{chart}")
        f.write(f"{highest_entropy_row['QCxMS_spec']}   "
                f"{highest_entropy_row['NIST23_spec']}   "
                f"{highest_entropy_row['Cos']:>6.2f}   "
                f"{highest_entropy_row['Wdot']:>6.2f}   "
                f"{highest_entropy_row['Entropy']:>6.2f}   "
                f"{int(highest_entropy_row['Match']):>4}\n")
        f.write(f"{emp_2}")
    
    highest_matching_peak_num = data_pd.loc[data_pd['Match'].idxmax()]    
    
    matching_peak_num = highest_matching_peak_num['Match']
    
    # Check if the value is greater than 1
    if matching_peak_num > 1:
        #print(f"The highest matching peak number ({matching_peak_num}) is greater than 1.")
        covered_filename = 'COVERED.out'
        
        with open(covered_filename, 'w') as f:
            f.write(f"{cid_id}")
            f.write(f"COVERED!")
    else:
        #covered_filename = 'notcovered.out'
        #
        #with open(covered_filename, 'w') as f:
        #    f.write(f"{cid_id}")
        #    f.write(f"nope!")
        pass
        print(f"The highest matching peak number ({matching_peak_num}) is not greater than 1.")




def parse_result_file( file_path):

    """used in : def lets_eval() for ### getting result summary files """
    
    data = []
    
    with open(file_path, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        match = re.match(r"(\S+)\s+(\S+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+(\d+)", line)
        
        if match:
            qcxms_spec, nist_spec, cos, wdot, entropy, match_value = match.groups()
            data.append({
                "QCxMS_spec": qcxms_spec,
                "NIST_spec": nist_spec,
                "Cos": float(cos),
                "Wdot": float(wdot),
                "Entropy": float(entropy),
                "Match": int(match_value)
            })
        else:
            continue
            #print("no match")
    
    return data





def save_sorted_result( data, output_file, theo_jdx_in):
    
    """used in : def lets_eval() for ### getting result summary files """
    
    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in)
    
    headers = ["QCxMS_spec", "NIST23_spec", "Cos", "Wdot", "Entropy", "Match"]
    
    with open(output_file, "w") as outfile:
        outfile.write("=" * 68 + "\n\n")
        outfile.write(cid_id+'\n\n')
        outfile.write(' ')
        outfile.write(f"{headers[0]:^21} {headers[1]:<11} {headers[2]:^8} {headers[3]:^8} {headers[4]:>8} {headers[5]:<4}\n")
        outfile.write("-" * 68 + "\n")

        for entry in data:
            outfile.write(f"{entry['QCxMS_spec']:<21} {entry['NIST_spec']:<11} {entry['Cos']:^8} {entry['Wdot']:^8} {entry['Entropy']:>8} {entry['Match']:>4}\n")
        
        outfile.write(f"\n")






def save_best_result(data, output_file):
    
    """used in : def lets_eval() for ### getting result summary files """
    
    bexp_id = os.getcwd().split('/')[-5]
    
    if data:
        
        with open(output_file, 'w') as file:
            #file.write("QCxMS_spec\tNIST_spec\tCos\tWdot\tEntropy\tMatch\n")
            file.write(f"{bexp_id} {data[0]['QCxMS_spec']}  {data[0]['NIST_spec']}\t{data[0]['Cos']:.2f}\t{data[0]['Wdot']:.2f}\t{data[0]['Entropy']:.2f}\t{data[0]['Match']}\n")






def lets_eval(theo_jdx_in, decimal_place_included, mol_info_file):
    
    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in)
    dp = decimal_place_included
   
    mol_name, mol_id, mol_cid, mol_inchi, molion, mol_xm, mol_mf, nistnumber = parse_mol_info(mol_info_file)
    mol_ion = process_molion(molion, dp)
    nist_jdx_files = glob.glob(f'{Bexp_nisthome}/NIST23_bexp_{mol_id}_*/bexp_*.JDX')
    nist_info_files =glob.glob(f'{Bexp_nisthome}/NIST23_bexp_{mol_id}_*/info_bexp_{mol_id}_*.txt')

    sorted_nist_info_files = sort_nistfiles_by_instrument(nist_info_files)
    matching_peak_num_list = []
   

    for nist_info_file in sorted_nist_info_files:
        
        nist_info_id = nist_info_file.split('/')[-1].split('.')[0].split('_')[-1]
        parent_dir = os.path.dirname(working_dir)
        file1 = f'{parent_dir}/info_spectrum.txt'
        file2 = nist_info_file

        spec_info1, spec_info2 = read_info_files(file1,file2)
        mol_ion_1, mol_ion_2, matching_peaks = find_matching_peaks(spec_info1, spec_info2)
        
        cid_id, nist_id, matching_peak_num = save_matching_peaks(mol_ion_1, mol_ion_2, matching_peaks, theo_jdx_in, nist_info_file) 
        
        matching_peak_num_list.append( matching_peak_num) 


    sorted_nist_jdx_files = sort_nistfiles_by_instrument(nist_jdx_files)
    
    for i in sorted_nist_jdx_files:
        nist_jdx_in = i
        
        nistfile = nist_jdx_in
        nist_id =nistfile.split('/')[-1].split('.')[0].split('_')[-1]
        nist_spec = parse_NISTjdx_df(nistfile)
        nistfile  = nist_jdx_in
        nist_spec = parse_NISTjdx_df(nistfile)
        qcxmsfile = theo_jdx_in
        theo_spec = parse_QCxMSjdx_df(qcxmsfile)
       
        theo_clean, nist_clean, tolerance = clean_dfspectra(theo_spec, nist_spec, decimal_place_included)
        
        merged_df = merge_with_tolerance(theo_clean, nist_clean, tolerance)
        filtered_df = merged_df[merged_df['mz'] != mol_ion]
        clean_merged_df = filtered_df.copy()
        
        merged_df = normalize_peaks(merged_df)
        clean_merged_df = normalize_peaks(clean_merged_df)
        
        cos_score_withmolion = calculate_cosine_similarity(merged_df)
        dot_score_withmolion = calculate_mass_weighted_dot_product(merged_df)
        ent_score_withmolion = calc_entropy_simscore(merged_df)
        
        save_simscores('raw.out',theo_jdx_in, nist_id, cos_score_withmolion, dot_score_withmolion, ent_score_withmolion )
        prettyprint_result_out('raw.out','result.out', theo_jdx_in, matching_peak_num_list)

        result = parse_result_file('result.out')
        sorted_result = sorted(result, key=lambda x: (-x["Match"], -x["Entropy"]))  # Sort by Match (desc), then Entropy (desc)
        
        save_sorted_result(sorted_result, "result_sorted.out",theo_jdx_in)
        get_best_entropy_out('result_sorted.out', 'result_best.out', theo_jdx_in)
       
        
        cos_score_nomolion = calculate_cosine_similarity(clean_merged_df)
        dot_score_nomolion = calculate_mass_weighted_dot_product(clean_merged_df)
        ent_score_nomolion = calc_entropy_simscore(clean_merged_df)
        
        #### for debugging ###
        #write_comparison_df(merged_df, clean_merged_df, theo_jdx_in, nist_id)
        
        save_simscores('clean.out',theo_jdx_in, nist_id, cos_score_withmolion, dot_score_withmolion, ent_score_withmolion )
        prettyprint_result_out('clean.out','clean_result.out', theo_jdx_in, matching_peak_num_list)
        result = parse_result_file('clean_result.out')
        sorted_result = sorted(result, key=lambda x: (-x["Match"], -x["Entropy"]))
        save_sorted_result(sorted_result, "clean_result_sorted.out", theo_jdx_in)
        get_best_entropy_out('clean_result_sorted.out', 'clean_result_best.out', theo_jdx_in)




        #### draw graphs ###
        generate_headtotail(theo_jdx_in, nist_jdx_in)
        generate_headtotail_slide(theo_jdx_in, nist_jdx_in)




def main():

    theo_jdx_in            = sys.argv[1]
    
    cid_id, bexp_id, mol_name, prot_id, param_id, short_id = parse_filename(theo_jdx_in)
    
    decimal_place_included = int(sys.argv[2])
    mol_info_file = str(sys.argv[3])
    #print('mol_info_file :' , mol_info_file)

    current_dir = os.getcwd()
    date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    print(" ") 
    print('-'*50)
    print(date, 'starting evaluation')
    print(theo_jdx_in)
    #print("cid_id  :", cid_id)
    print("bexp_id :", bexp_id)
    print("mol_name:", mol_name)
    print("prot_id :", prot_id)
    print("short_id:", short_id)

    lets_eval(theo_jdx_in, decimal_place_included, mol_info_file)

    date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print('\n*** ', date, 'done with evaluation\n')
    print('-'*50)

    




if __name__ == "__main__":
    main()

