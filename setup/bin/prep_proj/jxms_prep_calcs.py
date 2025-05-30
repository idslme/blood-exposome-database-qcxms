#!/usr/bin/env python3
# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                        — by Jesi Lee —                         ║
# ║ Project     : Blood Exposome CID-MS/MS Prediction via QCxMS    ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : jxms_prep_calcs.py                               ║
# ║ Description : automated calcs/ setup based on info/            ║ 
# ║               requires an input file:"info/jxms_pubchem_cid.in"║
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝


import sys, os
import glob
import shutil
import subprocess
import requests
from datetime import datetime

sys.path.append(os.path.abspath(".."))
from jxms_config import (
    jxms_projhome, jxms_infohome, jxms_nisthome,
    jxms_calcshome, jxms_codeBIN, jxms_prepBIN,
    jxms_runBIN, jxms_evalBIN, jxms_paramBIN
)



### The parameters used for 121 blood exposome compounds (100-200 Da):

A_parameter_options = [
    "n2-maxcol4-40", "n2-maxcol4-25", "n2-maxcol4-60",
    "n2-maxcol6-40", "n2-maxcol6-25", "n2-maxcol6-60",
    "n2-setcol4-40", "n2-setcol4-25", "n2-setcol4-60",
    "n2-setcol6-40", "n2-setcol6-25", "n2-setcol6-60",
    "n2-40", "n2-esi-maxcol4-40", "n2-esi-setcol4-40"
]

C_parameter_options = [
    "n2-maxcol4-60", "n2-maxcol4-80", "n2-maxcol4-45",
    "n2-maxcol6-60", "n2-maxcol6-80", "n2-maxcol6-45",
    "n2-setcol4-60", "n2-setcol4-80", "n2-setcol4-45",
    "n2-setcol6-60", "n2-setcol6-80", "n2-setcol6-45",
    "n2-60", "n2-esi-maxcol4-60", "n2-esi-setcol4-60"
]





def gather_molinfo(jxms_infohome):

    with open(jxms_infohome + "jxms_cid.txt", "r") as file:
        cidnumber_list = [line.strip() for line in file.readlines()]
    with open(jxms_infohome + "jxms_name.txt", "r") as file:
        chemname_list = [line.strip() for line in file.readlines()]
    with open(jxms_infohome + "jxms_mf.txt", "r") as file:
        molformula_list = [line.strip() for line in file.readlines()]
    with open(jxms_infohome + "jxms_monomass.txt", "r") as file:
        molmonoweight_list = [line.strip() for line in file.readlines()]
    with open(jxms_infohome + "jxms_inchikey.txt", "r") as file:
        molinchi_list = [line.strip() for line in file.readlines()]
    with open(jxms_infohome + "jxms_AorC.txt", "r") as file:
        molAorC_list = [line.strip() for line in file.readlines()]

    return cidnumber_list, chemname_list, molformula_list, molmonoweight_list,molinchi_list, molAorC_list




def write_mol_info(cidnumber_list, chemname_list, molformula_list, molmonoweight_list,molinchi_list, molseq, molAorC):

    # Index is n-1 because Python lists are 0-based
    index = int(molseq) - 1

    if index < len(cidnumber_list) and index < len(chemname_list) and index < len(molformula_list) and index < len(molmonoweight_list) and index < len(molinchi_list):
        cid = cidnumber_list[index]
        name = chemname_list[index]
        formula = molformula_list[index]
        weight = molmonoweight_list[index]
        inchi = molinchi_list[index]

        mol_id = f"{index + 1:03d}"
        mol_ion = float(weight) + 1.0079

        with open("mol_info.txt", "w") as outfile:
            outfile.write(f"###__ jxms project: mol_info __###\n")
            outfile.write(f"mol_name  = {name}\n")
            outfile.write(f"mol_id    = {mol_id}\n")
            outfile.write(f"mol_cid   = {cid}\n")
            outfile.write(f"mol_inchi = {inchi}\n")
            outfile.write(f"molion    = {mol_ion:.4f}\n")
            outfile.write(f"mol_xm    = {weight}\n")
            outfile.write(f"mol_mf    = {formula}\n")
            outfile.write(f"mol_AorC  = {molAorC}\n")
            #outfile.write(f"nistnumber= {nist_number_path}\n")
            outfile.write("\n")
    else:
        print("Error: The line number is out of range.")




def set_mol_cidseq(cidnumber_list):

    molseq_list=[]
    molcid_list=[]

    for molseq, cidnumber in enumerate(cidnumber_list, start=1):
        molcid = cidnumber.strip() 
        formatted_molseq = f"{molseq:03}"
        molseq_list.append(formatted_molseq)
        molcid_list.append(molcid)
    
    return molseq_list, molcid_list




def set_molstruc(molseq, molcid, cidnumber_list, chemname_list, molformula_list, molmonoweight_list,molinchi_list, molAofC):

    bexp_folder = f"bexp_{molseq}"
    cid_folder  = f"cid_{molcid}"
    cid_path    = os.path.join(bexp_folder, cid_folder)
    prep_path   = os.path.join(bexp_folder, "prep")
    qcrun_path  = os.path.join(bexp_folder, "qcrun")
    
    os.makedirs(bexp_folder, exist_ok=True)
    os.makedirs(os.path.join(bexp_folder, cid_folder), exist_ok=True)
    os.chdir(os.path.join(bexp_folder, cid_folder))
    
    print(f"Created {cid_path}")
    print("writing molinfo.txt file")
    write_mol_info(cidnumber_list, chemname_list, molformula_list, molmonoweight_list,molinchi_list, molseq, molAorC)

    #print(f"current folder: {os.getcwd()}")

    os.chdir("../")
    
    os.makedirs("prep" , exist_ok=True) 
    os.makedirs("qcrun", exist_ok=True) 
    print(f"Created {prep_path}")
    print(f"Created {qcrun_path}")
    #print(f"current folder: {os.getcwd()}")
    return None 




def set_mol_prep(molseq, molcid, chemname):
    
    os.chdir("prep")
    
    molecule_name = molseq 
    molecule_cidnumber = molcid
    
    write_start_prep_qcxms_log(molecule_name, molecule_cidnumber)
    setup_calc_mol(molecule_name,molecule_cidnumber)
    get_finalprotomers(chemname)
    
    protomer_files = [f for f in os.listdir() if f.endswith(".xyz") and f"{chemname}-protomer_" in f]
    write_done_prep_qcxms_log(molecule_name, molecule_cidnumber)
    #print(f"current folder: {os.getcwd()}")
    
    return protomer_files




def write_start_prep_qcxms_log(molecule_name, molecule_cidnumber):

    log_file = '_prep_qcxms.log'
    with open(log_file, 'a') as log:
        log.write('\n')
        log.write('#' * 69 + '\n')
        log.write(f"{datetime.now()}\n")
        log.write(f"Molecule name : {molecule_name}\n")
        log.write(f"The cidnumber for {molecule_name} : {molecule_cidnumber}.\n")
   
   return None




def setup_calc_mol(molecule_name,molecule_cidnumber):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{molecule_cidnumber}/SDF?record_type=3d"
    response = requests.get(url)

    with open(f"{molecule_name}.mol", 'wb') as file:
        file.write(response.content)
    
    subprocess.run(['obabel', '-imol', f'{molecule_name}.mol', '-oxyz', '-O', f'{molecule_name}.xyz'], check=True)
    
    with open('energy_txb.out', 'w') as outfile:
        subprocess.run(['xtb', f'{molecule_name}.xyz', '--opt', 'extreme'], stdout=outfile, stderr=subprocess.STDOUT, check=True)
    
    with open('protomer_crest.out', 'w') as outfile:
        subprocess.run(['crest', 'xtbopt.xyz', '-protonate', '-T', '80S'], stdout=outfile, stderr=subprocess.STDOUT, check=True)
    
    return None




def read_xyz(filename):

    frames = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        num_atoms = int(lines[0].strip())
        num_lines_per_frame = num_atoms + 2  
        num_frames = len(lines) // num_lines_per_frame
        
        for i in range(num_frames):
            frame_lines = lines[i * num_lines_per_frame:(i + 1) * num_lines_per_frame]
            frames.append(frame_lines)
    
    return frames




def write_xyz(frames, molname, base_filename):
    for i, frame in enumerate(frames):
        with open(f"{molname}-{base_filename}_{i}.xyz", 'w') as file:
            file.writelines(frame)




def get_finalprotomers(molname):

    infilename = ["protonated.xyz"]
    for filename in infilename:
        if os.path.isfile(filename):
            if filename.startswith("protonated"):
                output_base_filename = "protomer"
            #elif filename.startswith("deprotonated"):
            #    output_base_filename = "deprotomer"
            else:
                raise ValueError("Input filename should start with 'protonated' or 'deprotonated'")

            frames = read_xyz(filename)
            write_xyz(frames, molname, output_base_filename)
            print(f"Frames written for {filename} with base filename {output_base_filename}")
        else:
            print(f"{filename} not found in the folder.")
    
    return None




def write_done_prep_qcxms_log(molecule_name, molecule_cidnumber):

    log_file = '_prep_qcxms.log'
    with open(log_file, 'a') as log:
        log.write("prep is done\n")
        log.write(f"{datetime.now()}\n")
        log.write('#' * 69 + '\n')
        log.write('\n')
    
    return None




def get_protomers_qcrun(protomer_files):

    protomer_dir_list = []
    
    for protomer_file in protomer_files:
        protomer_dir_name = protomer_file.split(".")[0].rsplit("-", 2)[-1]
        protomer_number = protomer_dir_name.split("_")[-1]
        #dir_name = protomer_file.split("_", 2)[2].rsplit(".", 1)[0]
        target_dir = os.path.join("../qcrun", protomer_dir_name)
        os.makedirs(target_dir, exist_ok=True)
        shutil.copy(protomer_file, target_dir)
        print(f"created {target_dir}")
        protomer_dir_list.append(protomer_dir_name)
    os.chdir("../qcrun")
    
    return protomer_dir_list




def set_protomer_struc(protomer_dir_list, molAorC):

    for folder in protomer_dir_list:
        if os.path.isdir(folder):
            os.chdir(folder)
            print(f"Processing folder: {folder}")

            if molAorC == 'A': 
                for param in A_parameter_options:
                    param_folder = param
                    os.makedirs(param_folder, exist_ok=True)
                    shutil.copy(os.path.join(jxms_runBIN, "jxms_run_qcxms.sh"), param_folder)
                    shutil.copy(os.path.join(jxms_evalBIN, "jxms_run_eval.sh"), param_folder)
                    shutil.copy(os.path.join(jxms_evalBIN, "process_qcxms_jdx.py"), param_folder)
                    source_qcxmsin = os.path.join(jxms_paramBIN, param_folder, "qcxms.in")
                    target_qcxmsin = os.path.join(param_folder, "qcxms.in")
                    if os.path.exists(source_qcxmsin):
                        shutil.copy(source_qcxmsin, target_qcxmsin)
                    else:
                        print(f"Warning: {source_qcxmsin} not found. Skipping copy.")
                print(f"Current folder after setup: {os.getcwd()}")
                protomer_folder = os.getcwd() 
                source_run_protomer = os.path.join(jxms_paramBIN, "run_protomer_A.py")
                target_run_protomer = os.path.join(protomer_folder, "run_protomer_A.py")
                shutil.copy(source_run_protomer, target_run_protomer)
                os.chdir("..") 

                protomer_num = len(protomer_dir_list)
                write_run_molecule_script(protomer_num, 'A')

            if molAorC == 'C': 
                for param in C_parameter_options:
                    param_folder = param
                    os.makedirs(param_folder, exist_ok=True)
                    shutil.copy(os.path.join(jxms_runBIN, "jxms_run_qcxms.sh"), param_folder)
                    shutil.copy(os.path.join(jxms_evalBIN, "jxms_run_eval.sh"), param_folder)
                    shutil.copy(os.path.join(jxms_evalBIN, "process_qcxms_jdx.py"), param_folder)
                    source_qcxmsin = os.path.join(jxms_paramBIN, param_folder, "qcxms.in")
                    target_qcxmsin = os.path.join(param_folder, "qcxms.in")
                    if os.path.exists(source_qcxmsin):
                        shutil.copy(source_qcxmsin, target_qcxmsin)
                    else:
                        print(f"Warning: {source_qcxmsin} not found. Skipping copy.")
                print(f"Current folder after setup: {os.getcwd()}")
                protomer_folder = os.getcwd() 
                source_run_protomer = os.path.join(jxms_paramBIN, "run_protomer_C.py")
                target_run_protomer = os.path.join(protomer_folder, "run_protomer_C.py")
                shutil.copy(source_run_protomer, target_run_protomer)
                os.chdir("..") 
                
                protomer_num = len(protomer_dir_list)
                write_run_molecule_script(protomer_num, 'C')
            

    #protomer_num = len(protomer_dir_list)
    #write_run_molecule_script(protomer_num)
    print(f"set_prot_struc finished. Current folder: {os.getcwd()}")




def write_run_molecule_script(protomer_num, AorC ):

    script_content = "#!/bin/bash\n\n"
    #run_protomer_AorC = glob.glob(f"run_protomer_*.py")
    #AorC = run_protomer_AorC[0].split('.')[0].split('_')[-1]
    print(AorC)
    # Loop to generate commands for each protomer
    for i in range(protomer_num):
        script_content += f"cd protomer_{i}\n"
        #if AorC == 'A':
        #    script_content += "python3 run_protomer_A.py\n"
        #    script_content += "cd ..\n"
        script_content += f"python3 run_protomer_{AorC}.py\n"
        #script_content += f"python3 {run_protomer_AorC}\n"
        #script_content += "python3 run_protomer.py\n"
        script_content += "cd ..\n"


    with open("run_molecule.sh", "w") as f:
        f.write(script_content)
    os.chmod("run_molecule.sh", 0o755)
    print(f"run_molecule.sh generated for {protomer_num} protomers.")





def main():

    cidnumber_list, chemname_list, molformula_list, molmonoweight_list, molinchi_list, molAorC_list= gather_molinfo(jxms_infohome)
    
    molseq_list, molcid_list = set_mol_cidseq(cidnumber_list)
    
    
    for molseq, molcid, chemname, molAorC in zip(molseq_list, molcid_list, chemname_list, molAorC_list):
        print(f"Prepping for : bexp_{molseq}, cid_{molcid}  {chemname} {molAorC}")
    
        set_molstruc(molseq, molcid, cidnumber_list, chemname_list, molformula_list, molmonoweight_list, molinchi_list, molAorC)
        protomer_files = set_mol_prep(molseq, molcid, chemname)
        protomer_dir_list = get_protomers_qcrun(protomer_files)
        protomer_dir_list.sort()
        set_protomer_struc(protomer_dir_list, molAorC)
    
        print(f"Current folder: {os.getcwd()}")
        print(f"Next mol!")
        os.chdir(jxms_calcshome)




if __name__ == "__main__":
    main()

