#!/usr/bin/env python3
# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                     — Authored by Jesi Lee —                   ║
# ║ Project     : Blood Exposome MS/MS Prediction via QCxMS        ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : run_protomer_A.py                                ║
# ║ Description : Execute QCxMS production runs for A-protomers    ║
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝
#  Keep calm and simulate on 

import os, sys
import subprocess
import time 
import gc
from datetime import datetime

sys.path.append(os.path.abspath(".."))
from jxms_config import (
    jxms_projhome, jxms_infohome, jxms_nisthome,
    jxms_calcshome, jxms_codeBIN, jxms_prepBIN,
    jxms_runBIN, jxms_evalBIN, jxms_paramBIN
)

# Define parameter set (A protomer only)
A_parameter_options = [
    "n2-maxcol4-40", "n2-maxcol4-25", "n2-maxcol4-60",
    "n2-maxcol6-40", "n2-maxcol6-25", "n2-maxcol6-60",
    "n2-setcol4-40", "n2-setcol4-25", "n2-setcol4-60",
    "n2-setcol6-40", "n2-setcol6-25", "n2-setcol6-60",
    "n2-40", "n2-esi-maxcol4-40", "n2-esi-setcol4-40"
]

run_script = "./jxms_run_qcxms.sh"
eval_script = "./jxms_run_eval.sh"

full_path = os.getcwd()
bexp_id = os.path.basename(os.path.dirname(os.path.dirname(full_path)))
protomer_logfile = "../../_run_protomer.log"
molecule_logfile = "../../../_run_molecule.log"




def log(msg):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(molecule_logfile, 'a') as logf:
        logf.write(f"{timestamp} - {msg}\n")




def check_covered():
    return os.path.isfile("evaluate/COVERED.out")




for folder in A_parameter_options:
    if not os.path.isdir(folder):
        continue

    os.chdir(folder)
    protomer_id = os.path.basename(os.path.dirname(os.getcwd()))

    log(f"{'*' * 58}")
    log(f"Start  {bexp_id} {protomer_id} {folder}")
    
    try:
        log(f"Running QCxMS production: {run_script}")
        subprocess.run([run_script], check=True)

        gc.collect()
        time.sleep(5)

        log("Running evaluation script")
        subprocess.run([eval_script], check=True)

        if check_covered():
            log("*** COVERED ***")
        else:
            log("*** NOT COVERED ***")

    except subprocess.CalledProcessError as e:
        log(f"Subprocess error in {folder}: {e}")
    except Exception as e:
        log(f"Unexpected error in {folder}: {e}")
    finally:
        os.chdir("..")
        log(f"Done  {bexp_id} {protomer_id} {folder}")
        log(f"{'*' * 58}")

