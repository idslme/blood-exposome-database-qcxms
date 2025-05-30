#!/bin/bash
# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                        — by Jesi Lee —                         ║
# ║ Project     : Blood Exposome CID-MS/MS Prediction via QCxMS    ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : jxms_run_qcxms.sh                                ║
# ║ Description : automated qcxms MD simulation runs               ║ 
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝


date | tee -a _run_qcxms.log
echo "START run_qcxms.sh" | tee -a _run_qcxms.log

cp ../*.xyz .

full_path=$(pwd)
three_up=$(dirname "$(dirname "$(dirname "$full_path")")")
bexp_id=$(basename "$three_up")

xyzfile=$(ls *xyz)
mol_name=$(echo $xyzfile | cut -d "-" -f 1)
prot_id=$(echo $xyzfile | cut -d "-" -f 2 | cut -d "." -f 1)
param_id=$(basename $PWD)

echo "1_qcxms start" | tee -a _run_qcxms.log
date | tee -a 1_qcxms.log
qcxms -i ${xyzfile} | tee -a 1_qcxms.log
date | tee -a _run_qcxms.log
echo "1_qcxms done" | tee -a _run_qcxms.log

echo "2_qcxms start" | tee -a _run_qcxms.log
qcxms -i ${xyzfile} | tee -a 2_qcxms.log
echo "2_qcxms done" | tee -a _run_qcxms.log

echo "3_qcxms start" | tee -a _run_qcxms.log
date | tee -a 3_pqcxms.log
pqcxms -j 2 -t 10 | tee -a 3_pqcxms.log
date | tee -a 3_pqcxms.log
echo "3_qcxms done" | tee -a _run_qcxms.log
date | tee -a _run_qcxms.log

plotms
cp "result.jdx" "${mol_name}_${prot_id}-${param_id}.jdx" #Taurine-protomer_0_n2-maxcol4-25.jdx
cp "result.jdx" "${bexp_id}_${mol_name}_${prot_id}-${param_id}.jdx" #Taurine-protomer_0_n2-maxcol4-25.jdx
echo "qcxms calculation done" | tee -a 3_pqcxms.log

echo "DONE run_qcxms.sh" | tee -a _run_qcxms.log
date | tee -a _run_qcxms.log


