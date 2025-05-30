#!/bin/bash
# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                         — by Jesi Lee —                        ║
# ║ Project     : Blood Exposome CID-MS/MS Prediction via QCxMS    ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : jxms_run_eval.sh                                 ║
# ║ Description : automated evaluation of QCxMS result             ║
# ║                given nist/                                     ║ 
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝


rm -rf _run_eval.log
rm -rf evaluate

decimal_included="2"
jxms_evalBIN="$HOME/JxMS/setup/bin/proj_eval"

full_path=$(pwd)
bexp_id=$(basename "$(dirname "$(dirname "$(dirname "$full_path")")")")
xyzfile=$(ls *xyz)

mol_name=$(echo "$xyzfile" | awk -F'-protomer' '{print $1}')
prot_id=$(echo "$xyzfile" | awk -F'-' '{print $NF}' | sed 's/\.xyz//')
param_id=$(basename $PWD)
my_theo=${bexp_id}_${mol_name}_${prot_id}-${param_id}.jdx 

date | tee -a _run_eval.log
echo "***START analysis***************${my_theo}" | tee -a _run_eval.log
echo | tee -a _run_eval.log
echo "Running  jxms_run_eval.sh " | tee -a _run_eval.log


mv "result.jdx" "result_og.jdx"
python3 filter_result_jdx.py

cp "result.jdx" "${bexp_id}_${mol_name}_${prot_id}-${param_id}.jdx"
echo "Running  convert_jdx2msp.py" | tee -a _run_eval.log
python3 ${jxms_evalBIN}/convert_jdx2msp.py
echo "Running  process_qcxms_jdx.py" | tee -a _run_eval.log
python3 process_qcxms_jdx.py


mkdir evaluate
cd evaluate
cp ../${my_theo} .
mol_nistfolder="../../../../../../nist/NIST23_${bexp_id}_*"
ln -fs ${mol_nistfolder}/peaks_tohave_*.dat .
mol_infofile_pattern="../../../../cid_*/mol_info.txt"
mol_infofile=$(echo ${mol_infofile_pattern} | head -n 1)


echo "Running  run_eval_qcxms2nist.py" | tee -a ../_run_eval.log
python3  ${jxms_evalBIN}/run_eval_qcxms2nist.py ${my_theo} ${decimal_included} ${mol_infofile}
echo ">>> Best_result :" | tee -a ../_run_eval.log
grep n2 result_best.out | tee -a ../_run_eval.log
#cat result_best.out | tee -a ../_run_eval.log
echo "*********DONE analysis*********" | tee -a ../_run_eval.log
echo | tee -a ../_run_eval.log


