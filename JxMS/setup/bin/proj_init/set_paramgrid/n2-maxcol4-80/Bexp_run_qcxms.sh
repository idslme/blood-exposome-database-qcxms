
cp ../*.xyz .
xyzfile=$(ls *xyz)
mol_id=$(echo "$xyzfile" | cut -d "." -f 1)
current_directory=$(basename "$PWD")

date | tee -a 1_qcxms.log
qcxms -i ${xyzfile} | tee -a 1_qcxms.log
qcxms -i ${xyzfile} | tee -a 2_qcxms.log
date | tee -a 3_pqcxms.log
pqcxms | tee -a 3_pqcxms.log
date | tee -a 3_pqcxms.log
plotms

cp "result.jdx" "${mol_id}_${current_directory}.jdx"
echo "qcxms calculation done" | tee -a 3_pqcxms.log
date
