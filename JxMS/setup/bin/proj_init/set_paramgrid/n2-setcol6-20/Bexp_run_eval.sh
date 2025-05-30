#########################################################################
###  Last updated: Aug 14, 2024     Jesi Lee                          ###
###                                                                   ###
### under bexp_000/qcrun/protomer_0/n2-20/                            ###
###                                                                   ###
###                                                                   ###
#########################################################################

molecule_name="$1"
molecule_weight="$2"
molecule_formula="$3"

decimal_included=2
mol_ion=$(echo "$molecule_weight + 1.007825" | bc)
bexp_codeBIN="/home/jesi/projects/qcxms/BloodExposome_rand100/qcxms/bexp_${bexp_id}/cidnum_*"


-rf evaluate


xyzfile=$(ls *xyz)
mol_id=$(echo "$xyzfile" | cut -d "." -f 1)
bexp_id=$(echo $xyzfile | cut -d '_' -f 2)
prot=$(echo $mol | cut -d '_' -f 2)
prot_id=$(echo $prot | cut -d '.' -f 1)
param_id=$(basename $PWD)


mkdir evaluate
cd evaluate
my_theo=${mol_id}*.jdx
cp ../${my_theo} .
mol_nistfolder="../../../nist_bexp_*"
ln -fs ${mol_nistfolder}/*.JDX .
ln -fs ${mol_nistfolder}/peaks_tohave_*.dat .
#ln -fs $
nistfile=$(ls -1 *_CE*.JDX)

for i in $nistfile; do python3 ${bexp_codeBIN}/compare_qcxms_nist.py ${my_theo} ${i} ${decimal_included} ${mol_ion} | tee -a run_eval.log ; done

cat simscores.out



#for i in $nistfile; do python3 ${bexp_codeBIN}/Bexp_compare_qcxms_nist.py ${my_theo} ${i} | tee -a run_eval.log ; done
#echo         cid_fname,             nist_fname,      cos_score , dot_score , ent_score| tee -a simscores.out
#echo '* Filename: result_simscores.out is created.'






#CIDBIN="/home/jesilee/projects/CIDMD/pilot_4/scripts_bin"
#${CIDBIN}/make_mp4Ngif.sh



