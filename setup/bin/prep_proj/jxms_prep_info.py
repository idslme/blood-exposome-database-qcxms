# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                     — Authored by Jesi Lee —                   ║
# ║ Project     : Blood Exposome CID-MS/MS Prediction via QCxMS    ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : jxms_prep_info.py                                ║
# ║ Description : automated info/ setup via pubchem                ║ 
# ║               requires an input file:"info/jxms_pubchem_cid.in"║
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝


import sys, os 
import requests
import time
from rdkit import Chem
sys.path.append(os.path.abspath(".."))
from jxms_config import (
    jxms_projhome, jxms_infohome, jxms_nisthome,
    jxms_calcshome, jxms_codeBIN, jxms_prepBIN,
    jxms_runBIN, jxms_evalBIN, jxms_paramBIN
)


def read_cid_file(file_name):

    cids = []
    try:
        with open(file_name, 'r') as file:
            for line in file:
                cid = line.strip()
                if cid.isdigit():
                    cids.append(int(cid))
                else:
                    print(f"!! Invalid CID skipped: {cid}")
    
    except FileNotFoundError:
        print(f"!! File not found: {file_name}")
    except Exception as e:
        print(f"!! Error reading CID file: {e}")
    
    return cids




def get_cheminfo(cid, retries=3, delay=5):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey,MolecularFormula,MolecularWeight,ExactMass,IUPACName,CanonicalSMILES/JSON"

    for attempt in range(retries):
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            props = data['PropertyTable']['Properties'][0]
            return {
                'inchikey': props.get('InChIKey', ''),
                'mf': props.get('MolecularFormula', ''),
                'mw': props.get('MolecularWeight', ''),
                'monomass': props.get('ExactMass', ''),
                'name': props.get('IUPACName', ''),
                'smile': props.get('CanonicalSMILES', '')
            }
        except requests.exceptions.HTTPError as http_err:
            if response.status_code == 503:
                time.sleep(delay)
            else:
                return None
        except Exception as err:
            return None
   
    print(f"!! Failed after {retries} attempts: {cid}.")
    return None




def write_info_files(cid_info_map, prefix="jxms_"):
    
    fields = ['inchikey', 'mf', 'mw', 'monomass', 'name', 'smile']
    for field in fields:
        
        with open(f"{prefix}{field}.txt", 'w') as f:
            for cid, info in cid_info_map.items():
                f.write(f"{cid}\t{info.get(field, '')}\n")
        print(f"**{prefix}{field}.txt created")




def write_AorC_out(cid_info_map, outfilename="jxms_AorC.txt"):

    with open(outfilename, 'w') as f:
        for cid, info in cid_info_map.items():
            smile = info.get('smile', '')
            mol = Chem.MolFromSmiles(smile)
            
            if mol is None:
                print(f"!! invalid SMILES: RDKit failed")
                f.write(f"{cid}\tUnknown\n")
                continue
            is_cyclic = mol.GetRingInfo().NumRings() > 0
            f.write(f"{cid}\t{'C' if is_cyclic else 'A'}\n")

    print(f"**{outfilename} created")





def main():

    cid_infilename = os.path.join(jxms_infohome,"jxms_pubchem_cid.in")
    #Please change it to your filename: 
    #cid_infilename = "jxms_pubchem_cid.in" 
    cid_list = read_cid_file(cid_infilename)
    if not cid_list:
        return

    cheminfo_map = {}
    for cid in cid_list:
        info = get_cheminfo(cid)
        if info:
            cheminfo_map[cid] = info

    write_info_files(cheminfo_map)
    write_AorC_out(cheminfo_map)




if __name__ == "__main__":
    main()

