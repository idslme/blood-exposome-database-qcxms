# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                     — Authored by Jesi Lee —                   ║
# ║ Project     : Blood Exposome CID-MS/MS Prediction via QCxMS    ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : jxms_config.py                                   ║
# ║ Description : project folder structure setup                   ║ 
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝


import os

JXMS_HOME = os.path.expanduser("~/JxMS")
jxms_projhome = os.path.join(JXMS_HOME, "jxms_proj")

jxms_infohome = os.path.join(jxms_projhome, "info")
jxms_nisthome = os.path.join(jxms_projhome, "nist")
jxms_calcshome = os.path.join(jxms_projhome, "calcs")

jxms_codeBIN = os.path.join(jxms_projhome, "setup", "bin")

jxms_initBIN = os.path.join(jxms_codeBIN, "proj_init")
jxms_paramBIN = os.path.join(jxms_initBIN, "set_paramgrid")
jxms_evalBIN = os.path.join(jxms_codeBIN, "proj_eval")



