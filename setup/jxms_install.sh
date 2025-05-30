#!/bin/bash
# ╔════════════════════════════════════════════════════════════════╗
# ║                   JxMS Computational Workflow                  ║
# ║                     — Authored by Jesi Lee —                   ║
# ║ Project     : Blood Exposome CID-MS/MS Prediction via QCxMS    ║
# ╟────────────────────────────────────────────────────────────────╢
# ║ Script      : jxms_install.sh                                  ║
# ║ Description : automated installation of all required tools:    ║ 
# ║                QCxMS, PlotMS, xTB, and CREST.                  ║ 
# ║               requires sudo access.                            ║
# ║ Date        : 2025-05-28 (last updated)                        ║
# ╚════════════════════════════════════════════════════════════════╝


set -euo pipefail

QCXMS_VER="v.5.2.1"
PLOTMS_VER="v.6.2.0"

XTB_VER="6.7.0"
CREST_VER="v3.0.2"
CREST_BIN="crest-intel-2023.1.0-ubuntu-latest.tar.xz"

mkdir -p qcxms plotms xtb crest
echo ">>> Downloading..."

curl -L -o qcxms/QCxMS.${QCXMS_VER}.tar.xz https://github.com/qcxms/QCxMS/releases/download/${QCXMS_VER}/QCxMS.${QCXMS_VER}.tar.xz
curl -L -o plotms/PlotMS.${PLOTMS_VER}.tar.xz https://github.com/qcxms/PlotMS/releases/download/${PLOTMS_VER}/PlotMS.${PLOTMS_VER}.tar.xz

curl -L -o xtb/xtb-${XTB_VER}-linux-x86_64.tar.xz https://github.com/grimme-lab/xtb/releases/download/v${XTB_VER}/xtb-${XTB_VER}-linux-x86_64.tar.xz
curl -L -o crest/${CREST_BIN} https://github.com/crest-lab/crest/releases/download/${CREST_VER}/${CREST_BIN}


echo ">>> Extracting..."
tar -xf qcxms/QCxMS.${QCXMS_VER}.tar.xz -C qcxms
tar -xf plotms/PlotMS.${PLOTMS_VER}.tar.xz -C plotms
tar -xf xtb/xtb-${XTB_VER}-linux-x86_64.tar.xz -C xtb
tar -xf crest/${CREST_BIN} -C crest

echo ">>> Installing binaries to /usr/local/bin..."
sudo cp -f qcxms/qcxms /usr/local/bin/
sudo cp -f qcxms/pqcxms /usr/local/bin/
sudo cp -f plotms/PlotMS.${PLOTMS_VER}/plotms /usr/local/bin/
sudo cp -f xtb/xtb-dist/bin/xtb /usr/local/bin/
sudo cp -f crest/crest /usr/local/bin/

cp -f crest/crest xtb/xtb-dist/bin

echo ">>> Installed versions:"
qcxms --version || echo "QCxMS not found"
plotms --version || echo "PlotMS not found"
xtb --version || echo "xTB not found"
crest --version || echo "CREST not found"
echo ">>> QCxMS environment setup complete!"


