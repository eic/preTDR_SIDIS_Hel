#!/bin/bash
gr=$1

export MAINDIR=/direct/eic+u/cvhulse/epic/analysis/condor_scripts/18x275minQ1_${gr}
export OUTDIR=/gpfs/mnt/gpfs02/eic/cvhulse/epic/analysis/18x275minQ1_${gr}

/direct/eic+u/cvhulse/epic/eic-shell <<EOF
export LD_LIBRARY_PATH=/eic/u/cvhulse/epic/lhapdf/lib:/eic/u/cvhulse/epic/lhapdf/share/LHAPDF:/eic/u/cvhulse/epic/lhapdf/include:/eic/u/cvhulse/epic/local/lib:/opt/local/lib64:/opt/local/lib:$LD_LIBRARY_PATH
cd $MAINDIR
./event_eval.exe ${MAINDIR}/runlist.txt ${OUTDIR} 321
exit
EOF


