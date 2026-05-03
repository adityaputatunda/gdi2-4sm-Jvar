#!/bin/bash

# Runs all noncollinear 4-state VASP calculations.
# Execute with: nohup bash run_leo.sh &
#
# Walks: {m}x{n}/Gd*_Gd*/J*_bath*/{uu,dd,ud,du}/

ROOTDIR="/home/aditya/GdI2_Jzz/Jzz_var"
STATES="uu dd ud du"

cd "$ROOTDIR"

for scdir in [0-9]*x[0-9]*/; do
    scdir="${scdir%/}"
    echo "======== Supercell: $scdir ========"
    for pairdir in ${scdir}/Gd*_Gd*/; do
        pairdir="${pairdir%/}"
        for topdir in ${pairdir}/J*_bath*/; do
            topdir="${topdir%/}"
            for state in $STATES; do
                dir="${topdir}/${state}"
                if [ ! -d "$dir" ]; then
                    echo "WARNING: $dir not found, skipping."
                    continue
                fi
                echo "Running: $dir"
                cd "$ROOTDIR/$dir"
                cp "$ROOTDIR/job.sh" .
                sbatch job.sh
                cd "$ROOTDIR"
            done
        done
    done
done

echo "All calculations submitted."
