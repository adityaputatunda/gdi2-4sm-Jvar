#!/bin/bash

# Extracts TOTEN from all 4-state calculations and computes J values.
# Walks: {m}x{n}/Gd*_Gd*/J*_bath*/{uu,dd,ud,du}/

ulimit -s unlimited
export LC_NUMERIC=C

ROOTDIR="/home/aditya/GdI2_Jzz/Jzz_var"
STATES=(uu ud du dd)

cd "$ROOTDIR"

for scdir in [0-9]*x[0-9]*/; do
    scdir="${scdir%/}"
    echo "############## Supercell: $scdir ##############"

    for pairdir in ${scdir}/Gd*_Gd*/; do
        pairdir="${pairdir%/}"
        echo "======== $pairdir ========"

        for topdir in ${pairdir}/J*_bath*/; do
            topdir="${topdir%/}"
            enes=()

            for state in "${STATES[@]}"; do
                dir="${topdir}/${state}"
                if [ ! -d "$dir" ]; then
                    echo "WARNING: $dir not found"
                    enes+=("0")
                    continue
                fi
                cd "$ROOTDIR/$dir"
                ene=$(grep TOTEN OUTCAR | tail -1 | cut -c32-45)
                ene=$(echo "$ene" | bc)
                enes+=("$ene")
                echo -e "${topdir}\t${state}\t${ene}"
                cd "$ROOTDIR"
            done

            # J = (E_uu - E_ud - E_du + E_dd) / (4 * M^2)  [* conversion]
            J=$(bc -l <<< "(${enes[0]} - ${enes[1]} - ${enes[2]} + ${enes[3]}) * 250")
            echo -e "${topdir}\tJ = ${J} meV\n"
        done
    done
done
