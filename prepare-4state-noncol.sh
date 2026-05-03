#!/bin/bash

# Prepares directory tree for noncollinear 4-state exchange tensor calculations.
#
# Usage: bash prepare-4state-noncol.sh
#
# Prompts for supercell size (m n), then creates:
#   {m}x{n}/Gd1_Gd2/Jzz_bathx/{uu,ud,du,dd}/
#   {m}x{n}/Gd1_Gd2/Jzz_bathy/{uu,ud,du,dd}/
#   ...etc for all pairs and components.

ROOTDIR="/home/aditya/GdI2_Jzz/Jzz_var"
SCRIPT="GdI2_4state-noncol.py"
SUPERCELL_SCRIPT="supercell.py"
STATES="uu dd ud du"

# ---- Python path ----
PYTHON="/home/aditya/mamba_env/bin/python"
if [ ! -f "$PYTHON" ]; then
    echo "ERROR: Python not found at $PYTHON"
    exit 1
fi
echo "Using Python: $PYTHON"

# ================================================================
# STEP 1: Supercell size
# ================================================================
echo ""
echo "Enter supercell size (m n), e.g. '2 3' for 2x3x1:"
read -r m n
SCDIR="${m}x${n}"
echo "Building ${m}x${n}x1 supercell..."
cd "$ROOTDIR"
"$PYTHON" "$ROOTDIR/${SUPERCELL_SCRIPT}" "$m" "$n"
if [ $? -ne 0 ]; then
    echo "ERROR: supercell.py failed"; exit 1
fi

# Create supercell-level directory
mkdir -pv "$ROOTDIR/$SCDIR"

# ================================================================
# STEP 2: Define pairs and components
# ================================================================

# --- HARDCODED pairs (0-indexed atom indices from supercell.py) ---
# Gd1=0 (origin), Gd2=1 (+a), Gd3=2 (+b), Gd4=3 (+a+b)
# Pairs: Gd1-Gd2, Gd1-Gd3, Gd1-Gd4
PAIR_LABELS=("Gd1_Gd2" "Gd1_Gd3" "Gd1_Gd4")
PAIR_IDX1=(0 0 0)
PAIR_IDX2=(1 2 3)

# --- To make interactive later, uncomment below and comment above ---
# echo "Enter pairs as '0,1 0,2 0,3':"
# read -r pair_input
# ... parse into PAIR_LABELS, PAIR_IDX1, PAIR_IDX2

# --- HARDCODED components: Jzz with bath x, y, z ---
COMPONENTS=("zz")
BATHS=("x")

# --- To make interactive later, uncomment below and comment above ---
# echo "Enter tensor components (e.g. xx yz zz):"
# read -r comp_input
# COMPONENTS=($comp_input)
# for comp in "${COMPONENTS[@]}"; do
#     echo "  Bath direction for J${comp}? (x / y / z):"
#     read -r bdir
#     BATHS+=("$bdir")
# done

# ================================================================
# STEP 3: Create directories and generate VASP inputs
# ================================================================
echo ""
n_pairs=${#PAIR_LABELS[@]}

for p in $(seq 0 $((n_pairs-1))); do
    label="${PAIR_LABELS[$p]}"
    idx1="${PAIR_IDX1[$p]}"
    idx2="${PAIR_IDX2[$p]}"

    echo "=== Pair: $label  (atoms $idx1, $idx2) ==="
    mkdir -pv "$ROOTDIR/$SCDIR/$label"

    for comp in "${COMPONENTS[@]}"; do
        for bdir in "${BATHS[@]}"; do
            topdir="${SCDIR}/${label}/J${comp}_bath${bdir}"
            echo "  Creating: $topdir"
            mkdir -pv "$ROOTDIR/$topdir"

            for state in $STATES; do
                statedir="$ROOTDIR/${topdir}/${state}"
                mkdir -pv "$statedir"
                cd "$statedir"
                echo "    --> J${comp}  state=${state}  bath=${bdir}  atoms=${idx1},${idx2}"
                "$PYTHON" "$ROOTDIR/${SCRIPT}" "$comp" "$state" "$bdir" "$idx1" "$idx2" "$m" "$n"
                cd "$ROOTDIR"
            done
        done
    done
    echo ""
done

# ================================================================
# Summary
# ================================================================
echo "All done. Directory tree:"
for p in $(seq 0 $((n_pairs-1))); do
    label="${PAIR_LABELS[$p]}"
    for comp in "${COMPONENTS[@]}"; do
        for bdir in "${BATHS[@]}"; do
            echo "  ${SCDIR}/${label}/J${comp}_bath${bdir}/ : uu  dd  ud  du"
        done
    done
done
