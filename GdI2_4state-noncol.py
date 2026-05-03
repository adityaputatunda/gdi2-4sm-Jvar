from ase.io import read
from ase.calculators.vasp import Vasp

import numpy as np
import os, sys
from collections import OrderedDict

# ---------- USAGE ---------------
# python3 GdI2_4state-noncol.py <component> <state> <bath> <magatom1> <magatom2> <m> <n>
# e.g.:  python3 GdI2_4state-noncol.py zz ud x 0 2 2 3

component  = sys.argv[1]   # two-letter string: xx, xy, xz, yx, yy, yz, zx, zy, zz
state      = sys.argv[2]   # uu, ud, du, dd
bath       = sys.argv[3]   # x, y, or z
magatom1   = int(sys.argv[4])
magatom2   = int(sys.argv[5])
sc_m       = int(sys.argv[6])
sc_n       = int(sys.argv[7])

# ---------- PARAMETERS ---------------

ROOTDIR  = '/home/aditya/GdI2_Jzz/Jzz_var'

system   = 'GdI2'
relax    = False
ldau     = True
U, J     = 9.2, 1.2
spinorb  = True

magmom_constr = 1
LAMBDA        = 10

M = 8.0   # magnetic moment magnitude (mu_B)

import math

encut    = 500
# Unit-cell base k-mesh: 8, 8, 1 — scaled by supercell size
K1_base, K2_base, K3_base = 8, 8, 1
k1 = max(1, math.ceil(K1_base / sc_m))
k2 = max(1, math.ceil(K2_base / sc_n))
k3 = K3_base

# ---------- SPIN DIRECTIONS ---------------

dirs = {
    'x': np.array([1., 0., 0.]),
    'y': np.array([0., 1., 0.]),
    'z': np.array([0., 0., 1.]),
}

if component[0] not in dirs or component[1] not in dirs:
    print(f'ERROR: Invalid component "{component}". Each letter must be x, y, or z.')
    sys.exit(1)
if bath not in dirs:
    print(f'ERROR: Invalid bath direction "{bath}". Must be x, y, or z.')
    sys.exit(1)

dir1     = dirs[component[0]]
dir2     = dirs[component[1]]
bath_dir = dirs[bath]

sign = {'u': +1., 'd': -1.}
sign1 = sign[state[0]]
sign2 = sign[state[1]]

print(f'\nComponent: J{component}  |  State: {state}  |  Bath: {bath}')
print(f'  Gd1 (atom {magatom1}): {sign1 * M * dir1}')
print(f'  Gd2 (atom {magatom2}): {sign2 * M * dir2}')
print(f'  Bath Gd atoms:         {M * bath_dir}\n')

# ---------- BUILD STRUCTURE + MAGMOMS ---------------

# POSCAR is 3 levels up: Gd1_Gd3/Jzz_bathx/uu/ -> ../../../POSCAR
atoms = read(f'{ROOTDIR}/POSCAR', format='vasp')

magmoms = np.zeros((len(atoms), 3))
for i, atom in enumerate(atoms):
    if atom.symbol == 'Gd':
        magmoms[i] = M * bath_dir

magmoms[magatom1] = sign1 * M * dir1
magmoms[magatom2] = sign2 * M * dir2

atoms.set_initial_magnetic_moments(magmoms)

# ---------- CALCULATOR ---------------

calc = Vasp(
    system = system,

    # ---- PARALLELIZATION ----
    npar = 4,
    kpar = 1,

    # ---- GENERAL ----
    xc     = 'pbe',
    setups = {'Gd': '', 'I': ''},
    kpts   = (k1, k2, k3),
    gamma  = True,

    # ---- SCF ----
    encut   = encut,
    ediff   = 1e-6,
    istart  = 0,
    icharg  = 2,
    ismear  = 0, sigma = 0.05,
    ispin   = 2,
    prec    = 'Accurate',
    lorbit  = 10,
    lreal   = 'False',
    isym    = -1,
    nelm    = 500,
    lmaxmix = 6,
    algo    = 'All',

    # ---- OUTPUT ----
    lwave  = False,
    lcharg = False,

    # ---- NONCOLLINEAR ----
    lnoncollinear = True,
    gga_compat    = False,
    lasph         = True,
)

if spinorb:
    calc.set(
        lsorbit = True,
        lorbmom = True,
        saxis   = [0, 0, 1],
    )

if relax:
    calc.set(
        isif=2, ibrion=1, potim=0.8,
        ediffg=-0.005, nsw=200, nelm=100,
    )

if ldau:
    calc.set(
        ldautype  = 2,
        ldau_luj  = {'Gd': {'L': 3, 'U': U, 'J': J}},
        ldauprint = 4,
    )

calc.initialize(atoms)
calc.write_input(atoms)

# ---------- POST-PROCESS INCAR ---------------

with open('INCAR', 'r') as f:
    incar = f.readlines()
os.remove('INCAR')

species_order = list(OrderedDict.fromkeys(atoms.get_chemical_symbols()))

magmoms_reordered = np.zeros_like(magmoms)
idx = 0
for sp in species_order:
    for atom in atoms:
        if atom.symbol == sp:
            magmoms_reordered[idx] = magmoms[atom.index]
            idx += 1

def fmt_val(x):
    return '{:.1f}'.format(0.0 if x == 0.0 else x)

def fmt_triplet(v):
    return '{} {} {}'.format(fmt_val(v[0]), fmt_val(v[1]), fmt_val(v[2]))

n_Gd = sum(1 for a in atoms if a.symbol == 'Gd')
n_I  = sum(1 for a in atoms if a.symbol == 'I')

gd_block = '  '.join(fmt_triplet(magmoms_reordered[i]) for i in range(n_Gd))
i_block  = '{}*0.0'.format(n_I * 3)

magmom_values = gd_block + '  ' + i_block

magmom_str   = ' MAGMOM = '   + magmom_values + '\n'
m_constr_str = ' M_CONSTR = ' + magmom_values + '\n'

incar = [l for l in incar if 'MAGMOM' not in l]
incar.append('\n# --- Noncollinear block (added by GdI2_4state-noncol.py) ---\n')
incar.append(magmom_str)

if magmom_constr == 1:
    incar.append(' I_CONSTRAINED_M = 1\n')
    incar.append(' LAMBDA = {}\n'.format(LAMBDA))
    incar.append(m_constr_str)
    print('Magmom directions constrained (I_CONSTRAINED_M = 1)')
elif magmom_constr == 0:
    print('Magmoms free (no constraint)')

rwigs = []
with open('POTCAR', 'r') as f:
    for line in f:
        if 'RWIGS' in line:
            rwigs.append(float(line.split()[5]))
print('RWIGS:', rwigs)
incar.append(' RWIGS = ' + ' '.join(str(r) for r in rwigs) + '\n')

with open('INCAR', 'w') as f:
    f.writelines(incar)

print('Done. Written: INCAR, POSCAR, KPOINTS, POTCAR')
