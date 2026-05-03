import sys
import numpy as np
from ase.io import read, write

# --- Usage: python supercell.py m n
m, n = int(sys.argv[1]), int(sys.argv[2])

ROOTDIR = '/home/aditya/GdI2_Jzz/Jzz_var'

unit = read(f'{ROOTDIR}/POSCAR.3atom.transZ.cart.vasp')
a_vec = unit.cell[0]
b_vec = unit.cell[1]

supercell = unit * (m, n, 1)

# --- Identify Gd1(0,0) and its 3 NN: +a, +b, +a+b ---
symbols = np.array(supercell.get_chemical_symbols())
gd_mask = symbols == 'Gd'
gd_indices = np.where(gd_mask)[0]
i_indices  = np.where(~gd_mask)[0]

gd_pos = supercell.positions[gd_indices]

# Gd1 = the Gd closest to the original unit cell Gd position
gd1_ref = unit.positions[0]
dists = np.linalg.norm(gd_pos - gd1_ref, axis=1)
gd1_local = np.argmin(dists)

# 3 nearest neighbors of Gd1 by lattice vector matching
gd1_pos = gd_pos[gd1_local]
targets = {
    '+a':   gd1_pos + a_vec,
    '+b':   gd1_pos + b_vec,
    '+a+b': gd1_pos + a_vec + b_vec,
}

def find_nearest_gd(target, positions, cell):
    """Find Gd index closest to target under PBC."""
    diff = positions - target
    frac = diff @ np.linalg.inv(cell)
    frac -= np.round(frac)
    cart = frac @ cell
    dists = np.linalg.norm(cart, axis=1)
    return np.argmin(dists)

nn = {lbl: find_nearest_gd(t, gd_pos, supercell.cell[:]) for lbl, t in targets.items()}

# --- Build sorted index: [Gd1, 3 NN, remaining Gd, all I] ---
top4 = [gd1_local, nn['+a'], nn['+b'], nn['+a+b']]
remaining_gd = [i for i in range(len(gd_indices)) if i not in top4]

sorted_global = (
    [gd_indices[i] for i in top4] +
    [gd_indices[i] for i in remaining_gd] +
    list(i_indices)
)

sorted_supercell = supercell[sorted_global]

write(f'{ROOTDIR}/POSCAR', sorted_supercell,
      format='vasp', vasp5=True, sort=False, direct=True)

print(f'Written POSCAR ({m}x{n}x1)')
print(f'Gd1 + 3 NN are indices 0, 1, 2, 3 in the Gd block.')
print(f'  Gd1 (0,0):  {sorted_supercell.positions[0]}')
print(f'  Gd2 (+a):   {sorted_supercell.positions[1]}')
print(f'  Gd3 (+b):   {sorted_supercell.positions[2]}')
print(f'  Gd4 (+a+b): {sorted_supercell.positions[3]}')
