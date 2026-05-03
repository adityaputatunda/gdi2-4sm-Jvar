# 4-State Energy Mapping Workflow for 2D GdI₂

ASE-based workflow to study variation of nearest-neighbor Heisenberg magnetic exchange tensor components of 2D GdI₂ across:

- Different supercell sizes (SC size convergence check)
- Gd–Gd bonds
- Direction of moments in the surrounding bath

…for the 4 states using the **4-state energy mapping method**.

The workflow initializes VASP calculations via ASE and executes them with a bash script.

## Variations

| # | Variation | Type | Example |
|---|-----------|------|---------|
| 1 | Tensor component | Input `<string>` | `xy`, `zz` |
| 2 | One of 4 states | Optional input `<string>` | `ud`, `dd` |
| 3 | Direction of bath moments | Input `<string>` | `x`, `y` |
| 4 | 2 selected atoms (J1/J2/J3…) | Optional input `<2 ints>` | `3 1`, `5 2` |
| 5 | 2D supercell to create | Input `<2 ints>` | `4 5`, `3 4` |

## Files

- **`prepare-4state-noncol.sh`** — main script; run with the arguments above. Minimum requirement: two integers for SC size.
- **`POSCAR.3atom.transZ.cart.vasp`** — unit cell structure (POSCAR format), used to build supercells.
- **`supercell.py`** — supercell builder.
- **`run_leo.sh`** — submits jobs for execution.
- **`GdI2_4state-noncol.py`** — main ASE workflow that initializes the calculation.
- **`toten.sh`** — gathers energies and computes J values.

## Requirements

- ASE (Python)
- Environment variable: `VASP_PP_PATH`
