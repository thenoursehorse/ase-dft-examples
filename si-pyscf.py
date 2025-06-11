import matplotlib.pyplot as plt
import numpy as np
import pyscf.pbc.dft as pbcdft
import pyscf.pbc.gto as pbcgto
import pyscf.pbc.tools.pyscf_ase as pyscf_ase
from ase.build import bulk
from ase.dft.kpoints import get_bandpath
from ase.dft.kpoints import sc_special_points as special_points

# NOTE:
# gto : Gaussian type orbitals
# dft : density functional theory
# pbc : periodic boundary conditions

hartree_to_ev = 27.211

# Get Silicon structure
si = bulk("Si", "diamond", 5.43)

# Convert to pyscf object
cell = pbcgto.Cell()
cell.build(
    atom=pyscf_ase.ase_atoms_to_pyscf(si),
    a=si.cell,
    basis="gth-szv",
    pseudo="gth-pade",  # psuedopotential (effective core potentials)
)

# Do Gamma pt calculation to get ground state
mf_dft = pbcdft.RKS(cell)
print(mf_dft.kernel())

# Bandpath for calculating bands
points = special_points["fcc"]
L = points["L"]
G = points["G"]
X = points["X"]
K = points["K"]
path = get_bandpath([L, G, X, K, G], si.cell, npoints=50)
path_kpts = cell.get_abs_kpts(path.kpts)

# Show Brillouin zone and the kpath
path.plot()
plt.show()

# Without updating charge density, get the KS orbital energies along the band path
e_kn = mf_dft.get_bands(path_kpts)[0]
e_kn = np.array(e_kn)

nbands = cell.nao_nr()
for n in range(nbands):
    plt.plot(
        path.get_linear_kpoint_axis()[0], e_kn[:, n] * hartree_to_ev, color="tab:blue"
    )
plt.xticks(path.get_linear_kpoint_axis()[1], path.get_linear_kpoint_axis()[2])
plt.ylabel("eV")
plt.show()
