import matplotlib.pyplot as plt
import numpy as np
from pyscf import dft, gto, scf

# NOTE: 
# gto : Gaussian type orbitals
# scf : self-consistent field
# dft : density functional theory

hartree_to_ev = 27.211


def get_molecule(d, basis):
    """Get hydrogen molecule d Angstroms apart"""
    mol = gto.Mole()
    mol.build(
        atom=f"""H 0 0 0; H  {d} 0 0""",
        basis=basis,
        unit="Angstrom",
    )
    mol.build()
    return mol


# The basis set
basis = "sto-3g"

distances = np.linspace(0.5, 3, 50)

# Get energies using restricted Hartree-Fock
energies_hf = np.empty(shape=distances.shape)
for i, d in enumerate(distances):
    mol = get_molecule(d=d, basis=basis)
    mf_hf = scf.RKS(mol)
    mf_hf = mf_hf.newton()
    mf_hf.kernel()
    energies_hf[i] = mf_hf.energy_tot() * hartree_to_ev

# Get energies using restricted Kohn-Sham DFT
energies_dft = np.empty(shape=distances.shape)
for i, d in enumerate(distances):
    mol = get_molecule(d=d, basis=basis)
    mf_dft = dft.RKS(mol)
    mf_dft.xc = "lda"
    # mf_dft.xc = "lda,vwn" # With dispersion corrections
    # mf_dft.xc = 'b3lyp'
    mf_dft = mf_dft.newton()
    mf_dft.kernel()
    energies_dft[i] = mf_dft.energy_tot() * hartree_to_ev

# Minimum energy indices
idx_hf = np.argmin(energies_hf)
idx_dft = np.argmin(energies_dft)

print(f"Equilibrium distance Hartree-Fock = {distances[idx_hf]:.2f} Angstroms")
print(f"Equilibrium energy Hartree-Fock = {energies_hf[idx_hf]:.2f} eV")
print(f"Molecular orbital energies Hartree-Fock = {mf_hf.mo_energy}")
print()
print(f"Equilibrium distance DFT = {distances[idx_dft]:.2f} Angstroms")
print(f"Equilibrium energy DFT = {energies_dft[idx_dft]:.2f} eV")
print(f"Molecular orbital energies DFT = {mf_dft.mo_energy}")

plt.plot(distances, energies_hf, label="Hartree-Fock", color="tab:blue")
plt.plot(distances, energies_dft, label="LDA", color="tab:orange")
plt.plot(distances[idx_hf], energies_hf[idx_hf], "o", color="tab:blue")
plt.plot(distances[idx_dft], energies_dft[idx_dft], "o", color="tab:orange")
plt.xlabel(r"$d$ $\AA$")
plt.ylabel(r"eV")
plt.legend()
plt.show()
