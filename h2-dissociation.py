
from ase import Atoms
from ase.calculators.emt import EMT
import numpy as np
import matplotlib.pyplot as plt

def get_molecule(d):
    '''Get hydrogen molecule d Angstroms apart'''
    return Atoms('H2', [[0, 0, 0], [d, 0, 0]])

distances = np.linspace(0.5, 3, 50)
energies = np.empty(shape=distances.shape)

for i, d in enumerate(distances):
    mol = get_molecule(d=d)
    mol.calc = EMT()
    energies[i] = mol.get_potential_energy()

# Minimum energy indices
idx = np.argmin(energies)

print(f"Equilibrium distance = {distances[idx]:.2f} Angstroms")
print(f"Equilibrium energy = {energies[idx]:.2f} eV")

plt.plot(distances, energies)
plt.plot(distances[idx], energies[idx], 'o')
plt.xlabel(r"$d$ $\AA$")
plt.ylabel(r"eV")
plt.show()