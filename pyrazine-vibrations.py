import matplotlib.pyplot as plt
from ase.calculators.emt import EMT
from ase.calculators.loggingcalc import LoggingCalculator
from ase.data.pubchem import pubchem_atoms_search
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.visualize import view

# Get molecule
mol = pubchem_atoms_search(name="pyrazine")
# mol = pubchem_atoms_search(cid=9261)

# View molecule
view(mol)

# Set up logger
log_calc = LoggingCalculator(EMT())

# Set up Effective Medium Potential calculator
log_calc.label = "pyrazine"
mol.calc = log_calc

# Optmiize structure so that forces are fmax
BFGS(mol).run(fmax=0.01)

# Calculate vibrations
vib = Vibrations(mol)
vib.run()
vib.summary()

# Plot some things
log_calc.plot(energy=False)
plt.show()
