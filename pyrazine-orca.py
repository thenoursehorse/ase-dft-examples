from pathlib import Path

from ase.calculators.orca import ORCA, OrcaProfile
from ase.data.pubchem import pubchem_atoms_search
from ase.io import read

# Initial geometry
atoms = pubchem_atoms_search(name="pyrazine")
# atoms = pubchem_atoms_search(cid=9261)

# Directories
directory_ground = Path("ground")
directory_tddft = Path("tddft")
directory_cis = Path("cis")

# Chemistry software setup
profile = OrcaProfile(command='/Users/henry/Library/orca_6_0_1/orca')
#profile = OrcaProfile(command="/path/to/orca/")
# functional = 'B3LYP D4 def2-SVP'
functional = "B3LYP def2-SVP"

# Optimize groundstate geometry
calc = ORCA(
    profile=profile,
    directory=directory_ground,
    orcasimpleinput=functional + "OPT",
)
atoms.calc = calc
atoms.get_potential_energy()


# Optimize at the conical intersection of the second (i=2) and first (j=1) excited states
# Overlaps go like <J|Operator|I>
atoms = read(directory_ground / Path("orca.xyz"))
calc = ORCA(
    profile=profile,
    directory=directory_cis,
    orcasimpleinput=functional + "CI-OPT TIGHTSCF",
    orcablocks=("%TDDFT IROOT 2\n       JROOT 1\nEND"),
)
atoms.calc = calc
e0 = atoms.get_potential_energy()
