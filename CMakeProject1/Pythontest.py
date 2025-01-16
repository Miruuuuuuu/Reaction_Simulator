from rdkit import Chem
from rdkit.Chem import Descriptors

# Create a molecule from a SMILES string
smiles = "CCO"  # Ethanol
molecule = Chem.MolFromSmiles(smiles)

if molecule:
    print("Molecule successfully created!")
    
    # Print basic information about the molecule
    print("SMILES:", Chem.MolToSmiles(molecule))
    print("Number of atoms:", molecule.GetNumAtoms())
    print("Number of bonds:", molecule.GetNumBonds())
    
    # Calculate a molecular descriptor (e.g., molecular weight)
    mol_weight = Descriptors.MolWt(molecule)
    print("Molecular Weight:", mol_weight)
else:
    print("Failed to create a molecule. Check your SMILES string.")
