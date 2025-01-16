from rdkit import Chem
from rdkit.Chem import Draw

def generate_molecule_from_input(element_data):
    """
    Generate a molecule from user-provided element symbols and counts.
    Supports straight-chain organic molecules with hydrogens and other elements.
    """
    mol = Chem.RWMol()

    # Add carbon atoms first (primary chain)
    carbon_indices = []
    for _ in range(element_data.get("C", 0)):
        idx = mol.AddAtom(Chem.Atom("C"))
        carbon_indices.append(idx)

    # Connect carbon atoms in a straight chain
    for i in range(len(carbon_indices) - 1):
        mol.AddBond(carbon_indices[i], carbon_indices[i + 1], Chem.BondType.SINGLE)

    # Add hydrogens to saturate carbon valencies
    for c_idx in carbon_indices:
        atom = mol.GetAtomWithIdx(c_idx)
        needed_hydrogens = 4 - atom.GetDegree()  # Carbon needs 4 bonds
        for _ in range(needed_hydrogens):
            h_idx = mol.AddAtom(Chem.Atom("H"))
            mol.AddBond(c_idx, h_idx, Chem.BondType.SINGLE)

    # Add other elements (e.g., O, N, etc.)
    for symbol, count in element_data.items():
        if symbol == "C" or symbol == "H":
            continue  # Carbon and hydrogen already handled

        for _ in range(count):
            element_idx = mol.AddAtom(Chem.Atom(symbol))
            # For simplicity, attach to the first carbon (customize as needed)
            mol.AddBond(carbon_indices[0], element_idx, Chem.BondType.SINGLE)

    # Sanitize and generate the molecule
    Chem.SanitizeMol(mol)
    smiles = Chem.MolToSmiles(mol, canonical=True)
    smiles = smiles.replace("(", "").replace(")", "").replace("H", "")  
    return mol, smiles


if __name__ == "__main__":
    # User provides symbols and counts
    element_data = {}
    n = int(input("Enter the number of different elements in the molecule: "))
    for _ in range(n):
        symbol = input("Enter the element symbol (e.g., C, H, O): ")
        count = int(input(f"Enter the number of {symbol} atoms: "))
        element_data[symbol] = count

    try:
        mol, smiles = generate_molecule_from_input(element_data)
        print(f"Generated SMILES: {smiles}")
        # Optionally display or save the molecule image
        Draw.MolToImage(mol).show()
    except Exception as e:
        print(f"Error: {str(e)}")
