from rdkit import Chem
from rdkit.Chem import AllChem

# Define the double replacement reaction template
reaction_smarts = "[A:1][X:2].[B:3][Y:4]>>[A:1][Y:4].[B:3][X:2]"
reaction = AllChem.ReactionFromSmarts(reaction_smarts)

def predict_reaction(reactant1, reactant2):
    try:
        # Convert reactants to RDKit molecule objects
        mol1 = Chem.MolFromSmiles(reactant1)
        mol2 = Chem.MolFromSmiles(reactant2)

        if not mol1 or not mol2:
            raise ValueError("Invalid chemical structure for one or both reactants.")

        # Apply the reaction
        products = reaction.RunReactants((mol1, mol2))

        # Output the predicted products
        print("Predicted products:")
        for i, product_set in enumerate(products):
            product_smiles = [Chem.MolToSmiles(product) for product in product_set]
            print(f"Product set {i+1}: {' + '.join(product_smiles)}")

    except Exception as e:
        print(f"Error: {e}")

# Example usage
reactant1 = "NaCl"  # Sodium chloride
reactant2 = "AgNO3"  # Silver nitrate
predict_reaction(reactant1, reactant2)
