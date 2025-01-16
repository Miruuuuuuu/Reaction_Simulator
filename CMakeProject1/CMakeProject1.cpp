// CMakeProject1.cpp : Defines the entry point for the application.
//

#include "CMakeProject1.h"

using namespace std;
/*  int main()
{
	string name;
	cout << "Hello CMake." << endl;
	cin >> name;
	cout << "String display :" << name.length();
	return 0;
}
 */
int main() {
    try {
        // Parse a SMILES string to create a molecule
        string smiles = "CCO";  // Ethanol
        RDKit::ROMol* mol = RDKit::SmilesToMol(smiles);

        if (!mol) {
            cerr << "Error: Could not parse SMILES string!" << endl;
            return 1;
        }

        // Output basic molecule information
        cout << "Molecule: " << smiles << endl;
        cout << "Number of atoms: " << mol->getNumAtoms() << endl;
        cout << "Number of bonds: " << mol->getNumBonds() << endl;

        // Compute molecular weight
        double molWeight = RDKit::Descriptors::calcAMW(*mol);
        cout << "Molecular weight: " << molWeight << " g/mol" << endl;

        // Clean up
        delete mol;
    }
    catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }

    return 0;
}