#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

void generateSMILES(const string& molecule) {
    // Write molecule to input file
    ofstream inputFile("Organic_input.txt");
    if (!inputFile) {
        cerr << "Error creating input file!" << endl;
        return;
    }
    inputFile << molecule;
    inputFile.close();

    // Call Python script and rdkit environment
    int env= system("conda activate rdkit_env");
    int result = system("python Organic_Smile.py");
    if (result != 0 && env !=0) {
        cerr << "Error running Python script!" << endl;
        return;
    }

    // Read SMILES from output file
    ifstream outputFile("Organic_Smi_out.txt");
    if (!outputFile) {
        cerr << "Error reading output file!" << endl;
        return;
    }

    string smiles;
    getline(outputFile, smiles);
    outputFile.close();

    cout << "Generated SMILES: " << smiles << endl;
}

int main() {
    string molecule;
    cout << "Enter the molecule (e.g., C6H6 for benzene): ";
    cin >> molecule;

    generateSMILES(molecule);

    return 0;
}
