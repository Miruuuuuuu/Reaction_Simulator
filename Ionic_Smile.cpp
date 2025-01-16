#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <stdexcept>

using namespace std;

// Struct to hold element data
struct Element {
    int atomicNumber;
    string symbol;
    double atomicMass;
    double electronegativity;
    int valency;  // Valency (oxidation state)
    string nature;  // Positive, Negative, Neutral
    double ionizationEnergy;  // In eV
    double shieldingEnergy;  // In eV
};

// Function to safely convert string to double
double safeStod(const string& str) {
    try {
        return stod(str);
    } catch (const invalid_argument&) {
        // If the conversion fails, return a default value of 0.0
        return 0.0;
    } catch (const out_of_range&) {
        // If the value is out of range, return a default value of 0.0
        return 0.0;
    }
}

// Function to load a specific element from the CSV file
Element loadElement(const string& symbol, const string& filename) {
    ifstream file(filename);
    string line;

    // Skip header
    getline(file, line);

    while (getline(file, line)) {
        stringstream ss(line);
        Element element;
        string token;

        getline(ss, token, ','); // Atomic Number
        element.atomicNumber = stoi(token);

        getline(ss, token, ','); // Symbol
        element.symbol = token;

        getline(ss, token, ','); // Atomic Mass
        element.atomicMass = safeStod(token);

        getline(ss, token, ','); // Electronegativity
        element.electronegativity = safeStod(token);

        getline(ss, token, ','); // Valency
        element.valency = stoi(token);

        getline(ss, token, ','); // Nature (Positive, Negative, Neutral)
        element.nature = token;

        getline(ss, token, ','); // Ionization Energy (eV)
        element.ionizationEnergy = safeStod(token);

        getline(ss, token, ','); // Shielding Energy (eV)
        element.shieldingEnergy = safeStod(token);

        // If this is the element requested, return it
        if (element.symbol == symbol) {
            return element;
        }
    }

    // If the element is not found, return an empty Element with an error message
    return Element{-1, "", 0.0, 0.0, 0, "", 0.0, 0.0};
}

// Function to generate the SMILES notation for an element
string generateSMILES(const string& elementSymbol, const string& filename) {
    Element element = loadElement(elementSymbol, filename);

    // If the element is not found
    if (element.atomicNumber == -1) {
       return "Element not found!";
    }

    // Handle SMILES based on valency
    string smiles = element.symbol;

    if (element.valency != 0) {
        if (element.valency > 0) {
            // For positive valency, generate "+" signs equal to the valency
            smiles = "[" + element.symbol + string(element.valency, '+') + "]";
        } else if (element.valency < 0) {
            // For negative valency, generate "-" signs equal to the absolute value of the valency
            smiles = "[" + element.symbol + string(element.valency, '-') + "]";
        }
    }

    return smiles;
}

// Main program to take user input and output the SMILES
int main() {
    string filename = "Periodic Table.csv";

    string input;
    cout << "Enter an element symbol (e.g., Na, Cl, H): ";
    cin >> input;

    string smiles = generateSMILES(input, filename);
    cout << "The SMILES notation for " << input << " is: " << smiles << endl;

    return 0;
}
