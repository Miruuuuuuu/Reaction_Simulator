#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <windows.h>
#include <unordered_map>
using namespace std;

string toSubscript(int number) {
    string subscript = "";
    string normal = to_string(number);

    for (char c : normal) {
        switch (c) {
            case '0': subscript += "₀"; break;
            case '1': subscript += "₁"; break;
            case '2': subscript += "₂"; break;
            case '3': subscript += "₃"; break;
            case '4': subscript += "₄"; break;
            case '5': subscript += "₅"; break;
            case '6': subscript += "₆"; break;
            case '7': subscript += "₇"; break;
            case '8': subscript += "₈"; break;
            case '9': subscript += "₉"; break;
        }
    }
    return subscript;
}
unordered_map<string, int> loadReactivityLevels(const string& filename) {
    unordered_map<string, int> reactivityMap;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string symbol, reactivity;
        getline(ss, symbol, ',');     // Read element symbol
        getline(ss, reactivity, ','); // Read reactivity level
        reactivityMap[symbol] = stoi(reactivity);
    }

    return reactivityMap;
}
class Element {
public:
    int atomicNumber;
    string symbol;
    double atomicMass;
    double electronegativity;
    int valency;
    string nature;  // Positive, Negative, Neutral
    double ionizationEnergy;  // In eV
    double shieldingEnergy;   // In eV
    int atoms;  // Number of atoms in the molecule

    double safeStod(const string& str) {
        try {
            return stod(str);
        } catch (const invalid_argument&) {
            cerr << "Invalid argument for converting to double: " << str << endl;
            return 0.0;
        } catch (const out_of_range&) {
            cerr << "Out of range error for converting to double: " << str << endl;
            return 0.0;
        }
    }

    void loadElement(const string& searchSymbol, const string& filename) {
        ifstream file(filename);
        string line;

        // Skip header
        getline(file, line);

        while (getline(file, line)) {
            stringstream ss(line);
            string token;

            getline(ss, token, ','); // Read Atomic Number
            if (token.empty()) {
                cerr << "Error: Missing atomic number in the CSV file." << endl;
                continue;
            }
            this->atomicNumber = stoi(token);

            getline(ss, token, ','); // Read Symbol
            if (token.empty()) {
                cerr << "Error: Missing symbol in the CSV file." << endl;
                continue;
            }
            this->symbol = token;

            if (this->symbol == searchSymbol) {
                getline(ss, token, ','); // Atomic Mass
                this->atomicMass = token.empty() ? 0.0 : safeStod(token);

                getline(ss, token, ','); // Electronegativity
                this->electronegativity = token.empty() ? 0.0 : safeStod(token);

                getline(ss, token, ','); // Valency
                this->valency = token.empty() ? 0 : stoi(token);

                getline(ss, token, ','); // Nature
                this->nature = token.empty() ? "Unknown" : token;

                getline(ss, token, ','); // Ionization Energy
                this->ionizationEnergy = token.empty() ? 0.0 : safeStod(token);

                getline(ss, token, ','); // Shielding Energy
                this->shieldingEnergy = token.empty() ? 0.0 : safeStod(token);

                cout << "Element loaded: " << endl;
                cout << "Atomic Number: " << atomicNumber << endl;
                cout << "Symbol: " << symbol << endl;
                cout << "Atomic Mass: " << atomicMass << endl;
                cout << "Electronegativity: " << electronegativity << endl;
                cout << "Valency: " << valency << endl;
                cout << "Nature: " << nature << endl;
                cout << "Ionization Energy: " << ionizationEnergy << endl;
                cout << "Shielding Energy: " << shieldingEnergy << endl;
                return;
            }
        }
        cout << "Element " << searchSymbol << " not found in the file." << endl;
    }

    // Function to combine two elements and create a molecule
    string createMolecule(const Element& element1, const Element& element2) {
        // Generate molecular formula based on the new formula format
        string molecule = element1.symbol + to_string(element1.atoms) +
                          element2.symbol + to_string(element2.atoms);
        return molecule;
    }
};

class SingleReplacement {
private:
    unordered_map<string, int> reactivityMap;

public:
    SingleReplacement(const string& reactivityFile) {
        reactivityMap = loadReactivityLevels(reactivityFile);
    }

    void predictProduct(const string& mol1Symbol, const string& mol2Symbol, const string& freeElementSymbol) {
        // Check if elements exist in the reactivity map
        if (reactivityMap.find(mol1Symbol) == reactivityMap.end() ||
            reactivityMap.find(mol2Symbol) == reactivityMap.end() ||
            reactivityMap.find(freeElementSymbol) == reactivityMap.end()) {
            cout << "Error: Reactivity data for one or more elements not found." << endl;
            printReactivityMap(); // Debugging: Print available elements
            return;
        }

        // Retrieve reactivity levels
        int mol1Reactivity = reactivityMap[mol1Symbol];
        int mol2Reactivity = reactivityMap[mol2Symbol];
        int freeElementReactivity = reactivityMap[freeElementSymbol];

        // Debugging: Print reactivity levels
        cout << "Reactivity of " << mol1Symbol << ": " << mol1Reactivity << endl;
        cout << "Reactivity of " << mol2Symbol << ": " << mol2Reactivity << endl;
        cout << "Reactivity of " << freeElementSymbol << ": " << freeElementReactivity << endl;

        // Determine which element in the molecule will be replaced
        string prod1, prod2;
        if (freeElementReactivity > mol1Reactivity && freeElementReactivity > mol2Reactivity) {
            if (mol1Reactivity < mol2Reactivity) {
                // Replace mol1 with the free element
                prod1 = freeElementSymbol;
                prod2 = mol2Symbol;
            } else {
                // Replace mol2 with the free element
                prod1 = mol1Symbol;
                prod2 = freeElementSymbol;
            }

            cout << "Reaction occurs! The free element replaces the less reactive element in the compound." << endl;
            cout << "Predicted products: " << endl;
            cout << "Product 1: " << prod1 << endl;
            cout << "Product 2: " << prod2 << endl;
        } else {
            cout << "No reaction occurs: The free element is not reactive enough to replace any element in the compound." << endl;
        }
    }
};



int main() {
    SetConsoleOutputCP(CP_UTF8);
    string filename = "Periodic Table.csv";
    string reactivityFile = "Reactivity.csv";       // Reactivity data

    // Load elements and reactivity levels
    SingleReplacement sr("Reactivity_Series.csv");
     Element e1, e2,e3;
    string input1, input2, symbol1, symbol2;
    cout << "Enter an element 01 symbol (e.g., Na, Cl, H): ";
    cin >> input1;
    //cout<<"Enter the atoms of Element 01:";
    //cin>>e1.atoms;
    cout << "Enter an element 02 symbol (e.g., Na, Cl, H): ";
    cin >> input2;
    //cout<<"Enter the atoms of Element 02:";
    //cin>>e2.atoms;
  
  cout<<"Enter the Symbol of free Element :";
  cin>>e3.symbol;
  //cout<<"Enter the atoms of free Element :";
    //cin>>e3.atoms;
    // e1.loadElement(input1, filename);
    // symbol1 = e1.symbol;
    // e2.loadElement(input2, filename);
    // symbol2 = e2.symbol;

    // Create the molecule
    //string molecule = e1.createMolecule(e1, e2);

   // cout << endl << "Molecular Formula: " << molecule << endl;
        sr.predictProduct(e1.symbol,e2.symbol, e3.symbol);

    return 0;
}
