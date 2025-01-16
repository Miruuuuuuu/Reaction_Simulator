#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include<windows.h>
using namespace std;
#include <string>
using namespace std;

// Function to convert a number to its subscript representation
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
        }sxa
    }
    return subscript;
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

    // Function to safely convert string to double
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

    // Function to load a specific element from the CSV file
   void loadElement(const string& searchSymbol, const string& filename) {
    ifstream file(filename);
    string line;

    // Skip header
    getline(file, line);

    while (getline(file, line)) {
        stringstream ss(line);
        string token;

        // Read Atomic Number
        getline(ss, token, ',');
        if (token.empty()) {
            cerr << "Error: Missing atomic number in the CSV file." << endl;
            continue; // Skip this row
        }
        this->atomicNumber = stoi(token);

        // Read Symbol
        getline(ss, token, ',');
        if (token.empty()) {
            cerr << "Error: Missing symbol in the CSV file." << endl;
            continue; // Skip this row
        }
        this->symbol = token;

        // Check if it matches the searched symbol
        if (this->symbol == searchSymbol) {
            // Read Atomic Mass
            getline(ss, token, ',');
            this->atomicMass = token.empty() ? 0.0 : safeStod(token);

            // Read Electronegativity
            getline(ss, token, ',');
            this->electronegativity = token.empty() ? 0.0 : safeStod(token);

            // Read Valency
            getline(ss, token, ',');
            if (token.empty()) {
                cerr << "Error: Missing valency for element " << this->symbol << endl;
                this->valency = 0; // Default valency to 0 if missing
            } else {
                this->valency = stoi(token);
            }

            // Read Nature
            getline(ss, token, ',');
            this->nature = token.empty() ? "Unknown" : token;

            // Read Ionization Energy
            getline(ss, token, ',');
            this->ionizationEnergy = token.empty() ? 0.0 : safeStod(token);

            // Read Shielding Energy
            getline(ss, token, ',');
            this->shieldingEnergy = token.empty() ? 0.0 : safeStod(token);

            // Print loaded element details
            cout << "Element loaded: " << endl;
            cout << "Atomic Number: " << atomicNumber << endl;
            cout << "Symbol: " << symbol << endl;
            cout << "Atomic Mass: " << atomicMass << endl;
            cout << "Electronegativity: " << electronegativity << endl;
            cout << "Valency: " << valency << endl;
            cout << "Nature: " << nature << endl;
            cout << "Ionization Energy: " << ionizationEnergy << endl;
            cout << "Shielding Energy: " << shieldingEnergy << endl;
            return; // Exit once the element is found
        }
    }

    // If the element is not found
    cout << "Element " << searchSymbol << " not found in the file." << endl;
}


};
class Synthesis{
    public:
    Element reactant1;
    Element reactant2;
    string product;
    Synthesis(Element a,Element b)
    {
        reactant1=a;
        reactant2=b;
    }
   void createProduct() {
    // Debugging: Print values to check correctness
    cout << "Reactant 1 Symbol: " << reactant1.symbol << endl;
    cout << "Reactant 1 Valency: " << reactant1.valency << endl;
    cout << "Reactant 2 Symbol: " << reactant2.symbol << endl;
    cout << "Reactant 2 Valency: " << reactant2.valency << endl;

    // Forming the product string correctly
            product = reactant1.symbol + toSubscript(reactant2.valency) +
              reactant2.symbol + toSubscript(reactant1.valency);

    cout << endl << "Displaying Product: " << product;
}



};
int main() {
   // Element element;
    //element.loadElement("H", "elements.csv");  // Example usage
    SetConsoleOutputCP(CP_UTF8);
    string filename = "Periodic Table.csv";
   
    string input1,input2,symbol1,symbol2;
    cout << "Enter an element 01 symbol (e.g., Na, Cl, H): ";
    cin >> input1;
    cout << "Enter an element 02 symbol (e.g., Na, Cl, H): ";
    cin >> input2;
Element e1;
Element e2;
e1.loadElement(input1,filename);
 symbol1=e1.symbol;
 e2.loadElement(input2,filename);
 symbol2=e2.symbol;
 
Synthesis s1(e1,e2);
s1.createProduct();
   // string smiles = generateSMILES(input, filename);
   // cout << "The SMILES notation for " << input << " is: " << smiles << endl;
    return 0;
}
