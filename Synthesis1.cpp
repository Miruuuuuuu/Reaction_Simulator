
//#include"Synthesis1.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <windows.h>
#include <vector>
using namespace std;

// Function to set text color and background
void setColor(int textColor, int bgColor = 0) {
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, textColor | (bgColor << 4));
}

// Function to reset color to default
void resetColor() {
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, 7); // Default color
}

// Enum for colors
enum Colors {
    RED = FOREGROUND_RED | FOREGROUND_INTENSITY,
    GREEN = FOREGROUND_GREEN | FOREGROUND_INTENSITY,
    BLUE = FOREGROUND_BLUE | FOREGROUND_INTENSITY,
    WHITE_BG = BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE
};

// Converts numbers to subscript
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

class Element {
public:
    int atomicNumber;
    string symbol;
    double atomicMass;
    double electronegativity;
    int valency;
    string nature;
    double ionizationEnergy;
    double shieldingEnergy;

    Element() {
        atomicNumber = 0;
        symbol = "";
        atomicMass = 0.0;
        electronegativity = 0.0;
        valency = 0;
        nature = "";
        ionizationEnergy = 0.0;
        shieldingEnergy = 0.0;
    }

    double safeStod(const string& str) {
        try {
            return stod(str);
        } catch (const invalid_argument&) {
            setColor(RED, WHITE_BG);
            cerr << "Invalid argument for converting to double: " << str << endl;
            resetColor();
            return 0.0;
        } catch (const out_of_range&) {
            setColor(RED, WHITE_BG);
            cerr << "Out of range error for converting to double: " << str << endl;
            resetColor();
            return 0.0;
        }
    }

    void loadElement(const string& searchSymbol, const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            setColor(RED, WHITE_BG);
            cerr << "Error: Could not open file " << filename << endl;
            resetColor();
            return;
        }

        string line;
        getline(file, line);

        while (getline(file, line)) {
            stringstream ss(line);
            string token;
            vector<string> tokens;

            while (getline(ss, token, ',')) {
                tokens.push_back(token);
            }

            if (tokens.size() < 8) {
                continue;
            }

            symbol = tokens[1];
            if (symbol == searchSymbol) {
                atomicNumber = stoi(tokens[0]);
                atomicMass = safeStod(tokens[2]);
                electronegativity = safeStod(tokens[3]);
                valency = stoi(tokens[4]);
                nature = tokens[5];
                ionizationEnergy = safeStod(tokens[6]);
                shieldingEnergy = safeStod(tokens[7]);

                setColor(GREEN);
                cout << "\nElement Details:" << endl;
                setColor(BLUE);
                cout << "-------------------------------" << endl;
                cout << "Symbol:           " << symbol << endl;
                cout << "Nature:           " << nature << endl;
                cout << "Valency:          " << valency << endl;
                cout << "Atomic Number:    " << atomicNumber << endl;
                cout << "Atomic Mass:      " << atomicMass << endl;
                cout << "Ionization Energy: " << ionizationEnergy << endl;
                cout << "Shielding Energy: " << shieldingEnergy << endl;
                cout << "-------------------------------" << endl;
                resetColor();
                return;
            }
        }
        setColor(RED, WHITE_BG);
        cout << "Element " << searchSymbol << " not found in the file." << endl;
        resetColor();
    }
};

class Synthesis : public Element {
public:
    Element reactant1;
    Element reactant2;
    string product;

    Synthesis(string a, string b, string filename) {
        reactant1.loadElement(a, filename);
        reactant2.loadElement(b, filename);
    }

    bool canFormProduct() {
        cout << "\nChecking Compatibility:" << endl;
        cout << reactant1.symbol << " (Nature: " << reactant1.nature << ")" << endl;
        cout << reactant2.symbol << " (Nature: " << reactant2.nature << ")" << endl;

        if (reactant1.nature == "Neutral" || reactant2.nature == "Neutral") {
            setColor(RED, WHITE_BG);
            cout << "Error: Cannot form product with neutral elements!" << endl;
            resetColor();
            return false;
        }

        bool isCompatible = (reactant1.nature == "Positive" && reactant2.nature == "Negative") ||
                          (reactant1.nature == "Negative" && reactant2.nature == "Positive");

        if (!isCompatible) {
            setColor(RED, WHITE_BG);
            cout << "Error: Elements must have opposite charges (Positive-Negative) to form a product!" << endl;
            resetColor();
            return false;
        }

        return true;
    }

    void createProduct() {
        if (canFormProduct()) {
            product = reactant1.symbol + toSubscript(reactant2.valency) +
                     reactant2.symbol + toSubscript(reactant1.valency);

            setColor(GREEN);
            cout << "\nProduct Successfully Formed!" << endl;
            resetColor();
        } else {
            product = "No product formed";
            setColor(RED, WHITE_BG);
            cout << "Product formation failed due to incompatible element nature." << endl;
            resetColor();
        }
    }

    void displayReaction() {
        setColor(GREEN);
        cout << "\n-----------------------------------" << endl;
        cout << "            Reaction Box           " << endl;
        cout << "-----------------------------------" << endl;
        cout << "       " << reactant1.symbol << " + " << reactant2.symbol << " → " << product << endl;
        cout << "-----------------------------------" << endl;
        resetColor();
    }
};

int main() {
    SetConsoleOutputCP(CP_UTF8);
    string filename = "Periodic_Table2.csv";

    string input1, input2;
    cout << "\nWelcome to the Chemical Reaction Simulator!" << endl;
    setColor(BLUE);
    cout << "--------------------------------------------" << endl;
    cout << "Enter symbols for two elements (e.g., Na, Cl):" << endl;
    cout << "--------------------------------------------" << endl;
    resetColor();
    cout << "Enter Element 1: ";
    cin >> input1;
    cout << "Enter Element 2: ";
    cin >> input2;

    Synthesis s1(input1, input2, filename);
    s1.createProduct();
    s1.displayReaction();

    return 0;
}