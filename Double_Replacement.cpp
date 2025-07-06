#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cctype>
#include <regex>
#include <limits>
#include <windows.h>
#include <stack>
#include <vector>
#include <queue>

using namespace std;

// Global arrays for ion charge
vector<int> IonChargeArr;
queue<int> IonChargeQueue;
vector<string> positive;
vector<string> negative;
string symbols[4];
int atomCount[4];

// Console text colors
const int REACTANT_COLOR = 10;       // Light Green
const int PRODUCT_COLOR = 11;         // Light Cyan
const int ARROW_COLOR = 12;           // Light Red
const int BOX_COLOR = 15;             // White (Bright)
const int DEFAULT_COLOR = 7;          // Gray
const int INFO_COLOR = 14;            // Yellow
const int ERROR_COLOR = 12;           // Light Red
const int COMPOUND_COLOR = 13;        // Light Purple

struct ElementProperties {
    int valency;
    string nature;
};

class CompoundReaction {
public:
    // Function to set console text color
    void setConsoleColor(int color) {
        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
    }

    // Function to reset console to default color
    void resetConsoleColor() {
        setConsoleColor(DEFAULT_COLOR);
    }

    // Function to print styled header
    void printHeader(const string& text, int color = BOX_COLOR) {
        setConsoleColor(color);
        cout << "\n" << text << "\n";
        resetConsoleColor();
    }

    // Function to print styled text
    void printStyled(const string& text, int color) {
        setConsoleColor(color);
        cout << text;
        resetConsoleColor();
    }

    // Function to print the reaction in a formatted box
    void printReactionBox(const string& r1, const string& r2, const string& p1, const string& p2) {
        setConsoleColor(BOX_COLOR);
        cout << "\n╔══════════════════════════════════════════════════════════════╗\n";
        cout << "║                  DOUBLE REPLACEMENT REACTION                 ║\n";
        cout << "╠══════════════════════════════════════════════════════════════╣\n";
        cout << "║  ";
        
        printStyled(r1, REACTANT_COLOR);
        cout << "   +   ";
        printStyled(r2, REACTANT_COLOR);
        cout << "     ";
        printStyled("→", ARROW_COLOR);
        cout << "     ";
        printStyled(p1, PRODUCT_COLOR);
        cout << "   +   ";
        printStyled(p2, PRODUCT_COLOR);
        
        cout << "  ║\n";
        cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
        resetConsoleColor();
    }

    // Convert numbers to subscript format
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

    // Calculate compound reaction
    void CalCompReact() {
        int prev = IonChargeQueue.back();
        int x = 0 - prev;
        IonChargeQueue.push(x);
    }

    // Load element data from CSV file
    map<string, ElementProperties> loadElementData(const string& filename) {
        map<string, ElementProperties> elementMap;
        ifstream file(filename);
        
        if (!file.is_open()) {
            throw runtime_error("Unable to open file: " + filename);
        }

        string line;
        getline(file, line); // Skip header

        while (getline(file, line)) {
            stringstream ss(line);
            string atomicNum, symbol, atomicMass, electronegativity, valencyStr, nature, ionizationEnergy, shieldingEnergy;

            getline(ss, atomicNum, ',');
            getline(ss, symbol, ',');
            getline(ss, atomicMass, ',');
            getline(ss, electronegativity, ',');
            getline(ss, valencyStr, ',');
            getline(ss, nature, ',');
            getline(ss, ionizationEnergy, ',');
            getline(ss, shieldingEnergy, ',');

            try {
                int valency = 0;
                if (!valencyStr.empty()) {
                    valency = stoi(valencyStr);
                    if (nature == "Negative") {
                        valency = -valency;
                    } else if (nature == "Neutral") {
                        valency = 0;
                    }
                }
                elementMap[symbol] = {valency, nature};
            } catch (const exception& e) {
                cerr << "Error processing element " << symbol << ": " << e.what() << endl;
            }
        }

        return elementMap;
    }

    // Calculate ion charge for an element
    int calculateIonCharge(const string& element, int atomCount, const map<string, ElementProperties>& elementMap) {
        auto it = elementMap.find(element);
        if (it == elementMap.end()) {
            throw runtime_error("Element " + element + " not found in map");
        }
        return it->second.valency * atomCount;
    }

    // Calculate charge for a compound
    int calculateCompoundCharge(const string& element1, int atomCount1, 
                              const string& element2, int atomCount2,
                              const map<string, ElementProperties>& elementMap) {
        int charge1 = calculateIonCharge(element1, atomCount1, elementMap);
        int charge2 = calculateIonCharge(element2, atomCount2, elementMap);
        return charge1 + charge2;
    }

    // Get valid atom count from user
    int getValidAtomCount() {
        int atomCount;
        while (true) {
            if (cin >> atomCount) {
                if (atomCount > 0) {
                    return atomCount;
                }
                cout << "Please enter a positive number: ";
            } else {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Invalid input. Please enter a number: ";
            }
        }
    }

    // Parse compound into elements and atom counts
    void parseCompound(const string& compound, string& element1, int& atomCount1, string& element2, int& atomCount2) {
        regex elementRegex("([A-Z][a-z]?)(\\d*)");
        smatch match;
        
        element1 = "";
        element2 = "";
        atomCount1 = 1;
        atomCount2 = 1;

        if (regex_search(compound, match, elementRegex)) {
            element1 = match[1];
            string count = match[2];
            if (!count.empty()) {
                try {
                    atomCount1 = stoi(count);
                } catch (const exception&) {
                    atomCount1 = 1;
                }
            }

            string remainingCompound = compound.substr(match.position(0) + match.length(0));
            if (regex_search(remainingCompound, match, elementRegex)) {
                element2 = match[1];
                string count = match[2];
                if (!count.empty()) {
                    try {
                        atomCount2 = stoi(count);
                    } catch (const exception&) {
                        atomCount2 = 1;
                    }
                }
            }
        }
    }

    // Check if input is a compound
    bool isCompound(const string& input) {
        regex elementRegex("[A-Z][a-z]?[0-9]*[A-Z][a-z]?[0-9]*");
        return regex_match(input, elementRegex);
    }

    // Process user input (element or compound)
    bool processInput(const string& input, int atomCount, const map<string, ElementProperties>& elementMap) {
        bool proceed = true;
        try {
            if (isCompound(input)) {
                string element1, element2;
                int atomCount1, atomCount2;
                parseCompound(input, element1, atomCount1, element2, atomCount2);
                
                // Adjust atom counts based on the total number of compound molecules
                atomCount1 *= atomCount;
                atomCount2 *= atomCount;
                
                printHeader("Compound Analysis:", COMPOUND_COLOR);
                cout << "Compound: ";
                printStyled(input, COMPOUND_COLOR);
                cout << " contains:\n";
                cout << "- ";
                printStyled(element1, INFO_COLOR);
                cout << " with ";
                printStyled(to_string(atomCount1), INFO_COLOR);
                cout << " atoms\n";
                cout << "- ";
                printStyled(element2, INFO_COLOR);
                cout << " with ";
                printStyled(to_string(atomCount2), INFO_COLOR);
                cout << " atoms\n";
                
                // Calculate and display compound charge
                int totalCharge = calculateCompoundCharge(element1, atomCount1, element2, atomCount2, elementMap);
                cout << "Total ion charge: ";
                printStyled(to_string(totalCharge), totalCharge != 0 ? 
                         (totalCharge > 0 ? REACTANT_COLOR : PRODUCT_COLOR) : INFO_COLOR);
                cout << endl;
                
                if(totalCharge > 0) {
                    positive.push_back(input);
                } else if(totalCharge < 0) {
                    negative.push_back(input);
                }
                IonChargeArr.push_back(totalCharge);
                CalCompReact();
                
                // Display individual element charges
                int charge1 = calculateIonCharge(element1, atomCount1, elementMap);
                int charge2 = calculateIonCharge(element2, atomCount2, elementMap);
                cout << "Individual ion charges:\n";
                cout << "- " << element1 << ": ";
                printStyled(to_string(charge1), charge1 != 0 ? 
                          (charge1 > 0 ? REACTANT_COLOR : PRODUCT_COLOR) : INFO_COLOR);
                cout << "\n";
                cout << "- " << element2 << ": ";
                printStyled(to_string(charge2), charge2 != 0 ? 
                          (charge2 > 0 ? REACTANT_COLOR : PRODUCT_COLOR) : INFO_COLOR);
                cout << "\n";
            } else {
                if (elementMap.find(input) != elementMap.end()) {
                    int charge = calculateIonCharge(input, atomCount, elementMap);
                    printHeader("Element Analysis:", COMPOUND_COLOR);
                    cout << "Element: ";
                    printStyled(input, COMPOUND_COLOR);
                    cout << "\nIon charge: ";
                    printStyled(to_string(charge), charge != 0 ? 
                              (charge > 0 ? REACTANT_COLOR : PRODUCT_COLOR) : INFO_COLOR);
                    cout << endl;
                    
                    IonChargeArr.push_back(charge);
                    IonChargeQueue.push(charge);
                    if(charge > 0) {
                        positive.push_back(input);
                    } else if(charge < 0) {
                        negative.push_back(input);
                    } else if(charge == 0) {
                        printStyled("Error: Neutral Element " + input + " - Reaction can't proceed", ERROR_COLOR);
                        cout << endl;
                        proceed = false;
                    }
                } else {
                    printStyled("Unknown element or invalid compound: " + input, ERROR_COLOR);
                    cout << endl;
                    proceed = false;
                }
            }
        } catch (const exception& e) {
            printStyled("Error processing " + input + ": " + e.what(), ERROR_COLOR);
            cout << endl;
            proceed = false;
        }
        return proceed;
    }

    // Display compounds with styling
    void displayCompounds() {
        string r1 = symbols[0] + toSubscript(atomCount[0]) + symbols[1] + toSubscript(atomCount[1]);
        string r2 = symbols[2] + toSubscript(atomCount[2]) + symbols[3] + toSubscript(atomCount[3]);
        
        printHeader("\nREACTANTS SUMMARY:", BOX_COLOR);
        cout << "Reactant 1: ";
        printStyled(r1, REACTANT_COLOR);
        cout << "\nReactant 2: ";
        printStyled(r2, REACTANT_COLOR);
        cout << endl;
    }

    // Function to display ion charge queue
    void displayQueue(queue<int> arr) { 
        int i = 1;
        printHeader("ION CHARGES:", INFO_COLOR);
        while (!arr.empty()) {
            int x = arr.front();
            cout << "Element " << i << " has Ion Charge of ";
            printStyled(to_string(x), x != 0 ? 
                     (x > 0 ? REACTANT_COLOR : PRODUCT_COLOR) : INFO_COLOR);
            cout << endl;
            arr.pop();
            i++;
        }
    }

    // Function to calculate GCD
    int gcd(int a, int b) {
        if (b == 0)
            return a;
        return gcd(b, a % b);
    }

    // Function to calculate LCM
    int lcm(int a, int b) {
        return (a * b) / gcd(a, b);
    }

    // Function to balance atoms
    pair<int, int> BalanceAtoms(int pos, int neg, int ElemNo1, int ElemNo2) {
        int lcm_value = lcm(abs(pos), abs(neg));
        atomCount[ElemNo1] = lcm_value / abs(pos);
        atomCount[ElemNo2] = lcm_value / abs(neg);
        printHeader("BALANCED ATOMS:", BOX_COLOR);
        cout << "Balanced Atoms: ";
        printStyled(to_string(atomCount[ElemNo1]), REACTANT_COLOR);
        cout << ":";
        printStyled(to_string(atomCount[ElemNo2]), PRODUCT_COLOR);
        cout << endl;
        return {atomCount[ElemNo1], atomCount[ElemNo2]};
    }

    // Convert queue to array
    void QueuetoArr(queue<int> queue) {
        int i = 0;
        while (!queue.empty()) {
            IonChargeArr[i] = queue.front();
            queue.pop();
            i++;
        }
    }

    // Create products and display reaction
    void makeCompound() {
        QueuetoArr(IonChargeQueue);
        
        // Create first product
        pair<int, int> balanced1 = BalanceAtoms(IonChargeArr[0], IonChargeArr[3], 0, 3);
        if(balanced1.first == 0) balanced1.first = 1;
        if(balanced1.second == 0) balanced1.second = 1;
        string product1 = positive[0] + toSubscript(balanced1.first) + negative[1] + toSubscript(balanced1.second);
        
        // Create second product
        pair<int, int> balanced2 = BalanceAtoms(IonChargeArr[1], IonChargeArr[2], 1, 2);
        if(balanced2.first == 0) balanced2.first = 1;
        if(balanced2.second == 0) balanced2.second = 1;
        string product2 = positive[1] + toSubscript(balanced2.second) + negative[0] + toSubscript(balanced2.first);

        // Display products
        printHeader("\nPRODUCTS SUMMARY:", BOX_COLOR);
        cout << "Product 1: ";
        printStyled(product1, PRODUCT_COLOR);
        cout << "\nProduct 2: ";
        printStyled(product2, PRODUCT_COLOR);
        cout << endl;
        
        // Build reactant strings
        string r1 = symbols[0] + toSubscript(atomCount[0]) + symbols[1] + toSubscript(atomCount[1]);
        string r2 = symbols[2] + toSubscript(atomCount[2]) + symbols[3] + toSubscript(atomCount[3]);
        
        // Display final reaction in styled box
        printReactionBox(r1, r2, product1, product2);
    }
};

int main() {
    SetConsoleOutputCP(CP_UTF8);
    CompoundReaction reaction;
    
    try {
        map<string, ElementProperties> elementMap = reaction.loadElementData("Periodic_Table2.csv");
        
        reaction.printHeader("DOUBLE REPLACEMENT REACTION SIMULATOR", BOX_COLOR);
        cout << "Enter 4 chemical components (elements or compounds)\n";
        
        for (int i = 0; i < 4; i++) {
            reaction.printHeader("\nEnter chemical component " + to_string(i+1), COMPOUND_COLOR);
            cout << "Enter element/compound: ";
            cin >> symbols[i];
            cout << "Enter atom/molecule count: ";
            atomCount[i] = reaction.getValidAtomCount();

            if (!reaction.processInput(symbols[i], atomCount[i], elementMap)) {
                reaction.printStyled("\nReaction cannot proceed due to invalid input!", ERROR_COLOR);
                return 1;
            }
        }
        
        reaction.displayCompounds();
        reaction.displayQueue(IonChargeQueue);
        reaction.makeCompound();
        
    } catch (const exception& e) {
        reaction.printStyled("Fatal error: " + string(e.what()), ERROR_COLOR);
        return 1;
    }
    
    return 0;
}