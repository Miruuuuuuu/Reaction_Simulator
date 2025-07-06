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

// Keep your existing global variables and structures
vector<int> IonChargeArr;
queue<int> IonChargeQueue;
vector<string> positive;
vector<string> negative;
string symbols[4];
int atomCount[4];

// Keep all your existing functions (toSubscript, ElementProperties, etc.)

void writeReactionData() {
    ofstream file("reaction_data.txt");
    if (!file.is_open()) {
        throw runtime_error("Could not open reaction_data.txt for writing");
    }

    // Write symbols
    for (int i = 0; i < 4; i++) {
        file << symbols[i] << "\n";
    }

    // Write atom counts
    for (int i = 0; i < 4; i++) {
        file << atomCount[i] << "\n";
    }

    // Write ion charges
    file << IonChargeArr.size() << "\n";
    for (int charge : IonChargeArr) {
        file << charge << "\n";
    }

    // Write positive ions
    file << positive.size() << "\n";
    for (const string& ion : positive) {
        file << ion << "\n";
    }

    // Write negative ions
    file << negative.size() << "\n";
    for (const string& ion : negative) {
        file << ion << "\n";
    }

    file.close();
}

bool readReactionResults() {
    ifstream file("reaction_results.txt");
    if (!file.is_open()) {
        cout << "Could not open reaction results file" << endl;
        return false;
    }

    string line;
    // Read status
    getline(file, line);
    if (line != "success") {
        getline(file, line); // Read error message
        cout << "Error in reaction processing: " << line << endl;
        return false;
    }

    // Read balance status
    getline(file, line);
    bool is_balanced = (line == "true");

    if (!is_balanced) {
        cout << "\nReaction Imbalances Detected:" << endl;
        
        // Read number of suggestions
        getline(file, line);
        int num_suggestions = stoi(line);
        
        // Read suggestions
        for (int i = 0; i < num_suggestions; i++) {
            getline(file, line);
            cout << "- " << line << endl;
        }

        // Read updated atom counts
        for (int i = 0; i < 4; i++) {
            getline(file, line);
            atomCount[i] = stoi(line);
        }

        cout << "\nSuggested balanced equation:" << endl;
    } else {
        cout << "\nReaction is properly balanced!" << endl;
    }

    return true;
}

// Modified main function
int main() {
    SetConsoleOutputCP(CP_UTF8);
    string filename = "Periodic_Table2.csv";
    
    try {
        map<string, ElementProperties> elementMap = loadElementData(filename);
        
        // Get input compounds and process them
        for (int i = 0; i < 4; i++) {
            cout << "\nEnter element or compound " << i + 1 << " (e.g., Na, Cl, H, NaCl, H2O): ";
            cin >> symbols[i];
            
            cout << "Enter the number of atoms/molecules: ";
            atomCount[i] = getValidAtomCount();

            bool possible = processInput(symbols[i], atomCount[i], elementMap);
            if (!possible) {
                return 1;
            }
        }

        // Display initial compounds
        displayCompounds(symbols, atomCount);
        displayQueue(IonChargeQueue);

        // Write data for Python processing
        writeReactionData();

        // Call Python script
        cout << "\nVerifying and balancing reaction..." << endl;
        system("python reaction_balancer.py");

        // Read and process results
        if (readReactionResults()) {
            // Make final compound with possibly updated coefficients
            makeCompound();
        }

    } catch (const exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return 1;
    }

    return 0;
} 
