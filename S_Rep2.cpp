#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
using namespace std;

// Function to load reactivity data from CSV file
unordered_map<string, int> loadReactivityData(const string& filename) {
    unordered_map<string, int> reactivityMap;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return reactivityMap;
    }

    string line, symbol;
    int reactivity;

    // Skip the header
    getline(file, line);

    while (getline(file, line)) {
        stringstream ss(line);
        string token;

        getline(ss, symbol, ','); // Element symbol
        getline(ss, token, ','); // Reactivity level

        if (!symbol.empty() && !token.empty()) {
            reactivity = stoi(token);
            reactivityMap[symbol] = reactivity;
        }
    }

    file.close();
    return reactivityMap;
}

// Function to determine which element in the molecule is active
string getActiveElement(const string& element1, int reactivity1, const string& element2, int reactivity2) {
    if (reactivity1 > reactivity2) {
        return element1;
    } else {
        return element2;
    }
}

// Function to simulate the reaction and determine the products
void predictReaction(const string& element1, int reactivity1, const string& element2, int reactivity2, const string& freeElement, int freeReactivity) {
    cout << "Molecule: " << element1 << element2 << endl;
    cout << "Free Element: " << freeElement << endl;

    if (freeReactivity > reactivity1) {
        // Replace element1 with freeElement
        cout << "Reaction: " << freeElement << element2 << endl;
    } else if (freeReactivity > reactivity2) {
        // Replace element2 with freeElement
        cout << "Reaction: " << element1 << freeElement << endl;
    } else {
        // Reaction not possible under normal conditions
        cout << "Reaction not possible under normal conditions." << endl;

        // Check for specific conditions (e.g., high temperature)
        cout << "Under specific conditions:\n";
        if (reactivity1 > reactivity2) {
            cout << "Reaction: " << freeElement << element2 << endl;
        } else {
            cout << "Reaction: " << element1 << freeElement << endl;
        }
    }
}

int main() {
    // Load reactivity series from CSV file
    string filename = "Reactivity _Series.csv";
    unordered_map<string, int> reactivityMap = loadReactivityData(filename);

    // Input molecule and free element
    string element1, element2, freeElement;
    cout << "Enter first element of the molecule: ";
    cin >> element1;
    cout << "Enter second element of the molecule: ";
    cin >> element2;
    cout << "Enter the free element: ";
    cin >> freeElement;

    // Get reactivity levels for the elements
    int reactivity1 = reactivityMap.count(element1) ? reactivityMap[element1] : 0;
    int reactivity2 = reactivityMap.count(element2) ? reactivityMap[element2] : 0;
    int freeReactivity = reactivityMap.count(freeElement) ? reactivityMap[freeElement] : 0;

    // Display reactivity levels
    cout << "Reactivity of " << element1 << ": " << reactivity1 << endl;
    cout << "Reactivity of " << element2 << ": " << reactivity2 << endl;
    cout << "Reactivity of " << freeElement << ": " << freeReactivity << endl;

    // Predict the reaction
    predictReaction(element1, reactivity1, element2, reactivity2, freeElement, freeReactivity);

    return 0;
}
