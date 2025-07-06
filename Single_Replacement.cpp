#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
#include <string>
#include <windows.h>

using namespace std;

// Function to set text color and background
void setColor(int textColor, int bgColor = 0) {
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, textColor | (bgColor << 4));
}

// Function to reset color to default
void resetColor() {
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleTextAttribute(hConsole, 7);  // 7 is default (white text on black background)
}

// Enum for colors
enum Colors {
    RED = FOREGROUND_RED | FOREGROUND_INTENSITY,
    GREEN = FOREGROUND_GREEN | FOREGROUND_INTENSITY,
    BLUE = FOREGROUND_BLUE | FOREGROUND_INTENSITY,
    YELLOW = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY,
    WHITE_BG = BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE
};

// Function to load reactivity from file
int loadReactivity(string element, string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        setColor(RED, WHITE_BG);
        cerr << "Error: Unable to open the file: " << filename << endl;
        resetColor();
        return 0;
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string symbol;
        string reactivityStr;

        // Read the element symbol
        getline(ss, symbol, ',');
        if (symbol.empty()) {
            setColor(RED, WHITE_BG);
            cerr << "Error: Missing symbol in the CSV file." << endl;
            resetColor();
            continue; // Skip invalid rows
        }

        // Check if the symbol matches the input element
        if (symbol == element) {
            // Read the reactivity value
            getline(ss, reactivityStr, ',');
            if (reactivityStr.empty()) {
                setColor(RED, WHITE_BG);
                cerr << "Error: Missing reactivity in the CSV file." << endl;
                resetColor();
                return 0;
            }

            // Convert reactivity to integer and return
            try {
                return stoi(reactivityStr);
            } catch (const invalid_argument& e) {
                setColor(RED, WHITE_BG);
                cerr << "Error: Invalid reactivity value in the CSV file for element: " << symbol << endl;
                resetColor();
                return 0;
            }
        }
    }

    // If the element is not found, return a default reactivity of 0
    setColor(RED, WHITE_BG);
    cout << "Element not found: " << element << endl;
    resetColor();
    return 0;
}

// Struct to represent elements
struct Elements {
    string symbol;
    int Reactivity;
};

// Reaction class to form and display chemical reactions
class Reaction {
public:
    string Product1;
    string Product2;
    string Reactant1;
    string Reactant2;

    // Function to form the product based on reactivity
    void formProduct(Elements ElemArr[]) {
        int a1 = ElemArr[0].Reactivity;
        int a2 = ElemArr[1].Reactivity;
        int a3 = ElemArr[2].Reactivity;
        string a = ElemArr[0].symbol;
        string b = ElemArr[1].symbol;
        string c = ElemArr[2].symbol;

        Reactant1 = a + b;
        Reactant2 = c;

        if (a3 > a1) {
            Product1 = c + a;
            setColor(GREEN, WHITE_BG);
            cout << "Displaying Product: " << Product1 << endl;
            resetColor();
            Product2 = b;
            cout << "Product 2: " << Product2 << endl;
        } else if (a3 > a2) {
            Product1 = c + b;
            setColor(GREEN, WHITE_BG);
            cout << "Displaying Product: " << Product1 << endl;
            resetColor();
            Product2 = a;
            cout << "Product 2: " << Product2 << endl;
        } else if (a3 < a1 && a3 < a2) {
            cout << "Reaction not possible in Normal Conditions." << endl;
            cout << "But in Specific Conditions, The Chemical Reaction can be: ";
            if (a1 > a2) {
                Product1 = a + c;
                setColor(GREEN, WHITE_BG);
                cout << "Displaying Product: " << Product1 << endl;
                resetColor();
                Product2 = b;
                cout << "Product 2: " << Product2 << endl;
            } else if (a2 > a1) {
                Product1 = c + b;
                setColor(GREEN, WHITE_BG);
                cout << "Displaying Product: " << Product1 << endl;
                resetColor();
                Product2 = a;
                cout << "Product 2: " << Product2 << endl;
            }
        }
    }

    // Function to display the reaction
    void displayReaction() {
        setColor(BLUE, WHITE_BG);
        cout << "\nDisplaying the Reaction:" << endl;
        cout << Reactant1 << " + " << Reactant2 << " = " << Product1 << " + " << Product2 << endl;
        resetColor();
    }
};

int main() {
    // Setting the console output to UTF-8 for special characters
    SetConsoleOutputCP(CP_UTF8);

    // ASCII Art for "Virtual Chemistry Lab"
    

    string filename = "Reactivity _Series.csv";
    Elements elements[3];
    Reaction reaction;

    // Getting element symbols from user input
    for (int i = 0; i < 3; i++) {
        cout << "Enter the element " << i + 1 << " symbol (e.g., Na, Cl, H): ";
        cin >> elements[i].symbol;
    }

    cout << "\nDisplaying the elements and their reactivity:" << endl;
    for (int i = 0; i < 3; i++) {
        elements[i].Reactivity = loadReactivity(elements[i].symbol, filename);
        if (elements[i].Reactivity > 0) {
            setColor(GREEN);
            cout << "Reactivity of " << elements[i].symbol << " is: " << elements[i].Reactivity << endl;
            resetColor();
        }
    }

    // Form and display the product of the chemical reaction
    reaction.formProduct(elements);
    reaction.displayReaction();

    return 0;
}
