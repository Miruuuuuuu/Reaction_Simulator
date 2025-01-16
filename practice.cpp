#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
#include <string>

using namespace std;

void abc(string a, int a1, string b, int a2, string c, int a3) {
    if (a3 > a1) {
        string molecule = a + c;
        cout << "Displaying Product: " << molecule << endl;
        cout<<"Product 2 :"<<b<<endl;
    } else if (a3 > a2) {
        string molecule = b + c;
        cout << "Displaying Product: " << molecule << endl;
        cout<<"Product 2 : "<<a<<endl;
    } else if (a3 < a1 && a3 < a2) {
        cout << "Reaction not possible in normal conditions." << endl;
        cout << "But in specific conditions, the chemical reaction will be: ";
        if (a1 > a2) {
            string molecule = c + a;
            cout << "Displaying Product: " << molecule << endl;
            cout<<"Product 2 : "<<b<<endl;
        } else if (a2 > a1) {
            string molecule = b + c;
            cout << "Displaying Product: " << molecule << endl;
            cout<<"Product 2 : "<<a<<endl;
        }
    }
}

int loadReactivity(string element, string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filename << endl;
        return 0;
        
    }

    string line;
    // Skip header
   // getline(file, line);

    while (getline(file, line)) {
        stringstream ss(line);
        string symbol;
        string reactivityStr;

        // Read the element symbol
        getline(ss, symbol, ',');
        if (symbol.empty()) {
            cerr << "Error: Missing symbol in the CSV file." << endl;
            continue; // Skip invalid rows
        }

        // Check if the symbol matches the input element
        if (symbol == element) {
            // Read the reactivity value
            getline(ss, reactivityStr, ',');
            if (reactivityStr.empty()) {
                cerr << "Error: Missing reactivity in the CSV file." << endl;
                return 0; // Return default reactivity for invalid rows
            }

            // Convert reactivity to integer and return
            try {
                return stoi(reactivityStr);
            } catch (const invalid_argument& e) {
                cerr << "Error: Invalid reactivity value in the CSV file for element: " << symbol << endl;
                return 0;
            }
        }
    }

    // If the element is not found, return a default reactivity of 0
    return 0;
}

int main() {
    string filename = "Reactivity _Series.csv";
    int no_of_ele;
    cout<<"Enter the number of elements:";
    cin>>no_of_ele;
     string element[no_of_ele];
     int num[no_of_ele];
    for(int i=0;i<no_of_ele;i++)
    {
        cout<<"Enter the element "<<i+1<<":";
        cin>>element[i];
    }
    cout<<"Displaying the elements:"<<endl;
    for(int i=0;i<no_of_ele;i++)
    {
         num[i] = loadReactivity(element[i], filename);
        cout<<"Reactivity of "<<element[i]<<" is :"<<num[i]<<endl;
    }
    // int num1 = loadReactivity("Na", filename);
    // cout << "Reactivity of Na: " << num1 << endl;

    // int num2 = loadReactivity("Cl", filename);
    // cout << "Reactivity of Cl: " << num2 << endl;

    // int num3 = loadReactivity("K", filename);
    // cout << "Reactivity of K: " << num3 << endl;

    // Call the `abc` function with sample data
  abc(element[0], num[0], element[1], num[1], element[2], num[2]);

    return 0;
}
