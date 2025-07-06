#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cctype>
#include <regex>
#include <limits>
#include<windows.h>
#include<stack>
#include<vector>
#include<queue>

using namespace std;


// global array for ion cahrge

vector <int> IonChargeArr;
//stack<int >IonChargeStack;
queue<int>IonChargeQueue;
vector<string>positive;
vector<string>negative;
string symbols[4];
int atomCount[4];
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

struct ElementProperties {
    int valency;
    string nature;
};

void CalCompReact()
{
    int prev=IonChargeQueue.back();
    int x= 0-prev;
    IonChargeQueue.push(x);
}

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

int calculateIonCharge(const string& element, int atomCount, const map<string, ElementProperties>& elementMap) {
    auto it = elementMap.find(element);
    if (it == elementMap.end()) {
        throw runtime_error("Element " + element + " not found in map");
    }
    return it->second.valency * atomCount;
}
/*int calculateIonCharge(const string& element, int atomCount) {
    if (elementMap.find(element) == elementMap.end()) {
        throw runtime_error("Error: Element '" + element + "' not found in the periodic table data.");
    }

    if (atomCount <= 0) {
        throw runtime_error("Error: Atom count must be greater than zero.");
    }

    int valency = elementMap[element].valency;
    string nature = elementMap[element].nature;

    if (nature == "neutral") {
        return 0;  // No charge for neutral elements
    }

    long long charge = 0;  // Use long long for large values

    if (nature == "cation") {
        charge = valency * atomCount;  // Positive charge for cations
    } else if (nature == "anion") {
        charge = -(valency * atomCount);  // Negative charge for anions
    } else {
        throw runtime_error("Error: Invalid ion nature for element '" + element + "'.");
    }

    return charge;
}
*/
int calculateCompoundCharge(const string& element1, int atomCount1, 
                          const string& element2, int atomCount2,
                          const map<string, ElementProperties>& elementMap) {
    int charge1 = calculateIonCharge(element1, atomCount1, elementMap);
    int charge2 = calculateIonCharge(element2, atomCount2, elementMap);
    return charge1 + charge2;
}

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

bool isCompound(const string& input) {
    regex elementRegex("[A-Z][a-z]?[0-9]*[A-Z][a-z]?[0-9]*");
    return regex_match(input, elementRegex);
}

bool processInput(const string& input, int atomCount, const map<string, ElementProperties>& elementMap) {
     bool proceed=true;
    try {
       
        if (isCompound(input)) {
            string element1, element2;
            int atomCount1, atomCount2;
            parseCompound(input, element1, atomCount1, element2, atomCount2);
            
            // Adjust atom counts based on the total number of compound molecules
            atomCount1 *= atomCount;
            atomCount2 *= atomCount;
            
            cout << "Compound " << input << " contains:" << endl;
            cout << "- " << element1 << " with " << atomCount1 << " atoms" << endl;
            cout << "- " << element2 << " with " << atomCount2 << " atoms" << endl;
            
            // Calculate and display compound charge
            int totalCharge = calculateCompoundCharge(element1, atomCount1, element2, atomCount2, elementMap);
            cout << "Total ion charge for compound: " << totalCharge << endl;
            if(totalCharge>0)
            {
              positive.push_back(input);
            }
            else if(totalCharge<0)
            {
              negative.push_back(input);
            }
            IonChargeArr.push_back(totalCharge);
            
                  CalCompReact();
            // Display individual element charges
            int charge1 = calculateIonCharge(element1, atomCount1, elementMap);
            int charge2 = calculateIonCharge(element2, atomCount2, elementMap);
            cout << "Individual ion charges:" << endl;
            cout << "- " << element1 << ": " << charge1 << endl;
            cout << "- " << element2 << ": " << charge2 << endl;
            // return totalCharge;
        } else {
            if (elementMap.find(input) != elementMap.end()) {
                int charge = calculateIonCharge(input, atomCount, elementMap);
                cout << "Ion charge for element " << input << ": " << charge << endl;
                IonChargeArr.push_back(charge);
                IonChargeQueue.push(charge);
                if(charge>0)
                {
                    positive.push_back(input);
                }
                else if(charge<0)
                {
                    negative.push_back(input);
                }
                else if(charge==0)
                {
                    cout<<"Error: Neutral Element "<<input<<" Reaction don`t proceed ";
                    proceed=false;
                    
                }
               // return charge;
            } else {
                cout << "Unknown element or invalid compound: " << input << endl;
               proceed= false;
            }
        }
    } catch (const exception& e) {
        cerr << "Error processing " << input << ": " << e.what() << endl;
        proceed=false;
    }
    return proceed;
}
// int CompoundCharge(int )
// {

// }
void displayIonCharge(vector <int>arr)
{
    //int size = sizeof(arr)/sizeof(arr[0]);
    for(int i=0;i<arr.size();i++)
    {
        cout<<"Ion charge for "<<i+1<<" element is: "<<arr[i]<<endl;
    }
}
void displayCompounds(string arr[4], int arr2[4]) {
   // cout << "Reactant 01 is " << arr[0] + to_string(arr2[0]) + arr[1] + to_string(arr2[1]) << endl;
   cout << "Reactant 01 is " << arr[0] + toSubscript(arr2[0]) + arr[1] + toSubscript(arr2[1]) << endl;
   // cout << "Reactant 02 is " << arr[2] + to_string(arr2[2]) + arr[3] + to_string(arr2[3]) << endl;
        cout << "Reactant 02 is " << arr[2] + toSubscript(arr2[2]) + arr[3] + toSubscript(arr2[3]) << endl;

}
// Function to calculate GCD
int gcd(int a, int b) {
    if (b == 0)
        return a;
    return gcd(b, a % b);
}
void displayQueue(queue<int>arr)
{ int i=1;
    // while (!arr.empty())
    // {
    //     int x= arr.top();
    //     cout<<"Element "<<i<<" has Ion Charge of "<<x<<endl;
    //     arr.pop();
    //    i++;
    // }
    while (!arr.empty())
    {
        int x= arr.front();
        cout<<"Element "<<i<<" has Ion Charge of "<<x<<endl;
        arr.pop();
       i++;
    }
    
}
// Function to calculate LCM
int lcm(int a, int b) {
    return (a * b) / gcd(a, b);
}

// Function to balance the ratio
pair<int, int> BalanceAtoms(int pos,int neg,int ElemNo1,int ElemNo2)
{
        int lcm_value = lcm(pos,neg);
        atomCount[ElemNo1]=lcm_value / pos;
        atomCount[ElemNo2]=lcm_value / neg;
        cout<<"Balanced Atoms: "<<atomCount[ElemNo1]<<":"<<atomCount[ElemNo2]<<endl;
         return {atomCount[ElemNo1],atomCount[ElemNo2]};

}
void QueuetoArr(queue <int > queue)
{
    int i=0;
    while (!queue.empty())
    {
        IonChargeArr[i]=queue.front();
        queue.pop();
        i++;
    }
}
void makeCompound()
{
    QueuetoArr(IonChargeQueue);
          pair<int, int> balanced1 = BalanceAtoms(IonChargeArr[0],IonChargeArr[3],0,3);
          if(balanced1.first==0)
    {
        balanced1.first=1;
    }
    else if(balanced1.second==0)
    {
        balanced1.second=1;
    }
    string product1= positive[0]+toSubscript(balanced1.first)+negative[1]+toSubscript(balanced1.second);
    pair<int, int> balanced2 = BalanceAtoms(IonChargeArr[1],IonChargeArr[2],1,2);
    if(balanced2.first==0)
    {
        balanced2.first=1;
    }
    else if(balanced2.second==0)
    {
        balanced2.second=1;
    }
    string product2= positive[1]+ toSubscript(balanced2.second)+negative[0]+toSubscript(balanced2.first);

    cout<<"Displaying Products"<<endl;
    cout<<"Prod1 : "<<product1<<"+"<<"  "<<"Prod2: "<<product2;
}

int main() {


       SetConsoleOutputCP(CP_UTF8);


    string filename = "Periodic_Table2.csv";
  //  int ionchargearr[4];
   
    // vector <string>Compound1;
    // vector <string>Compound2;
   
    try {
        map<string, ElementProperties> elementMap = loadElementData(filename);
        
        for (int i = 0; i < 4; i++) {
            cout << "\nEnter element or compound " << i << " (e.g., Na, Cl, H, NaCl, H2O): ";
            //string input;
            cin >> symbols[i];
            
            
            cout << "Enter the number of atoms/molecules: ";
             atomCount[i] = getValidAtomCount();

       bool possible = processInput(symbols[i], atomCount[i], elementMap);
        if(!possible)
        {
            return 1;
        }
 // ionchargearr[i]=processInput(symbols[i], atomCount, elementMap);

        }
    } catch (const exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return 1;
    }
    
     
   
 displayCompounds(symbols,atomCount);

 displayQueue(IonChargeQueue);
 //displayIonCharge(IonChargeArr);
  //displayQueue(IonChargeQueue);
 makeCompound();
 //BalanceAtoms(IonChargeArr[0],IonChargeArr[2],0,2);
    return 0;
}