#ifndef SYNTHESIS1_H
#define SYNTHESIS1_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <windows.h>
#include <vector>

using namespace std;

// Function declarations and class definitions from Synthesis1.cpp
void setColor(int textColor, int bgColor = 0);
void resetColor();
string toSubscript(int number);

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

    Element();
    double safeStod(const string& str);
    void loadElement(const string& searchSymbol, const string& filename);
};

class Synthesis : public Element {
public:
    Element reactant1;
    Element reactant2;
    string product;

    Synthesis(string a, string b, string filename);
    bool canFormProduct();
    void createProduct();
    void displayReaction();
};

#endif // SYNTHESIS1_H