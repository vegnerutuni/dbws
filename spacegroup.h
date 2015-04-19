#ifndef SPACEGROUP_H
#define SPACEGROUP_H

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>


using namespace std;

class SpaceGroup
{
public:
    SpaceGroup();
    void set(string);
    string get(void);
    SpaceGroup& operator=(string);
    operator string();
    friend ostream& operator<<(ostream&, SpaceGroup&);

    //ofstream* outputfile;
    int NSPGRP;
    int NAXIS;
    int NCONT[10+1];
    int NC;
    int NPOL;
    int MULTP; // multpommon

    void SPGP(string SPG);
    void GOTOER(void);
private:
    string value;
};

#endif // SPACEGROUP_H

// TODO: talvez mudar LAU para ca
/*
// Laue classes
const string LAU[14] = {
    "1BAR",			// Triclini
    "2/M",			// Monoclinic
    "MMM",			// Orthorhombic
    "4/M",			// Tetragonal
    "4/MMM",		// "
    "3BAR   R",		// Trigonal, romboedric
    "3BAR M R",		// Trigonal
    "3BAR",			// Trigonal
    "3BAR M 1",		// Trigonal
    "3BAR 1 M",
    "6/M",			// Hexagonal
    "6/MMM",		// Hexagonal
    "M3",			// Cubico
    "M3M"			// Cubico
};
*/
