#ifndef REFLECTIONS_H
#define REFLECTIONS_H


#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cmath>

#include "param_inc.h"
#include "dbwsexception.h"

using namespace std;

class Reflections
{
public:
    Reflections();

    double REFS[3+1];
    int IREFS;
    double FMGNTD;
    double HALFL;
    double HALFG;
    double GAM;
    double FWHM[2+1];

    double APHASE;
    double TAVIX;
    double SRIX;
};

#endif // REFLECTIONS_H
