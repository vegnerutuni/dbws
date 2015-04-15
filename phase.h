#ifndef PHASE_H
#define PHASE_H

#include <string>
#include "parameter.h"


using namespace std;

class Phase
{
public:
    Phase();

    string SYMB;
    string PHSNM;

    Parameter PAR_[30];
};

#endif // PHASE_H
