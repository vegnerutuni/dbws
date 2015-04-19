#ifndef PHASE_H
#define PHASE_H

#include <string>
#include "parameter.h"
#include "spacegroup.h"


using namespace std;

class Phase
{
public:
    Phase();

    SpaceGroup SYMB;   //string SYMB;
    string name;       // Old PHSNM

    Parameter PAR[30];  // INDICE JA ESTA EM -1

    ParameterInteger AtomCount;
    ParameterInteger NMOL;
    double PREF[3+1];

    ParticleAbsorption SAQF;
    double WTIS;
};

#endif // PHASE_H
