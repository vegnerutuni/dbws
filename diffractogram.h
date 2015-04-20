#ifndef DIFRATOGRAM_H
#define DIFRATOGRAM_H

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "param_inc.h"
#include "dbwsexception.h"

using namespace std;

class Diffractogram
{
public:
    Diffractogram();

    ifstream file4;
    ofstream file4b;

    string file4name;

    ofstream* file6;


    double THMIN;
    double STEP;
    double THMAX;
    string DATAID;

    double Y[IDSZ+1];
    double VAR[IDSZ+1];
    double YC[IDSZ+1];
    int KR[IDSZ+1];  // phs
    double BK[IDSZ+1];
    int NPTS;
    double AMORPHOUS[IDSZ+1];

    double BKCOM[IDSZ+1];
    double BKDIS[IDSZ+1];
    double BKAM[IDSZ+1];
    double BKPOL[IDSZ+1];

    double ALOW[100+1];
    double AHIGH[100+1];
    double POS[100+1];
    double BCK[100+1];
    int NBCKGD;
    int NEXCRG;
    int NSCAT; //nscatx

    double XMAXINT;

    void ReadDBWS(void);
    void ReadDBWSF(void);
    void GSASREAD(void);           // SUBROUTINE GSASREAD * READ GSAS FORMATTED DATA FILE
    void PHILIPSREAD(void);        // SUBROUTINE TO READ PHILIPS UDF DATA FILE
    void scintag(void);            // subroutine to read SCINTAG TXT file
    void SIEMENSREAD(void);        // SUBROUTINE TO READ SIEMENS UXD DATA FILE
    void rigakuread(void);         // SUBROUTINE TO READ RIGAKU DATA FILE
private:
    int LB1_,LB2_,LB3_;
    string TEST_;
    void readasc(void);
};

#endif // DIFRATOGRAM_H
