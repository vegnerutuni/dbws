#ifndef DIFRATOGRAM_H
#define DIFRATOGRAM_H

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

// Adding the code from Ian Madsen (10 feb 98)  (Au,Cr,Fe,Co,Cu,
//                                               Mo,Ag,Ta,W ,Au)
const double XRYZ[10+1] = { 0.0,
    2.748510,2.289620, 1.935970, 1.788965, 1.540520 ,
    0.709260,0.559360, 0.215947, 0.209010, 0.180195 };

class Diffractogram
{
public:
    Diffractogram();

    ifstream file3;
    string file3name;

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

    double S1_,S2_,SS2_,S3_,SS4_,D1_,D2_,D4_,R1_,R2_,R3_,R2NOBK_,R3NOBK_;


    double LAMDA[2+1];
    double LAMDAM;
    double RATIO[2+1];

    void CalcLamdaM(void);
    int SearchWavelength(void);

    void ReadDBWS(void);
    void ReadDBWSF(void);
    void GSASREAD(void);                // SUBROUTINE GSASREAD * READ GSAS FORMATTED DATA FILE
    void PHILIPSREAD(void);             // SUBROUTINE TO READ PHILIPS UDF DATA FILE
    void scintag(void);                 // subroutine to read SCINTAG TXT file
    void SIEMENSREAD(void);             // SUBROUTINE TO READ SIEMENS UXD DATA FILE
    void rigakuread(void);              // SUBROUTINE TO READ RIGAKU DATA FILE
    void ReadSynchrotronData(void);     // SYNCHROTRON X-RAY DATA
    void ReadNeutronData(void);

    void ReadFileBackground();
    void InitPolynomialBackground();
private:
    int LB1_,LB2_,LB3_;
    string TEST_;
    void readasc(void);
};

#endif // DIFRATOGRAM_H
