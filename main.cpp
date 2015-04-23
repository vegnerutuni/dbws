
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <stdlib.h>

#include "param_inc.h"
#include "tabelas.h"
#include "parameter.h"
#include "phase.h"
#include "diffractogram.h"
#include "dbwsexception.h"
#include "reflections.h"

using namespace std;


// NPROF Profile selection
const int _Gaussian			= 0;
const int _Lorentzian		= 1;
const int _Mod1				= 2;
const int _Mod2				= 3;
const int _SplitPearsonVII	= 4;
const int _pseudoVoigt		= 5;
const int _PearsonVII		= 6;
const int _TCHZ				= 7;






class COEFF
{
public:
    double AC[10+1][16+1];
    double POSI[30+1];
    double SCAT[30+1];
    double DFP[16+1];
    double DFPP[16+1];
    double XMAS[16+1];

    int N;
    double SINLAM[30+1];
    double F[30+1];

    void CV1(double AB[],double V[]);
    void CV2(double AB[],double* G);
    void STEEP(double X[], int N, int M);
    void COEF(int* J, int* K);
};










class DBWS
{
public:
    DBWS(void);
    virtual ~DBWS(void);
    void run(void);

public:
    string title;

    int JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,IPLCOM,IPLDIS,IPLAM,IPBIG,MCYCLE;
    double EPS;
    int MAXS, INSTRM,MAXSX,MCYCLX,ICYRUN,IPREF,IABSR;
    int FONDO;
    int IAS,IASYM;
    double SW;
    int IBGD,IDATA;
    int ISPHASE; // Internal Standard Phase

    Parameter GLB_[20];  // GLOBAL PARAMETERS

    Parameter XL_[NATS][11 +1];

    int PTR_[NATS +1];
    double RELAX[4 +1];


    double SAVE[99+1][6+1];
    int NSAVE;



    double TMASSA[99+1];
    int MLTPHASE;


    Phase phases[99+1];


    int MURT[NATS+1];

    string ATEXT[NATS+1];
    string NTYP[NATS+1];


    double U_;
    double V_;
    double W_;
    double TMV_;
    double CTHM_;
    double RLIM_;
    double BKPOS_;
    double WDT_;
    double ULOR_;
    double VLOR_;
    double ZZZ_;
    double UC_;

    double DELT_;
    double TL_;
    double GAM1_;
    double GAM2_;
    double PRFDER_;
    int IPH_;
    double DELTA_;

    double TF1_,TF2_,TF4_,TF6_,TF8_,TF9_,C4_;

    double RL_,DA1L_,DA2L_,DA3L_,DA4L_,DA5L_,DA6L_,DA7L_,
        RH_,DA1H_,DA2H_,DA3H_,DA4H_,DA5H_,DA6H_,DA7H_;

    Reflections refs[IRS+1];

    Diffractogram diffractogram;


    double FINAL_[NFINAL+1][2+1];
    int ILOC_;

    string NAM_[16+1];

    COEFF coeff;

    double SIZEG[15+1];
    double STRAING[15+1];
    double SIZEL[15+1];
    double STRAINL[15+1];
    double SIZ[15+1];
    double STRAIN[15+1];
    int NSIZESTRAIN;

    double RJAC[MSZ+1][MSZ+1]; // SMM
    double VX[MSZ+1]; // V1, V

    double _COND;
    int _IORD1;
    int _IORD2;
    double _TH;
    int _NUM;

    int IHKL[3+1][48+1];
    double AZ[48+1];
    int NC1[10+1][7+1][99+1];
    int ICHKL[99+1];
    int N1HKL[99+1];
    int IER;


    double cell_A;
    double cell_B;
    double cell_C;
    double cell_ALPHA;
    double cell_BETA;
    double cell_GAMMA;
    double cell_AL[3+1][3+1];

    double DCSM[6+1][6+1];
    double DCV[6+1];

    int IBCKCODE;

    double TANN_[NOV+1];
    double DERSTO_[NOV+1][MSZ+1];

    double CC_[4+1][16+1];
    double ZEFF_[16+1];
    int PTC_[NATS+1];


    double FI2_[2+1];










    ifstream file5;
    ofstream file5b;

    string file5name;

    ofstream file6;
    string file6name;


    ofstream file8o;
    ifstream file8i;
    string file8name;

    ifstream file11;
    string file11name;

    ofstream file20;


public:
    int sign(int i);
    double sign(double i);
    void GOTOER(void);


    void SORT(int IPHASE);
    void ASSIGN_(void);
    void CHISQ(void);
    void COMPTON(int K, double STH, double* CISK);
    void DISORDER(int K, double STH, int IDERIV, double* SDK, double* DYC, int FONDO, double DERISO[]);
    void SMTRY2(int* IPHASE);
    void CALCUL(int NN);
    void PRSVII(double T);
    double gamma(double X);
    void mspvii(double A, double W);
    double PROFIL(int N, double X);
    void SUMMAT(int IPM, double ISODER[], double TOTCS);
    double DPINV(double A1[][MSZ+1], double B1[], int* N);
    double ERROR(double A[][6+1], double B[], double* OMEGA);
    void ESD(double SM[][6+1], double V[], double SUM);
    void DIRECT(double SM[][6+1], double V[], int* IPH);
    void OUTSCR(int I, double R2, double R3, double X);
    void OUTPTR(int ICYCLE);
    void ITER(void);
    void CEL000(int* MULTX_, double Y_[][3+1]);
    void OPERTR(int* MULTX_, double Y_[][3+1], int* IPHASE, int* I, int* L2);
    void SYMOPR(int* MULTX_, double Y_[][3+1], double* XLT_, int* IXB_, int NCTR_[], int* IPHASE);
    void RTMT(int* MULTX_, double Y_[][3+1], int* IPRT, int NCTR_[], int* ISIMOP_, int* IPHASE);
    void OP1(int* IPHASE, int NCTR_[]);
    void LOOKUP(int K, int N, int NSCAT, int IXRAY);
    void CELL2(int NPHASE);
    double MULT(int IPHASE, int IH, int IK, int IL, int KXIS);
    void REWRIT(int ISCALE, int IDIF);
    void size(int K);
    void EXPUT(void);
    void REFGEN(int IPHASE, double ZERO, double DIS, double TRANS, double PREFOR, int NCTR_[]);
    void FINDC(int K, int NSCAT);
    void ABSORP(double MU, double SW, double TH, double* ABC);

    void qpainit(void);



    void INPTR(void);
    void inpam(void);

};




void COEFF::CV1(double AB[],double V[])
{
    int I,J,K,MAXX,NN;

    MAXX=min(9,N);
    for (K=1; K <= MAXX; K = K + 2)
    {
        V[K]=0;
        for (I=1; I <= N; ++I)
        {
            V[K]=V[K]+F[I];
            for (J=1; J <= MAXX; J = J + 2)
            {
                V[K]=V[K]-AB[J]*exp(-AB[J+1]*SINLAM[I]);
            }
        }
        V[K]=-2.0*V[K]*exp(-AB[K+1]*SINLAM[K]);
        V[K+1]=-V[K]*SINLAM[K]*AB[K];
    }
    NN=MAXX+1;
    for (I=NN; I <= 10; ++I)
    {
        V[I]=0;
    }
}

void COEFF::CV2(double AB[],double* G)
{
    int I,J,MAXX;
    double T;

    MAXX=min(9,N);
    *G=0;
    for (I=1; I <= N; ++I)
    {
        T=F[I];
        for (J=1; J <= MAXX; J = J + 2)
        {
            T=T-AB[J]*exp(-AB[J+1]*SINLAM[I]);
        }
        *G=*G+T*T;
    }
}

void COEFF::STEEP(double X[], int N, int M)
{
    int I,J,K;
    double Z[16+1],T[16+1],ZNT[16+1],SIGT[16+1],SIG[16+1],Y[16+1][16+1],XI[16+1],ZN[16+1];
    double ALPHA,D,G,G1,U,V,W,T1,T2,ETA;

    if (N == 0)
    {
        CV1(X,Z);
        CV2(X,&G);
        for (I=1; I <= M; ++I)
        {
            for (J=1; J <= M; ++J)
            {
                Y[I][J]=0.0;
                if (I == J) Y[I][J]=1.0;
            }
        }
        for (K=1; K <= N; ++K)
        {
            U=0;
            for (I=1; J <= M; ++J)
            {
                T1=0;
                for (J=1; J <= M; ++J) T1=T1+Z[J]*Y[I][J];
                SIG[I]=T1;
                U=U-SIG[I]*Z[I];
            }
            ETA=min(1.0,abs(2.0*G/U));
            for (I=1; I <= M; ++I) T[I]=X[I]-ETA*SIG[I];
            CV1(T,ZN);
            V=0.0;
            for (I=1; I <= M; ++I) V=V-ZN[I]*SIG[I];
            CV2(T,&G1);
            D=3.0*abs(G-G1)/ETA+U+V;
            W=sqrt(D*D-U*V);
            ALPHA=ETA*(1.0-(V+W-D)/(V-U+2.0*W));
            for (I=1; I <= M; ++I)
            {
                SIG[I]=-1.0*ALPHA*SIG[I];
                X[I]=X[I]+SIG[I];
                ZN[I]=Z[I];
            }
            CV1(X,Z);
            CV2(X,&G);
            T1=0;
            for (I=1; I <= M; ++I)
            {
                XI[I]=Z[I]-ZN[I];
                T1=T1+XI[I]*SIG[I];
            }
            for (I=1; I <= M; ++I) SIGT[I]=SIG[I]/T1;
            T1=0;
            for (I=1; I <= M; ++I)
            {
                T2=0;
                for (J=1; J <= M; ++J) T2=T2+XI[J]*Y[I][J];
                ZN[I]=T2;
                T1=T1+ZN[I]*XI[I];
            }
            for (I=1; I <= M; ++I)
            {
                T2=0;
                for (J=1; J <= M; ++J) T2=T2+XI[J]*Y[J][I];
                ZNT[I]=T2/T1;
            }
            for (I=1; I <= M; ++I)
            {
                for (J=1; J <= M; ++J) Y[I][J]=Y[I][J]+SIG[I]*SIGT[J]-ZN[I]*ZNT[J];
            }
        }
    }
}

void  COEFF::COEF(int* J, int* K)
{
    int I,NN,NN1;
    double AB[10+1];

    for(I=1; I <= *K; ++I)
    {
        SINLAM[I]=POSI[I]*POSI[I];
        F[I]=SCAT[I];
    }
    N=*K;
    NN = min(N,9);
    for (I=1; I <= NN; ++I)
    {
        AB[I]=I;
    }
    NN1=NN+1;
    for (I=NN1; I <= 10; ++I)
    {
        AB[I]=0.0;
    }
    STEEP(AB,3,10);
    if (AB[3] == 0.0) AB[3]=-1E-6;
    for (I=1; I <= 9; ++I)
    {
        AC[I][*J]=AB[I];
    }
}




DBWS::DBWS(void)
{
    diffractogram.file6 = &file6;

}

DBWS::~DBWS(void)
{
}

void DBWS::run(void)
{
    try
    {
        INPTR();
        // Canton et all code starts here !cp may 03 97
        //-----OPEN,if NECESSARY, FILE CONTAINING AMORPHOUS SCATTERING
        if (GLB_[20-1] != 0.0 || GLB_[20-1].codeword != 0.0 ) inpam();
        // Canton et all code stops here

        ITER();
        if ( MAXS > 0 )
        {
            MCYCLX =MCYCLE;
            MAXSX  = MAXS;
            MCYCLE = 1;
            MAXS   = 0;
            //-------LAST CALL TO ITER if MCYCLE = 1 & MAXS = 0
            ITER();
            MCYCLE = MCYCLX;
            MAXS   = MAXSX;
        }
        EXPUT();
    }
    catch (DBWSException& e)
    {
        cout << e.Show() << endl;
    }
}

int  DBWS::sign(int i)
{
    if (i < 0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

double DBWS::sign(double i)
{
    if (i < 0)
    {
        return -1.0;
    }
    else
    {
        return 1.0;
    }
}

void DBWS::GOTOER(void)
{
    file6 << "subroutine gotoer was called" << endl;
    throw DBWSException("subroutine gotoer was called");
}

void DBWS::SORT(int IPHASE)
{
    int I,J,K,M,IC, LL, IZ, LX, IOF, INP, INP1,IIPHAS;
    int L[IRS+1];
    int ITEMP[IRS+1];
    double TEMP[IRS+1];
    double R;

    IOF=0;
    if (IPHASE >= 2 && IPHASE != 1000)
    {
        for (IIPHAS=2; IIPHAS <= IPHASE; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].ICR;
    }
    if (IPHASE == 1000)
    {
        IC=0;
        for (IIPHAS=1; IIPHAS <= 99; ++IIPHAS) IC  = IC  + phases[IIPHAS].ICR;
    }
    else
    {
        IC=phases[IPHASE].ICR;
    }
    for (I=1; I <= IC; ++I) L[I]=I;
    M=IC-1;

    // TODO: possivel erro aqui!!!!
    //////////////////////////////
    for (LL=1; LL <= 5000; ++LL)
    {
        K=0;
        for (J=1; J <= 2; ++J)
        {
            for (I=J; I <= M; I = I+2)
            {
                INP=L[I]+IOF;
                INP1=L[I+1]+IOF;
                R=refs[INP1].REFS[2]-refs[INP].REFS[2];
                if (R < 0) goto L20; else goto L10;
L20:
                IZ=L[I];
                L[I]=L[I+1];
                L[I+1]=IZ;
                K=K+1;
L10:;
            }
        }
        if (K == 0) goto L5;
    }
    ////////////////////////////////


L5:
    for (I=1; I <= IC; ++I) ITEMP[I]=refs[I+IOF].IREFS;
    for (I=1; I <= IC; ++I)
    {
        LX=L[I];
        refs[I+IOF].IREFS=ITEMP[LX];
    }
    for (J=1; J <= 3; ++J)
    {
        for (I=1; I <= IC; ++I) TEMP[I]=refs[I+IOF].REFS[J];
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refs[I+IOF].REFS[J]=TEMP[LX];
        }
    }
    for (I=1; I <= IC; ++I) TEMP[I]=refs[I+IOF].FMGNTD;
    for (I=1; I <= IC; ++I)
    {
        LX=L[I];
        refs[I+IOF].FMGNTD=TEMP[LX];
    }
    if (NPROF == _TCHZ)
    {
        for (I=1; I <= IC; ++I) TEMP[I]=refs[I+IOF].HALFL;
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refs[I+IOF].HALFL=TEMP[LX];
        }
        for (I=1; I <= IC; ++I) TEMP[I]=refs[I+IOF].HALFG;
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refs[I+IOF].HALFG=TEMP[LX];
        }
        for (I=1; I <= IC; ++I) TEMP[I]=refs[I+IOF].GAM;
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refs[I+IOF].GAM=TEMP[LX];
        }
    }
}

void DBWS::ASSIGN_(void)
{
    int I,J,K,IN1,IN2,IRC,ICX,IRK,IRL,IRH,MIN,MAX,IPHAS;
    double PX, WT,YUN, RMIN, RMAX;
    int KRR[2+1][IDSZ+1];

    ICX=0;
    for (J=1; J <= NPHASE; ++J)
    {
        ICX = ICX + phases[J].ICR;
    }
    if (NPHASE > 1) SORT(1000);
    for (J=1; J <= 2; ++J)
    {
        for (K=1; K <= IDSZ; ++K) KRR[J][K]=0;
    }
    if(LST3 == 0) goto L7;
    for (I=1; I <= ICX; ++I)
    {
        IRL=( refs[I].IREFS % 256)-128;
        IRK=((refs[I].IREFS/256) % 256)-128;
        IRH=((refs[I].IREFS/(256*256)) % 256)-128;
        IRC=(refs[I].IREFS/(256*256*256)) % 8;
        IPHAS=refs[I].IREFS/(256*256*256*8);
        if ( (I-1 % 60) == 0) file6 << "NO.  CODE    H   K   L  PHASE  HW       POSN" << endl;
        file6 << setw(4) << I
            << setw(4) << IRC << "   "
            << setw(4) << IRH
            << setw(4) << IRK
            << setw(4) << IRL
            << setw(6) << IPHAS
            << setw(8) << setprecision(3) << refs[I].REFS[1]
            << setw(8) << setprecision(3) << refs[I].REFS[2] << endl;
    }
    L7:
    for(I=1; I <= ICX; ++I)
    {
        if (NPROF == _SplitPearsonVII)
        {
            RMIN=refs[I].REFS[2]-WDT_*refs[I].FWHM[1];
            RMAX=refs[I].REFS[2]+WDT_*refs[I].FWHM[2];
        }
        else
        {
            RMIN=refs[I].REFS[2]-WDT_*refs[I].REFS[1];
            RMAX=refs[I].REFS[2]+WDT_*refs[I].REFS[1];
        }
        MIN= static_cast<int>( (RMIN-diffractogram.THMIN)/diffractogram.STEP+1.5 );
        MAX= static_cast<int>( (RMAX-diffractogram.THMIN)/diffractogram.STEP+1.5 );
        MIN=max(MIN,1);
        MIN=min(MIN,diffractogram.NPTS);
        MAX=min(MAX,diffractogram.NPTS);
        for (K=MIN; K <= MAX; ++K)
        {
            KRR[2][K]=I;
            if(KRR[1][K] == 0)KRR[1][K]=I;
        }
    }
    for (J=1; J <= diffractogram.NEXCRG; ++J)
    {
        if (diffractogram.AHIGH[J] <= diffractogram.THMIN) goto L482;
        IN1= static_cast<int>( (diffractogram.ALOW[J]-diffractogram.THMIN)/diffractogram.STEP+1.5 );
        IN1=max(IN1,1);
        IN2= static_cast<int>( (diffractogram.AHIGH[J]-diffractogram.THMIN)/diffractogram.STEP+1.5 );
        IN2=min(IN2,diffractogram.NPTS);
        for(I=IN1; I <= IN2; ++I)
        {
            KRR[2][I]=1;
            KRR[1][I]=0;
        }
        if(IN2 == diffractogram.NPTS)goto L484;
        L482:;
    }
    L484:
    if(LST2 == 0)goto L530;
    file6 << "PATTERN FROM"
          << setw(8) << setprecision(4) << diffractogram.THMIN
          << " TO" << setw(8) << setprecision(4) << diffractogram.THMAX
          << " IN STEPS OF"  << setw(8) << setprecision(4) << diffractogram.STEP << "DEGREES" << endl;
    file6 << "POSN      I+B     B     I       100*W   K11  K21" << endl;
    for(J=1; J <= diffractogram.NPTS; ++J)
    {
        PX=diffractogram.THMIN + static_cast<double>(J-1) * diffractogram.STEP;
        YUN=diffractogram.Y[J]-diffractogram.BK[J];
        if(JOBTYP < 3) WT=1.0/diffractogram.VAR[J];
        file6 << setw(8) << setprecision(4) << PX
              << setw(7) << setprecision(0) << diffractogram.Y[J]
              << setw(7) << setprecision(0) << diffractogram.BK[J]
              << setw(7) << setprecision(0) << YUN
              << setw(9) << setprecision(4) << WT
              << setw(5) << KRR[1][J]
              << setw(5) << KRR[1][J] << endl;
    }
    L530:;
    for (I=1; I <= diffractogram.NPTS; ++I)
    {
        if(KRR[2][I]-KRR[1][I] > NOV)
        {
            file6 << "EXCESSIVE PEAK OVERLAP" << endl
                  << "     AT THE" << I << "TH STEP THERE ARE "
                  << KRR[2][I]-KRR[1][I] << " REFLECTIONS" << endl
                  << "     INCREASE THE VALUE OF *NOV* WHICH IS NOW " << NOV << endl;
            DBWSException("EXCESSIVE PEAK OVERLAP");
        }
    }
    for (I=1; I <= diffractogram.NPTS; ++I) diffractogram.KR[I]=KRR[1][I]+IRS *KRR[2][I];
    for (I=1; I <= diffractogram.NPTS; ++I) if(diffractogram.KR[I] != 0 && diffractogram.KR[I] != IRS )goto L603;
    DBWSException("NO REFLECTIONS FOUND");
    // test for detecting the asymmetry model required    !cp ap 16 97
    // new code included in line 2 of ICF
    // iasym = 0 (usual Rietveld asymmetry)
    // iasym = 1 (new asymmetry. Riello, Canton & Fagherazzi.PD 10,3,204-206,1997)
    L603:
    if (IASYM == 0)
    {
        // THE FOLLOWING TEST IS REALLY ONLY VALID FOR THE SINGLE PHASE CASE
        if(refs[ KRR[1][I] ].REFS[2] >= RLIM_ && phases[1].PAR[13].L != 0)
        {
            file6 << "ASYMMETRY PARAMETER USAGE INVALID" << endl;
            DBWSException("ASYMMETRY PARAMETER USAGE INVALID");
        }
    }
    else
    {
        if(refs[ KRR[1][I] ].REFS[2] >= (90.0-RLIM_) && phases[1].PAR[13].L != 0)
        {
            file6 << "ASYMMETRY PARAMETER USAGE INVALID" << endl;
            DBWSException("ASYMMETRY PARAMETER USAGE INVALID");
        }
    }
    for (I=1; I <= NPHASE; ++I)
    {
        if(phases[I].PAR[11] == 0.0 && phases[I].PAR[12] == 0.0 && phases[I].PAR[12].L != 0)
        {
            file6 << "PREFERRED ORIENTATION USAGE INVALID" << endl;
            DBWSException("PREFERRED ORIENTATION USAGE INVALID");
        }
    }
}

void DBWS::CHISQ(void)
{
    int I,IEXC,NPTHI,NPTLOW;
    double DEL,SDELW,D1NOBK,D2NOBK,S1NOBK,S2NOBK,SDELWP;

    diffractogram.S1_=0.0;
    S1NOBK=0.0;
    diffractogram.S2_=0.0;
    S2NOBK = 0.0;
    diffractogram.SS2_=0.0;
    diffractogram.S3_=0.0;
    diffractogram.SS4_=0.0;
    diffractogram.D1_=0.0;
    diffractogram.D2_=0.0;
    diffractogram.D4_=0.0;
    SDELW = 1.0E+25;
    SDELWP= 1.0E+25;
    D1NOBK = 0.0;
    D2NOBK = 0.0;
    for (I=1; I <= diffractogram.NPTS; ++I)
    {
        if((diffractogram.NBCKGD != 0 && diffractogram.KR[I] == 0) || diffractogram.KR[I] == IRS )goto L10;
        if(diffractogram.NEXCRG > 0)
        {
            for (IEXC=1; IEXC <= diffractogram.NEXCRG; ++IEXC)
            {
                NPTLOW = static_cast<int>( (diffractogram.ALOW[IEXC]-diffractogram.THMIN)/diffractogram.STEP + 1.5 );
                NPTHI  = static_cast<int>( (diffractogram.AHIGH[IEXC]-diffractogram.THMIN)/diffractogram.STEP + 1.5 );
                if (I >= NPTLOW && I <= NPTHI) goto L10;
            }
        }
        if (SDELW < 1.0e+24) SDELWP = SDELW;
        DEL=diffractogram.Y[I]-diffractogram.BK[I]-diffractogram.YC[I];
        SDELW = DEL;
        diffractogram.S1_ = diffractogram.S1_ + abs(DEL);
        S1NOBK = S1NOBK + DEL*(diffractogram.Y[I]-diffractogram.BK[I])/diffractogram.Y[I];                              //numerador para calcular Rp(-bcg)
        diffractogram.S2_=diffractogram.S2_+DEL*DEL/diffractogram.VAR[I];
        S2NOBK=S2NOBK+DEL*DEL*(diffractogram.Y[I]-diffractogram.BK[I])*(diffractogram.Y[I]-diffractogram.BK[I])/(diffractogram.Y[I]*diffractogram.Y[I])/diffractogram.VAR[I];   //numerador para calcular Rwp(-bcg)
        diffractogram.SS2_ = diffractogram.SS2_ + DEL*DEL;
        if (SDELWP < 1.0e+24) diffractogram.SS4_=diffractogram.SS4_+(SDELW-SDELWP)*(SDELW-SDELWP);
        L10:;
    }
    for (I=1; I <= diffractogram.NPTS; ++I)
    {
        if((diffractogram.NBCKGD != 0 && diffractogram.KR[I] == 0) || diffractogram.KR[I] == IRS )goto L15;
        if(diffractogram.NEXCRG > 0)
        {
            for (IEXC=1; IEXC <= diffractogram.NEXCRG; ++IEXC)
            {
                NPTLOW = static_cast<int>( (diffractogram.ALOW[IEXC]-diffractogram.THMIN)/diffractogram.STEP + 1.5 );
                NPTHI  = static_cast<int>( (diffractogram.AHIGH[IEXC]-diffractogram.THMIN)/diffractogram.STEP + 1.5 );
                if (I >= NPTLOW && I <= NPTHI) goto L15;
            }
        }
        diffractogram.D1_=diffractogram.D1_+diffractogram.Y[I];
        D1NOBK = D1NOBK+(diffractogram.Y[I]-diffractogram.BK[I]);                                //denominador para calcular Rp(-bcg)
        diffractogram.D2_=diffractogram.D2_+diffractogram.Y[I]*diffractogram.Y[I]/diffractogram.VAR[I];
        D2NOBK = D2NOBK + (diffractogram.Y[I]-diffractogram.BK[I])*(diffractogram.Y[I]-diffractogram.BK[I])/diffractogram.VAR[I];          //denominador para calcular Rwp(-bcg)
        diffractogram.S3_=diffractogram.S3_+diffractogram.BK[I]+diffractogram.YC[I];
        L15:;
    }
    diffractogram.R1_=0.0;
    diffractogram.R2_=100.0*diffractogram.S1_/diffractogram.D1_;
    diffractogram.R2NOBK_ = 100.0*S1NOBK/D1NOBK;                                //Rp without background 03nov00 (no writable)
    diffractogram.R3_=100.0*sqrt(diffractogram.S2_/diffractogram.D2_);
    diffractogram.R3NOBK_ = 100.0*sqrt(S2NOBK/D2NOBK);                           //Rwp without background 03nov00 (no writable)
}

void DBWS::COMPTON(int K, double STH, double* CISK)
{
    //-----THIS SUBROUTINE IS USED TO COMPUTE:
    //     CISK = COMPTON INTENSITY SCATTERED BY THE K-TH CRYSTALLINE PHASE
    //            AT THE IPM-TH POINT OF THE X-RAY SPECTRUM.
    //     THE METHOD FOLLOWED IS THAT REPORTED IN :
    //*****A NEW ANALYTIC APPROXIMATION TO ATOMIC INCOHERENT X-RAY SCATTERING
    //     INTENSITIES.
    //     BY VEDENE H. SMITH JR. AJIT J. THAKKAR AND DOUGLAS C. CHAPMAN
    //     ACTA CRYST. (1975). A31, 391-392.
    //-----COMPTON SCATTERING OF THE AMORPHOUS PHASE IS CONTAINED INTO AM(S)
    //     MOREOVER CISK IS MULTIPLIED BY THE FOLLOWING GEOMETRICAL FACTOR:
    //     SEE:
    //     RULAND, ACTA CRYST. 1961, 14, 1180.
    //     ASS1 = CORRECTION FOR ABSORPTION OF COMPTON SCATTERING BY THE SAMPLE
    //     BDF  = BREIT-DIRAC RECOIL FACTOR
    //     ASS2 = CORRECTION FOR ABSORPTION BY AIR CONTAINED BETWEEN THE SAMPLE
    //            AND THE COUNTER
    //     ASS3 = CORRECTION FOR THE CUTTING UP OF COMPTON SCATTERING BY THE
    //            MONOCHROMATOR ( ALSO COMPTON SCATTERING IS DifFRACTING )

    const double HMC = {0.024263935};			// HMC = H/MC

    int I,J,NK,IOF, ICX,IRL;
    double S, S2, S4,BDF,ASS1, ASS2, ASS3, FDEN, FNUM;
    double CSP[2 +1];

    //-----NK = ATOMS IN THE ASSIMETRIC UNIT OF K-TH PHASE
    NK = phases[K].AtomCount;
    //-----IRL = NUMBER OF EQUIVALENT POSITION - ALSO THE IDENTITY
    //           AND if IT IS PRESENT THE SIMMETRY CENTRE
    IRL = phases[K].MLTPHS;

    //-----CALCULATE IOF = ALL ATOMS OF THE K-1 PHASES
    IOF = 0;
    if (K > 1)
    {
        for (I = 2; I <= K; ++I) IOF = IOF + phases[I-1].AtomCount;
    }
    for (ICX = 1; ICX <= 2; ++ICX)
    {
        S  =  STH / diffractogram.LAMDA[ICX];
        S2 =  S * S;
        S4 = S2 * S2;
        CSP[ICX] = 0.0;
        for (I = 1; I <= NK; ++I)
        {
            J = PTC_[I+IOF];
            FNUM    =   1.0 + CC_[1][J] * S2 + CC_[2][J] * S4;
            FDEN    = pow( ( 1.0 + CC_[3][J] * S2 + CC_[4][J] * S4 ) ,2);
            CSP[ICX]=CSP[ICX] + ZEFF_[J] * XL_[I+IOF][5] * (1.0-(FNUM/FDEN));
        }
        CSP[ICX] = CSP[ICX] * static_cast<double>(IRL);
        //-----UPDATE S
        S  = 2.0 * S;
        S2 =   S * S;
        //-----COMPUTE ASS1
        ASS1 =   1.0 + 0.75 * HMC * diffractogram.LAMDA[ICX] * S2;
        //-----COMPUTE BDF
        BDF  = pow(( 1.0 + 0.50 * HMC * diffractogram.LAMDA[ICX] * S2 ) , 2);
        //-----COMPUTE ASS2 = EXP (-DELTA(MU)*D)
        //      WHERE:
        //      MU       = ABSORPTION COEFFICIENT OF AIR
        //      DELTA(MU)= VARIATION OF MU WITH LAMDA, WAS CALCULATED BY RIELLO
        //                 USING REGRESSION ANALYSIS ON ESTIMATED MU BY ASSUMING
        //                 AIR COMPOSITION 20% O2 AND 80% N2 AT 300�K AND 1 ATM.:
        //                 MU = (E ** 3.40089)*(1.0E-04)*( LAMDA ** 2.79287)
        //      D = 17.3 CM. DISTANCE SPECIMENT-COUNTER OR RADIUS OF THE CAMERA
        ASS2 = exp(-1.5*HMC*17.3*(29.99E-04)* pow(diffractogram.LAMDA[ICX],3.79287)*S2);
        //-----COMPUTE ASS3
        //     FORMULA MUST BE CHANGED FOR OTHER RADIATION AND MONOCHROMATOR
        //-----  ASS3 IS A LORENTZIAN FUNCTION
        ASS3=1/(1+GLB_[18-1]* pow(S,GLB_[19-1]));
        //----- ASS3 IS A GAUSSIAN FUNCTION
        //       ASS3=EXP(-GLB_[18-1]*S2)
        //     ASS3=0.38*EXP(-4.0*S2)+0.62*EXP(-0.15*S2)
        CSP[ICX] = CSP[ICX] * ASS2 * ASS3  / ( ASS1 * BDF );
    }
    //     COMPUTE CISK = COMPTON INTENSITY SCATTERED BY THE K-TH PHASE
    *CISK = diffractogram.RATIO[1] * CSP[1] + diffractogram.RATIO[2] * CSP[2];
}

void DBWS::DISORDER(int K, double STH, int IDERIV, double* SDK, double* DYC, int FONDO, double DERISO[])
{
    //-----THIS SUBROUTINE COMPUTES SDK = THE SCATTERING DUE TO THE THERMAL
    //     OR LATTICE DISORDER IN THE SAMPLE FOR K-TH PHASE AT THE IPM POINT,
    //     (THEN ALSO IN STEP POINTS THAT DO NOT CONTAIN BRAGG REFLECTIONS),
    //     P.RIELLO, G. FAGHERAZZI, D. CLEMENTE AND P.CANTON IN
    //     J. APPL. CRYST. (1995) 28,115-120

    int I,II, NI, NK, LK,IOF,ICX, IRL;
    double FI;
    double DER[2+1];
    double DIS[2+1];
    double STHL2[2+1];
    double DERDIS[NATS+1][2+1];

    //-----COMPUTE IOF = ATOM NUMBER IN THE K-1 PHASE
    IOF = 0;
    if (K > 1)
    {
        for (I = 2; I <= K; ++I) IOF = IOF + phases[I-1].AtomCount;
    }
    for (I=1; I <= NATS; ++I)
    {
        DERDIS[I][1]=0.0;
        DERDIS[I][2]=0.0;
    }
    //-----IRL = NUMBER OF EQUIVALENT POSITIONS- ALSO THE IDENTITY POSITION
    //           AND THE SIMMETRY CENTRE if PRESENT
    IRL = phases[K].MLTPHS;
    //-----NK = NUMBER OF ATOMS IN THE K-TH PHASE
    NK  = phases[K].AtomCount;
    for (ICX = 1; ICX <= 2; ++ICX)
    {
        STHL2[ICX] = pow(( STH / diffractogram.LAMDA[ICX] ) , 2);
        if(FONDO == 1)
        {
            //-----COMPUTE FI2 = SUM OF SQUARE SCATTERING FACTORS DUE TO ALL ATOMS
            //                   IN THE CELL AT LAMDA(ICX) CORRECTED FOR THE ISOTROPIC
            //                    THERMAL FACTORS.
            FI2_[ICX] = 0.0;
            for (I = 1; I <= NK; ++I)
            {
                NI = PTR_[I+IOF];
                FI = 0.0;
                coeff.AC[10][NI] = 0.0;
                for (II = 1; II <= 9; II = II + 2) FI = FI + coeff.AC[II][NI]*exp(-coeff.AC[II+1][NI]*STHL2[ICX]);
                FI = FI + coeff.DFP[NI];
                FI2_[ICX] += XL_[I+IOF][5] *
                    static_cast<double>(IRL)*( 1.0 - exp(-XL_[I+IOF][4]*2.0*STHL2[ICX]) )* (pow(FI,2)+pow(coeff.DFPP[NI],2));
                //-----NEXT LINES EVALUATE THE DERIVATES
                if(XL_[I+IOF][4].L != 0 && IDERIV == 2)
                {
                    DERDIS[I+IOF][ICX] = XL_[I+IOF][5] * 2.0 * STHL2[ICX] *
                        exp(-XL_[I+IOF][4]*2.0*STHL2[ICX])*static_cast<double>(IRL)*(pow(FI,2)+pow(coeff.DFPP[NI],2));
                }
                else
                {
                    DERDIS[I+IOF][ICX]=0.0;
                }
            }
            DIS[ICX] = FI2_[ICX];
        }
        if(FONDO == 2)
        {
            //-----COMPUTE FI2 = SUM OF SQUARE SCATTERING FACTORS DUE TO ALL ATOMS
            //                   IN THE CELL AT LAMDA(ICX) CORRECTED FOR THE OVERALL
            //                    THERMAL FACTORS.
            FI2_[ICX] = 0.0;
            for (I = 1; I <= NK; ++I)
            {
                NI = PTR_[I+IOF];
                FI = 0.0;
                coeff.AC[10][NI] = 0.0;
                for (II = 1; II <= 9; II = II + 2) FI += coeff.AC[II][NI]*exp(-coeff.AC[II+1][NI]*STHL2[ICX]);
                FI += coeff.DFP[NI];
                FI2_[ICX] += XL_[I+IOF][5] * (pow(FI,2)+pow(coeff.DFPP[NI],2));
            }
            FI2_[ICX] += static_cast<double>(IRL);
            DIS[ICX] = ( 1.0 - exp(-phases[K].PAR[1]*2.0*STHL2[ICX]) ) * FI2_[ICX];
            LK =phases[K].PAR[1].L;
            if(LK != 0 && IDERIV == 2)
            {
                DER[ICX] = 2.0 * STHL2[ICX] * FI2_[ICX] * exp(-phases[K].PAR[1]*2.0*STHL2[ICX]);
            }
            else
            {
                DER[ICX] = 0.0;
            }
        }
    }
    //-----COMPUTE SDK = SCATTERING DISORDER DUE TO THE K-TH PHASE
    *SDK = diffractogram.RATIO[1] * DIS[1] + diffractogram.RATIO[2] * DIS[2];
    //-----COMPUTE  DERIVATIVE OF YC RESPECT TO ISOTROPIC THERMAL FACTORS
    //                   XL(I+IOF,4) IN THE K-TH PHASE
    if (FONDO == 1)
    {
        for (I=1; I <= NK; ++I)
        {
            if(XL_[I+IOF][4].L != 0 && IDERIV == 2)
            {
                DERISO[I+IOF]=diffractogram.RATIO[1] * DERDIS[I+IOF][1] + diffractogram.RATIO[2] * DERDIS[I+IOF][2];
            }
            else
            {
                DERISO[I+IOF]=0.0;
            }
        }
    }
    //-----COMPUTE  DERIVATIVE OF YC RESPECT TO OVERALL THERMAL FACTOR
    //                   PAR[K][2] IN THE K-TH PHASE
    if(FONDO == 2)
    {
        LK =phases[K].PAR[1].L;				// TODO: Coloquei talez esteja incorreto, verificar depois de retirar os gotos
        if(LK != 0 && IDERIV == 2)
        {
            *DYC = diffractogram.RATIO[1] * DER[1] + diffractogram.RATIO[2] * DER[2];
        }
        else
        {
            *DYC = 0.0;
        }
    }
}

void DBWS::SMTRY2(int* IPHASE)
{
    //                             THIS SR GENERATES THE FULL SET OF
    //                           EQUIVALENT HKLS FROM ANY MEMBER OF THE SET.
    //                           IT ALSO DETERMINES THE PHASE SHIFTS RELATIVE
    //                           TO THE INPUT HKL PHASE.

    const int KD[6+1] = { 0,   0,6,3,9,4,8};
    const int KE[3+1] = { 0,   0,4,8};
    const int KF[3+1] = { 0,   0,8,4};

    int I,J,L,M,N,LA,JD = 0,MS,NCO,NCX,K1JD,IXIT;
    double CI;
    int LD[3+1][48+1];
    int NCH[48+1];
    int NCK[48+1];
    int NCL[48+1];
    int IPH[48+1];

    if ( phases[*IPHASE].SYMB.NC > 0 ) goto L14;
    if ( IHKL[1][1] < 0 )
    {
        goto L11;
    }
    else if ( IHKL[1][1] == 0 )
    {
        goto L9;
    }
    else
    {
        goto L14;
    }
L9:
    if ( IHKL[2][1] < 0 )
    {
        goto L12;
    }
    else if ( IHKL[2][1] == 0 )
    {
        goto L10;
    }
    else
    {
        goto L14;
    }
L10:
    if ( IHKL[3][1] < 0 )
    {
        goto L13;
    }
    else
    {
        goto L14;
    }
L11:
    IHKL[1][1] = -IHKL[1][1];
L12:
    IHKL[2][1] = -IHKL[2][1];
L13:
    IHKL[3][1] = -IHKL[3][1];
L14:
    IPH[1] = IHKL[3][1]+512*(IHKL[2][1]+512*IHKL[1][1]);
    CI = 1.0;
    ICHKL[*IPHASE] = 1;
    AZ[1] = 0.0;
    LD[1][1] = 0;
    LD[2][1] = 0;
    LD[3][1] = 0;
    NCH[1] = 0;
    NCK[1] = 0;
    NCL[1] = 1;
    if ( N1HKL[*IPHASE] == 0 ) goto L1000;
    L = 2;
    for (I=1; I <= N1HKL[*IPHASE]; ++I)
    {
        CI = CI*2.0;
        if ( NC1[I][1][*IPHASE] > 0 ) CI = CI*1.5;
        for (J=1; J <= ICHKL[*IPHASE]; ++J)
        {
            NCO = NC1[I][1][*IPHASE]+1;
            switch (NCO) {
            case 1:
                goto L20;
                break;
            case 2:
                goto L700;
                break;
            case 3:
                goto L800;
                break;
            case 4:
                goto L15;
                break;
            }
            GOTOER();
L15:
            if( ((IHKL[1][1]-IHKL[2][1]-IHKL[3][1]) % 3) != 0 ) goto L1002;
            goto L900;
L20:
            NCH[L] =  NCH[J] ^ NC1[I][2][*IPHASE];
            if ( (NC1[I][2][*IPHASE] % 2) != 0 ) NCH[L] = NCK[J] ^ NC1[I][2][*IPHASE];
            M = 1+(NCH[L] % 2);
            MS = 1-2*((NCH[L]/2) % 2);
            N = abs(NCL[J]);
            IHKL[M][L] = IHKL[1][N]*MS;
            M = 1+(NC1[I][2][*IPHASE] % 2);
            MS = 1-2*((NC1[I][2][*IPHASE]/2) % 2);
            NCX = NC1[I][3][*IPHASE]+1;
            LD[1][L] = KD[NCX]+LD[M][J]*MS;
            NCK[L] =  NCK[J] ^ NC1[I][4][*IPHASE];
            if ( (NC1[I][4][*IPHASE] % 2) != 0 ) NCK[L]= NCH[J] ^ NC1[I][4][*IPHASE];
            M = 2-(NCK[L] % 2);
            MS = 1-2*((NCK[L]/2) % 2);
            IHKL[M][L] = IHKL[2][N]*MS;
            M = 2-(NC1[I][4][*IPHASE] % 2);
            MS = 1-2*((NC1[I][4][*IPHASE]/2) % 2);
            NCX = NC1[I][5][*IPHASE]+1;
            LD[2][L] = KD[NCX]+LD[M][J]*MS;
            MS = (1-2*NC1[I][6][*IPHASE]) *  sign(NCL[J]);
            NCL[L] = MS*NCL[N];
            IHKL[3][L] = IHKL[3][N]*MS;
            //L70:
            NCO = NC1[I][7][*IPHASE]+1;
            if ( NCO > 5 ) NCO=6;
            MS = 1-2*NC1[I][6][*IPHASE];
            if ( N > 1 && NCK[J] >= 4 && M == 1 ) MS=-MS;
            LD[3][L] = KD[NCO]+LD[3][J]*MS;
            //L80:
            IXIT = 0;
L81:
            IPH[L] = IHKL[3][L]+512*(IHKL[2][L]+512*IHKL[1][L]);
            //L87:
            LA = LD[1][L]*IHKL[1][1]+LD[2][L]*IHKL[2][1]+LD[3][L]*IHKL[3][1];
            AZ[L] = static_cast<double>(LA)/12.0;
            if ( phases[*IPHASE].SYMB.NC != 0 ) goto L89;
            if ( phases[*IPHASE].SYMB.NSPGRP < 12 ) goto L89;
            if ( NC1[I][1][*IPHASE] > 0 ) goto L89;
            if ( (NC1[I][2][*IPHASE] % 2) == 1 ) goto L89;
            if ( IPH[L] <= 0 ) goto L680;
L89:
            for (M=2; M <= L; ++M)
            {
                if ( IPH[M-1] == IPH[L] ) goto L690;
                if ( phases[*IPHASE].SYMB.NC != 0 ) goto L90;
                if ( IPH[M-1] == -IPH[L] ) goto L690;
L90:;
            }
            L = L+1;
L680:
            if ( IXIT-1 < 0 )
            {
                goto L890;
            }
            else if ( IXIT-1 == 0 )
            {
                goto L750;
            }
            else
            {
                goto L850;
            }
L690:
            if ( fmod(abs(AZ[L]-AZ[M-1]),1.0) != 0.0 ) goto L1002;
            goto L680;
L700:
            IXIT = 2;
            JD = J;
L750:
            IXIT = IXIT-1;
            LD[1][L] = LD[1][JD];
            LD[2][L] = LD[2][JD];
            LD[3][L] = LD[3][JD];
            NCH[L] = 0;
            NCK[L] = 0;
            NCL[L] = L;
            K1JD = IHKL[1][JD];
            IHKL[1][L] = IHKL[2][JD];
            IHKL[2][L] = IHKL[3][JD];
            IHKL[3][L] = K1JD;
            JD = L;
            goto L81;
L800:
            IXIT = 4;
            JD = J;
            NCO = NC1[I][6][*IPHASE];
L850:
            IXIT = IXIT-2;
            NCH[L] = 0;
            NCK[L] = 4;
            NCL[L] = L;
            LD[1][L] = 0;
            LD[2][L] = 0;
            IHKL[1][L] = -IHKL[1][JD]-IHKL[2][JD];
            IHKL[2][L] = IHKL[1][JD];
            IHKL[3][L] = IHKL[3][JD];
            LD[3][L] = KF[NCO]+LD[3][JD];
            if ( (NCH[J] % 2) != 0 ) LD[3][L]=KE[NCO]+LD[3][JD];
            JD = L;
            goto L81;
L890:;
        }
L900:
        ICHKL[*IPHASE] = L-1;
    }
L1000:
    IER = 0;
    if ( phases[*IPHASE].SYMB.NC != 0 ) goto L1001;
    for (M=2; M <= ICHKL[*IPHASE]; ++M)
    {
        if ( IPH[M] > 0 ) goto L1003;
        IHKL[1][M] = -IHKL[1][M];
        IHKL[2][M] = -IHKL[2][M];
        IHKL[3][M] = -IHKL[3][M];
L1003:;
    }
L1001:
    return;
L1002:
    IER = 1;
    goto L1001;
}

void DBWS::CALCUL(int NN)
{
    double DPRECORX,X, B1, B2, BB, AH, BH, CH, DH, EH = 0.0, AV, PI, CV, DV, BV,TR, SS, TT, SR = 0.0,
        YY,DA1, DA3, DA4, DA5, AH2 = 0.0, BH2 = 0.0,BNI, ARG,PAK, TAV, TLR, SAI,
        SBI, ARG2, FFX, FNN, DER, SRD,TLG, TLL, COSA, SINA,
        SSNN, SNXI, SUMA, SUMB, DHDHG, DHDHL,PAKNN,SINTH, COSTH,TANTH,PREXP = 0.0,
        SNEXI, SLABDA,TANTHE, EXPARG, PRECOR, PREXPX, DPRECOR = 0.0, PRECORX,ISITH = 0.0;
    int I, J, K, N,II, IJ,JJ,NX,NM,IR, IV, NI,KL,IOF, IRL, ICX,KKL,ICENT,
        IIPHAS;
    bool VERT,PAC;

    int HNN[3+1];
    double XI[14+1];
    double DERIV[MSZ+1];
    double SNX[NATS+1];
    double SA[NATS+1];
    double SB[NATS+1];
    double TEMP[NATS+1];
    double SNEX[NATS+1];
    double T[3+1];
    double SUMAX[NATS+1][9+1];
    double SUMBX[NATS+1][9+1];
    double H[3+1];
    double AL[3+1][3+1];
    double SM[3+1][3+1];

    PI = 90.0/atan(1.0);              // 360./3.14159265359
    IPH_ = refs[NN].IREFS/(256*256*256*8);  // 256/256/256/8
    IOF = 0;
    if ( IPH_ > 1 )
    {
        for (IIPHAS=2; IIPHAS <= IPH_; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
    }
    IRL = phases[IPH_].MLTPHS;
    N = phases[IPH_].AtomCount;
    ICENT = phases[IPH_].ICNTPHS;
    //-----ZEROIZE THE DERIVATIVES OF THIS REFLECTION W.R.T. TO PARAMETERS
    NX=0;
    for (IIPHAS=1; IIPHAS <= NPHASE; ++IIPHAS) NX = NX+phases[IIPHAS].AtomCount;
    for (I=1; I <= MSZ; ++I) DERIV[I]=0.0;
    for (I=1; I <= NX; ++I)
    {
        for (J=1; J <= 9; ++J)
        {
            SUMBX[I][J] = 0.0;
            SUMAX[I][J] = 0.0;
        }
    }
    CV=0.0;
    DV=0.0;
    AV=0.0;
    BV=0.0;
    PAC = phases[IPH_].PAR[11] != 0.0  ||  phases[IPH_].PAR[11].L != 0  ||
          phases[IPH_].PAR[12] != 0.0  ||  phases[IPH_].PAR[12].L != 0;
    for (I=1; I <= 3; ++I) AL[I][I]=phases[IPH_].PAR[I+5-1];
    AL[3][2]=phases[IPH_].PAR[8];
    AL[2][3]=AL[3][2];
    AL[3][1]=phases[IPH_].PAR[9];
    AL[1][3]=AL[3][1];
    AL[2][1]=phases[IPH_].PAR[10];
    AL[1][2]=AL[2][1];
    AH=phases[IPH_].PAR[2];
    BH=phases[IPH_].PAR[3];
    CH=phases[IPH_].PAR[4];
    DH=phases[IPH_].PAR[19];
    if (NPROF == _pseudoVoigt) EH=phases[IPH_].PAR[20];
    if (NPROF == _TCHZ)
    {
        AH2=phases[IPH_].PAR[14];
        BH2=phases[IPH_].PAR[15];
    }
    B1=phases[IPH_].PAR[11];
    B2=phases[IPH_].PAR[12];
    NM=(NN % NOV)+1;
    ICX=(refs[NN].IREFS/(256*256*256)) % 8;
    SLABDA=diffractogram.LAMDA[ICX]*diffractogram.LAMDA[ICX]/4.0;
    N=phases[IPH_].AtomCount;
    refs[NN].FMGNTD=0.0;
    PAKNN=0.0;
    PRECOR=1.0;
    TR=0.0;
    HNN[3]=(refs[NN].IREFS % 256)-128;
    HNN[2]=((refs[NN].IREFS/256) % 256)-128;
    HNN[1]=((refs[NN].IREFS/(256*256)) % 256)-128;

    for (I=1; I <= 3; ++I) H[I]=HNN[I];
    //-----CALCULATION OF TEMP.FACTOR,POSITION AND FWHM
    SS=0.0;
    for (I=1; I <= 3; ++I)
    {
        IHKL[I][1] = HNN[I];
        for (J=I; J <= 3; ++J) SS = HNN[I]*AL[I][J]*HNN[J]+SS;                         // GSAS DH2
    }
    if ( PAC )
    {
        TT = 0.0;
        for (I=1; I <= 3; ++I)
        {
            for (J=I; J <= 3; ++J) TT = TT + phases[IPH_].PREF[I]*AL[I][J] * phases[IPH_].PREF[J];             // GSAS DP2
        }
        SMTRY2(&IPH_);
        PRECOR = 0.0;
        PREXP = 0.0;
        DPRECOR = 0.0;
        for (IJ=1; IJ <= ICHKL[IPH_]; ++IJ)
        {
            PAK = 0.0;
            for (I=1; I <= 3; ++I)
            {
                for (J=I; J <= 3; ++J) PAK = phases[IPH_].PREF[I]*AL[I][J]*static_cast<double>(IHKL[J][IJ])+PAK;     // GSAS CA
            }
            PAK = PAK*PAK/(TT*SS);
            PAKNN = pow((PI/2.0),2);
            if (PAK != 0) PAKNN=pow(atan(sqrt(abs((1.0-PAK)/PAK))),2);
            if (IPREF == 0)
            {
                PREXPX = exp(B1*PAKNN);
                PRECORX = B2+(1.0-B2)*PREXPX;
                DPRECORX = PAKNN*(1.0-B2)*PREXPX/PRECORX;
            }
            else
            {
                PREXPX = B1*B1*PAK+(1.0-PAK)/B1;
                PRECORX = 1.0/pow(PREXPX,1.5);
                DPRECORX = -1.5*(2.0*B1*PAK-(1.0-PAK)/B1/B1)/ (pow(PREXPX,2.5)*PRECORX);
            }
            PREXP = PREXP+PREXPX;
            PRECOR = PRECOR+PRECORX;
            DPRECOR = DPRECOR+DPRECORX;
        }
        PREXP = PREXP/static_cast<double>(ICHKL[IPH_]);
        PRECOR = PRECOR/static_cast<double>(ICHKL[IPH_]);
        DPRECOR = DPRECOR/static_cast<double>(ICHKL[IPH_]);
    }
    //L13:
    //SINTL[NM] = sqrt(SS);
    SSNN = 0.25*SS;
    TAV = exp(-2.0*phases[IPH_].PAR[1]*SSNN)*PRECOR;
    SINTH = SLABDA*SS;
    COSTH = 1.0-SINTH;
    TANTH = sqrt(SINTH/COSTH);
    TANN_[NM] = TANTH;
    //   Correction of microabsorption
    if (IABSR == 1)
    {
        ISITH=(1.0)/(sqrt(SINTH));
        SR = GLB_[12-1]*(1.0-GLB_[8-1]*exp(-GLB_[9-1])+GLB_[8-1]*exp(-GLB_[9-1]/sqrt(SINTH)))+(1.0-GLB_[12-1])*(1+GLB_[13-1]*(asin(sqrt(SINTH)))-1.5707963268);
    }
    else if (IABSR == 2)
    {
        ISITH = sqrt(SINTH);
        SR = 1.0-GLB_[13-1]*(asin(ISITH)-1.5707963268);
    }
    else if (IABSR == 3)
    {
        ISITH=(1.0)/(sqrt(SINTH));
        SR = 1.0-GLB_[8-1]*exp(-GLB_[9-1])+GLB_[8-1]*exp(-GLB_[9-1]*ISITH);
    }
    else if (IABSR == 4)
    {
        ISITH=(1.0)/(sqrt(SINTH));
        SR = 1.0-GLB_[8-1]*GLB_[9-1]*(1.0-GLB_[9-1])-ISITH*GLB_[8-1]*GLB_[9-1]*(1.0-GLB_[9-1]*ISITH);
    }

    refs[NN].REFS[2]=atan(TANTH)*PI;
    refs[NN].HALFG=(AH*TANTH*TANTH+BH*TANTH+CH+DH/COSTH+EH/(TANTH*TANTH));
    if (refs[NN].HALFG > 0.0)
    {
        refs[NN].HALFG = sqrt(refs[NN].HALFG);
    }
    else
    {
        cout << "1 = " << HNN[1] << " " << HNN[2] << " " << HNN[3]  << endl;

        file6 << "   SQUARE OF FWHM NEGATIVE AT TWO-THETA=" << setw(8) << setprecision(3) << refs[NN].REFS[2] << " FOR PHASE NO. " << setw(4) << IPH_ << endl;
        cout << "   SQUARE OF FWHM NEGATIVE AT TWO-THETA=" << setw(8) << setprecision(3) << refs[NN].REFS[2] << " FOR PHASE NO. " << setw(4) << IPH_ << endl;
        DBWSException("SQUARE OF FWHM IS NEGATIVE");
    }
    if (NPROF == _TCHZ)
    {
        refs[NN].HALFL = AH2*TANTH+BH2*sqrt(1.0+TANTH*TANTH);
        BB = pow((pow(refs[NN].HALFG,5.0)+2.69269*pow(refs[NN].HALFG,4.0)*refs[NN].HALFL+2.42843*pow(refs[NN].HALFG,3.0)*pow(refs[NN].HALFL,2.0)+4.47163*pow(refs[NN].HALFG,2.0)*pow(refs[NN].HALFL,3.0) +0.07842*refs[NN].HALFG*pow(refs[NN].HALFL,4.0) + pow(refs[NN].HALFL,5.0)),0.2);
        TLR = refs[NN].HALFL/BB;
        refs[NN].GAM = 1.36603*TLR-0.47719*TLR*TLR+0.11116*pow(TLR,3.0);
    }
    else
    {
        BB = refs[NN].HALFG;
    }
    TL_=BB;
    refs[NN].REFS[1]=BB;
    BB=BB*BB;
    //-----VERT=.TRUE. if ASYMMETRY CORRECTION IS TO BE CALCULATED
    VERT=false;
    if(IASYM == 0)
    {
        if(refs[NN].REFS[2] <= RLIM_ && NPROF != _SplitPearsonVII) VERT=true;
    }
    else
    {
        if (abs(refs[NN].REFS[2]-90.0) >= RLIM_) VERT=true;
    }

    //-----CALCULATION OF COS(H.X),SIN(H.X) AND TEMP. FACTOR FOR EACH ATOM
    for (I=1; I <= N; ++I)
    {
        SNXI=0.0;
        SAI=0.0;
        SBI=0.0;
        for (J=1; J <= 11; ++J) XI[J]=XL_[I+IOF][J];
        for (IR=1; IR <= IRL; ++IR)
        {
            for (J=1; J <= 3; ++J)
            {
                IV=phases[IPH_].IVEC_[IR] /32768 / static_cast<int>(pow(32,3-J));
                IV=(IV % 32);
                SM[J][1]=IV/9-1;
                SM[J][2]=((IV/3) % 3)-1;
                SM[J][3]=(IV % 3)-1;
                T[J] = static_cast<double>( ((phases[IPH_].IVEC_[IR]/
                    static_cast<int>(pow(32,3-J))
                    ) % 32) - 16) / 12.0;
            }
            X=0.0;
            for (II=1; II <= 3; ++II)
            {
                YY=0.0;
                X=T[II]*HNN[II]+X;
                for (J=1; J <= 3; ++J) YY=SM[J][II]*HNN[J]+YY;
                H[II]=YY;
            }
            TR=X;
            ARG=TR;
            for (J=1; J <= 3; ++J)
            {
                ARG=H[J]*XI[J]+ARG;
            }
            ARG=6.28318530718*ARG;
            ARG2=H[1]*H[1]*XI[6]+H[2]*H[2]*XI[7]+H[3]*H[3]*XI[8]+2.0*H[1]*H[2]*XI[9]+2.0*H[1]*H[3]*XI[10]+2.0*H[2]*H[3]*XI[11];
            EXPARG=exp(-ARG2);
            COSA=cos(ARG)*EXPARG;
            SINA=sin(ARG)*EXPARG;
            SAI=SAI+COSA;
            if(ICENT == 1)SBI=SINA+SBI;
            for (JJ=1; JJ <= 3; ++JJ)
            {
                SUMAX[I][JJ]=SUMAX[I][JJ]+H[JJ]*SINA;
                if(ICENT == 1)
                {
                    SUMBX[I][JJ]=SUMBX[I][JJ]+H[JJ]*COSA;
                }
            }
            SUMAX[I][4]=SUMAX[I][4]+H[1]*H[1]*COSA;
            SUMAX[I][5]=SUMAX[I][5]+H[2]*H[2]*COSA;
            SUMAX[I][6]=SUMAX[I][6]+H[3]*H[3]*COSA;
            SUMAX[I][7]=SUMAX[I][7]+H[1]*H[2]*COSA;
            SUMAX[I][8]=SUMAX[I][8]+H[1]*H[3]*COSA;
            SUMAX[I][9]=SUMAX[I][9]+H[2]*H[3]*COSA;
            if(ICENT == 1)
            {
                SUMBX[I][4]=SUMBX[I][4]+H[1]*H[1]*SINA;
                SUMBX[I][5]=SUMBX[I][5]+H[2]*H[2]*SINA;
                SUMBX[I][6]=SUMBX[I][6]+H[3]*H[3]*SINA;
                SUMBX[I][7]=SUMBX[I][7]+H[1]*H[2]*SINA;
                SUMBX[I][8]=SUMBX[I][8]+H[1]*H[3]*SINA;
                SUMBX[I][9]=SUMBX[I][9]+H[2]*H[3]*SINA;
            }
        }
        TEMP[I]=exp(-XL_[I+IOF][4]*SSNN);
        SA[I]=SAI;
        SB[I]=SBI;
        NI=PTR_[I+IOF];
        BNI=coeff.DFPP[NI];
        FFX=coeff.DFP[NI];
        coeff.AC[10][NI]=0.0;
        for (II=1; II <= 9; II = II + 2) FFX=FFX+coeff.AC[II][NI]*exp(-coeff.AC[II+1][NI]*SSNN);
        SNEXI=FFX*XL_[I+IOF][5]*TEMP[I];
        SNXI=BNI*XL_[I+IOF][5]*TEMP[I];
        //-----CALCULATE A AND B OF F
        AV=SNEXI*SAI+AV;
        BV=SNEXI*SBI+BV;
        SNEX[I]=2.0*SNEXI*TAV*diffractogram.RATIO[ICX];
        SNX[I]=2.0*SNXI*TAV*diffractogram.RATIO[ICX];
        CV=CV+SNXI*SAI;
        DV=DV+SNXI*SBI;
    }
    FNN=diffractogram.RATIO[ICX]*(CV*CV+AV*AV+DV*DV+BV*BV)*TAV*SR;

    // PREPARING PHASE and struc fact TO BE PRINTED  !cp june 98
    refs[NN].TAVIX=TAV;
    refs[NN].SRIX=SR;
    if(AV == 0)AV=1E-6;
    refs[NN].APHASE=atan(BV/AV);
    if(AV < 0.0) refs[NN].APHASE = refs[NN].APHASE+M_PI;  // PI
    if(BV < 0 && AV == 0) refs[NN].APHASE =1.5*M_PI;     // 2/3 PI
    if(BV > 0 && AV == 0) refs[NN].APHASE =M_PI_2;  // PI/2
    //L120:
    refs[NN].FMGNTD=FNN;
    if(MAXS == 0) return;

    //-----CALCULATE DERIVATIVES
    for (I=1; I <= N; ++I)
    {
        SNEXI=SNEX[I];
        SNXI=SNX[I];
        SAI=SA[I];
        SBI=SB[I];
        for (J=1; J <= 11; ++J)
        {
            K=XL_[I+IOF][J].L;
            if(K != 0)
            {
                if(J > 5) goto L221;
                if(J <= 3)
                {
                    SUMA=SUMAX[I][J];
                    SUMB=SUMBX[I][J];
                    DER=-((AV*SUMA-BV*SUMB)*SNEXI+(CV*SUMA-DV*SUMB)*SNXI)*6.2831853071;
                    goto L26;
                }

                if(J > 4)
                {
                    DER=((SAI*AV+SBI*BV)*SNEXI+(SAI*CV+SBI*DV)*SNXI)/XL_[I+IOF][5];
                    goto L26;
                }

                DER=-((SAI*AV+SBI*BV)*SNEXI+(SAI*CV+SBI*DV)*SNXI)*SSNN;
                goto L26;

    L221:
                SUMA=SUMAX[I][J-2];
                SUMB=SUMBX[I][J-2];
                DER=-((AV*SUMA+BV*SUMB)*SNEXI+(CV*SUMA+DV*SUMB)*SNXI);
                if(J >= 9) DER = 2.0*DER;
    L26:
                DERIV[K] = sign(XL_[I+IOF][J].codeword)*DER+DERIV[K];
            }
        }
    }
    //-----CALCULATE DERIVATIVES
    //-----Preferred Orientation Derivatives
    K = phases[IPH_].PAR[11].L;
    if ( K != 0 )
    {
        DERIV[K] = DERIV[K]+FNN*DPRECOR;
        //         print '(3x,i3,3f10.5)',k,fnn,dprecor,deriv(k)         ! ********
    }

    //////////////////////////////////////////////////////////
    // TODO: remover este isso. Aideia é dividir a rotina INPTR em duas e fazer estes teste por lá

    K = phases[IPH_].PAR[12].L;
    if (IPREF == 0)
    {
        if(K != 0) DERIV[K]=DERIV[K]+(1.0-PREXP)*FNN/PRECOR;
    }
    else
    {
        if (K != 0)
        {
            file6 << "G2 IS NOT A REFINABLE PARAMETER FOR IPREF = 1" << endl;
            DBWSException("");
        }
    }
    ////////////////////////////////////////////////////

    //-----Derivatives for microabsorption parameter
    if (IABSR == 1)
    {
        K = GLB_[13-1].L;
        SRD = (1.0-GLB_[12-1])*(asin(sqrt(SINTH))-1.5707963268);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = GLB_[12-1].L;
        SRD = 1.0 -GLB_[8-1]*exp(-GLB_[9-1])+GLB_[8-1]*exp(-GLB_[9-1]/sqrt(SINTH))-1.0-GLB_[13-1]*(asin(sqrt(SINTH))-1.5707963268);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = GLB_[9-1].L;
        SRD = GLB_[8-1]*GLB_[12-1]*(exp(-GLB_[9-1])-exp(-GLB_[9-1]/sqrt(SINTH))/sqrt(SINTH));
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = GLB_[8-1].L;
        SRD = -GLB_[12-1]*(exp(-GLB_[9-1])+exp(-GLB_[9-1]/sqrt(SINTH)));
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
    }
    else if (IABSR == 2)
    {
        ///////////////////////////////
        // TODO: Eliminar este if
        KKL = GLB_[12-1].L;
        K  = GLB_[9-1].L;
        KL = GLB_[8-1].L;
        if (K != 0 || KL != 0 || KKL != 0)
        {
            file6 << "P AND/OR Q AND/OR R ARE NOT REFINABLE PARAMETERS FOR IABSR=2" << endl;
            cout << "P AND/OR Q AND/OR R ARE NOT REFINABLE PARAMETERS FOR IABSR=2" << endl;
            DBWSException("");
        }
        ///////////////////////////

        K = GLB_[13-1].L;
        SRD = 1.5707963268 - asin(sqrt(SINTH));
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
    }
    else if (IABSR == 3)
    {
        /////////////////////////////////////////////////////////
        // TODO: Eliminar este if!!!
        KKL = GLB_[13-1].L;
        K = GLB_[12-1].L;
        if (K != 0 || KKL != 0)
        {
            file6 << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            cout << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            DBWSException("");
        }
        ////////////////////////////////////////////////////////////

        K = GLB_[9-1].L;
        SRD = GLB_[8-1]*exp(-GLB_[9-1])-GLB_[8-1]*ISITH*exp(-GLB_[9-1]*ISITH);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN;
        K = GLB_[8-1].L;
        SRD = -exp(-GLB_[9-1])+exp(-GLB_[9-1]*ISITH);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN;
    }
    else if (IABSR == 4)
    {
        //////////////////////////////////////////
        // TODO: Eliminar este if!!!!
        KKL = GLB_[13-1].L;
        K = GLB_[12-1].L;
        if (K != 0 || KKL != 0)
        {
            file6 << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            cout << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            DBWSException("");
        }
        /////////////////////////////

        K = GLB_[9-1].L;
        SRD = GLB_[8-1]*(2*GLB_[9-1]-1)+GLB_[8-1]*ISITH*(2*GLB_[9-1]*ISITH-1);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = GLB_[8-1].L;
        SRD = GLB_[9-1]*(GLB_[9-1]-1)-GLB_[9-1]*ISITH*(1-GLB_[9-1]*ISITH);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
    }
    //----Overall Temperature and Scale Factor
    K=phases[IPH_].PAR[1].L;
    if(K != 0) DERIV[K]=DERIV[K]-2.0*SSNN*FNN;
    K=phases[IPH_].PAR[0].L;
    if(K != 0) DERIV[K]=DERIV[K]+FNN/phases[IPH_].PAR[0];
    SINTH=FNN*PI*SLABDA/(sqrt(SINTH*COSTH)*BB);
    SS=FNN/TL_;
    X=TANTH*TANTH;
    //-----Broadening Derivatives
    if (NPROF == _SplitPearsonVII || NPROF == _TCHZ) goto L9212;
    for (J=3; J <= 5; ++J)
    {
        K=phases[IPH_].PAR[J-1].L;
        if(K == 0) goto L78;
        DERIV[K]=X*SS+DERIV[K];
L78:
        X=X/TANTH;
    }
    K=phases[IPH_].PAR[20].L;
    if(K == 0) goto L9212;
    // for cot^2 case                                    !cp nov 29 96
    DERIV[K]=SS/(TANTH*TANTH)+DERIV[K];
    //-----Split Pearson VII Broadening Derivatives
L9212:
    if (NPROF == _SplitPearsonVII)
    {
        if (DELTA_ < 0.0)
        {
            DA3 = DA3L_;
            DA1 = DA1L_;
        }
        else
        {
            DA3 = DA3H_;
            DA1 = DA1H_;
        }
        for (J=3; J <= 5; ++J)
        {
            K=phases[IPH_].PAR[J-1].L;
            if(K == 0) goto L780;
            DERIV[K]=DERIV[K]+X*SS*((DA3*DELT_/(1.0+DA1*DELT_))-(1.0/TL_));
L780:
            X=X/TANTH;
        }
    }
    //L9211:

    //-----TCHZ Broadening Derivatives
    if (NPROF == _TCHZ)
    {
        TL_ = refs[NN].REFS[1];
        TLG = refs[NN].HALFG;
        TLL = refs[NN].HALFL;
        DHDHG = 0.2/pow(TL_,4.0)*(5.*pow(TLG,4.0)+10.77076*pow(TLG,3.0)*TLL+ 7.28529*TLG*TLG*TLL*TLL+8.94326*TLG*pow(TLL,3.0) + 0.07842*pow(TLL,4.0));
        DHDHL = 0.2/pow(TL_,4.0)*(2.69269*pow(TLG,4.0)+ 4.85686*pow(TLG,3.0)*TLL +13.41489*TLG*TLG*TLL*TLL + 0.31368*TLG*pow(TLL,3.0)+5.*pow(TLL,4.0));
        for (J=3; J <= 5; ++J)
        {
            K=phases[IPH_].PAR[J-1].L;
            if(K == 0)goto L9078;
            DERIV[K]=DHDHG*X*SS+DERIV[K];
L9078:
            X=X/TANTH;
        }
        K=phases[IPH_].PAR[19].L;
        if  (K != 0) DERIV[K]=DHDHG*SS/COSTH+DERIV[K];
        K = phases[IPH_].PAR[14].L;
        if (K == 0) goto L9213;
        DERIV[K] = 2.0*FNN*DHDHL*TANTH+DERIV[K];
L9213:
        K = phases[IPH_].PAR[15].L;
        if (K == 0) goto L9214;
        DERIV[K] = 2.0*FNN*DHDHL/sqrt(COSTH) + DERIV[K];
L9214:;
    }
    //-----Profile Shape Derivatives
    K = phases[IPH_].PAR[16].L;
    if(K != 0 && (NPROF == _pseudoVoigt || NPROF == _PearsonVII)) DERIV[K]=DERIV[K]+ FNN;
    K = phases[IPH_].PAR[17].L;
    if(K != 0 && NPROF == _pseudoVoigt) DERIV[K]=DERIV[K]+ FNN * refs[NN].REFS[2];
    if(K != 0 && NPROF == _PearsonVII) DERIV[K]=DERIV[K]+ FNN / refs[NN].REFS[2];
    K = phases[IPH_].PAR[18].L;
    if(K != 0 && NPROF == _PearsonVII)DERIV[K]=DERIV[K]+FNN/refs[NN].REFS[2] /refs[NN].REFS[2];
    //-----Split Pearson VII Shape Derivative
    if (NPROF == _SplitPearsonVII)
    {
        K = phases[IPH_].PAR[16].L;
        if (K != 0.0)
        {
            if (DELTA_ < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+DA1L_*DELT_/BB)+ DA7L_*DELT_/BB/(1.0+DA1L_*DELT_/BB));
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*DA6L_;
            }
        }
        K = phases[IPH_].PAR[17].L;
        if (K != 0.0)
        {
            if (DELTA_ < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+DA1L_*DELT_/BB)+ DA7L_*DELT_/BB/(1.0+DA1L_*DELT_/BB))/refs[NN].REFS[2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*DA6L_/refs[NN].REFS[2];
            }
        }
        K = phases[IPH_].PAR[18].L;
        if (K != 0.0)
        {
            if (DELTA_ < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+DA1L_*DELT_/BB)+DA7L_*DELT_/BB/(1.0+DA1L_*DELT_/BB))/refs[NN].REFS[2]/refs[NN].REFS[2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*DA6L_/refs[NN].REFS[2]/refs[NN].REFS[2];
            }
        }
        K = phases[IPH_].PAR[23].L;
        if (K != 0.0)
        {
            if (DELTA_ < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*DA6H_;
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+DA1H_*DELT_/BB)+DA7H_*DELT_/BB/(1.0+DA1H_*DELT_/BB));
            }
        }
        K = phases[IPH_].PAR[24].L;
        if (K != 0.0)
        {
            if (DELTA_ < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*DA6H_/refs[NN].REFS[2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+DA1H_*DELT_/BB)+DA7H_*DELT_/BB/(1.0+DA1H_*DELT_/BB))/refs[NN].REFS[2];
            }
        }
        K = phases[IPH_].PAR[25].L;
        if (K != 0.0)
        {
            if (DELTA_ < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*DA6H_/refs[NN].REFS[2]/refs[NN].REFS[2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+DA1H_*DELT_/BB)+DA7H_*DELT_/BB/(1.0+DA1H_*DELT_/BB))/refs[NN].REFS[2]/refs[NN].REFS[2];
            }
        }
    }
    //-----Split Pearson VII Asymmetry Derivative
    K = phases[IPH_].PAR[26].L;
    if (DELTA_ < 0.0)
    {
        DA1 = DA1L_;
        DA4 = DA4L_;
        DA5 = DA5L_;
    }
    else
    {
        DA1 = DA1H_;
        DA4 = DA4H_;
        DA5 = DA5H_;
    }
    if (K != 0 && NPROF == _SplitPearsonVII)
    {
        DERIV[K]=DERIV[K]+FNN*(DA4+DA5*DELT_/BB/(1.0+DA1*DELT_/BB));
    }
    //-----Zero, Displacement, and Transparancy Derivatives
    K=GLB_[1-1].L;
    if(K != 0)DERIV[K]=DERIV[K]+2.0*FNN/BB;
    K=GLB_[10-1].L;
    if (K != 0) DERIV[K]=DERIV[K]+2.0*FNN/BB*sqrt(COSTH);
    K=GLB_[11-1].L;
    if (K != 0) DERIV[K]=DERIV[K]+2.0*FNN/BB*sin(refs[NN].REFS[2]/57.2958);
    //c-----Lattice Parameter Derivatives
    for (J=1; J <= 6; ++J)
    {
        K=phases[IPH_].PAR[J+5-1].L;
        if(K == 0) goto L79;
        if(J < 4) X=HNN[J]*HNN[J];
        if(J == 4) X=HNN[2]*HNN[3];
        if(J == 5) X=HNN[1]*HNN[3];
        if(J == 6) X=HNN[1]*HNN[2];
        DERIV[K]= X*SINTH + DERIV[K];
L79:;
    }
    //-Asymmetry Derivative.  Test for asymmetry model included !cp may 97
    K=phases[IPH_].PAR[13].L;
    if((K != 0) && VERT)
    {
        if (IASYM == 0)
        {
            DERIV[K]=-FNN/TANTH+DERIV[K];
        }
        else
        {
            TANTHE=TANTH;
            if (TANTHE >= 1.0) TANTHE=tan(atan(TANTHE-M_PI_2));
            DERIV[K]=-FNN/TANTH+DERIV[K];
        }
    }
    //-----STORE DERIVATIVES FOR LIMO REFLECTIONS AT A TIME
    for (I=1; I <= MAXS; ++I) DERSTO_[NM][I]=DERIV[I];
    //L48:;
}

void DBWS::PRSVII(double T)
{
    //     CALCULATES THE COEFFICIENT C4 AND ITS DERIVATIVE WRT T FOR PHASE K
    //     THIS SUBROUTINE OBTAINED FROM IMMIRZI'S CODE
    //     RGD = GAMMA(T)/GAMMA(T-0.5)
    //     DGD = D(RGD)/D(T)

    const double RGD[] = { 0.0,
        0.156535,  0.183771,  0.209926,  0.235089,  0.259339,
        0.282749,  0.305380,  0.327289, 0.348527,  0.369141,  0.389170,
        0.408654,  0.427625,  0.446115,  0.464153,  0.481764,  0.498972,
        0.515799,  0.532265,  0.548390,  0.564190,  0.579680,  0.594877,
        0.609794,  0.624443,  0.638837,  0.652986,  0.666902,  0.680594,
        0.694071,  0.707342,  0.720415,  0.733297,  0.754996,  0.758519,
        0.770871,  0.783059,  0.795089,  0.806966,  0.818696,  0.830282,
        0.841731,  0.853045,  0.864230,  0.875290,  0.886227,  0.897046,
        0.907750,  0.918344,  0.928828,  0.939207,  0.949484};
    const double DGD[] = { 0.0,
        1.390576,  1.334036,  1.282282,  1.234747,  1.190975,
        1.150537,  1.113079,  1.078304,  1.045950,  1.015831,  0.987556,
        0.961144,  0.936314,  0.913013,  0.891052,  0.870228,  0.850726,
        0.832193,  0.814684,  0.798032,  0.782274,  0.767112,  0.752807,
        0.739098,  0.726022,  0.713542,  0.701621,  0.690110,  0.679232,
        0.668690,  0.658669,  0.648871,  0.639521,  0.630580,  0.621825,
        0.613555,  0.605434,  0.597760,  0.590086,  0.582896,  0.575967,
        0.569075,  0.562482,  0.555962,  0.549890,  0.543743,  0.538044,
        0.532381,  0.526905,  0.521578,  0.516400,  0.511296};
    const double AL2 = 0.6931472;

    int I,N1,IT;
    double DG, RG, FT,FL1;


    if (T < 0.6)
    {
        file6 << endl
            << "UNACCEPTABLE M VALUE IN CALCULATING PEARSON VII FUNCTION "
            << setw(6) << setprecision(2) << T << endl;
        DBWSException("");
    }
    IT=static_cast<int>(T-0.6);
    FT=T-static_cast<double>(IT);
    N1=static_cast<int>((FT-0.6)/0.02+1.0001);
    DG=DGD[N1]+(DGD[N1+1]-DGD[N1])*(FT-0.58-0.02*static_cast<double>(N1))/0.02;
    RG=RGD[N1]+(RGD[N1+1]-RGD[N1])*(FT-0.58-0.02*static_cast<double>(N1))/0.02;
    if (IT > 0)
    {
        for (I=1; I <= IT; ++I)
        {
            FL1=FT+static_cast<double>(I-1);
            DG=DG*(FL1)/(FL1-0.5)-0.5*RG/pow((FL1-0.5),2);
            RG=RG*FL1/(FL1-0.5);
        }
    }
    TF1_=sqrt( pow(2.0,(1.0/T))-1.0);
    C4_=2.0*RG*TF1_/M_SQRT2;
    TF2_=AL2* pow(2.0,(1.0/T))/(2.0*pow(TF1_,2));
    TF4_=4.0*AL2* pow(2.0,(1.0/T))/T;
    TF6_=4.0*T*pow(TF1_,2);
    TF8_=DG/RG-TF2_/pow(T,2);
    TF9_=4.0*pow(TF1_,2);
}

double DBWS::gamma(double X)
{
    double r;
    int I,IN;
    double GG,GX,DU1, DU2, DU4;


    if (X-1.0 < 0.0) goto L11; else goto L13;
L11:
    if (X-0.1 < 0.0) goto L97; else goto L12;
L12:
    GX=X+1.0;
    goto L21;
L13:
    if (X-2.0 <= 0.0) goto L14; else goto L15;
L14:
    GX=X;
    goto L21;
L15:
    if (X-8.0 <= 0.0) goto L16; else goto L98;
L16:
    for (I=1; I <= 8; ++I)
    {
        IN=I;
        GX=X-static_cast<double>(I);
        if (GX >= 1.0 && GX <= 2.0) goto L21;
    }
L21:
    DU1=GX-1.0;
    DU2=DU1*DU1;
    DU4=DU2*DU2;
    GG=1.0-0.5771917*DU1+0.9882059*DU2-0.8970569*DU1*DU2 +0.9182069*DU4-0.7567041*DU4*DU1+0.4821994*DU4*DU2 -0.1935278*DU1*DU2*DU4+0.0358683*DU4*DU4;
    if (X-1.0 < 0.0) goto L31; else goto L32;
L31:
    r=GG/X;
    goto L99;
L32:
    if (X-2.0 <= 0.0) goto L33; else goto L34;
L33:
    r=GG;
    goto L99;
L34:
    r=X-1.0;
    if (IN-1 <= 0.0) goto L35; else goto L36;
L35:
    r=r*GG;
    goto L99;
L36:
    for (I=2; I <= IN; ++I) r=r*(X-static_cast<double>(I));
    r=r*GG;
    goto L99;
L97:
    cout << "X IN GAMMA(X) IS LESS THAN 0.1   X =" << setw(8) << setprecision(4) << X << endl;
    r=9.513508;
    goto L99;
L98:
    cout << "X IN GAMMA(X) IS GREATER THAN 8.0   X =" << setw(8) << setprecision(4) << X << endl;
    r=5040.0;
L99:
    return r;
}

//    Split Pearson VII coding for function and derivatives
//    Based on Toraya's Code from Profit, Ver 1.22N
void DBWS::mspvii(double A, double W)
{
    double DH, DL, DHD, DLD;

    DL = gamma(RL_-0.5)/(sqrt( pow(2.0,(1.0/RL_))-1.0)*gamma(RL_));
    DH = gamma(RH_-0.5)/(sqrt( pow(2.0,(1.0/RH_))-1.0)*gamma(RH_));
    DLD = gamma(RL_-0.5+0.001)/(sqrt(pow(2.0,(1.0/(RL_+0.001)))-1.0)*gamma(RL_+0.001));
    DHD = gamma(RH_-0.5+0.001)/(sqrt( pow(2.0,(1.0/(RH_+0.001))) -1.0)*gamma(RH_+0.001));
    DA1L_ = ( pow(((1.0+A)/A),2.0))*( pow(2.0,(1.0/RL_))-1.0);
    DA1H_ = (pow((1.0+A),2.0))*(pow(2.0,(1/RH_))-1.0);
    DA2L_ = 1.128379167*(1.0+A)*(1.0/(A*DL+DH))/W;
    DA2H_ = DA2L_;
    DA3L_ = 2.0*RL_*DA1L_/W;
    DA3H_ = 2.0*RH_*DA1H_/W;
    DA4L_ = (1.0/(1.0+A))-(DL/(A*DL+DH));
    DA4H_ = DA4L_;
    DA5L_ = 2.0*RL_*pow(((1+A)/A),2.0)*(1.0+A)/(pow(A,3));
    DA5H_ = -2.0*RH_*pow((1+A),2.0)*(1+A);
    DA6L_ = -A*1000.0*(DLD-DL)/(A*DL+DH);
    DA6H_ = -1000.0*(DHD-DH)/(A*DL+DH);
    DA7L_ = log(2.0)*pow(((1.0+A)/A),2.0)*(pow(2.0,(1.0/RL_))-1.0+1.0)/RL_;
    DA7H_ = log(2.0)*pow((1.0+A),2.0)*(pow(2.0,(1.0/RH_))-1.0+1.0)/RH_;
}

double DBWS::PROFIL(int N, double X)
{
    double r;

    switch (N) {
    case _Gaussian:
        r=0.939437279*exp(-2.772588722*X);
        PRFDER_=2.772588722;
        break;
    case _Lorentzian:
        r=0.636619772/(1.0+4.0*X);
        PRFDER_=4.0/(1.+4.0*X);
        break;
    case _Mod1:
        r=0.819449653/ pow((1.0+1.656854248*X),2.0);
        PRFDER_=3.313708496/(1.0+1.656854248*X);
        break;
    case _Mod2:
        r=0.766420937/pow((1.0+2.349604208*X),1.5);
        PRFDER_=3.5244063/(1.0+2.349604208*X);
        break;
    case _SplitPearsonVII:
        if (DELTA_ < 0.0)
        {
            r = pow((1.0+DA1L_*X),(-RL_)) *DA2L_;
            PRFDER_ = DA1L_*RL_/(1.0+DA1L_*X);
        }
        else
        {
            r = pow((1.0+DA1H_*X),(-RH_)) *DA2H_;
            PRFDER_ = DA1H_*RH_/(1.0+DA1H_*X);
        }
        break;
    case _pseudoVoigt:
        r=GAM1_*0.636619772/(1.0+4.0*X)+(1-GAM1_)*0.939437279*exp(-2.772588722*X);
        PRFDER_=(GAM1_*2.546479088/ pow((1.0+4.0*X),2) +(1-GAM1_)*2.6046732048*exp(-2.772588722*X))/r;
        break;
    case _PearsonVII:
        r=C4_/pow((1.0+4.0*( pow(2,(1.0/GAM1_)) -1)*X) , GAM1_);
        PRFDER_=TF6_*GAM1_/(1.0+4.0*( pow(2,(1.0/GAM1_))-1)*X);
        break;
    case _TCHZ:
        r=GAM1_*0.636619772/(1.0+4.0*X)+(1-GAM1_)*0.939437279*exp(-2.772588722*X);
        PRFDER_=(GAM1_*2.546479088/ pow((1.0+4.0*X),2) +(1-GAM1_)*2.6046732048*exp(-2.772588722*X))/r;
        break;
    default:
        DBWSException("ILLEGAL PROFILE FUNCTION REQUEST");
    }
    return r;
}

void DBWS::SUMMAT(int IPM, double ISODER[], double TOTCS)
{
    int I, J, K, M,II, IL, KK, LK, KM,LK1, LK2,IOF,KRP1, IISO;
    double X, Z, X1, BB,YX,DER, TLL, TLR,ASS5,ESSE,OMEGA, YCALC,
        TANNJ, SHIFT, OMEGA8,DERMON,PRTEMP;
    bool VERT;
    double DERIV[MSZ+1];

    SHIFT = 0.0;
    for (J=1; J <= MSZ; ++J) DERIV[J]=0.0;
    YCALC=0.0;
    IL=0;
    //-----CALCULATE THE CONTRIBUTION OF THE REFLECTIONS ORD1 TO ORD2 TO THE
    //-----DERIVATIVES W.R.T. THE PROFILE INTENSITY YOBS
    if( _IORD1 == 0) goto L12;
    for (I=_IORD1; I <= _IORD2; ++I)
    {
        IPH_=refs[I].IREFS/(256*256*256*8);
        IL=IL+1;
        J=(I % NOV)+1;
        // test for asymmetry function !cp ap 12 97
        if (IASYM == 0)
        {
            VERT = refs[I].REFS[2] <= RLIM_;
        }
        else
        {
            VERT = abs(refs[I].REFS[2]-90.0) >= RLIM_;
        }
        SHIFT = GLB_[10-1] * cos(refs[I].REFS[2]/2.0/57.2958) + GLB_[11-1] * sin(refs[I].REFS[2]/57.2958);
        DELTA_=_TH-refs[I].REFS[2]-GLB_[1-1]-SHIFT;
        //TANTH=tan(_TH*3.14159265359/360.0);
        DELT_=DELTA_*DELTA_;
        TL_=refs[I].REFS[1];
        if (NPROF == _pseudoVoigt) GAM1_ = phases[IPH_].PAR[16] + phases[IPH_].PAR[17] * refs[I].REFS[2];
        if (NPROF == _PearsonVII) GAM1_ = phases[IPH_].PAR[16] + phases[IPH_].PAR[17] / refs[I].REFS[2] + phases[IPH_].PAR[18]/refs[I].REFS[2]/refs[I].REFS[2];
        if (NPROF == _PearsonVII) PRSVII(GAM1_);
        if (NPROF == _SplitPearsonVII)
        {
            RL_=phases[IPH_].PAR[16]+(phases[IPH_].PAR[17]+phases[IPH_].PAR[18]/refs[I].REFS[2])/refs[I].REFS[2];
            RH_=phases[IPH_].PAR[23]+(phases[IPH_].PAR[24]+phases[IPH_].PAR[25]/refs[I].REFS[2])/refs[I].REFS[2];
            mspvii(phases[IPH_].PAR[26],TL_);
        }
        if(NPROF == _TCHZ)
        {
            TLL = refs[I].HALFL;
            GAM1_ = refs[I].GAM;
            TLR = TLL/TL_;
        }
        BB=TL_*TL_;
        //-----NEXT LINE IS NECESSEARY FOR 2 PHASES WITH VERY DifFERENT FWHM.
        if (DELT_/BB > WDT_*WDT_) goto L33;
        if (VERT)
        {
            //       test for asymmetry model               !cp may 01 97
            if (IASYM == 0)
            {
                YX=DELT_*sign(DELTA_);
                Z=1.0-phases[IPH_].PAR[13]*YX/TANN_[J];
                if ( Z <= 0.0 ) Z=0.0001;
            }
            else
            {
                YX=sign(DELTA_)*DELTA_/(2*TL_);
                TANNJ=TANN_[J];
                if (TANNJ >= 1.0) TANNJ=tan(atan(TANNJ)-M_PI_2);
                Z=(phases[IPH_].PAR[13]/TANNJ) * (2.0*(DELTA_/(2*TL_))*exp(-YX));
                Z=1+Z;
                if ( Z <= 0.0 ) Z=0.0001;
            }
        }
        else
        {
            Z=1.0;
        }
        //L5:
        PRTEMP = PROFIL(NPROF,DELT_/BB);
        if (NPROF == _SplitPearsonVII)
        {
            OMEGA = refs[I].REFS[3]*Z*PRTEMP*phases[IPH_].PAR[0];
        }
        else
        {
            OMEGA = refs[I].REFS[3]*Z*PRTEMP*phases[IPH_].PAR[0]/TL_;
        }
        YCALC = YCALC+OMEGA*refs[I].FMGNTD;
        if ( JOBTYP > 2 ) goto L33;
        X = PRFDER_*2.0* DELT_/BB-1.0;
        for (K=1; K <= MAXS; ++K)
        {
            DER=1.0;
            //-----Broadening Coeficients Derivatives
            if(NPROF != _SplitPearsonVII)
            {
                for (M=3; M <= 5; ++M) if(phases[IPH_].PAR[M-1].L == K) DER=X/TL_/2.0;
            }
            if (phases[IPH_].PAR[19].L == K) DER=X/TL_/2.0;
            X1=0.0;
            //-----Asymmetry Derivative
            if (VERT)
            {
                if (IASYM == 0)
                {
                    X1=phases[IPH_].PAR[13]*sign(DELTA_)*BB/TANN_[J]/Z;
                }
                else
                {
                    X1=-phases[IPH_].PAR[13]*exp(-YX)*(TL_/(2*DELTA_)-sign(DELTA_)*1.0/4)/TANNJ/Z;
                }
            }
            //-----Zero, Displacement, and Transparancy Derivative
            if ( GLB_[1-1].L == K ) DER=DELTA_*(PRFDER_+X1);
            if ( GLB_[10-1].L == K ) DER=DELTA_*(PRFDER_+X1);
            if ( GLB_[11-1].L == K ) DER=DELTA_*(PRFDER_+X1);
            if ( (phases[IPH_].PAR[13].L == K)  &&  VERT )
            {
                if ( IASYM == 0 )
                {
                    DER = YX/Z;
                }
                else
                {
                    DER = -2.0*(DELTA_/(2.0*TL_))*exp(-YX)/Z;
                }
            }
            if ( NPROF == _TCHZ ) goto L8;
            //-----Pseudo-Voigt Shape Derivatives
            if(NPROF == _pseudoVoigt)
            {
                KRP1=phases[IPH_].PAR[16].L;
                if (K == KRP1) DER=(0.636619772/(1.0+4.0*DELT_/BB)-0.939437279*exp(-2.772588722*DELT_/BB))/PRTEMP;
                KRP1=phases[IPH_].PAR[17].L;
                if (K == KRP1) DER=(0.636619772/(1.0+4.0*DELT_/BB)-0.939437279*exp(-2.772588722*DELT_/BB))/PRTEMP;
            }
            //-----Pearson VII Shape Derivatives
            if (NPROF == _PearsonVII)
            {
                KRP1=phases[IPH_].PAR[16].L;
                if(K == KRP1) DER=-log(1.0+TF9_*DELT_/BB)+TF4_*(DELT_/BB)/(1.0+TF9_*DELT_/BB)+TF8_;
                KRP1=phases[IPH_].PAR[17].L;
                if(K == KRP1)   DER=-log(1.0+TF9_*DELT_/BB)+TF4_*(DELT_/BB)/(1.0+TF9_*DELT_/BB)+TF8_;
                KRP1=phases[IPH_].PAR[18].L;
                if(K == KRP1) DER=-log(1.0+TF9_*DELT_/BB)+TF4_*(DELT_/BB)/(1.0+TF9_*DELT_/BB)+TF8_;
            }
            //-----Lattice Parameter Derivatives
L8:
            for (M=6; M <= 11; ++M) if(phases[IPH_].PAR[M-1].L == K) DER=(PRFDER_+X1)*DELTA_;
            DERIV[K]=DERSTO_[J][K]*DER*OMEGA+DERIV[K];
            //L3:;
        }
        //----TCHZ Profile Derivatives
        if(NPROF == _TCHZ)
        {
            OMEGA8 = Z*phases[IPH_].PAR[0]*refs[I].REFS[3]/TL_;
            for (K = 1; K <= MAXS; ++K)
            {
                for (M=3; M <= 5; ++M)
                {
                    if (phases[IPH_].PAR[M-1].L != K) goto L1001;
                    DERIV[K] = DERIV[K]+ OMEGA8*DERSTO_[J][K]/2.0*(0.939437279*exp(-2.772588722*DELT_/BB) - 0.636619772/(1.0+4.0*DELT_/BB)) * (1.36603*TLR/TL_-0.95438*TLR*TLR/TL_+0.33348 * pow(TLR,3.0)/TL_);
L1001:;
                }
                if (phases[IPH_].PAR[19].L == K) DERIV[K] = DERIV[K]+ OMEGA8*DERSTO_[J][K]/2.0*(0.939437279*exp(-2.772588722*DELT_/BB) - 0.636619772/(1.0+4.0*DELT_/BB)) * (1.36603*TLR/TL_-0.95438*TLR*TLR/TL_+0.33348*pow(TLR,3.0)/TL_);
                if (phases[IPH_].PAR[14].L != K) goto L1002;
                DERIV[K] = DERIV[K]+ OMEGA8*(0.939437279*exp(-2.772588722*DELT_/BB) -0.636619772/(1.0+4.0*DELT_/BB)) *((1.36603*TLR/TL_-0.95438*TLR*TLR/TL_+0.33348*pow(TLR,3.0)/TL_)* DERSTO_[J][K]/2.0 - refs[I].FMGNTD*TANN_[J]*(1.36603/TL_-0.95438*TLR/TL_+0.33348*TLR*TLR/TL_));
L1002:
                if (phases[IPH_].PAR[15].L != K) goto L1003;
                DERIV[K] = DERIV[K]+ OMEGA8*(0.939437279*exp(-2.772588722*DELT_/BB) -0.636619772/(1.0+4.0*DELT_/BB)) *((1.36603*TLR/TL_-0.95438*TLR*TLR/TL_+0.33348*pow(TLR,3.0)/TL_)* DERSTO_[J][K]/2.0 - refs[I].FMGNTD*sqrt(1+TANN_[J]*TANN_[J])*(1.36603/TL_-0.95438*TLR/TL_+0.33348*TLR*TLR/TL_));
L1003:;
            }
        }
L33:;
    }
    //-----FORM SUMS
L12:
    if(diffractogram.NBCKGD != 0)goto L11;
    for (II=2; II <= 7; ++II)
    {
        if(GLB_[II-1].L == 0)goto L10;
        KM=GLB_[II-1].L;
        if(II == 2)DERIV[KM]=DERIV[KM]+1.0;
        if(II == 2)goto L10;
        DERIV[KM]=DERIV[KM]+ pow( ((diffractogram.THMIN+static_cast<double>(IPM-1)*diffractogram.STEP)/BKPOS_-1.0) , (II-2));
L10:;
    }
L11:
    diffractogram.YC[IPM]=YCALC;
    if (JOBTYP > 2) goto L20;
    for (K = 1; K <= NPHASE; ++K)
    {
        LK1  = phases[K].PAR[0].L;
        LK2  = phases[K].PAR[1].L;
        //
        //-----UPDATING GLOBAL SCALE DERIVATE FOR BKG CONTRIBUTE
        //
        if(FONDO == 1 || FONDO == 2)
        {
            if(LK1 != 0) DERIV[LK1]=DERIV[LK1]+phases[K].GCOM*phases[K].CSK+phases[K].GCOM*phases[K].DISK;
        }
        //
        //-----UPDATING DERIVATE OF Q OVERALL FOR BKG CONTRIBUTE
        //
        if(FONDO == 2)
        {
            if(LK2 != 0) DERIV[LK2] = DERIV[LK2] +  phases[K].DYCDD;
        }
        //
        //-----UPDATING DERIVATE OF ISOTROPIC THERMAL PARAMETERS FOR BKG CONTRIBUTE
        //
        if(FONDO == 1)
        {
            IOF = 0;
            if(K > 1)
            {
                for (I = 2; I <= K; ++I) IOF = IOF + phases[I-1].AtomCount;
            }
            for (I = 1; I <= phases[K].AtomCount; ++I)
            {
                IISO=XL_[I+IOF][4].L;
                if(IISO != 0) DERIV[IISO] = DERIV[IISO] + ISODER[I+IOF];
                ISODER[I+IOF]=0.0;
            }
        }
    }
    //-----DYC RESPECT TO AMORPHOUS SCALE FACTOR
    LK = GLB_[20-1].L;
    if(LK != 0) DERIV[LK] = DERIV[LK] + diffractogram.AMORPHOUS[IPM];
    //-----MONOCHROMATOR PARAMETERS DERIVATIVES

    diffractogram.CalcLamdaM();
    ESSE=2*sin((_TH-GLB_[1-1]-SHIFT)*0.00872665)/diffractogram.LAMDAM;
    LK=GLB_[18-1].L;
    if ( LK != 0 )
    {
        //-------NEXT LINE FOR A LORENTZIAN MONOCHROMATOR BASS-BAND  FUNCTION
        ASS5=1/(1+GLB_[18-1]*  pow(ESSE,GLB_[19-1]) );
        DERMON=-( pow(ESSE,GLB_[19-1]) )/ ( pow((1+GLB_[18-1]*  pow(ESSE,GLB_[19-1]) ),2) );
        //--------------------------------------------------------------------
        DERIV[LK]=DERIV[LK]+TOTCS/ASS5*DERMON;
    }
    //   !cp ap 20 97
    LK=GLB_[19-1].L;
    ASS5=1/(1+GLB_[18-1]*  pow(ESSE,GLB_[19-1]) );
    if ( LK != 0 )
    {
        DERMON=-GLB_[18-1]*log(ESSE)*( pow(ESSE,GLB_[19-1]) )/ pow( (1+GLB_[18-1]*(  pow(ESSE,GLB_[19-1]) )) , 2);
        DERIV[LK]=DERIV[LK]+TOTCS/ASS5*DERMON;
    }
    //-----FORM THE UPPER TRIANGULAR OF "RJAC" = MATRIX OF NORMAL EQUATIONS
    DELTA_ = diffractogram.Y[IPM]-diffractogram.BK[IPM]-YCALC;
    for (J=1; J <= MAXS; ++J)
    {
        X = DERIV[J]/(diffractogram.VAR[IPM]);
        VX[J] = VX[J]+DELTA_*X;
        for (KK=J; KK <= MAXS; ++KK) RJAC[J][KK] = RJAC[J][KK]+X*DERIV[KK];
    }
L20:;
}

double DBWS::DPINV(double A1[][MSZ+1], double B1[], int* N)
{
    int I,J,K,I1,MI,MJ,MI1,MK1;
    double r,C,D,S,T;
    double A[MSZ+1][MSZ+1];
    double B[MSZ+1];
    double U[MSZ+1];
    int M[MSZ+1];
    //DIMENSION ,A1(MSZ,MSZ),,B1(MSZ),

    if ( *N != 1 ) goto L50;
    A1[1][1] = 1.0/A1[1][1];
    B1[1] = B1[1]*A1[1][1];
    r = 1;
    goto L99;
L50:
    r = 0;
    for (I=1; I <= *N; ++I)
    {
        B[I]=B1[I];
        for (J=1; J <= *N; ++J) A[I][J]=A1[I][J];
    }
    for (I=1; I <= *N; ++I)
    {
        U[I]=abs(A[I][1]);
        C=U[I];
        for (J=2; J <= *N; ++J)
        {
            U[I]=max(U[I],abs(A[I][J]));
            C=C+abs(A[I][J]);
        }
        r=max(r,C);
        M[I]=I;
        B[I]=B[I]/U[I];
        for (J=1; J <= *N; ++J) A[I][J]=A[I][J]/U[I];
    }
    for (I=1; I <= *N; ++I)
    {
        if(I == *N)goto L15;
        J=I;
        I1=I+1;
        MI=M[I];
        S=abs(A[MI][I]);
        for (K=I1; K <= *N; ++K)
        {
            MI=M[K];
            if(abs(A[MI][I]) <= S)goto L4;
            S=abs(A[MI][I]);
            J=K;
L4:;
        }
        MK1=M[J];
        M[J]=M[I];
        M[I]=MK1;
L15:
        MI1=M[I];
        A[MI1][I]=1./A[MI1][I];
        B[MI1]=B[MI1]*A[MI1][I];
        for (J=1; J <= *N; ++J)
        {
            if(I == J)goto L20;
            A[MI1][J]=A[MI1][J]*A[MI1][I];
L20:;
        }
        for (J=1; J <= *N; ++J)
        {
            if(J == I)goto L30;
            MJ=M[J];
            T=A[MJ][I];
            MI1=M[I];
            B[MJ]=B[MJ]-T*B[MI1];
            A[MJ][I]=-T*A[MI1][I];
            for (K=1; K <= *N; ++K)
            {
                if(K == I)goto L40;
                A[MJ][K]=A[MJ][K]-T*A[MI1][K];
L40:;
            }
L30:;
        }
        //L10:;
    }
    for (I=1; I <= *N; ++I)
    {
        MI1=M[I];
        B1[I]=B[MI1];
        for (J=1; J <= *N; ++J)
        {
            MJ=M[J];
            A1[I][MJ]=A[MI1][J]/U[MJ];
        }
    }
    C=0;
    for (I=1; I <= *N; ++I)
    {
        D=abs(A1[I][1]);
        for (J=2; J <= *N; ++J) D=D+abs(A1[I][J]);
        C=max(C,D);
    }
    r=r*C;
L99:;
    return r;
}

double DBWS::ERROR(double A[][6+1], double B[], double* OMEGA)
{
    double r,X,SUM;
    int I,J;
    //DIMENSION A(6,6),B(6)

    SUM=0.0;
    for (I=1; I <= 6; ++I)
    {
        for (J=I; J <= 6; ++J)
        {
            X=2.0;
            if(I == J) X=1.0;
            SUM=SUM+A[I][J]*B[I]*B[J]*X;
        }
    }
    if (SUM < 0.0) SUM=0.0;
    r=sqrt(SUM * *OMEGA);
    return r;
}

void DBWS::ESD(double SM[][6+1], double V[], double SUM)
{
    int I,IA,IB,IC,ID,IE,IP;
    double DEN,DENL,FNUM,FNUML;
    double S[6+1],R[6+1],E[6+1];
    //DIMENSION SM(6,6),V(6)

    I=0;
    for (IA=1; IA <= 3; ++IA)
    {
        I=I+1;
        IB=IA+1;
        if(IB > 3)IB=IB-3;
        IC=IB+1;
        if(IC > 3)IC=IC-3;
        ID=IA+3;
        IE=IB+3;
        IP=IC+3;
        FNUM=4.0*V[IB]*V[IC]-V[ID]*V[ID];
        DEN=4.0*V[IA]*V[IB]*V[IC]-V[IA]*V[ID]*V[ID]-V[IB]*V[IE]*V[IE]-V[IC]*V[IP]*V[IP]+V[ID]*V[IE]*V[IP];
        S[I]=sqrt(FNUM/DEN);
        R[IA]= pow(-S[I],4) ;
        DENL=1./(DEN*DEN);
        R[IB]=DENL*(4.0*V[IC]*DEN-(4.0*V[IA]*V[IC]-V[IE]*V[IE])*FNUM);
        R[IC]=DENL*(4.0*V[IB]*DEN-(4.0*V[IA]*V[IB]-V[IP]*V[IP])*FNUM);
        R[ID]=-DENL*(2.0*V[ID]*DEN-(2.0*V[IA]*V[ID]-V[IE]*V[IP])*FNUM);
        R[IE]=DENL*FNUM*(2.0*V[IB]*V[IE]-V[ID]*V[IP]);
        R[IP]=DENL*FNUM*(2.0*V[IC]*V[IP]-V[ID]*V[IE]);
        E[I]=ERROR(SM,R,&SUM);
        E[I]=E[I]/(2.0*S[I]);
        FNUM=V[IE]*V[IP]-2.0*V[IA]*V[ID];
        DEN=16.0*V[IA]*V[IA]*V[IB]*V[IC]+V[IE]*V[IE]*V[IP]*V[IP]-4.0*(V[IA]*V[IB]*V[IE]*V[IE]+V[IC]*V[IA]*V[IP]*V[IP]);
        DENL=1.0/(DEN*DEN);
        FNUML=FNUM*FNUM;
        S[I+3]=FNUML/DEN;
        R[IA]=-4.0*FNUM*DENL*(V[ID]*DEN+FNUM*(8.0*V[IA]*V[IB]*V[IC]-V[IB]*V[IE]*V[IE]-V[IC]*V[IP]*V[IP]));
        R[IB]=-4.0*FNUML*DENL*(4.0*V[IA]*V[IA]*V[IC]-V[IA]*V[IE]*V[IE]);
        R[IC]=-4.0*FNUML*DENL*(4.0*V[IA]*V[IA]*V[IB]-V[IA]*V[IP]*V[IP]);
        R[ID]=-4.0*FNUM/DEN*V[IA];
        R[IE]=2.0*FNUM*DENL*(V[IP]*DEN-FNUM*(V[IE]*V[IP]*V[IP]-4.0*V[IA]*V[IB]*V[IE]));
        R[IP]=2.0*FNUM*DENL*(V[IE]*DEN-FNUM*(V[IE]*V[IE]*V[IP]-4.0*V[IC]*V[IA]*V[IP]));
        E[I+3]=ERROR(SM,R,&SUM);
        if(S[I+3] == 0.0) goto L3;
        E[I+3]=E[I+3]/(2.0*sqrt(S[I+3]*(1.0-S[I+3])))*180.0/M_PI;
        S[I+3]=atan(sqrt((1.0-S[I+3])/S[I+3]))*180.0/M_PI;
        goto L4;
L3:
        S[I+3]=90.0;
        E[I+3]=0.0;
L4:
        if(FNUM < 0.0)S[I+3]=180.0-S[I+3];
        if( SM[ID][ID] == 0.0) E[ID]=0.0;
    }
    if(abs(S[1]-S[2]) < 0.00008) E[1]=E[2];
    for (I=1; I <= 6; ++I)
    {
        SM[I][I]=E[I];
        V[I]=S[I];
    }
}

void DBWS::DIRECT(double SM[][6+1], double V[], int* IPH)
{
    int J,K,L,M;
    double X;

    for (J=1; J <= 6; ++J)
    {
        V[J]=phases[*IPH].PAR[J+5-1];
        K=phases[*IPH].PAR[J+5-1].L;
        X=phases[*IPH].PAR[J+5-1].codeword;
        if(K == 0)X=0.0;
        for (L=J; L <= 6; ++L)
        {
            M=phases[*IPH].PAR[L+5-1].L;
            SM[L][J]=0.0;
            if((M == 0) || (K == 0))goto L1;
            SM[L][J]=RJAC[M][K]*X*phases[*IPH].PAR[L+5-1].codeword;
            L1:
            SM[J][L]=SM[L][J];
        }
    }
    ESD(SM,V,1.0);
}

void DBWS::OUTSCR(int I, double R2, double R3, double X)
{
    cout << "CYCLE NUMBER - " << I << endl;
    cout << "R-P =" << setw(8) << setprecision(2) << R2 << "%"
        << "    R-WP =" << setw(8) << setprecision(2) << R3 << "%"
        << "       R-EXPECTED =" << setw(8) << setprecision(2) << X << "%"
        << "  S = " << setw(8) << setprecision(2) << R3/X << endl;
}

void DBWS::OUTPTR(int ICYCLE)
{
    int I, J, N, IP, KM,IOF, MMM,ILOC1, ILOC2, MATCH, ICOCO,
        IIPHAS,IINNMOL;
    double RATIODCV,X, V0, SRATIODCV,FT,ARG1, ARG2, ARG3, XFAC,SFIS,VOSQ,
        ARGCOS, WTOTAL;

    double SY[30+1];
    double SZ[30+1];
    double DUMMY[2*MSZ+4  +1];
    double VOL[99+1];
    double DVOL[99+1];
    double W[99+1];
    double DW[99+1];
    double XMASS[99+1];
    double DMASS[99+1];
    double FR[99+1];
    double FRP[99+1];
    double DFR[99+1];
    double STMASSA[99+1];
    double SMASS[99+1];

    for (I=1; I <= 2*MSZ+4; ++I) DUMMY[I] = 0.0;
    file6 << endl
        << "R-P        = " << setw(8) << setprecision(2) << diffractogram.R2_ << "%" << endl
        << "R-WP       = " << setw(8) << setprecision(2) << diffractogram.R3_ << "%     "
        << "R-WP(Background Removed) = " << setw(8) << setprecision(2) << diffractogram.R3NOBK_ << "%" << endl;
    if ( MAXS == 0  &&  MAXSX != 0 )
    {
        X=100.0*sqrt((static_cast<double>(_NUM)-static_cast<double>(MAXSX))*1.0/diffractogram.D2_);
    }
    else
    {
        X=100.0*sqrt((static_cast<double>(_NUM)-static_cast<double>(MAXS))*1.0/diffractogram.D2_);
    }
    if ( MAXS == 0  &&  MCYCLE == 1 )
    {
        OUTSCR(MCYCLX,diffractogram.R2_,diffractogram.R3_,X);
    }
    else
    {
        OUTSCR(ICYCLE-1,diffractogram.R2_,diffractogram.R3_,X);
    }
    file6
        << "R-EXPECTED = " << setw(8) << setprecision(2) << X << "%" << endl
        << "S          = " << setw(8) << setprecision(2) << sqrt(diffractogram.S2_/static_cast<double>(_NUM-MAXS))
        << "     SQRT(RESIDUAL/N-P)GOODNESS OF FIT" << endl
        << "D - W D    = " << setw(8) << setprecision(2) << diffractogram.SS4_/diffractogram.SS2_ << "     UNWEIGHTED DURBIN-WATSON STATISTIC D" << endl;
    I=_NUM-MAXS;
    file6
        << "N-P        = " << setw(8) << I << endl

        << endl
        << scientific
        << "SUMYDIF     = " << setw(14) << setprecision(8) << diffractogram.S1_ << endl
        << "SUMYOBS     = " << setw(14) << setprecision(8) << diffractogram.D1_ << endl
        << "SUMYCALC    = " << setw(14) << setprecision(8) << diffractogram.S3_ << endl
        << "SUMWYOBSSQ  = " << setw(14) << setprecision(8) << diffractogram.D2_ << endl
        << "RESIDUAL    = " << setw(14) << setprecision(8) << diffractogram.S2_ << endl
        << "CONDITION   = " << setw(14) << setprecision(8) << _COND << endl
        << fixed << endl;

    DUMMY[2*MSZ+1] = diffractogram.R2_;
    DUMMY[2*MSZ+2] = diffractogram.R3_;
    DUMMY[2*MSZ+3] = sqrt(diffractogram.S2_/static_cast<double>(_NUM));
    DUMMY[2*MSZ+4] = diffractogram.SS4_/diffractogram.SS2_;

    //     FINAL PARAMETERS AND R-FACTORS
    if (IPLST != 0 && MAXS == 0)
    {
        for (IP=1; IP <= NPHASE; ++IP)
        {
            IOF=0;
            if (IP > 1)
            {
                for (IIPHAS=2; IIPHAS <= IP; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
            }
            N=phases[IP].AtomCount;
            for (I=1; I <= N; ++I)
            {
                for (J=1; J <= 5; ++J)
                {
                    KM=XL_[I+IOF][J].L;
                    if(KM != 0)
                    {
                        //  !cp jun 96 start
                        if (J ==  5)
                        {
                            DUMMY[KM] = XL_[I+IOF][J]*phases[IP].XMLTP/MURT[I+IOF];
                        }
                        else
                        {
                            // !cp jun 96 stop
                            DUMMY[KM] = XL_[I+IOF][J];
                        }
                        // !cp sept 96 start
                        // !cp sept 96 stop
                    }
                }
            }
            for (I=1; I <= N; ++I)
            {
                for (J=6; J <= 11; ++J)
                {
                    KM=XL_[I+IOF][J].L;
                    if(KM != 0)
                    {
                        //  !cp jun 96 start
                        if (J ==  5)
                        {
                            DUMMY[KM] = XL_[I+IOF][J]*phases[IP].XMLTP/MURT[I+IOF];
                        }
                        // !cp jun 96 stop
                        DUMMY[KM]  =XL_[I+IOF][J];
                    }
                }
            }
            for (J=1; J <= 27; ++J)
            {
                KM=phases[IP].PAR[J-1].L;
                if (KM != 0) DUMMY[KM] = phases[IP].PAR[J-1];
            }
            DIRECT(DCSM,DCV,&IP);
            for (I=1; I <= 6; ++I)
            {
                KM = phases[IP].PAR[I+5-1].L;
                MATCH = 0;
                if (KM != 0)
                {
                    if (I >= 4)
                    {
                        for (MMM=1; MMM <= 3; ++MMM) if (KM == phases[IP].PAR[MMM+5-1].L) MATCH=1;
                    }
                    if (MATCH == 0) DUMMY[KM] = SAVE[IP][I];
                }
            }
        }
        for (J=1; J <= 20; ++J)
        {
            if(NPROF == _TCHZ && (J == 14 || J == 15 || J == 16)) goto L6019;
            if(NPROF == _PearsonVII && (J == 14 || J == 15 || J == 16)) goto L6019;
            if(NPROF == _pseudoVoigt && (J == 14 || J == 15 || J == 16)) goto L6019;
            if(NPROF == _SplitPearsonVII && (J == 14 || J == 15 || J == 16)) goto L6019;
            KM=GLB_[J-1].L;
            if(KM != 0) DUMMY[KM] = GLB_[J-1];
L6019:;
        }
        for (I=1; I <= MAXSX; ++I) file8o << DUMMY[I];
        for (I=MSZ+1; I <= MSZ+MAXSX; ++I) file8o << DUMMY[I];
        for (I=2*MSZ+1; I <= 2*MSZ+4; ++I) file8o << DUMMY[I];
        for (I=1; I <= 4; ++I) if ((I % 2) == 1) ILOC_ = ILOC_ + 1;
        //L9487:
        FINAL_[ILOC_][2-(I % 2)] = DUMMY[2*MSZ+I];
    }
    if (MAXS == 0 && MCYCLE == 1) return;
    ILOC_ = 0;
    for (I=1; I <= NFINAL; ++I)
    {
        for (J=1; J <= 2; ++J) FINAL_[I][J] = 0.0;
    }
    file6 << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CYCLE NUMBER=" << setw(4) << ICYCLE << endl;
    for (IP=1; IP <= NPHASE; ++IP)
    {
        IOF=0;
        if (IP > 1)
        {
            for (IIPHAS=2; IIPHAS <= IP; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
        }
        file6
            << "PHASE " << setw(2) << IP << ": " << setw(50) << phases[IP].name << endl
            << "NEW PARAMETERS, SHIFTS, AND STANDARD DEVIATIONS=" << endl << endl
            << "ATOM     X        DX       SX        Y        DY       SY        Z        DZ       SZ        B      DB     SB      So       DSo      SSo" << endl;
        N=phases[IP].AtomCount;
        TMASSA[IP]=0.0;
        STMASSA[IP]=0.0;
        for (I=1; I <= N; ++I)
        {
            for (J=1; J <= 5; ++J)
            {
                ILOC_ = ILOC_ + 1;
                KM=XL_[I+IOF][J].L;
                if(KM != 0)goto L6;
                SY[J]=0.0;
                SZ[J]=0.0;
                //  !cp jun 97 start
                if(J == 5)
                {
                    FINAL_[ILOC_][1] = XL_[I+IOF][J]*phases[IP].XMLTP/MURT[I+IOF];
                }
                else
                {
                    FINAL_[ILOC_][1] = XL_[I+IOF][J];
                }
                //  !cp jun 97 stop
                goto L5;
L6:
                SZ[J]=sqrt(abs(RJAC[KM][KM]));
                SY[J]=VX[KM]*XL_[I+IOF][J].codeword*RELAX[1];
                //  !cp jun 96 start
                if(J == 5)
                {
                    DUMMY[KM] = XL_[I+IOF][J]*phases[IP].XMLTP/MURT[I+IOF];
                    XL_[I+IOF][J]=XL_[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ]  = SY[J]*phases[IP].XMLTP/MURT[I+IOF];
                    FINAL_[ILOC_][1] = XL_[I+IOF][J]*phases[IP].XMLTP/MURT[I+IOF];
                    FINAL_[ILOC_][2] = SZ[J]*phases[IP].XMLTP/MURT[I+IOF];
                }
                else
                {
                    // !cp jun 96 stop
                    DUMMY[KM]  =XL_[I+IOF][J];
                    XL_[I+IOF][J]=XL_[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ] = SY[J];
                    FINAL_[ILOC_][1] = XL_[I+IOF][J];
                    FINAL_[ILOC_][2] = SZ[J];
                    SMASS[IP]=SZ[5];
                }
L5:;
            }
            // !cp jun 96 start
            //L4:
            file6
                << setw(4) << ATEXT[I+IOF] << "  "
                << setw(9) << setprecision(6) << XL_[I+IOF][1]
                << setw(9) << setprecision(6) << SY[1]
                << setw(9) << setprecision(6) << SZ[1]
                << setw(9) << setprecision(6) << XL_[I+IOF][2]
                << setw(9) << setprecision(6) << SY[2]
                << setw(9) << setprecision(6) << SZ[2]
                << setw(9) << setprecision(6) << XL_[I+IOF][3]
                << setw(9) << setprecision(6) << SY[3]
                << setw(9) << setprecision(6) << SZ[3]
                << setw(9) << setprecision(6) << XL_[I+IOF][4]
                << setw(9) << setprecision(6) << SY[4]
                << setw(9) << setprecision(6) << SZ[4]
                << setw(9) << setprecision(6) << XL_[I+IOF][5]*phases[IP].XMLTP/MURT[I+IOF]
                << setw(9) << setprecision(6) << SY[5]*phases[IP].XMLTP/MURT[I+IOF]
                << setw(9) << setprecision(6) << SZ[5]*phases[IP].XMLTP/MURT[I+IOF] << endl;
        }
        file6
            << "ATOM       B11      DB11      SB11       B22      DB22      SB22       B33      DB33      SB33" << endl
            << "           B12      DB12      SB12       B13      DB13      SB13       B23      DB23      SB23" << endl;
        for (I=1; I <= N; ++I)
        {
            for (J=6; J <= 11; ++J)
            {
                ILOC_ = ILOC_ + 1;
                KM=XL_[I+IOF][J].L;
                if(KM != 0)goto L9;
                SY[J]=0.0;
                SZ[J]=0.0;
                FINAL_[ILOC_][1] = XL_[I+IOF][J];
                goto L8;
L9:
                SZ[J]=sqrt(abs(RJAC[KM][KM]));
                SY[J]=VX[KM]*XL_[I+IOF][J].codeword*RELAX[2];
                //  !cp jun 96 start
                if(J ==  5)
                {
                    DUMMY[KM] = XL_[I+IOF][J]*phases[IP].XMLTP/MURT[I+IOF];
                    XL_[I+IOF][J]=XL_[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ] = SY[J]*phases[IP].XMLTP/MURT[I+IOF];
                    FINAL_[ILOC_][1] = XL_[I+IOF][J]*phases[IP].XMLTP/MURT[I+IOF];
                    FINAL_[ILOC_][2] = SZ[J]*phases[IP].XMLTP/MURT[I+IOF];
                }
                else
                {
                    DUMMY[KM]  =XL_[I+IOF][J];
                    XL_[I+IOF][J]=XL_[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ] = SY[J];
                    FINAL_[ILOC_][1] = XL_[I+IOF][J];
                    FINAL_[ILOC_][2] = SZ[J];
                }
                // !cp jun 96 stop
L8:;
            }
            //L7:
            file6 << setw(4) << ATEXT[I+IOF]
                << setw(10) << setprecision(6) << XL_[I+IOF][6]
                << setw(10) << setprecision(6) << SY[6]
                << setw(10) << setprecision(6) << SZ[6]
                << setw(10) << setprecision(6) << XL_[I+IOF][7]
                << setw(10) << setprecision(6) << SY[7]
                << setw(10) << setprecision(6) << SZ[7]
                << setw(10) << setprecision(6) << XL_[I+IOF][8]
                << setw(10) << setprecision(6) << SY[8]
                << setw(10) << setprecision(6) << SZ[8] << endl
                << "    "
                << setw(10) << setprecision(6) << XL_[I+IOF][9]
                << setw(10) << setprecision(6) << SY[9]
                << setw(10) << setprecision(6) << SZ[9]
                << setw(10) << setprecision(6) << XL_[I+IOF][10]
                << setw(10) << setprecision(6) << SY[10]
                << setw(10) << setprecision(6) << SZ[10]
                << setw(10) << setprecision(6) << XL_[I+IOF][11]
                << setw(10) << setprecision(6) << SY[11]
                << setw(10) << setprecision(6) << SZ[11] << endl;
        }
        for (J=1; J <= 27; ++J)
        {
            if ( J == 12 )
            {
                ILOC_ = ILOC_ + 1;
                FINAL_[ILOC_][1] = phases[IP].PREF[1];
                FINAL_[ILOC_][2] = phases[IP].PREF[2];
                ILOC_ = ILOC_ + 1;
                FINAL_[ILOC_][1] = phases[IP].PREF[3];
            }
            ILOC_ = ILOC_ + 1;
            if ( J == 6 ) ILOC1 = ILOC_;
            KM = phases[IP].PAR[J-1].L;
            if ( KM == 0 )
            {
                SY[J] = 0.0;
                SZ[J] = 0.0;
                FINAL_[ILOC_][1] = phases[IP].PAR[J-1];
            }
            else
            {
                SZ[J] = sqrt(abs(RJAC[KM][KM]));
                SY[J] = VX[KM]*phases[IP].PAR[J-1].codeword*RELAX[3];
                DUMMY[KM] = phases[IP].PAR[J-1];
                phases[IP].PAR[J-1] = phases[IP].PAR[J-1]+SY[J];
                DUMMY[KM+MSZ] = SY[J];
                FINAL_[ILOC_][1] = phases[IP].PAR[J-1];
                FINAL_[ILOC_][2] = SZ[J];
            }
        }
        DIRECT(DCSM,DCV,&IP);
        ILOC2 = ILOC_;
        ILOC_ = ILOC1-1;
        for (I=1; I <= 6; ++I)
        {
            ILOC_ = ILOC_ + 1;
            KM = phases[IP].PAR[I+5-1].L;
            MATCH=0;
            SY[I+5]=DCV[I]-SAVE[IP][I];
            if (KM != 0)
            {
                if (I >= 4)
                {
                    for (MMM=1; MMM <= 3; ++MMM) if (KM == phases[IP].PAR[MMM+5-1].L) MATCH=1;
                }
                if (MATCH == 0)
                {
                    DUMMY[KM] = SAVE[IP][I];
                    DUMMY[KM+MSZ] = SY[I+5];
                }
            }
            FINAL_[ILOC_][1] = DCV[I];
            FINAL_[ILOC_][2] = DCSM[I][I];
            SAVE[IP][I]=DCV[I];
        }
        ILOC_ = ILOC2;
        if(DCV[1] == DCV[2] && DCV[2] != DCV[3])
        {
            if(DCV[4] == DCV[5] && DCV[5] != 90.0)
            {
                if(DCV[6] == 90.0 || DCV[6] == 120.0)
                {
                    RATIODCV=DCV[3]/DCV[1];
                    SRATIODCV=(DCV[1]*DCSM[3][3]+DCV[3]*DCSM[1][1])/(DCV[3] * DCV[3]);
                    file6
                        << "OVERALL SCALE FACTOR=" << scientific
                        << setw(9) << setprecision(3) << phases[IP].PAR[0]
                        << setw(9) << setprecision(3) << SY[1]
                        << setw(9) << setprecision(3) << SZ[1] << endl << fixed
                        << "OVERALL TEMP. FACTOR="
                        << setw(9) << setprecision(4) << phases[IP].PAR[1]
                        << setw(9) << setprecision(4) << SY[2]
                        << setw(9) << setprecision(4) << SZ[2] << endl
                        << "CELL PARAMETERS="
                        << setw(11) << setprecision(6) << DCV[1]
                        << setw(11) << setprecision(6) << SY[1+5]
                        << setw(11) << setprecision(6) << DCSM[1][1] << endl
                        << setw(11) << setprecision(6) << DCV[2]
                        << setw(11) << setprecision(6) << SY[2+5]
                        << setw(11) << setprecision(6) << DCSM[2][2] << endl
                        << setw(11) << setprecision(6) << DCV[3]
                        << setw(11) << setprecision(6) << SY[3+5]
                        << setw(11) << setprecision(6) << DCSM[3][3] << endl
                        << setw(11) << setprecision(4) << DCV[4]
                        << setw(11) << setprecision(4) << SY[4+5]
                        << setw(11) << setprecision(4) << DCSM[4][4] << endl
                        << setw(11) << setprecision(4) << DCV[5]
                        << setw(11) << setprecision(4) << SY[5+5]
                        << setw(11) << setprecision(4) << DCSM[5][5] << endl
                        << setw(11) << setprecision(4) << DCV[6]
                        << setw(11) << setprecision(4) << SY[6+5]
                        << setw(11) << setprecision(4) << DCSM[6][6] << endl << endl
                        << "            c/a= "  << setw(11) << setprecision(6) << RATIODCV << "  +/-" << setw(11) << setprecision(6) << SRATIODCV << endl << endl
                        << "PREFERRED ORIENTATION PARAMETERS="
                        << " "
                        << setw(8) << setprecision(5) << phases[IP].PAR[11]
                        << setw(8) << setprecision(5) << SY[12]
                        << setw(8) << setprecision(5) << SZ[12]
                        << " "
                        << setw(8) << setprecision(5) << phases[IP].PAR[12]
                        << setw(8) << setprecision(5) << SY[13]
                        << setw(8) << setprecision(5) << SZ[13] << endl
                        << "ASYMMETRY PARAMETER="
                        << " "
                        << setw(8) << setprecision(4) << phases[IP].PAR[13]
                        << setw(8) << setprecision(4) << SY[14]
                        << setw(8) << setprecision(4) << SZ[14] << endl
                        << "LORENTZIAN HALF WIDTH PARAMS (X AND Y) "
                        << " "
                        << setw(8) << setprecision(5) << phases[IP].PAR[14]
                        << setw(8) << setprecision(5) << SY[15]
                        << setw(8) << setprecision(5) << SZ[15] << endl
                        << " "
                        << setw(8) << setprecision(5) << phases[IP].PAR[15]
                        << setw(8) << setprecision(5) << SY[16]
                        << setw(8) << setprecision(5) << SZ[16] << endl;
                    goto L1032;
                }
            }
        }
        file6 << "OVERALL SCALE FACTOR=" << scientific
            << setw(9) << setprecision(3) << phases[IP].PAR[0]
        << setw(9) << setprecision(3) << SY[1]
        << setw(9) << setprecision(3) << SZ[1] << fixed << endl
            << "OVERALL TEMP. FACTOR="
            << setw(9) << setprecision(4) << phases[IP].PAR[1]
        << setw(9) << setprecision(4) << SY[2]
        << setw(9) << setprecision(4) << SZ[2] << endl
            << "CELL PARAMETERS=" << endl
            << setw(11) << setprecision(6) << DCV[1]
        << setw(11) << setprecision(6) << SY[1+5]
        << setw(11) << setprecision(6) << DCSM[1][1] << endl
            << setw(11) << setprecision(6) << DCV[2]
        << setw(11) << setprecision(6) << SY[2+5]
        << setw(11) << setprecision(6) << DCSM[2][2] << endl
            << setw(11) << setprecision(6) << DCV[3]
        << setw(11) << setprecision(6) << SY[3+5]
        << setw(11) << setprecision(6) << DCSM[3][3] << endl
            << setw(11) << setprecision(4) << DCV[4]
        << setw(11) << setprecision(4) << SY[4+5]
        << setw(11) << setprecision(4) << DCSM[4][4] << endl
            << setw(11) << setprecision(4) << DCV[5]
        << setw(11) << setprecision(4) << SY[5+5]
        << setw(11) << setprecision(4) << DCSM[5][5] << endl
            << setw(11) << setprecision(4) << DCV[6]
        << setw(11) << setprecision(4) << SY[6+5]
        << setw(11) << setprecision(4) << DCSM[6][6] << endl
            << "PREFERRED ORIENTATION PARAMETERS=" << endl
            << " "
            << setw(8) << setprecision(5) << DCV[1]
        << setw(8) << setprecision(5) << SY[1+5]
        << setw(8) << setprecision(5) << DCSM[1][1] << endl
            << " "
            << setw(8) << setprecision(5) << DCV[2]
        << setw(8) << setprecision(5) << SY[2+5]
        << setw(8) << setprecision(5) << DCSM[2][2] << endl
            << " "
            << setw(8) << setprecision(5) << DCV[3]
        << setw(8) << setprecision(5) << SY[3+5]
        << setw(8) << setprecision(5) << DCSM[3][3] << endl
            << "ASYMMETRY PARAMETER=" << endl
            << " "
            << setw(8) << setprecision(4) << DCV[4]
        << setw(8) << setprecision(4) << SY[4+5]
        << setw(8) << setprecision(4) << DCSM[4][4] << endl
            << "LORENTZIAN HALF WIDTH PARAMS (X AND Y)" << endl
            << " "
            << setw(8) << setprecision(5) << DCV[5]
        << setw(8) << setprecision(5) << SY[5+5]
        << setw(8) << setprecision(5) << DCSM[5][5] << endl
            << " "
            << setw(8) << setprecision(5) << DCV[6]
        << setw(8) << setprecision(5) << SY[6+5]
        << setw(8) << setprecision(5) << DCSM[6][6] << endl
            << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP].PAR[11]
        << setw(8) << setprecision(5) << SY[12]
        << setw(8) << setprecision(5) << SZ[12] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP].PAR[12]
        << setw(8) << setprecision(5) << SY[13]
        << setw(8) << setprecision(5) << SZ[13] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP].PAR[13]
        << setw(8) << setprecision(5) << SY[14]
        << setw(8) << setprecision(5) << SZ[14] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP].PAR[14]
        << setw(8) << setprecision(5) << SY[15]
        << setw(8) << setprecision(5) << SZ[15] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP].PAR[15]
        << setw(8) << setprecision(5) << SY[16]
        << setw(8) << setprecision(5) << SZ[16] << endl;
L1032:
        //-----CHECK FOR THE SPLIT PEARSON VII PROFILE
        if(NPROF == _SplitPearsonVII)
        {
            file6 << "LOW SIDE EXPONENT PARAMETERS (NA, NB, NC)=" << endl
                //<< scientific
                << " "
                << setw(10) << setprecision(4) << phases[IP].PAR[16]
                << setw(10) << setprecision(4) << SY[17]
                << setw(10) << setprecision(4) << SZ[17] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP].PAR[17]
                << setw(10) << setprecision(4) << SY[18]
                << setw(10) << setprecision(4) << SZ[18] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP].PAR[18]
                << setw(10) << setprecision(4) << SY[19]
                << setw(10) << setprecision(4) << SZ[19] << endl
                << "HIGH SIDE EXPONENT PARAMETERS (NA, NB, NC)=" << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP].PAR[23]
                << setw(10) << setprecision(4) << SY[24]
                << setw(10) << setprecision(4) << SZ[24] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP].PAR[24]
                << setw(10) << setprecision(4) << SY[25]
                << setw(10) << setprecision(4) << SZ[25] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP].PAR[25]
                << setw(10) << setprecision(4) << SY[26]
                << setw(10) << setprecision(4) << SZ[26] << endl
                << fixed
                << "SPLIT PEARSON VII ASSYMETRY PARAMETER="
                << " "
                << setw(10) << setprecision(6) << phases[IP].PAR[26]
                << setw(10) << setprecision(6) << SY[27]
                << setw(10) << setprecision(6) << SZ[27] << endl;
        }
        else
        {
            //-----if NOT A SPLIT PEARSON VII PROFILE
            file6 << "MIXING PARAMETERS" << endl
                << "NA= "
                //<< scientific
                << setw(10) << setprecision(3) << phases[IP].PAR[16]
                << setw(10) << setprecision(3) << SY[17]
                << setw(10) << setprecision(3) << SZ[17] << endl
                << "NB= "
                << setw(10) << setprecision(3) << phases[IP].PAR[17]
                << setw(10) << setprecision(3) << SY[18]
                << setw(10) << setprecision(3) << SZ[18] << endl
                << "NC= "
                << setw(10) << setprecision(3) << phases[IP].PAR[18]
                << setw(10) << setprecision(3) << SY[19]
                << setw(10) << setprecision(3) << SZ[19] << fixed << endl;
        }
        file6 << "FWHM PARAMETERS=" << endl
            << "U = "
            << setw(10) << setprecision(6) << phases[IP].PAR[2]
            << setw(10) << setprecision(6) << SY[3]
            << setw(10) << setprecision(6) << SZ[3] << endl
            << "V = "
            << setw(10) << setprecision(6) << phases[IP].PAR[3]
            << setw(10) << setprecision(6) << SY[4]
            << setw(10) << setprecision(6) << SZ[4] << endl
            << "W = "
            << setw(10) << setprecision(6) << phases[IP].PAR[4]
            << setw(10) << setprecision(6) << SY[5]
            << setw(10) << setprecision(6) << SZ[5] << endl
            << "CT= "
            << setw(10) << setprecision(6) << phases[IP].PAR[20]
            << setw(10) << setprecision(6) << SY[21]
            << setw(10) << setprecision(6) << SZ[21] << endl
            << "Z = "
            << setw(10) << setprecision(6) << phases[IP].PAR[19]
            << setw(10) << setprecision(6) << SY[20]
            << setw(10) << setprecision(6) << SZ[20] << endl;

        //-----Modification introduced by Carlos O. Paiva-Santos to perform
        //-----Quantitative phase analysis, 03/94. Added 05/94, T.S. Moss
        //-----CHANGES TO INCORPORATE THE REFINED OCCUPANCY. Paiva-Santos (Feb-Mar/95)
        for (I=1; I <= N; ++I)
        {
            ICOCO=PTR_[I+IOF];
            TMASSA[IP] = TMASSA[IP] + XL_[I+IOF][5]*coeff.XMAS[ICOCO]*phases[IP].XMLTP;
            STMASSA[IP] = STMASSA[IP] + SMASS[IP]*coeff.XMAS[ICOCO]*phases[IP].XMLTP;
        }
        XFAC = M_PI / 180.000000;
        DCV[4] = XFAC * DCV[4];
        DCSM[4][4] = DCSM[4][4] * XFAC;
        DCV[5] = XFAC * DCV[5];
        DCSM[5][5] = DCSM[5][5] * XFAC;
        DCV[6] = DCV[6] * XFAC;
        DCSM[6][6] = DCSM[6][6] * XFAC;
        //-----Calculations of VOLUME and SVZM (=W) for each phase
        //-----and respectives standard deviations
        //-----New standard deviation code introduced in nov 96 !cp
        ARGCOS= 1-pow((cos(DCV[4])),2)-pow((cos(DCV[5])),2)-pow((cos(DCV[6])),2) + 2 * (cos(DCV[4])) * (cos(DCV[5])) * (cos(DCV[6]));
        V0 = DCV[1] * DCV[2] * DCV[3];
        VOL[IP] = V0 * sqrt(ARGCOS);
        VOSQ = 0.5*VOL[IP]/ARGCOS;
        ARG1 = VOSQ*(2 * cos(DCV[4]) * sin(DCV[4]) - 2*sin(DCV[4]) *cos(DCV[5]) *cos(DCV[6])) * DCSM[4][4];
        ARG2 = VOSQ*(2 * cos(DCV[5]) * sin(DCV[5]) - 2*sin(DCV[5]) *cos(DCV[4]) *cos(DCV[6])) * DCSM[5][5];
        ARG3 = VOSQ*(2 * cos(DCV[6]) * sin(DCV[6]) - 2*sin(DCV[6]) *cos(DCV[4]) *cos(DCV[5])) * DCSM[6][6];
        DVOL[IP] = sqrt(pow((VOL[IP] * DCSM[1][1] / DCV[1]),2) + pow((VOL[IP] * DCSM[2][2] / DCV[2]),2) + pow((VOL[IP] * DCSM[3][3] / DCV[3]),2) + pow(ARG1,2) + pow(ARG2,2) + pow(ARG3,2));
        // standard deviations are calculed below                      !cp nov 96
        W[IP] = phases[IP].PAR[0] * TMASSA[IP] * VOL[IP]/phases[IP].NMOL;
        DW[IP] = (SZ[1]/phases[IP].PAR[0]) + (DVOL[IP]/VOL[IP]) + (STMASSA[IP]/TMASSA[IP])/phases[IP].SAQF;
        //   end of std
        file6 << "Volume= "
            << setw(9) << setprecision(3) << VOL[IP]
        << "(+/-)"
            << setw(7) << setprecision(3) << DVOL[IP]
        << " UCW= "
            << setw(7) << setprecision(2) << TMASSA[IP]
        << " U.C.Density = "
            << setw(7) << setprecision(3) << 1.66113*TMASSA[IP]/VOL[IP]
        << " gr/cm^3" << endl
            << "                     _____________________________" << endl;
        //L10:;
    }
    // ****** QUANTITATIVE ANALYSIS ***************
    WTOTAL = 0.000000;
    //DWTOTAL = 0.000000;
    for (I = 1; I <= NPHASE; ++I) WTOTAL = WTOTAL + W[I];
    for (I = 1; I <= NPHASE; ++I)
    {
        //      !cp nov 10 96
        if (NPHASE == 1)
        {
            XMASS[I] = 100.0 * W[I] / WTOTAL;
            DMASS[1]= 0.0;
        }
        else
        {
            XMASS[I] = 100.0 * W[I] / WTOTAL;
            DMASS[I] = 100.0 * DW[I];
        }
    }
    IINNMOL = 0;
    for (I = 1; I <= NPHASE; ++I)
    {
        if (IINNMOL == 1) goto L2713;
        if (phases[I].NMOL == 0) IINNMOL=1;
L2713:;
    }
    if (IINNMOL == 1)
    {
        for (I = 1; I <= NPHASE; ++I)
        {
            // ** printing results
            file6 << "PHASE = "
                << setw(2) << I
                << " => %MASS = "
                << setw(6) << setprecision(2) << XMASS[I]
            << "(+/-)"
                << setw(6) << setprecision(2) << DMASS[I]
            << "  %MOLAR = NOT COMPUTED" << endl;
        }
    }
    else
    {
        // ****    CALCULATION OF MOLAR FRACTION  ****
        FT = 0.0000000;
        for (I = 1; I <= NPHASE; ++I)
        {
            FRP[I] = XMASS[I] * phases[I].NMOL / TMASSA[I];
            FT = FT + FRP[I];
        }
        for (I = 1; I <= NPHASE; ++I)
        {
            FR[I] = 100.0 * FRP[I] / FT;
            DFR[I] = FR[I] * DMASS[I]/XMASS[I];
            // ** printing results
            file6 << "PHASE = "
                << setw(2) << I
                << " => %MASS = "
                << setw(6) << setprecision(2) << XMASS[I]
            << "(+/-)"
                << setw(6) << setprecision(2) << DMASS[I]
            << "  %MOLAR = "
                << setw(6) << setprecision(2) << FR[I]
            << "(+/-)"
                << setw(6) << setprecision(2) << DFR[I] << endl;
        }
    }
    if (ISPHASE != 0)
    {
        file6 << "Considering Amorphous Content:" << endl;
        SFIS=phases[ISPHASE].WTIS/XMASS[ISPHASE];
        if(SFIS > 1.0)
        {
            file6 << "PROBLEM:Amount of Internal Standard (Phase #"
                << setw(2) << ISPHASE
                << ") is less than the specified "
                << setw(6) << setprecision(2) << phases[ISPHASE].WTIS
            << "%." << endl
                << "Amorphous content not computed. Check ISWT in line 11.2 for this phase" << endl;
            goto L2720;
        }
        for (I=1; I <= NPHASE; ++I)
        {
            file6 << "PHASE = " << setw(2) << I
                << " => %MASS = " << setw(6) << setprecision(2) << XMASS[I]*SFIS << endl;
        }
        file6 << "AMORPHOUS  => %MASS = " << setw(6) << setprecision(2) << 100*(1.0-SFIS) << endl;
    }
L2720:
    file6 << endl;
    //L291:
    for (J=1; J <= 20; ++J)
    {
        ILOC_ = ILOC_ +1;
        KM=GLB_[J-1].L;
        if(KM != 0)goto L20;
        SY[J]=0.0;
        SZ[J]=0.0;
        FINAL_[ILOC_][1] = GLB_[J-1];
        goto L19;
L20:
        SZ[J]=sqrt(abs(RJAC[KM][KM]));
        SY[J]=VX[KM]*GLB_[J-1].codeword*RELAX[4];
        DUMMY[KM] = GLB_[J-1];
        GLB_[J-1]=GLB_[J-1]+SY[J];
        DUMMY[KM+MSZ]  = SY[J];
        FINAL_[ILOC_][1] = GLB_[J-1];
        FINAL_[ILOC_][2] = SZ[J];
L19:;
    }
    if (ILOC_  >  NFINAL) DBWSException("Parameter NFINAL in PARAM.INC file too small");
    file6 << "GLOBAL PARAMETERS" << endl
        << "ZEROPOINT (ZER)             :"
        << setw(8) << setprecision(4) << GLB_[1-1]
    << setw(8) << setprecision(4) << SY[1]
    << setw(8) << setprecision(4) << SZ[1] << endl;
    file6 << "SAMPLE DISPLACEMENT (DISP)  :"
        << setw(8) << setprecision(4) << GLB_[10-1]
    << setw(8) << setprecision(4) << SY[10]
    << setw(8) << setprecision(4) << SZ[10] << endl
        << "SAMPLE TRANSPARENCY (TRANSP):"
        << setw(8) << setprecision(4) << GLB_[11-1]
    << setw(8) << setprecision(4) << SY[11]
    << setw(8) << setprecision(4) << SZ[11] << endl
        << "ROUGHNESS PARAMETERS        :"
        << "             P              :"
        << setw(8) << setprecision(4) << GLB_[8-1]
    << setw(8) << setprecision(4) << SY[8]
    << setw(8) << setprecision(4) << SZ[8] << endl
        << "             Q              :"
        << setw(8) << setprecision(4) << GLB_[9-1]
    << setw(8) << setprecision(4) << SY[9]
    << setw(8) << setprecision(4) << SZ[9] << endl
        << "             R              :"
        << setw(8) << setprecision(4) << GLB_[12-1]
    << setw(8) << setprecision(4) << SY[12]
    << setw(8) << setprecision(4) << SZ[12] << endl
        << "             T              :"
        << setw(8) << setprecision(4) << GLB_[13-1]
    << setw(8) << setprecision(4) << SY[13]
    << setw(8) << setprecision(4) << SZ[13] << endl;
    // !cp ap 20 97  !from It. codes
    file6 << "AMORPHOUS SCALE (SCAM):"
        << setw(11) << setprecision(4) << GLB_[20-1]
    << setw(11) << setprecision(4) << SY[20]
    << setw(11) << setprecision(4) << SZ[20] << endl;
    // !cp ap 20 97   !from It. codes
    file6 << "MONOCROMATOR BANDPASS PARAMETERS (PMONI)" << endl
        << setw(8) << setprecision(4) << GLB_[18-1]
    << setw(8) << setprecision(4) << SY[18]
    << setw(8) << setprecision(4) << SZ[18]
    << setw(8) << setprecision(4) << GLB_[19-1]
    << setw(8) << setprecision(4) << SY[19]
    << setw(8) << setprecision(4) << SZ[19] << endl;
    if (diffractogram.NBCKGD != 0) goto L90;
    file6 << "BACKGROUND PARAMETERS" << endl //<< scientific
        << setw(12) << setprecision(6) << GLB_[2-1]
    << setw(12) << setprecision(6) << SY[2]
    << setw(12) << setprecision(6) << SZ[2] << endl
        << setw(12) << setprecision(6) << GLB_[3-1]
    << setw(12) << setprecision(6) << SY[3]
    << setw(12) << setprecision(6) << SZ[3] << endl
        << setw(12) << setprecision(6) << GLB_[4-1]
    << setw(12) << setprecision(6) << SY[4]
    << setw(12) << setprecision(6) << SZ[4] << endl
        << setw(12) << setprecision(6) << GLB_[5-1]
    << setw(12) << setprecision(6) << SY[5]
    << setw(12) << setprecision(6) << SZ[5] << endl
        << setw(12) << setprecision(6) << GLB_[6-1]
    << setw(12) << setprecision(6) << SY[6]
    << setw(12) << setprecision(6) << SZ[6] << endl
        << setw(12) << setprecision(6) << GLB_[7-1]
    << setw(12) << setprecision(6) << SY[7]
    << setw(12) << setprecision(6) << SZ[7] << endl;
L90:
    for (I=1; I <= MAXS; ++I) file8o << DUMMY[I];
    for (I=MSZ+1; I <= MSZ+MAXS; ++I) file8o << DUMMY[I];
    for (I=2*MSZ+1; I <= 2*MSZ+4; ++I) file8o << DUMMY[I];
    return;
    //L99990:
    //	DBWSException("ERROR IN WRITING TO UNIT 8 IN OUTPTR");
}

void DBWS::ITER(void)
{
    int I,J,K,L,IA,IX,IOF,IPM,LIM,IPREV,IDERIV,IORDLIM2;
    int KS[50+1];
    double X,X1,YY,SDK,DYC,OOM,STH,THX,CISK,FPOL,TOTAS,TOTCS,TOTDS;
    //double CSK[99+1];
    //double DISK[99+1];
    //double DYCDD[99+1];
    double VLAST[MSZ+1];
    double DERISO[NATS+1];
    double ISODER[NATS+1];


    for (I=1; I <= MSZ; ++I) VLAST[I]=0.0;

    //     **** start cycles ****
    for (IX=1; IX <= MCYCLE; ++ IX)
    {
        ASSIGN_();
        for (I=1; I <= MSZ; ++I)
        {
            VX[I]=0.0;
            for (J=1; J <= MSZ; ++J) RJAC[I][J]=0.0;
        }
        for (I=1; I <= NATS; ++I)
        {
            if (IBGD == 1)
            {
                ISODER[I] = 1;
            }
            else
            {
                ISODER[I]=0;
            }
        }
        _NUM=0;
        IPREV=1;
        if(diffractogram.NBCKGD == 0)
        {
            _TH=diffractogram.THMIN-diffractogram.STEP;
            for (I=1; I <= diffractogram.NPTS; ++I)
            {
                _TH=_TH+diffractogram.STEP;
                THX=_TH/BKPOS_-1.0;
                diffractogram.BK[I]=GLB_[2-1];
                for (J=2; J <= 6; ++J) diffractogram.BK[I]=diffractogram.BK[I]+GLB_[J+1-1] * pow(THX,J-1);
                if(MCYCLE == 1 && MAXS == 0 && IPLPOL == 1) diffractogram.BKPOL[I]=diffractogram.BK[I];
            }
        }

        _TH=diffractogram.THMIN-diffractogram.STEP;

        //           ****** START GREAT LOOP ON NPTS STEP-POINTS ******
        for (IPM=1; IPM <= diffractogram.NPTS; ++IPM)
        {
            if ( (IPM % (diffractogram.NPTS/80+1)) == 0) cout << ".";
            _TH=_TH+diffractogram.STEP;
            // !cp jun 97 test for bck option
            if (IBGD == 1)
            {
                phases[1].CSK   = 1;
                phases[1].DISK  = 1;
                phases[1].DYCDD = 1;
                TOTCS  = 1;
                goto L4334;
            }
            //-----COMPUTE FPOL = POLARIZATION FACTOR FOR THIS IPM POINT
            //     TH  = 2.0 * THETA
            //     STH = SIN(THETA)
            STH = sin(_TH * 0.008726646);
            FPOL = (1.0 + pow(1.0 - 2.0 * STH * STH,2) * CTHM_ ) / 2.0;
            if(diffractogram.KR[IPM] == IRS)
            {
                // ----IDERIV = 1 MEANS ONLY BACKGROUND (FOR THE STEP-POINTS THAT ARE IN EXCLUDED REGIONS)
                IDERIV = 1;
            }
            else
            {
                //-----IDERIV = 2 MEANS BACKGROUND + DERIVATIVES (FOR STEP-POINTS THAT ARE IN INCLUDED REGIONS)
                IDERIV = 2;
            }
            //-----if NECESSARY COMPUTE TOTCS :
            //     TOTCS = TOTAL COMPTON INTENSITY SCATTERED BY ALL CRYSTALLINE PHASES AT THE IPM-TH POINT OF THE X-RAY PATTERN.
            TOTCS = 0.0;
            if ( FONDO == 1  ||  FONDO == 2 )
            {
                for (K = 1; K <= NPHASE; ++K)
                {
                    COMPTON(K,STH,&CISK);
                    phases[K].CSK = CISK * FPOL;
                    TOTCS = TOTCS + phases[K].SCABKG * phases[K].CSK;
                }
                //-----COMPTON UPDATE
                if ( MCYCLE == 1  &&  MAXS == 0  &&  IPLCOM == 1 )   diffractogram.BKCOM[IPM] = TOTCS;
                diffractogram.BK[IPM] = diffractogram.BK[IPM] + TOTCS;
            }
            //-----if NECESSARY CALL DIS = SUBROUTINE TO COMPUTE THERMAL AND LATTICE
            //                             DISORDER SCATTERING SDK  AND  DERIVATIVES
            //                             DYC IN THE K-TH PHASE
            //     TOTDS = TOTAL DISORDER SCATTERING
            TOTDS = 0.0;
            if (FONDO == 1 || FONDO == 2)
            {
                for (K = 1; K <= NPHASE; ++K)
                {
                    IOF = 0;
                    if(K > 1)
                    {
                        for (I = 2; I <= K; ++I) IOF = IOF + phases[I-1].AtomCount;
                    }
                    for (I=1; I <= NATS; ++I) DERISO[I]=0;
                    DISORDER(K,STH,IDERIV,&SDK,&DYC,FONDO,DERISO);
                    //int K, double STH, int IDERIV, double* SDK, double* DYC, int FONDO, double DERISO[])
                    phases[K].DISK  =   SDK * FPOL;
                    TOTDS    = TOTDS + phases[K].SCABKG * phases[K].DISK;
                    //     UPDATING DERIVATE OF ISOTROPIC THERMAL FACTORS
                    if(FONDO == 1)
                    {
                        for (I=1; I <= phases[K].AtomCount; ++I) ISODER[IOF+I]= ISODER[IOF+I]+phases[K].SCABKG*FPOL*DERISO[IOF+I];
                    }
                    //     UPDATING DERIVATE OF OVERALL THERMAL FACTOR
                    if(FONDO == 2) phases[K].DYCDD=DYC*phases[K].SCABKG*FPOL;
                }
                //-----DISORDER UPDATE
                if(MCYCLE == 1  &&  MAXS == 0  &&  IPLDIS == 1 ) diffractogram.BKDIS[IPM] = TOTDS;
                diffractogram.BK[IPM] = diffractogram.BK[IPM] + TOTDS;
            }
            //------AMORPHOUS EVALUATIONS
            TOTAS =0.0;
            TOTAS = GLB_[20-1] * diffractogram.AMORPHOUS[IPM];
            //-----AMORPHOUS UPDATE
            if ( MCYCLE == 1  &&  MAXS == 0  &&  IPLAM == 1 ) diffractogram.BKAM[IPM] = TOTAS;
            diffractogram.BK[IPM] = diffractogram.BK[IPM] + TOTAS;
L4334:
            if(diffractogram.KR[IPM] != IRS )
            {
                _IORD1=diffractogram.KR[IPM] % IRS;
                _IORD2=(diffractogram.KR[IPM]/IRS) % IRS;
                if ( !((_IORD2 == 0  ||  _IORD1 == 0)  &&  diffractogram.NBCKGD != 0) )
                {
                    _NUM=_NUM+1;
                    IORDLIM2=_IORD2;
                    if(IPREV <= IORDLIM2)
                    {
                        for (J=IPREV; J <= IORDLIM2; ++J) CALCUL(J);
                    }
                    IPREV=max(IPREV,IORDLIM2+1);
                    SUMMAT(IPM,ISODER,TOTCS);
                }
            }
        }
        cout << ":" << endl;
        if ( JOBTYP > 2 ) return;
        for (I=1; I <= MAXS; ++I)
        {
            for (J=I; J <= MAXS; ++J) RJAC[J][I] = RJAC[I][J];
        }
        _COND = DPINV(RJAC,VX,&MAXS);
        CHISQ();
        OOM=diffractogram.S2_ / static_cast<double>(_NUM-MAXS);
        for (I=1; I <= MAXS; ++I)
        {
            for (J=1; J <= MAXS; ++J) RJAC[I][J] = RJAC[I][J]*OOM;
        }
        //              CODE TO ATTEMPT TO STABILIZE OSCILLATIONS
        for (I=1; I <= MAXS; ++I)
        {
            if( sign(VX[I]) == sign(VLAST[I])) goto L60;
            if(abs(VX[I]) > 1.2*abs(VLAST[I])) goto L60;
            if(abs(VX[I]) < 0.8*abs(VLAST[I])) goto L60;
            VX[I]=VX[I]/2.0;
L60:
            VLAST[I]=VX[I];
        }
        OUTPTR(IX);
        for (I=1; I <= MAXS; ++I)
        {
            X1=sqrt(abs(RJAC[I][I]))*EPS;
            if(abs(VX[I]) > X1)goto L10;
        }
        if (MAXS > 0)
        {
            file6 << endl
                << "          ***** EPSED OUT *****" << endl;
        }
        if (MAXS > 0) ICYRUN = IX;
        goto L20;
L10:;
    }
    //                 ***** END CYCLES *****
L20:
    //     CODE FOR PRINTING CORRELATION MATRIX MOVED FROM EXPUT (65-85)
    if(MAT == 0 || MAXS == 0) return;
    file6 << endl << endl << "CORRELATION MATRIX=" << endl;
    IA=1;
    LIM=19;
L38:
    LIM=min(MAXS,LIM);
    for (I=IA; I <= LIM; ++I) file6 << setw(6) << I;
    file6 << endl;
    for (I=1; I <= MAXS; ++I)
    {
        L=0;
        X=RJAC[I][I];
        for (J=IA; J <= LIM; ++J)
        {
            L=L+1;
            YY=RJAC[J][J]*X;
            KS[L]= static_cast<int>( 100.0*RJAC[I][J]/(sqrt(abs(YY)) * sign(YY) )+0.5 );
        }
        file6 << setw(5) << I;
        for (J=1; J <= L; ++J) file6 << setw(6) << KS[J];
        file6 << endl;
    }
    if(LIM >= MAXS) goto L37;
    IA=LIM+1;
    LIM=LIM+19;
    goto L38;
L37:;
}

void DBWS::CEL000(int* MULTX_, double Y_[][3+1])
{
    //
    //   THIS SR IS USED TO FORCE A POINT XYZ TO LIE IN THE UNIT CELL AT 000.
    //

    Y_[*MULTX_-1][1] = fmod(Y_[*MULTX_-1][1]+8.0,1.0);
    if ( Y_[*MULTX_-1][1]-0.99999 >= 0 )
    {
        Y_[*MULTX_-1][1] = 0.0;
    }
    Y_[*MULTX_-1][2] = fmod(Y_[*MULTX_-1][2]+8.0,1.0);
    if ( Y_[*MULTX_-1][2]-0.99999 >= 0 )
    {
        Y_[*MULTX_-1][2] = 0.0;
    }
    Y_[*MULTX_-1][3] = fmod(Y_[*MULTX_-1][3]+8.0,1.0);
    if ( Y_[*MULTX_-1][3]-0.99999 >= 0 )
    {
        Y_[*MULTX_-1][3] = 0.0;
    }
}

void DBWS::OPERTR(int* MULTX_, double Y_[][3+1], int* IPHASE, int* I, int* L2)
{
    //
    //                             THIS SR GENERATES THE SYMMETRY EQUIVALENT
    //                           POSITIONS. IT IS CALLED BY SYMOPR.
    //

    int I1,I2,I3,I4,I5,I6,I7,I8;

    I1 = (phases[*IPHASE].SYMB.NCONT[*I] % 4)+1;
    I2 = ((phases[*IPHASE].SYMB.NCONT[*I]/4) % 4)+1;
    I3 = ((phases[*IPHASE].SYMB.NCONT[*I]/16) % 4)+1;
    I4 = ((phases[*IPHASE].SYMB.NCONT[*I]/64) % 4)+1;
    I5 = ((phases[*IPHASE].SYMB.NCONT[*I]/256) % 2);
    I6 = ((phases[*IPHASE].SYMB.NCONT[*I]/512) % 4)+1;
    I7 = ((phases[*IPHASE].SYMB.NCONT[*I]/2048) % 4)+1;
    I8 = phases[*IPHASE].SYMB.NCONT[*I]/8192+1;

    switch (I8) {
    case 1:
        switch (I1) {
        case 1:
            Y_[*MULTX_][1] = Y_[*L2][1];
            break;
        case 2:
            Y_[*MULTX_][1] = Y_[*L2][2];
            break;
        case 3:
            Y_[*MULTX_][1] = -Y_[*L2][1];
            break;
        case 4:
            Y_[*MULTX_][1] = -Y_[*L2][2];
            break;
        default:
            GOTOER();
            break;
        }
        switch (I2) {
        case 1:
            break;
        case 2:
            Y_[*MULTX_][1] = Y_[*MULTX_][1]+0.5;
            break;
        case 3:
            Y_[*MULTX_][1] = Y_[*MULTX_][1]+0.25;
            break;
        case 4:
            Y_[*MULTX_][1] = Y_[*MULTX_][1]+0.75;
            break;
        default:
            GOTOER();
            break;
        }
        switch (I3) {
        case 1:
            Y_[*MULTX_][2] = Y_[*L2][2];
            break;
        case 2:
            Y_[*MULTX_][2] = Y_[*L2][1];
            break;
        case 3:
            Y_[*MULTX_][2] = -Y_[*L2][2];
            break;
        case 4:
            Y_[*MULTX_][2] = -Y_[*L2][1];
            break;
        default:
            GOTOER();
            break;
        }
        switch (I4) {
        case 1:
            break;
        case 2:
            Y_[*MULTX_][2] = Y_[*MULTX_][2]+0.5;
            break;
        case 3:
            Y_[*MULTX_][2] = Y_[*MULTX_][2]+0.25;
            break;
        case 4:
            Y_[*MULTX_][2] = Y_[*MULTX_][2]+0.75;
            break;
        default:
            GOTOER();
            break;
        }
        Y_[*MULTX_][3] = Y_[*L2][3] * static_cast<double>(1-2*I5);
        switch (I6) {
        case 1:
            break;
        case 2:
            Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.50;
            break;
        case 3:
            Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.25;
            break;
        case 4:
            Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.75;
            break;
        default:
            GOTOER();
            break;
        }
        switch (I7) {
        case 1:
            break;
        case 2:
            Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.333333333;
            break;
        case 3:
            Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.666666666;
            break;
        default:
            GOTOER();
            break;
        }
        *MULTX_ = *MULTX_+1;
        CEL000(MULTX_,Y_);
        return;
        break;
    case 2:
        Y_[*MULTX_][1] = -Y_[*L2][2];
        Y_[*MULTX_][2] = Y_[*L2][1]-Y_[*L2][2];
        Y_[*MULTX_][3] = Y_[*L2][3]+0.333333333 * static_cast<double>(I2-2);
        *MULTX_ = *MULTX_+1;
        CEL000(MULTX_,Y_);
        Y_[*MULTX_][1] = -Y_[*MULTX_-1][2];
        Y_[*MULTX_][2] = -Y_[*L2][1];
        Y_[*MULTX_][3] = Y_[*L2][3]-0.333333333 * static_cast<double>(I2-2);
        *MULTX_ = *MULTX_+1;
        CEL000(MULTX_,Y_);
        return;
        break;
    case 3:
        I1 = *L2;
        L301:
        Y_[*MULTX_][1] = Y_[I1][3];
        Y_[*MULTX_][2] = Y_[I1][1];
        Y_[*MULTX_][3] = Y_[I1][2];
        if ( I1-*L2 <= 0 )
        {
            I1 = *MULTX_;
            *MULTX_ = *MULTX_+1;
            CEL000(MULTX_,Y_);
            goto L301;
        }
        else
        {
            *MULTX_ = *MULTX_+1;
            CEL000(MULTX_,Y_);
            return;
        }
        break;
    default:
        GOTOER();
        break;
    }
}

void DBWS::SYMOPR(int* MULTX_, double Y_[][3+1], double* XLT_, int* IXB_, int NCTR_[], int* IPHASE)
{
    //
    //                             THIS SR TAKES A POINT XYZ, DETERMINES THE
    //                           SITE MULTIPLICITY, THE SITE SYMMETRY, AND
    //                           THE SYMMETRY CONSTRAINTS ON THE POSITIONAL
    //                           AND THERMAL PARAMETERS.
    //

    int I,J,L1,L2,L3,L4,L5,L6,KX,KY,KZ,IOP,ISCK,ISRH;
    double X3;

    ISCK = 0;

L103:
    *MULTX_ = 2;
    ISRH = 0;
    *IXB_ = 0;
    *XLT_ = 1.0;
    J = 1;
    CEL000(MULTX_,Y_);
    L1=J;
    for (I=1; I <= 8; ++I)
    {
        IOP = 0;
        if ( phases[*IPHASE].SYMB.NCONT[I] == 0 ) goto L104;
        for (L2=J; L2 <= L1; ++L2)
        {
            L3=*MULTX_-1;
            if ( phases[*IPHASE].SYMB.NCONT[I]-8192 == 0 )
            {
                Y_[*MULTX_][1]=Y_[L2][1]+1.0/3.0;
                Y_[*MULTX_][2]=Y_[L2][2]+2.0/3.0;
                Y_[*MULTX_][3]=Y_[L2][3]+2.0/3.0;
                Y_[*MULTX_+1][1]=Y_[L2][1]+2.0/3.0;
                Y_[*MULTX_+1][2]=Y_[L2][2]+1.0/3.0;
                Y_[*MULTX_+1][3]=Y_[L2][3]+1.0/3.0;
                *MULTX_ = *MULTX_+1;
                CEL000(MULTX_,Y_);
                *MULTX_ = *MULTX_+1;
                CEL000(MULTX_,Y_);
            }
            else
            {
                OPERTR(MULTX_,Y_,IPHASE,&I,&L2);
            }
            L5=*MULTX_;
            L6=phases[*IPHASE].SYMB.NCONT[I];
            *MULTX_ = L3+1;
L8:
            L4 = 1;
            if ( abs(Y_[*MULTX_][1]-Y_[L4][1])-0.0001 <= 0 )
            {
                if ( abs(Y_[*MULTX_][2]-Y_[L4][2])-0.0001 <= 0 )
                {
                    if ( abs(Y_[*MULTX_][3]-Y_[L4][3])-0.0001 <= 0 )
                    {
                        if ( NCTR_[*MULTX_]-8192 < 0)
                        {
                            goto L92;
                        }
                        else if ( NCTR_[*MULTX_]-8192 == 0)
                        {
                            *IXB_ |= 2916928;
                            *XLT_ = *XLT_/3.0;
                            goto L99;
                        }
                        else
                        {
                            if ( (NCTR_[*MULTX_] % 16384) == 0 )
                            {
                                *IXB_ |= 17924240;
                                *XLT_ = *XLT_/3.0;
                                goto L99;
                            }
                            ISRH = 1;
                            if ( (307 & (NCTR_[*MULTX_]-phases[*IPHASE].SYMB.NCONT[I])) == 0 )
                            {
                                goto L92;
                            }
                            goto L90;
                        }
                    }
                    else
                    {
                        goto L90;
                    }
                }
                else
                {
                    goto L90;
                }
            }
            else
            {
                goto L90;
            }



L92:
            KX = (NCTR_[*MULTX_] % 4)+1;
            KY = ((NCTR_[*MULTX_]/16) % 4)+1;
            KZ = (NCTR_[*MULTX_]/256) % 2;
            if ( IOP == 0 ) *XLT_=*XLT_/2.0;
            IOP = 1;
            switch (KX) {
            case 1:
                switch (KY) {
                case 1:
                    if ( KZ == 0 )
                    {
                        file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                        DBWSException("7701");
                    }
                    *IXB_ |= 2621450;
                    goto L90;
                    break;
                case 2:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                case 3:
                    if ( KZ == 0 ) *IXB_ |= 2228292;
                    if ( KZ == 1 ) *IXB_ |= 655432;
                    goto L90;
                    break;
                case 4:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                }
                GOTOER();
                break;
            case 2:
                *IXB_ |= 32768;
                switch (KY) {
                case 1:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                case 2:
                    if ( KZ == 0 ) *IXB_ |= 4194464;
                    if ( KZ == 1 ) *IXB_ |= 12583048;
                    goto L90;
                    break;
                case 3:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                case 4:
                    *IXB_ |= 2753088;
                    if ( KZ == 1 ) *IXB_ |= 8;
                    goto L90;
                    break;
                }
                GOTOER();
                break;
            case 3:
                *IXB_ |= 512;
                switch (KY) {
                case 1:
                    *IXB_ |= 131072;
                    if ( KZ == 1 )
                    {
                        *IXB_ |= 2097152;
                        *IXB_ |= 8;
                        goto L90;
                    }
                    *IXB_ |= 525312;
                    goto L90;
                    break;
                case 2:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                case 3:
                    *IXB_ |= 64;
                    if ( KZ == 1 )
                    {
                        *IXB_ |= 8;
                        goto L90;
                    }
                    *IXB_ |= 2623488;
                    goto L90;
                    break;
                case 4:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                }
                GOTOER();
                break;
            case 4:
                *IXB_ |= 32768;
                switch (KY) {
                case 1:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                case 2:
                    *IXB_ |= 2753088;
                    if ( KZ == 1 ) *IXB_ |= 8;
                    goto L90;
                    break;
                case 3:
                    file6 << "THE SYMMETRY OPERATOR " << NCTR_[*MULTX_] << " IS WRONG" << endl;
                    DBWSException("7701");
                    break;
                case 4:
                    if ( KZ == 0 ) *IXB_ |= 8388864;
                    if ( KZ == 1 ) *IXB_ |= 4194568;
                    goto L90;
                    break;
                }
                GOTOER();
                break;
            }
            GOTOER();




L90:
            if ( L6-8192 < 0 )
            {
                goto L99;
            }
            else
            {
                L6=1;
                ++*MULTX_;
                goto L8;
            }

L99:
            *MULTX_ = L5;
        }
        L1=*MULTX_-1;
    }
L104:
    if ( ISRH == 0 ) goto L200;
    if ( ISCK > 1 ) goto L2000;
    ++ISCK;
    for (J=1; J <= I; ++J) if ( phases[*IPHASE].SYMB.NCONT[J] == 16384 )
    {
        if ( abs(Y_[1][1]-Y_[1][2]) < 0.0001 ) goto L120;
        file6 << "THE ATOM AT "
            << setw(8) << setprecision(5) << Y_[1][1]
        << setw(8) << setprecision(5) << Y_[1][2]
        << setw(8) << setprecision(5) << Y_[1][3]
        << " WAS PUT AT"
            << setw(8) << setprecision(5) << Y_[1][2]
        << setw(8) << setprecision(5) << Y_[1][3]
        << setw(8) << setprecision(5) << Y_[1][1] << endl;
        X3 = Y_[1][1];
        Y_[1][1] = Y_[1][2];
        Y_[1][2] = Y_[1][3];
        Y_[1][3] = X3;
        goto L103;
    }
    for (L2 = 1; L2 <= L1; ++L2)
    {
        if ( abs(Y_[L2][1]-Y_[L2][2]) <= 0.0001 ) goto L102;
        if ( abs(Y_[L2][1]+Y_[L2][2]-1.0) <= 0.0002 ) goto L102;
    }
    ISCK = 2;
    file6 << "NO X,X,Z OR X,-X,Z SET WAS FOUND FOR"
        << setw(8) << setprecision(5) << Y_[1][1]
    << setw(8) << setprecision(5) << Y_[1][2]
    << setw(8) << setprecision(5) << Y_[1][3]
    << setw(8) << setprecision(5) << *XLT_ << endl;
    goto L200;



L102:
    if ( L2 == 1 ) goto L120;
    file6 << "THE ATOM AT "  << setw(8) << setprecision(5) << Y_[1][1]  << setw(8) << setprecision(5) << Y_[1][2]  << setw(8) << setprecision(5) << Y_[1][3]
          << " WAS PUT AT"   << setw(8) << setprecision(5) << Y_[L2][1] << setw(8) << setprecision(5) << Y_[L2][2] << setw(8) << setprecision(5) << Y_[L2][3] << endl;
    Y_[1][1] = Y_[L2][1];
    Y_[1][2] = Y_[L2][2];
    Y_[1][3] = Y_[L2][3];
    goto L103;

L2000:
    file6 << "THERE ARE UNRESOLVABLE PROBLEMS WITH THE POSITION SET"
        << setw(8) << setprecision(5) << Y_[1][1]
    << setw(8) << setprecision(5) << Y_[1][2]
    << setw(8) << setprecision(5) << Y_[1][3] << endl
        << "ATTEMPTS TO DEFINE THE TRUE SITE SYMMETRY WERE ABANDONED." << endl;
L120:
    if ( abs(Y_[1][1]-Y_[1][2])+abs(Y_[1][1]-Y_[1][3]) <= 0.0002 ) goto L200;
    for (L2=2; L2 <= L1; ++L2) if ( abs(Y_[L2][1]-Y_[L2][2])+abs(Y_[L2][1]-Y_[L2][3]) <= 0.0002 ) goto L102;
L200:
    *IXB_ = max( (*IXB_ & 18874367), (*IXB_ & 16777215));
    if (  (*IXB_ & 384) == 384 ) *IXB_ |= 512;
}

void DBWS::RTMT(int* MULTX_, double Y_[][3+1], int* IPRT, int NCTR_[], int* ISIMOP_, int* IPHASE)
{
    //    THE MATRICES ARE IN THE FIRST THREE ROWS OF VCTR AND THE VECTORS
    //    ARE IN THE FOURTH ROW OF VCTR.
    //    THE MATRICES ARE PACKED INTO IVEC.

    int I,K,IX,IY,IZ;
    double VCTR[3+1][4+1];
    double XLT_;
    int IXB_;



    Y_[1][1] = 2.0/108.0;
    Y_[1][2] = 3.0/108.0;
    Y_[1][3] = 4.0/108.0;
    SYMOPR(MULTX_, Y_,&XLT_, &IXB_,NCTR_, IPHASE);
    phases[*IPHASE].SYMB.MULTP = *MULTX_-1;
    MLTPHASE=phases[*IPHASE].SYMB.MULTP;
    for (I=1; I <= phases[*IPHASE].SYMB.MULTP; ++I)
    {
        phases[*IPHASE].IVEC_[I] = 0;
        for (K=1; K <= 3; ++K)
        {
            IX = static_cast<int>(Y_[I][K]*108.0+0.1);
            VCTR[K][4] = static_cast<double>((IX+4)/9) /12.0;
            VCTR[K][4] = fmod(VCTR[K][4],1.0);
            IY = ((IX+4) % 9)-4;
            IZ = abs(IY);
            switch (IZ) {
            case 1:
                VCTR[K][1] = -1.0 * sign(IY);
                VCTR[K][2] =  1.0 * sign(IY);
                VCTR[K][3] =  0.0;
                break;
            case 2:
                VCTR[K][1] = 1.0 * sign(IY);
                VCTR[K][2] = 0.0;
                VCTR[K][3] = 0.0;
                break;
            case 3:
                VCTR[K][1] = 0.0;
                VCTR[K][2] = 1.0 * sign(IY);
                VCTR[K][3] = 0.00;
                break;
            case 4:
                VCTR[K][1] = 0.0;
                VCTR[K][2] = 0.0;
                VCTR[K][3] = 1.0 * sign(IY);
                break;
            default:
                GOTOER();
                break;
            }
            phases[*IPHASE].IVEC_[I] = phases[*IPHASE].IVEC_[I]+(
              9*(static_cast<int>(VCTR[K][1]+1.5))+
              3*(static_cast<int>(VCTR[K][2]+1.5))+
              static_cast<int>(VCTR[K][3]+1.5))*32768 * static_cast<int>(pow(32,3-K))+
              (static_cast<int>(VCTR[K][4]*12.0+.5)+16) * static_cast<int>(pow(32,3-K));
        }
        if ( *IPRT > 0 )
        {
            if(*ISIMOP_ != 1)
            {
                if ( (I-1 % 15) == 0) file6 << "The operations of the space group are" << endl;
                file6
                    << " ("	<< setw(3) << setprecision(0) << VCTR[1][1]	<< setw(3) << setprecision(0) << VCTR[1][2]	<< setw(3) << setprecision(0) << VCTR[1][3]	<< " )   ( X )   ("	<< setw(6) << setprecision(3) << VCTR[1][4] << " )   ( X2 )" << endl
                    << " ("	<< setw(3) << setprecision(0) << VCTR[2][1]	<< setw(3) << setprecision(0) << VCTR[2][2]	<< setw(3) << setprecision(0) << VCTR[2][3]	<< " ) * ( Y ) + (" << setw(6) << setprecision(3) << VCTR[2][4] << " ) = ( Y2 )" << endl
                    << " ("	<< setw(3) << setprecision(0) << VCTR[3][1]	<< setw(3) << setprecision(0) << VCTR[3][2]	<< setw(3) << setprecision(0) << VCTR[3][3]	<< " )   ( Z )   (" << setw(6) << setprecision(3) << VCTR[3][4] << " )   ( Z2 )"
                    << "          VEC(" << setw(3) << I << ") =" << setw(15) << phases[*IPHASE].IVEC_[I] << endl << endl;
            }

        }
    }
    phases[*IPHASE].MLTPHS=phases[*IPHASE].SYMB.MULTP;
    phases[*IPHASE].ICNTPHS=1;
    if ( phases[*IPHASE].SYMB.NC == 0 ) phases[*IPHASE].ICNTPHS=2;
}

void DBWS::OP1(int* IPHASE, int NCTR_[])
{
    //
    //                            THIS SR GENERATES THE OPERATORS FOR
    //                           SMTRY2, DETERMINES THE POLARITY OF THE SPACE
    //                           GROUP, NOTES THE PRESENCE OF A CENTER(1BAR),
    //                           AND GENERATES A SET OF PSEUDO-OPERATORS FOR
    //                           USE IN SYMOPR TO DEFINE SITE SYMMETRY.

    int I,J,K,L,N2;

    K = 2;
    phases[*IPHASE].SYMB.NPOL = 7;
    for (I=1; I <= 7; ++I) NC1[1][I][*IPHASE] = 0;
    NCTR_[1] = 0;
    for (I=1; I <= 8; ++I)
    {
        if ( phases[*IPHASE].SYMB.NCONT[I] <= 0 ) goto L101;
        if ( phases[*IPHASE].SYMB.NCONT[I] < 8192 ) goto L80;
        if ( phases[*IPHASE].SYMB.NCONT[I] <= 8192 ) goto L60;		// TODO: Possivel erro aqui estava 08192
        if ( phases[*IPHASE].SYMB.NCONT[I] == 16384 ) goto L70;
        phases[*IPHASE].SYMB.NPOL = 4;
        NC1[I][1][*IPHASE] = 2;
L50:
        for (J=2; J <= 7; ++J) NC1[I][J][*IPHASE] = 0;
        NC1[I][6][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/4) % 32;
        goto L90;
L60:
        NC1[I][1][*IPHASE] = 3;
        phases[*IPHASE].SYMB.NPOL = 1;
        goto L50;
L70:
        NC1[I][1][*IPHASE] = 1;
        goto L50;

L80:
        NC1[I][1][*IPHASE] = 0;
        NC1[I][2][*IPHASE] = phases[*IPHASE].SYMB.NCONT[I] % 4;
        NC1[I][3][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/4) % 4;
        NC1[I][4][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/16) % 4;
        NC1[I][5][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/64) % 4;
        NC1[I][6][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/256) % 2;
        NC1[I][7][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/512) % 16;
L90:
        L = K-1;
        for (J=1; J <= L; ++J)
        {
            if ( phases[*IPHASE].SYMB.NCONT[I] == 8192 ) goto L95;
            NCTR_[K] = (phases[*IPHASE].SYMB.NCONT[I] ^ NCTR_[J]) & 57651;
            if( (phases[*IPHASE].SYMB.NCONT[I] & 17) == 0) goto L92;
            if( (NCTR_[J] & 34) == 0) goto L92;
            if( (NCTR_[J] & 34) == 34) goto L92;
            NCTR_[K]= NCTR_[K] ^ 34;
L92:
            ++K;
            if ( phases[*IPHASE].SYMB.NCONT[I] < 8192 ) goto L100;
            NCTR_[K] = NCTR_[K-1]+32768;
            K = K+1;
            goto L100;
L95:
            NCTR_[K] = 0;
            NCTR_[K+1] = 0;
            K = K+2;
L100:;
        }
    }
L101:
    phases[*IPHASE].SYMB.MULTP = K-1;
    N1HKL[*IPHASE] = I-1;
    if ( N1HKL[*IPHASE] < 2 ) goto L105;
    for (N2=2; N2 <= N1HKL[*IPHASE]; ++N2) if ( NC1[N2-1][1][*IPHASE] == 1 ) goto L103;
    goto L105;
L103:
    for (I=N2; I <= N1HKL[*IPHASE]; ++I)
    {
        for (J=1; J <= 7; ++J) NC1[I-1][J][*IPHASE] = NC1[I][J][*IPHASE];
    }
    NC1[N1HKL[*IPHASE]][1][*IPHASE] = 1;
L105:
    phases[*IPHASE].SYMB.NC = 1;
    for (I=1; I <= phases[*IPHASE].SYMB.MULTP; ++I)
    {
        if ( (NCTR_[I] ^ 290) == 0 ) phases[*IPHASE].SYMB.NC = 0;
        if ( (NCTR_[I] & 3) > 0 ) phases[*IPHASE].SYMB.NPOL = (phases[*IPHASE].SYMB.NPOL & 6);
        if ( (NCTR_[I] & 48) > 0 ) phases[*IPHASE].SYMB.NPOL = (phases[*IPHASE].SYMB.NPOL & 5);
        if ( (NCTR_[I] & 256) > 0 ) phases[*IPHASE].SYMB.NPOL = (phases[*IPHASE].SYMB.NPOL & 3);
    }
}

void DBWS::LOOKUP(int K, int N, int NSCAT, int IXRAY)
{
    int I,J,L,NS,IOF,NSL,IIPHAS,TBXPTR;

    IOF=0;
    NS=0;
    if (K > 1)
    {
        for (IIPHAS=2; IIPHAS <= K; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
    }
    if (K > 1) NS=NSAVE;
    if (JOBTYP == 1 || JOBTYP == 3)
    {
        for (I=1; I <= N; ++I)
        {
            NS=max(NSCAT,NS);
            for (J=1; J <= NS; ++J)
            {
                if (NTYP[I+IOF] == NAM_[J])
                {
                    PTR_[I+IOF]=J;
                    goto L10;
                }
            }
            for (J=1; J <= 212; ++J)
            {
                if (NTYP[I+IOF] == TBXC[J])
                {
                    NS=NS+1;
                    PTR_[I+IOF]=NS;
                    for (L=1; L <= 9; ++ L) coeff.AC[L][NS]=TBX[J][L];
                    NAM_[NS]=TBXC[J];
                    TBXPTR=static_cast<int>(TBX[J][10]);			//TBXPTR=static_cast<int>(TBX[J][10]+0.5);
                    coeff.DFP[NS]=TBD[TBXPTR][IXRAY];
                    coeff.DFPP[NS]=TBD[TBXPTR][IXRAY+10];
                    coeff.XMAS[NS]=TBM[TBXPTR];
                    goto L10;
                }
            }
            file6 << " SCATTERING COEFFICIENTS NOT FOUND FOR " << NTYP[I+IOF] << endl;
            DBWSException("SCATTERING DATA MISSING");
            L10:;
        }
        NSAVE=NS;
        return;
    }


    for (I=1; I <= N; ++I)
    {
        NS=max(NSCAT,NS);
        for (J=1; J <= NS; ++J)
        {
            if(NTYP[I+IOF] == NAM_[J])
            {
                PTR_[I+IOF]=J;
                goto L80;
            }
        }
        for (J=1; J <= 85; ++J)
        {
            if(NTYP[I+IOF] == TABNC[J])
            {
                NS=NS+1;
                PTR_[I+IOF]=NS;
                coeff.DFP[NS]=TABN[J];
                NAM_[NS]=TABNC[J];
                if (J  >=  61  &&  J  <=  81)
                {
                    NSL = J+2;
                    goto L180;
                }
                if (J  ==  82)
                {
                    NSL = 90;
                    goto L180;
                }
                if (J  >=  83  &&  J  <=  85)
                {
                    NSL = J+9;
                    goto L180;
                }
                NSL = J;
                goto L180;
            }
        }
        file6 << " SCATTERING LENGTHS NOT FOUND FOR " << NTYP[I+IOF] << endl;
        //111     FORMAT(34H SCATTERING LENGTHS NOT FOUND FOR ,A4)
        DBWSException("SCATTERING DATA MISSING");


L180:
        coeff.XMAS[NS]=TBM[NSL];
L80:;
    }
    NSAVE=NS;
}

void DBWS::CELL2(int NPHASE)
{
    const string LAU[14+1] = {"",  "1BAR","2/M","MMM","4/M","4/MMM","3BAR   R","3BAR M R","3BAR","3BAR M 1","3BAR 1 M","6/M","6/MMM","M3","M3M"};
    double COSA,COSB,COSC,SINA,SINB,SINC,SINASR,SINBSR,COSASR,COSBSR,COSCSR;


    file6 << "THE LAUE SYMMETRY IS " << LAU[phases[NPHASE].SYMB.NSPGRP] << endl;
    if ( phases[NPHASE].SYMB.NAXIS > 3 ) DBWSException("5001");
    switch (phases[NPHASE].SYMB.NSPGRP) {
    case 1:
        break;
    case 2:
        if ( phases[NPHASE].SYMB.NAXIS != 1 ) cell_ALPHA=90.0;
        if ( phases[NPHASE].SYMB.NAXIS != 2 ) cell_BETA=90.0;
        if ( phases[NPHASE].SYMB.NAXIS != 3) cell_GAMMA=90.0;
        break;
    case 3:
        cell_ALPHA = 90.0;
        cell_BETA = 90.0;
        cell_GAMMA = 90.0;
        break;
    case 4:
    case 5:
        cell_B = cell_A;
        cell_ALPHA = 90.0;
        cell_BETA = 90.0;
        cell_GAMMA = 90.0;
        break;
    case 6:
    case 7:
        cell_B = cell_A;
        cell_C = cell_A;
        cell_BETA = cell_ALPHA;
        cell_GAMMA = cell_ALPHA;
        break;
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        cell_B = cell_A;
        cell_ALPHA = 90.;
        cell_BETA = 90.0;
        cell_GAMMA = 120.0;
        break;
    case 13:
    case 14:
        cell_B = cell_A;
        cell_C = cell_A;
        cell_ALPHA = 90.0;
        cell_BETA = 90.0;
        cell_GAMMA = 90.0;
        break;
    default:
        GOTOER();
        break;
    }

    COSA = cos(6.28318531*cell_ALPHA/360.0);
    //L105:
    COSB = cos(6.28318531*cell_BETA/360.0);
    //L1052:
    COSC = cos(6.28318531*cell_GAMMA/360.0);
    //L1053:
    SINA = sin(6.28318531*cell_ALPHA/360.0);		// SINA=sqrt(1.0-pow(COSA,2));
    SINB = sin(6.28318531*cell_BETA/360.0);		// SINB=sqrt(1.0-pow(COSB,2));
    SINC = sin(6.28318531*cell_GAMMA/360.0);		// SINC=sqrt(1.0-pow(COSC,2));




    //
    //-----THE VOLUME OF THE CELL IN THE DIRECT SPACE OF THE REAL CELL OF
    //-----THE K-TH PHASE=NPHASE IS CALCULATED
    //-----AND THE RAY GMAX OF THE SPHERE WITH THE VOLUME EQUIVALENT TO THE
    //-----BRILLOUIN CELL 4/3*PI*GMAX**3=1/VOL

    // Volume of cell in direct space
    phases[NPHASE].VOLI = cell_A * cell_B * cell_C * sqrt(1.0 - pow(COSA,2) - pow(COSB,2) - pow(COSC,2) + 2.0 * COSA * COSB * COSC );

    phases[NPHASE].GCOM = 0.877298169 * phases[NPHASE].VOLI / pow(diffractogram.LAMDAM,3);

    // cos alpha in space reciprocal. Giacovazzo, p67 table 2.1
    COSASR=(COSB*COSC-COSA)/(SINB*SINC);

    // cos beta in space reciprocal. Giacovazzo, p67 table 2.1
    COSBSR=(COSA*COSC-COSB)/(SINA*SINC);

    // cos gamma in space reciprocal. Giacovazzo, p67 table 2.1
    COSCSR=(COSA*COSB-COSC)/(SINA*SINB);

    // TODO: Verifica se isso esta funcionando em double
    SINASR=sqrt(1.0-pow(COSASR,2));
    SINBSR=sqrt(1.0-pow(COSBSR,2));

    // a*, a in space reciprocal. Giacovazzo, p67 table 2.1
    //
    // sin(betha*) = V/(a*b*c*sin(alpha)*sin(gamma))
    //
    // a* = b*c*sin(alpha)/V
    //
    // a* = 1/(sin(beta*)*sin(gamma)*a)
    //
    cell_AL[1][1]=1.0/(cell_A * SINC * SINBSR);

    // b*, a in space reciprocal. Giacovazzo, p67 table 2.1
    cell_AL[2][2]=1.0/(cell_B * SINC * SINASR);

    // c*, a in space reciprocal. Giacovazzo, p67 table 2.1
    cell_AL[3][3]=1.0/(cell_C * SINA * SINBSR);

    cell_AL[2][3]=cell_AL[2][2]*cell_AL[3][3]*COSASR*2.0;
    cell_AL[1][3]=cell_AL[3][3]*cell_AL[1][1]*COSBSR*2.0;
    cell_AL[1][2]=cell_AL[2][2]*cell_AL[1][1]*COSCSR*2.0;

    cell_AL[1][1]=cell_AL[1][1]*cell_AL[1][1];
    cell_AL[2][2]=cell_AL[2][2]*cell_AL[2][2];
    cell_AL[3][3]=cell_AL[3][3]*cell_AL[3][3];
}

double DBWS::MULT(int IPHASE, int IH, int IK, int IL, int KXIS)
{
    double P;

    P = 2.0;
    switch (phases[IPHASE].SYMB.NSPGRP)
    {
    case 1:
        break;
    case 2:
        if ( IH > 0 )
        {
            if ( KXIS == 3 )
            {
                if ( IL > 0  &&  IH+abs(IK) > 0 ) P=2.0*P;
            }
            else
            {
                if ( IK > 0  &&  IH+abs(IL) > 0 ) P=2.0*P;
            }
        }
        break;
    case 3:
        P = 1.0;
        if ( IH > 0 )
        {
            P = 2.0*P;
        }
        if ( IK > 0)
        {
            P = 2.0*P;
        }
        if ( IL > 0 )
        {
            P = 2.0*P;
        }
        break;
    case 4:
        if ( IK > 0 )
        {
            P = 2.0*P;
            if ( IL > 0 )
            {
                P = 2.0*P;
            }
        }
        break;
    case 5:
        if ( IK > 0)
        {
            P = 2.0*P;
            if ( IH > 0 )
            {
                if ( IH-IK < 0 )
                {
                    P = 2.0*P;
                }
            }
            if ( IL > 0 )
            {
                P = 2.0*P;
            }
        }
        break;
    case 6:
        if ( IH+IK+IL != 0 )
        {
            P = 2.0*P;
        }
        if ( IH-IK == 0 )
        {
            if ( IH-IL != 0)
            {
                P = 3.0*P;
            }
        }
        else
        {
            P = 3.0*P;
        }
        break;
    case 7:
        if ( IH-IK == 0 )
        {
            if ( IH-IL != 0 )
            {
                P = P/2.0;
                P = 3.0*P;
                if ( IH+IK+IL != 0 )
                {
                    P = 2.0*P;
                }
            }
        }
        else
        {
            if ( IK-IL == 0 )
            {
                P = P/2.0;
            }
            P = 3.0*P;
            if ( IH+IK+IL != 0 )
            {
                P = 2.0*P;
            }
        }
        break;
    case 8:
        if ( IH+abs(IK) > 0 )
        {
            P = 3.0*P;
        }
        if ( IL > 0 )
        {
            P = 2.0*P;
        }
        break;
    case 9:
        if ( IL > 0 )
        {
            P = 2.0*P;
        }
        if ( IH <= 0 )
        {
            if ( IK > 0 )
            {
                P = 3.0*P;
            }
        }
        else
        {
            P = 3.0*P;
            if ( IK > 0 )
            {
                P = 2.0*P;
            }
        }
        break;
    case 10:
        if ( IK > 0 )
        {
            P = 3.0*P;
            if ( IK-IH > 0 )
            {
                if ( IH > 0 )
                {
                    P = 2.0*P;
                }
            }
        }
        break;
    case 11:
        if ( IL > 0 )
        {
            P = 2.0*P;
        }
        if ( IH <= 0 )
        {
            if ( IK > 0 )
            {
                P = 3.0*P;
            }
        }
        else
        {
            P = 3.0*P;
            if ( IK > 0 )
            {
                P = 2.0*P;
            }
        }
        break;
    case 12:
        if ( IK > 0 )
        {
            P = 6.0*P;
            if ( IH != 0 )
            {
                if ( IK-IH != 0 )
                {
                    P = 2.0*P;
                }
            }
        }
        if ( IL != 0 )
        {
            P = 2.0*P;
        }
        break;
    case 13:
        if ( IH-IK == 0 )
        {
            if ( IK-IL == 0 )
            {
                P = 4.0*P;
            }
            else
            {
                P = 3.0*P;
                if ( IK > 0 )
                {
                    P = 2.0*P;
                }
                if ( IH > 0 )
                {
                    P = 2.0*P;
                }
            }
        }
        else
        {
            P = 3.0*P;
            if ( IK > 0 )
            {
                P = 2.0*P;
            }
            if ( IH > 0 )
            {
                P = 2.0*P;
            }
        }
        break;
    case 14:
        if ( IH == 0 )
        {
            P = 3.0*P;
            if ( IK != 0 )
            {
                P = 2.0*P;
                if ( IK-IL != 0 )
                {
                    P = 2.0*P;
                }
            }
        }
        else
        {
            P = 4.0*P;
            if ( IH-IL != 0 )
            {
                P = 3.0*P;
                if ( IH-IK != 0 )
                {
                    if ( IK-IL != 0 )
                    {
                        P = 2.0*P;
                    }
                }
            }
        }
        break;
    default:
        GOTOER();
        break;
    }


    //-----THIS LINE ADDED TO CORRECT THE PROBLEMS WITH THE -3, -3M1, 6/M,
    //-----AND 6/MMM SPACE GROUPS THAT WERE OFF BY A FACTOR OF TWO
    if (phases[IPHASE].SYMB.NSPGRP == 8  ||  phases[IPHASE].SYMB.NSPGRP == 9  ||  phases[IPHASE].SYMB.NSPGRP == 11  ||  phases[IPHASE].SYMB.NSPGRP == 12) P=0.5*P;
    return P;
}

void DBWS::REWRIT(int ISCALE, int IDIF)
{
    int I, J, K, N, IOF, ISOF, NPLOF, IIPHAS;

    file5.close();
    file5b.open(file5name.data(),ios::trunc);			//	REWIND 5;
    file5b << fixed;


    // line 1
    file5b << title << "          OpenDBWS" <<endl;
    JOBTYP=JOBTYP-1;
    //NPROF=NPROF-1;
    NPLOF = NPROF;
    if (NSIZESTRAIN == 9) NPLOF = 9;
    INSTRM = INSTRM-1;

    // line 2.1
    // line 2.1 changed due to size-strain calculation (NsizeStrain)
    file5b << setw(4) << JOBTYP
        << setw(4) << NPLOF
        << setw(4) << NPHASE
        << setw(4) << IBCKCODE
        << setw(4) << diffractogram.NEXCRG
        << setw(4) << diffractogram.NSCAT
        << setw(4) << INSTRM
        << setw(4) << IPREF
        << setw(4) << IASYM
        << setw(4) << IABSR
        << setw(4) << IDATA
        << setw(4) << ISPHASE
        << "         LINE 2.1" << endl;

    // line 2.2
    if(IBCKCODE == -1)
    {
        file5b << setw(4) << IAS
            << setw(4) << FONDO
            << "                                                 LINE 2.2" << endl;
    }

    // line 3
    file5b
        << setw(1) << IOT
        << setw(1) << IPL
        << setw(1) << IPC
        << setw(1) << MAT
        << setw(1) << NXT
        << " "
        << setw(1) << LST1
        << setw(1) << LST2
        << setw(1) << LST3
        << setw(1) << IPL1
        << setw(1) << IPL2
        << " "
        << setw(1) << IPLST
        << setw(1) << IPLOSS
        << setw(1) << IPLCAL
        << setw(1) << IPLPOL
        << setw(1) << IPLCOM
        << " "
        << setw(1) << IPLDIS
        << setw(1) << IPLAM
        << setw(1) << IPBIG
        << "                                    LINE 3" << endl;

    // line 4
    file5b
        << setw(8) << setprecision(5) << diffractogram.LAMDA[1]
        << setw(8) << setprecision(5) << diffractogram.LAMDA[2]
        << setw(8) << setprecision(5) << diffractogram.RATIO[2]
        << setw(8) << setprecision(4) << BKPOS_
        << setw(8) << setprecision(4) << WDT_
        << setw(8) << setprecision(4) << CTHM_
        << setw(8) << setprecision(4) << TMV_
        << setw(8) << setprecision(4) << RLIM_
        << setw(8) << setprecision(4) << SW
        << endl;

    // line 5
    file5b
        << setw(4) << MCYCLE
        << setw(4) << setprecision(2) << EPS
        << setw(4) << setprecision(2) << RELAX[1]
        << setw(4) << setprecision(2) << RELAX[2]
        << setw(4) << setprecision(2) << RELAX[3]
        << setw(4) << setprecision(2) << RELAX[4]
        << "                                 CYCLS EPS RELAX P_CALC" << endl;
    if(diffractogram.NBCKGD < 2)goto L120;

    // line 6(*)
    for (I=1; I <= diffractogram.NBCKGD; ++I) file5b << setw(8) << setprecision(2) << diffractogram.POS[I] << diffractogram.BCK[I] << endl;
L120:
    if(diffractogram.NEXCRG <= 0)goto L122;

    // line 7(*)
    for (I=1; I <= diffractogram.NEXCRG; ++I)
        file5b
            << setw(8) << setprecision(2) << diffractogram.ALOW[I]
            << setw(8) << setprecision(2) << diffractogram.AHIGH[I]
            << "                                         EXCLUDED REGION" << endl;
L122:
    if (diffractogram.NSCAT <= 0) goto L124;
    for (I=1; I <= diffractogram.NSCAT; ++I)
    {
        if (JOBTYP  ==  1 || JOBTYP == 3) goto L1228;
        //C line 8.1 XRD (*)
        file5b
            << setw(4) << NAM_[I]
            << setw(8) << setprecision(4) << coeff.DFP[I]
            << setw(8) << setprecision(4) << coeff.DFPP[I]
            << setw(8) << setprecision(4) << coeff.XMAS[I]
            << "                             SCATTERING SET" << setw(2) << diffractogram.NSCAT << endl;
        goto L126;
        // line 8.1 ND(*)
L1228:
        file5b
            << setw(4) << NAM_[I]
            << setw(8) << setprecision(4) << coeff.DFP[I]
            << setw(8) << setprecision(4) << coeff.XMAS[I]
            << "                                     SCATTERING SET " << setw(2) << diffractogram.NSCAT << endl;
        // line 8.2 XRD(*)
L126:
        if (JOBTYP == 0 || JOBTYP == 2)
        {
            for (J=1; J <= 9; ++J) file5b << setw(8) << setprecision(5) << coeff.AC[J][I];
            file5b << endl;
        }
    }
L124:;

    // line 9
    file5b << setw(8) << MAXS << "                                                 PARAMS REFINED" << endl;
    N=0;
    for (IIPHAS=1; IIPHAS <= NPHASE; ++IIPHAS) N=N+phases[IIPHAS].AtomCount;
    for (I=1; I <= N; ++I)
    {
        for (J=1; J <= 11; ++J) XL_[I][J].codeword=sign(XL_[I][J].codeword)*(static_cast<double>(10*XL_[I][J].L)+abs(XL_[I][J].codeword));
    }
    for (I=1; I <= NPHASE; ++I)
    {
        for (J=1; J <= 6; ++J)
        {
            phases[I].PAR[J+5-1]=SAVE[I][J];
        }
        for (J=1; J <= 27; ++J)
        {
            phases[I].PAR[J-1].codeword=sign(phases[I].PAR[J-1].codeword)*(static_cast<double>(10*phases[I].PAR[J-1].L)+abs(phases[I].PAR[J-1].codeword));
        }
    }
    for (J=1; J <= 20; ++J)
    {
        GLB_[J-1].codeword=sign(GLB_[J-1].codeword)*(static_cast<double>(10*GLB_[J-1].L)+abs(GLB_[J-1].codeword));
    }

    // line 10.1
    file5b
        << setw(8) << setprecision(4) << GLB_[1-1]
        << setw(8) << setprecision(4) << GLB_[10-1]
        << setw(8) << setprecision(4) << GLB_[11-1]
        << setw(8) << setprecision(4) << GLB_[8-1]
        << setw(8) << setprecision(4) << GLB_[9-1]
        << setw(8) << setprecision(4) << GLB_[12-1]
        << setw(8) << setprecision(4) << GLB_[13-1] << " ZER DISP TRANS p q r t" << endl;

    // line 10.11
    file5b
        << setw(8) << setprecision(4) << GLB_[1-1]
        << setw(8) << setprecision(4) << GLB_[10-1].codeword
        << setw(8) << setprecision(4) << GLB_[11-1].codeword
        << setw(8) << setprecision(4) << GLB_[8-1].codeword
        << setw(8) << setprecision(4) << GLB_[9-1].codeword
        << setw(8) << setprecision(4) << GLB_[12-1].codeword
        << setw(8) << setprecision(4) << GLB_[13-1].codeword << " CODEWORDS" << endl;
    if (IBGD == 1) goto L4600;
    // line 10.2
    file5b
        << setw(8) << setprecision(4) << GLB_[20-1]
        << setw(8) << setprecision(4) << GLB_[18-1]
        << setw(8) << setprecision(4) << GLB_[19-1]
        << "                                 AM MON1 MON2" << endl;

    // line 10.21
    file5b
        << setw(8) << setprecision(4) << GLB_[20-1].codeword
        << setw(8) << setprecision(4) << GLB_[18-1].codeword
        << setw(8) << setprecision(4) << GLB_[19-1].codeword
        << "                                 CODEWORS" << endl;
L4600:
    if(diffractogram.NBCKGD == 0)
    {
        file5b
            << setw(9) << setprecision(2) << GLB_[2-1]
            << setw(9) << setprecision(2) << GLB_[3-1]
            << setw(9) << setprecision(2) << GLB_[4-1]
            << setw(9) << setprecision(2) << GLB_[5-1]
            << setw(9) << setprecision(2) << GLB_[6-1]
            << setw(9) << setprecision(2) << GLB_[7-1]
            << "   BACKGROUND" << endl;
        file5b
            << setw(9) << setprecision(2) << GLB_[2-1].codeword
            << setw(9) << setprecision(2) << GLB_[3-1].codeword
            << setw(9) << setprecision(2) << GLB_[4-1].codeword
            << setw(9) << setprecision(2) << GLB_[5-1].codeword
            << setw(9) << setprecision(2) << GLB_[6-1].codeword
            << setw(9) << setprecision(2) << GLB_[7-1].codeword
            << "   CODEWORDS" << endl;
    }
    //L477:
    for (K=1; K <= NPHASE; ++K)
    {
        IOF=0;
        if (K > 1)
        {
            for (IIPHAS=2; IIPHAS <= K; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
        }
        // line 11.1
        file5b << setw(50) << phases[K].name << "       PHASE NUMBER " << setw(2) << K << endl;

        // line 11.2
        file5b
            << setw(4) << phases[K].AtomCount
            << setw(4) << phases[K].NMOL
            << setw(7) << setprecision(4) << phases[K].SAQF
            << " "
            << setw(4) << setprecision(1) << phases[K].PREF[1]
            << setw(4) << setprecision(1) << phases[K].PREF[2]
            << setw(4) << setprecision(1) << phases[K].PREF[3]
            << setw(7) << setprecision(2) << phases[K].WTIS
            << "                      #ATMS #FU AFQPA PREFDIR ISWT" << endl;
        N=phases[K].AtomCount;

        // line 11.3
        file5b << setw(20) << phases[K].SYMB << "                                     SPACE GROUP" << endl;
        // Changing N to 'so' !cp oct 96
        for (ISOF=1; ISOF <= N; ++ISOF) XL_[ISOF+IOF][5]=XL_[ISOF+IOF][5]*phases[K].XMLTP/MURT[ISOF+IOF];
        // !cp oct 96 #6 murt parametrs included below. FORMAT modified...
        for (I=1; I <= N; ++I)
        {
            file5b
                << setw(4) << ATEXT[I+IOF] << " "
                << setw(4) << MURT[I+IOF] << " "
                << setw(4) << NTYP[I+IOF] << "  "
                << setw(8) << setprecision(5) << XL_[I+IOF][1]
                << setw(8) << setprecision(5) << XL_[I+IOF][2]
                << setw(8) << setprecision(5) << XL_[I+IOF][3]
                << setw(8) << setprecision(5) << XL_[I+IOF][4]
                << setw(8) << setprecision(5) << XL_[I+IOF][5]
                << "  LBL M NTYP x y z B So" << endl
                << "                "
                << setw(8) << setprecision(5) << XL_[I+IOF][1].codeword
                << setw(8) << setprecision(5) << XL_[I+IOF][2].codeword
                << setw(8) << setprecision(5) << XL_[I+IOF][3].codeword
                << setw(8) << setprecision(5) << XL_[I+IOF][4].codeword
                << setw(8) << setprecision(5) << XL_[I+IOF][5].codeword
                << "  CODEWORDS" << endl
                << setw(8) << setprecision(5) << XL_[I+IOF][6]
                << setw(8) << setprecision(5) << XL_[I+IOF][7]
                << setw(8) << setprecision(5) << XL_[I+IOF][8]
                << setw(8) << setprecision(5) << XL_[I+IOF][9]
                << setw(8) << setprecision(5) << XL_[I+IOF][10]
                << setw(8) << setprecision(5) << XL_[I+IOF][11]
                << "          BETAS" << endl
                << setw(8) << setprecision(2) << XL_[I+IOF][6].codeword
                << setw(8) << setprecision(2) << XL_[I+IOF][7].codeword
                << setw(8) << setprecision(2) << XL_[I+IOF][8].codeword
                << setw(8) << setprecision(2) << XL_[I+IOF][9].codeword
                << setw(8) << setprecision(2) << XL_[I+IOF][10].codeword
                << setw(8) << setprecision(2) << XL_[I+IOF][11].codeword
            << "          CODEWORDS" << endl;
        }
        file5b
            << scientific << setw(8) << setprecision(2) << phases[K].PAR[0]
            << fixed
            << setw(8) << setprecision(4) << phases[K].PAR[1]
            << "                                         SCALE Bo(OVERALL)" << endl
            << setw(8) << setprecision(2) << phases[K].PAR[0].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[1].codeword << endl
            << setw(8) << setprecision(5) << phases[K].PAR[2]
            << setw(8) << setprecision(5) << phases[K].PAR[3]
            << setw(8) << setprecision(5) << phases[K].PAR[4]
            << setw(8) << setprecision(5) << phases[K].PAR[20]
            << setw(8) << setprecision(5) << phases[K].PAR[19]
            << setw(8) << setprecision(5) << phases[K].PAR[14]
            << setw(8) << setprecision(5) << phases[K].PAR[15]
            << " U V W CT Z X Y" << endl

            << setw(8) << setprecision(2) << phases[K].PAR[2].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[3].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[4].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[20].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[19].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[14].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[15].codeword << endl

            << setw(8) << setprecision(4) << phases[K].PAR[5]
            << setw(8) << setprecision(4) << phases[K].PAR[6]
            << setw(8) << setprecision(4) << phases[K].PAR[7]
            << setw(8) << setprecision(4) << phases[K].PAR[8]
            << setw(8) << setprecision(4) << phases[K].PAR[9]
            << setw(8) << setprecision(4) << phases[K].PAR[10]
            << "         CELL PARAMETERS" << endl

            << setw(8) << setprecision(2) << phases[K].PAR[5].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[6].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[7].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[8].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[9].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[10].codeword << endl

            << setw(8) << setprecision(5) << phases[K].PAR[11]
            << setw(8) << setprecision(5) << phases[K].PAR[12]
            << setw(8) << setprecision(5) << phases[K].PAR[13]
            << "                                 PREF1 PREF2 R/RCF_ASYM" << endl

            << setw(8) << setprecision(2) << phases[K].PAR[11].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[12].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[13].codeword << endl

            << setw(8) << setprecision(4) << phases[K].PAR[16]
            << setw(8) << setprecision(4) << phases[K].PAR[17]
            << setw(8) << setprecision(4) << phases[K].PAR[18]
            << "                                 NA NB NC (MIX_PARAMS)" << endl

            << setw(8) << setprecision(2) << phases[K].PAR[16].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[17].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[18].codeword << endl

            << setw(8) << setprecision(4) << phases[K].PAR[23]
            << setw(8) << setprecision(4) << phases[K].PAR[24]
            << setw(8) << setprecision(4) << phases[K].PAR[25]
            << "                                 NA NB NC (HIGH SIDE)" << endl

            << setw(8) << setprecision(2) << phases[K].PAR[23].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[24].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[25].codeword << endl

            << setw(8) << setprecision(4) << phases[K].PAR[26] << "" << endl

            << setw(8) << setprecision(2) << phases[K].PAR[26].codeword << endl;
    }
    if (IPL != 0) file5b << setw(8) << ISCALE << IDIF << "                                         LINE PRINTER INFO" << endl;
    if (file5b.is_open()) file5b.close();
    //L151:
    JOBTYP=JOBTYP+1;
}

// subroutine to compute size&strain (NsizeStrain)
void DBWS::size(int K)
{
    double A7, B7, C7, D7,HD, UA, XA, YA, ZA, HS, UI, XI, YI, ZI, HDG, HDL, HSG, HSL,RADO;
    int ISIZETEST,ISTRAINTEST;
    ifstream file17;
    string s;

    RADO=M_PI/180.0;
    A7 = 2.69269;
    B7 = 2.42843;
    C7 = 4.47163;
    D7 = 0.07842;
    file17.open("instr.dat");
    getline(file17,s);				// na primeira linha de padrao.dat sera apenas para comentario
    // a segunda linha ira ler U,Z,X,Y da FWHM do padrao
    file17 >> UI >> ZI >> XI >> YI;
    file17.close();
    UA = phases[K].PAR[2];
    ZA = phases[K].PAR[19];
    XA = phases[K].PAR[14];
    YA = phases[K].PAR[15];
    if (ZI  !=  0.0)
    {
        file6 << "SIZE NOT COMPUTED * bad standard? - check Z for standard." << endl;
        cout << "SIZE NOT COMPUTED * bad standard? - check Z for standard." << endl;
        goto L20;
    }
    if(ZA  <  0.0  ||  YA  <  YI)
    {
        file6 << "Phase " << setw(2) << K << ": Size-strain not Computed * check Z and',' Y for sample and standard" << endl;
        cout << "Phase " << setw(2) << K << ": Size-strain not Computed * check Z and',' Y for sample and standard" << endl;
        goto L20;
    }
    if(YA  ==  YI  &&  ZA  == 0.0)
    {
        file6 << "Phase " << setw(2) << K << ": Infinite Size * check Z and Y values',' for sample and standard" << endl;
        cout << "Phase " << setw(2) << K << ": Infinite Size * check Z and Y values',' for sample and standard" << endl;
        ISIZETEST = 1;
    }
    else
    {
        ISIZETEST = 0;
    }
    //
    // COMPUTE SIZE
    //
    //        COMPUTE SIZE BASED ON Gauss-COS Z parameter [sizeG(k)]
    //
    if(ZA > 0.0)
    {
        HSG = sqrt(ZA)*RADO;
        SIZEG[K] =  diffractogram.LAMDA[1] / HSG;
    }
    else
    {
        HSG = 0.0;
        SIZEG[K] = 99999.0;
    }

    //
    //        COMPUTE SIZE BASED ON Loren-COS Y parameter [sizeL(k)]
    //  PAR[K][16] = Y (Lorentz)
    //
    if(YA  >  YI)
    {
        HSL = (YA - YI) * RADO;
        SIZEL[K] =  diffractogram.LAMDA[1] / HSL;
    }
    else
    {
        HSL = 0.0;
        SIZEL[K] = 99999.0;
    }
    //
    //        COMPUTE WEIGHTED SIZE BASED on Gauss & Lorentz [siz(k)]
    //
    if(ISIZETEST == 0)
    {
        HS  = pow((pow(HSG,5) + A7 * pow(HSG,4) * HSL + B7 * pow(HSG,3) * pow(HSL,2) + C7 * pow(HSG,2) * pow(HSL,3) + D7 * HSG * pow(HSL,4) + pow(HSL,5)),0.2);
        SIZ[K] = diffractogram.LAMDA[1] / HS;
    }
    else
    {
        SIZ[K] = 0.0;
    }
    //
    //       STRAIN COMPUTATIONS
    //L10:
    if(UA  <  UI  ||  XA  <  XI)
    {
        file6 << "Phase " << setw(2) << K << ": Size-Strain not computed * check U and X for sample and standard" << endl;
        cout << "Phase " << setw(2) << K << ": Size-Strain not computed * check U and X for sample and standard" << endl;
        goto L20;
    }
    if(XA  ==  XI  &&  UA  ==  UI)
    {
        file6 << "Phase " << setw(2) << K << ": Strain not Present * check U and X values for sample and standard" << endl;
        cout << "Phase " << setw(2) << K << ": Strain not Present * check U and X values for sample and standard" << endl;
        ISTRAINTEST = 1;
    }
    else
    {
        ISTRAINTEST = 0;
    }

    //
    // COMPUTE STRAIN
    //
    // compute strain based on Gauss-U * [strainG(k)]
    //
    if(UA  >  UI)
    {
        HDG = sqrt( (UA - UI) )*RADO;
        STRAING[K] =  HDG / 2;
    }
    else
    {
        HDG = 0.0;
        STRAING[K] = 0.0;
    }
    //
    // compute strain based on Lorentz-X * [strainL(k)]
    //
    if(XA  >  XI)
    {
        HDL = (XA - XI) * RADO;
        STRAINL[K] =  HDL / 2;
    }
    else
    {
        HDL = 0.0;
        STRAINL[K] = 0.0;
    }
    // compute weighted strain: based on Gauss-U and Lorentz-X * [strain(k)]
    //
    if (ISTRAINTEST == 0)
    {
        HD  = pow(( pow(HDG,5) + A7 * pow(HDG,4) * HDL + B7 * pow(HDG,3) * pow(HDL,2) + C7 * pow(HDG,2) * pow(HDL,3) + D7 * HDG * pow(HDL,4) + pow(HDL,5)),0.2);
        STRAIN[K] = HD / 2;
    }
    else
    {
        STRAIN[K] = 0.0;
    }
    file6 << "Gauss(phase " << setw(2) << K << "): Size =" //<< scientific
        << setw(9) << setprecision(4) << SIZEG[K]
    << " Strain ="
        << setw(9) << setprecision(3) << STRAING[K] << endl
        << "Lorentz(phase " << setw(2) << K << "): Size ="
        << setw(9) << setprecision(4) << SIZEL[K]
    << " Strain ="
        << setw(9) << setprecision(3) << STRAINL[K] << endl
        << "Weighted(phase " << setw(K) << "): Size ="
        << setw(9) << setprecision(4) << SIZ[K]
    << " Strain ="
        << setw(9) << setprecision(3) << STRAIN[K] << endl;
L20:
    return;
}



void DBWS::EXPUT(void)
{
    const char ISPACE		= ' ';
    const char IDOT			= '.';
    const char ISTAR		= '*';
    const char IPLUS		= '+';
    const char IMINUS		= '0';
    const char IBAR			= '-';

    int I,J,K,L,N,J20,JJ,IX,IS,NP,IY,IRH,IRK,IRL,ICZ,IRC,MIN,MAX,IXX,
        IDIF,J20IK,IEXC,NCOL,IXXX,NCOL1,NPTS2,ISCALE,
        IIPHAS,NPAGES,IFINIS,NLINES,ISTART,INUMB;
    double Z,BB,TH,TH2,HW,YX,RAD,TLG,TLL,TIC,RMIN,RMAX,XXXX,
        T2OBS,AFCAL,BFCAL,OMEGA,TFCAL,AFOBS,BFOBS,TDOBS,RFDNR,SHIFT,TIOBS,
        TFOBS,TLABEL,RFNMR;
    string IOUT,TITLE1,TITLE2,TITLE3,TITLE4,TITLE6,TITLE7,TITLE8,OUTFILE;
    bool VERT;
    int KI[10+1];
    double DEL[10+1];
    int LABEL[12+1];
    double DUMMY[2*MSZ+1 +1];
    string s;
    ofstream file69,file690,file9,file10,file31,file32,file33,file34,file37,file36,file38;



    ISCALE = 0;
    IDIF = 0;
    TITLE1 = "OSS.COR\"";
    TITLE2 = "TOTALE \"";
    TITLE3 = "COM.TOT\"";
    TITLE4 = "DIS.TOT\"";
    TITLE6 = "AM.TOT \"";
    TITLE7 = "POL.TOT\"";
    TITLE8 = "ALL.TOT\"";
    RAD = 45.0/atan(1.0);
    if (IPL != 0)
    {
        getline(file5,s);
        stringstream(s.substr(0,8)) >> ISCALE;
        stringstream(s.substr(8,8)) >> IDIF;
    }
    if (NXT != 0 && JOBTYP < 3) REWRIT(ISCALE,IDIF);
    if (IPL2 != 0)
    {
        file69.open("bragg.dat");
        file690.open("xy-int.dat");
        file690 << "  2theta_i      y_o       y_c        yo-yc  w(yo-yc)^2" << endl;
    }
    if (JOBTYP <= 2)
    {
        file6 << endl
            << "AVERAGE INTENSITY DIFFERENCE FOR PATTERN" << endl
            << "GIVEN FOR BLOCKS OF 20 OBSERVATIONS." << endl;
        for (I=1; I <= diffractogram.NPTS; I = I + 200)
        {
            for (J=1; J <= 10; ++J)
            {
                J20=20*(J-1);
                DEL[J]=0;
                for (K=1; K <= 20; ++K)
                {
                    J20IK=J20+I+K-1;
                    DEL[J]=DEL[J]+diffractogram.Y[J20IK]-diffractogram.BK[J20IK]-diffractogram.YC[J20IK];
                }
                KI[J]=J+I/20;
                DEL[J]=DEL[J]/20.0;
                if (KI[J]*20 >= diffractogram.NPTS) goto L15;
            }
            J=10;
    L15:
            file6 << "     ";
            for (JJ=1; JJ <= J; ++JJ)
            {
                file6 << setw(4) << KI[JJ] << "  " << setw(5) << setprecision(1) << DEL[JJ];
            }
            file6 << endl;
        }
    }


    for (I=1; I <= diffractogram.NPTS; ++I)
    {
        diffractogram.YC[I]=diffractogram.YC[I]+diffractogram.BK[I];
        if (JOBTYP > 2) diffractogram.Y[I]=diffractogram.YC[I];
    }
    if (IPC != 0)
    {
        file9.open("plotinfo");
        file10.open("plotinfo.bin");
        for (I = 1; I <= diffractogram.NPTS; ++I) file10 << diffractogram.BK[I];
        file9 << title << endl;
        file9 << "NO. OF PHASESZ  " << NPHASE << endl
              << "NO. OF REFLECTIONS IN EACH PHASEQ  ";
        for (IIPHAS=1; IIPHAS <= NPHASE; ++IIPHAS) file9 << setw(4) << phases[IIPHAS].ICR;
        file9 << endl;
        file9 << "BRAGG POSITIONSZ" << endl;
        for (K=1; K <= NPHASE; ++K)
        {
            for (I=1; I <= diffractogram.NPTS; ++I) diffractogram.KR[I] = 0;
            if ((K % 2) == 1) ILOC_ = ILOC_ + 1;
            file6 << "PHASE NO. = " << K << "     PHASE NAME " << phases[K].name << endl;
            file9 << phases[K].name;
            T2OBS=0.0;
            TDOBS=0.0;
            RFNMR=0.0;
            RFDNR=0.0;
            ICZ=0;
            for (IIPHAS=1; IIPHAS <= NPHASE; ++IIPHAS) ICZ = ICZ + phases[IIPHAS].ICR;
            IXX=0;
            IXXX=0;
            for (IX=1; IX <= ICZ; ++IX)
            {
                if( K == refs[IX].IREFS/(256*256*256*8) )
                {
                    SHIFT = GLB_[10-1] * cos(refs[IX].REFS[2]/2.0/57.2958) + GLB_[11-1] * sin(refs[IX].REFS[2]/57.2958);
                    IXX=IXX+1;
                    for (J=1; J <= diffractogram.NEXCRG; ++J)
                    {
                        //-----CHECK FOR SPLIT PEARSON VII PROFILE
                        if (NPROF == _SplitPearsonVII)
                        {
                            if ((refs[IX].REFS[2]+SHIFT+GLB_[1-1]) >= (diffractogram.ALOW[J]-WDT_*refs[IX].FWHM[1]) && (refs[IX].REFS[2]+SHIFT+GLB_[1-1]) <= (diffractogram.AHIGH[J]+WDT_*refs[IX].FWHM[2]))  goto L481;
                            //-----FOR ALL OTHER PROFILES
                        }
                        else
                        {
                            if ((refs[IX].REFS[2]+SHIFT+GLB_[1-1]) >= (diffractogram.ALOW[J]-WDT_*refs[IX].REFS[1]) && (refs[IX].REFS[2]+SHIFT+GLB_[1-1]) <= (diffractogram.AHIGH[J]+WDT_*refs[IX].REFS[1])) goto L481;
                        }
                    }
                    IXXX=IXXX+1;                           //test !cp 29 jun 98
                    IRL=(refs[IX].IREFS % 256)-128;
                    IRK=((refs[IX].IREFS/256) % 256)-128;
                    IRH=((refs[IX].IREFS/(256*256)) % 256)-128;
                    IRC=(refs[IX].IREFS/(256*256*256)) % 8;
                    //-----CHECK FOR THE SPLIT PEARSON VII PROFILE
                    //-----if SO CHANGE THE PROFILE LIMITS
                    if (NPROF == _SplitPearsonVII)
                    {
                        RMIN=refs[IX].REFS[2]+GLB_[1-1]+SHIFT-WDT_*refs[IX].FWHM[1];
                        RMAX=refs[IX].REFS[2]+GLB_[1-1]+SHIFT+WDT_*refs[IX].FWHM[2];
                    }
                    else
                    {
                        //-----FOR ALL OTHER PROFILES
                        RMIN=refs[IX].REFS[2]+GLB_[1-1]+SHIFT-WDT_*refs[IX].REFS[1];
                        RMAX=refs[IX].REFS[2]+GLB_[1-1]+SHIFT+WDT_*refs[IX].REFS[1];
                    }
                    MIN=static_cast<int>((RMIN-diffractogram.THMIN)/diffractogram.STEP+1.5);
                    MAX=static_cast<int>((RMAX-diffractogram.THMIN)/diffractogram.STEP+1.5);
                    MIN=max(MIN,1);
                    MIN=min(MIN,diffractogram.NPTS);
                    MAX=min(MAX,diffractogram.NPTS);
                    //-----PATCH TO CALCULATE R-BRAGG
                    TL_=refs[IX].REFS[1];
                    VERT=refs[IX].REFS[2] <= RLIM_;
                    if (NPROF == _pseudoVoigt)
                    {
                        GAM1_=phases[K].PAR[16] + phases[K].PAR[17] * refs[IX].REFS[2];
                    }
                    else if (NPROF == _PearsonVII)
                    {
                        GAM1_=phases[K].PAR[16]+(phases[K].PAR[17]+phases[K].PAR[18]/refs[IX].REFS[2])/refs[IX].REFS[2];
                        PRSVII(GAM1_);
                    }
                    else if (NPROF == _TCHZ)
                    {
                        GAM1_=refs[IX].GAM;
                        TLG=refs[IX].HALFG;
                        TLL=refs[IX].HALFL;
                    }
                    else if (NPROF == _SplitPearsonVII)
                    {
                        RL_=phases[K].PAR[16]+(phases[K].PAR[17]+phases[K].PAR[18]/refs[IX].REFS[2])/refs[IX].REFS[2];
                        RH_=phases[K].PAR[23]+(phases[K].PAR[24]+phases[K].PAR[25]/refs[IX].REFS[2])/refs[IX].REFS[2];
                        mspvii(phases[K].PAR[26],TL_);
                    }
                    BB=TL_*TL_;
                    TIOBS=0.0;
                    TIC=0.0;
                    refs[IX].FMGNTD = refs[IX].FMGNTD*refs[IX].REFS[3]*phases[K].PAR[0];
                    for (IS=MIN; IS <= MAX; ++IS)
                    {
                        TH=diffractogram.THMIN+static_cast<double>(IS-1)*diffractogram.STEP;
                        if(diffractogram.NEXCRG > 0)
                        {
                            for (IEXC=1; IEXC <= diffractogram.NEXCRG; ++IEXC)
                            {
                                if (TH >= diffractogram.ALOW[IEXC] && TH <= diffractogram.AHIGH[IEXC]) goto L410;
                            }
                        }
                        TH2 = 0.5*TH/RAD;
                        DELTA_=TH-refs[IX].REFS[2]-GLB_[1-1]-SHIFT;
                        DELT_=DELTA_*DELTA_;
                        //     NEXT LINE IS NECESSEARY FOR 2 PHASES WITH VERY DIFFERENT FWHM.
                        if (DELT_/BB > WDT_*WDT_) goto L410;
                        if( !VERT )
                        {
                            Z=1.0;
                        }
                        else
                        {
                            if (IASYM == 0)
                            {
                                YX=DELT_ * sign(DELTA_);
                                Z=1.0-phases[K].PAR[13]*YX/tan(TH2);
                            }
                            else
                            {
                                YX=sign(DELTA_)*DELTA_/(2*TL_);
                                if (TH2 > (45.0/RAD)) TH2 = TH2-(90.0/RAD);
                                Z=(phases[K].PAR[13]/tan(TH2)) * (2.0*(DELTA_/(2*TL_))*exp(-YX));
                                Z=1+Z;
                            }
                            if(Z <= 0.0)Z=0.0001;
                        }
                        OMEGA= Z*PROFIL(NPROF,DELT_/BB);
                        if (NPROF == _SplitPearsonVII)
                        {
                            diffractogram.KR[IS] = diffractogram.KR[IS] + static_cast<int>(OMEGA*refs[IX].FMGNTD);
                        }
                        else
                        {
                            diffractogram.KR[IS] = diffractogram.KR[IS] + static_cast<int>(OMEGA*refs[IX].FMGNTD/TL_);
                        }
                        TIC   = TIC + OMEGA;
                        TIOBS=TIOBS+OMEGA*(diffractogram.Y[IS]-diffractogram.BK[IS])/ max((diffractogram.YC[IS]-diffractogram.BK[IS]),1.0);
        L410:;
                    }
                    TIC   = TIC   * refs[IX].FMGNTD/TL_;
                    TIOBS = TIOBS * refs[IX].FMGNTD/TL_;
                    // ! FROM ITALIAN CODE  !cp ap 97
                    //        ***************************************************************
                    //        NEXT LINE IS FOR NOT EVALUATING R-BRAGG FOR REFLECTIONS
                    //        WHICH ARE OUTSIDE THE MEASUREMENTS RANGE BUT WHOSE TAILS ARE
                    //        IN THE PATTERN
                    //        ***************************************************************
                    if((refs[IX].REFS[2]+SHIFT+GLB_[1-1]) <= diffractogram.THMAX)
                    {
                        T2OBS = T2OBS+TIOBS;
                        TDOBS = TDOBS + abs(TIOBS-TIC);
                        if (IPC == 2 || IPC == 3)
                        {
                            TFCAL= sqrt(abs(TIC/refs[IX].REFS[3]/phases[K].PAR[0]/refs[IX].TAVIX/refs[IX].SRIX));
                            TFOBS= sqrt(abs(TIOBS/refs[IX].REFS[3]/phases[K].PAR[0]/refs[IX].TAVIX/refs[IX].SRIX));
                            AFCAL=TFCAL*cos(refs[IX].APHASE);
                            BFCAL=TFCAL*sin(refs[IX].APHASE);
                            AFOBS=TFOBS*cos(refs[IX].APHASE);
                            BFOBS=TFOBS*sin(refs[IX].APHASE);
                            //        Line below also from italian code   !cp ap 97
                            if(refs[IX].REFS[2] <= diffractogram.THMAX)
                            {
                                RFNMR = RFNMR + abs(TFOBS-TFCAL);
                                RFDNR = RFDNR + abs(TFOBS);
                            }
                        }
                    }
                    if (NPROF == _TCHZ)
                    {
                        if((IXXX-1 % 60) == 0)
                        {
                            if (IPC == 1)
                            {
                                file6 << "NO. CODE     H   K   L     HW     POSN     ICALC    IOBS       HG      HL      ETA   PHASE_C" << endl;
                                if(IPL2 != 0) file69 << "  H   K   L     HW      POSN     ICALC    IOBS    HGHL      ETA   PHASE_C" << endl;
                            }
                            if (IPC == 2)
                            {
                                file6 << "NO. CODE     H   K   L     HW     POSN    ICALC      IOBS      HG      HL      ETA     FCALC     FOBS   PHASE_C" << endl;
                                if(IPL2 != 0) file69 << "  H   K   L     HW       POSN    ICALC    IOBS    HG     HL      ETA   FCALC     FOBS   PHASE_C" << endl;
                            }
                            if (IPC == 3)
                            {
                                file6 << "NO. CODE     H   K   L     HW     POSN    ICALC      IOBS      HG      HL      ETA    A_CALC    B_CALC     A_OBS     B_OBS" << endl;
                                if(IPL2 != 0) file69 << "  H   K   L     HW       POSN    ICALC    IOBS    HG     HL      ETA   A_CALC    B_CALC     A_OBS     B_OBS" << endl;
                            }
                        }
                        HW=refs[IX].REFS[1];
                        if (IPC == 1)
                        {
                            file6 << setw(4) << IXX
                                << setw(4) << IRC << "   "
                                << setw(4) << IRH
                                << setw(4) << IRK
                                << setw(4) << IRL
                                << setw(8) << setprecision(3) << HW
                                << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                << setw(10) << setprecision(0) << TIC
                                << setw(10) << setprecision(0) << TIOBS
                                << setw(8) << setprecision(3) << TLG
                                << setw(8) << setprecision(3) << TLL
                                << setw(8) << setprecision(3) << GAM1_
                                << setw(10) << setprecision(3) << RAD*refs[IX].APHASE << endl;
                            if(IPL2 != 0)
                            {
                                file69 << setw(3) << IRH << " "
                                    << setw(3) << IRK << " "
                                    << setw(3) << IRL << " "
                                    << setw(8) << setprecision(4) << HW << " "
                                    << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                    << setw(8) << setprecision(0) << TIC << " "
                                    << setw(8) << setprecision(4) << TIOBS << " "
                                    << setw(8) << setprecision(4) << TLG << " "
                                    << setw(8) << setprecision(4) << TLL << " "
                                    << setw(8) << setprecision(4) << GAM1_ << " "
                                    << setw(9) << setprecision(4) << RAD*refs[IX].APHASE
                                    << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                            }
                        }
                        else if(IPC == 2)
                        {
                            if (IRC == 1)
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(8) << setprecision(3) << TLG
                                    << setw(8) << setprecision(3) << TLL
                                    << setw(8) << setprecision(3) << GAM1_
                                    << setw(10) << setprecision(3) << TFCAL
                                    << setw(10) << setprecision(3) << TFOBS << " "
                                    << setw(8) << setprecision(2) << refs[IX].APHASE << endl;
                                if(IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << HW << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(0) << TIC << " "
                                        << setw(8) << setprecision(0) << TIOBS << " "
                                        << setw(8) << setprecision(4) << TLG << " "
                                        << setw(8) << setprecision(4) << TLL << " "
                                        << setw(8) << setprecision(4) << GAM1_ << " "
                                        << setw(8) << setprecision(4) << TFCAL << " "
                                        << setw(8) << setprecision(4) << TFOBS << " "
                                        << setw(9) << setprecision(4) << refs[IX].APHASE << " "
                                        << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                            else
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(8) << setprecision(3) << TLG
                                    << setw(8) << setprecision(3) << TLL
                                    << setw(8) << setprecision(3) << GAM1_ << endl;
                            }
                        }
                        else if(IPC == 3)
                        {
                            if (IRC == 1)
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(8) << setprecision(3) << TLG
                                    << setw(8) << setprecision(3) << TLL
                                    << setw(8) << setprecision(3) << GAM1_
                                    << setw(10) << setprecision(4) << AFCAL
                                    << setw(10) << setprecision(4) << BFCAL
                                    << setw(10) << setprecision(4) << AFOBS
                                    << setw(10) << setprecision(4) << BFOBS << endl;
                                if(IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << HW << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(0) << TIC << " "
                                        << setw(8) << setprecision(0) << TIOBS << " "
                                        << setw(8) << setprecision(4) << TLG << " "
                                        << setw(8) << setprecision(4) << TLL << " "
                                        << setw(8) << setprecision(4) << GAM1_ << " "
                                        << setw(8) << setprecision(4) << AFCAL << " "
                                        << setw(8) << setprecision(4) << BFCAL << " "
                                        << setw(8) << setprecision(4) << AFOBS << " "
                                        << setw(8) << setprecision(4) << BFOBS << " "
                                        << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                            else
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(8) << setprecision(3) << TLG
                                    << setw(8) << setprecision(3) << TLL
                                    << setw(8) << setprecision(3) << GAM1_ << endl;
                            }
                        }
                    }
                    else if (NPROF == _SplitPearsonVII)
                    {
                        if( (IXXX-1 % 60) == 0)
                        {
                            if(IPC == 1)
                            {
                                file6 << "NO.  CODE    H   K   L      HWL     HWH  POSN       ICALC     IOBS  PHASE" << endl;
                                if(IPL2 != 0) file69 << "  H   K   L      HWL       HWH      POSN     ICALC      IOBS    PHASE  Bragg_Pos" << endl;
                            }
                            if(IPC == 2)
                            {
                                file6 << "NO.  CODE    H   K   L      HWL     HWH  POSN      ICALC     IOBS        FCALC    FOBS       PHASE_C" << endl;
                                if(IPL2 != 0) file69 << "  H   K   L      HWL       HWH      POSN     ICALC      IOBS    FCALC     FOBS     PHASE_C  Bragg_Pos" << endl;
                            }
                            if(IPC == 3)
                            {
                                file6 << "NO.  CODE    H   K   L      HWL     HWH  POSN      ICALC     IOBS        A_CALC     B_CALC     A_OBS     B_OBS" << endl;
                                if(IPL2 != 0) file69 << "  H   K   L      HWL      HWH      POSN     ICALC      IOBS   A_CALC   B_CALC   A_OBS    B_OBS  Bragg_Pos" << endl;
                            }
                        }
                        if(IPC == 1)
                        {
                            file6 << setw(4) << IXX
                                << setw(4) << IRC << "   "
                                << setw(4) << IRH
                                << setw(4) << IRK
                                << setw(4) << IRL
                                << setw(8) << setprecision(3) << refs[IX].FWHM[1]
                                << setw(8) << setprecision(3) << refs[IX].FWHM[2]
                                << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                << setw(10) << setprecision(3) << TIC
                                << setw(10) << setprecision(3) << TIOBS << endl;
                            if(IRC == 1)
                            {
                                if (IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << refs[IX].FWHM[1] << " "
                                        << setw(8) << setprecision(4) << refs[IX].FWHM[2] << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(0) << TIC << " "
                                        << setw(8) << setprecision(0) << TIOBS << " "
                                        << setw(8) << setprecision(4) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                        }
                        else if(IPC == 2)
                        {
                            if (IRC == 1)
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[1]
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[2]
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(10) << setprecision(3) << TFCAL
                                    << setw(10) << setprecision(3) << TFOBS << " "
                                    << setw(8) << setprecision(2) << RAD*refs[IX].APHASE << endl;
                                if (IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << refs[IX].FWHM[1] << " "
                                        << setw(8) << setprecision(4) << refs[IX].FWHM[2] << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(0) << TIC << " "
                                        << setw(8) << setprecision(0) << TIOBS << " "
                                        << setw(8) << setprecision(4) << TFCAL << " "
                                        << setw(8) << setprecision(4) << TFOBS << " "
                                        << setw(9) << setprecision(4) << RAD*refs[IX].APHASE
                                        << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                            else
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[1]
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[2]
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC << endl;
                            }
                        }
                        else if(IPC == 3)
                        {
                            if (IRC == 1)
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[1]
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[2]
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(10) << setprecision(3) << AFCAL
                                    << setw(10) << setprecision(3) << BFCAL
                                    << setw(10) << setprecision(3) << AFOBS
                                    << setw(10) << setprecision(3) << BFOBS << endl;
                                if (IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << refs[IX].FWHM[1] << " "
                                        << setw(8) << setprecision(4) << refs[IX].FWHM[2] << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(0) << TIC << " "
                                        << setw(8) << setprecision(0) << TIOBS << " "
                                        << setw(8) << setprecision(2) << AFCAL << " "
                                        << setw(8) << setprecision(2) << BFCAL << " "
                                        << setw(8) << setprecision(2) << AFOBS
                                        << setw(8) << setprecision(2) << BFOBS
                                        << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                            else
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC << "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[1]
                                    << setw(8) << setprecision(3) << refs[IX].FWHM[2]
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS << endl;
                            }
                        }
                        file9 << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " K"
                            << setw(1) << IRC << " "
                            << setw(4) << IRH
                            << setw(4) << IRK
                            << setw(4) << IRL
                            << setw(10) << setprecision(4) << refs[IX].FMGNTD << endl;
                    }
                    else
                    {
                        if( (IXXX-1 % 60) == 0)
                        {
                            // alterei aqui para criar arquivo com dados para determinacao de estruturas e graficos ! cp 22abril01
                            if(IPC == 1)
                            {
                                file6 << "NO.  CODE    H   K   L   HW     POSN      ICALC     IOBS" << endl;
                                if(IPL2 != 0) file69 << "   H   K   L   HW     POSN      ICALC       IOBS   Bragg_Pos" << endl;
                            }
                            if(IPC == 2)
                            {
                                file6 << "NO.  CODE     H   K   L    FWHM    POSN    ICALC      IOBS      FCALC     FOBS   PHASE_C" << endl;
                                if(IPL2 != 0) file69 << "   H   K   L     FWHM    POSN    ICALC    IOBS      FCALC      FOBS   PHASE_C  Bragg_Pos" << endl;
                            }
                            if(IPC == 3)
                            {
                                file6 << "NO.  CODE     H   K   L     FWHM    POSN    ICALC      IOBS     A_CALC    B_CALC     A_OBS     B_OBS" << endl;
                                if(IPL2 != 0) file69 << "   H   K   L     FWHM    POSN    ICALC       IOBS     A_CAL  B_CALC      A_OBS     B_OBS  Bragg_Pos" << endl;
                            }
                        }
                        HW=refs[IX].REFS[1];
                        if (IPC == 1)
                        {
                            file6 << setw(4) << IXX
                                << setw(4) << IRC	<< "   "
                                << setw(4) << IRH
                                << setw(4) << IRK
                                << setw(4) << IRL
                                << setw(8) << setprecision(3) << HW
                                << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                << setw(10) << setprecision(0) << TIC
                                << setw(10) << setprecision(0) << TIOBS << endl;
                            if(IRC == 1)
                            {
                                if (IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << HW << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(0) << TIC << " "
                                        << setw(8) << setprecision(0) << TIOBS << " "
                                        << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                        }
                        else if (IPC == 2)
                        {
                            if (IRC == 1)
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC	<< "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(10) << setprecision(3) << TFCAL
                                    << setw(10) << setprecision(3) << TFOBS
                                    << setw(8) << setprecision(2) << RAD*refs[IX].APHASE;
                                if (IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << HW << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(1) << TIC << " "
                                        << setw(8) << setprecision(1) << TIOBS << " "
                                        << setw(8) << setprecision(1) << TFCAL << " "
                                        << setw(8) << setprecision(1) << TFOBS << " "
                                        << setw(9) << setprecision(4) << RAD*refs[IX].APHASE << " "
                                        << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                            else
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC	<< "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS << endl;
                            }
                        }
                        else if (IPC == 3)
                        {
                            if (IRC == 1)
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC	<< "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS
                                    << setw(10) << setprecision(3) << AFCAL
                                    << setw(10) << setprecision(3) << BFCAL
                                    << setw(10) << setprecision(3) << AFOBS
                                    << setw(10) << setprecision(3) << BFOBS << endl;
                                if (IPL2 != 0)
                                {
                                    file69 << setw(3) << IRH << " "
                                        << setw(3) << IRK << " "
                                        << setw(3) << IRL << " "
                                        << setw(8) << setprecision(4) << HW << " "
                                        << setw(8) << setprecision(4) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT << " "
                                        << setw(8) << setprecision(1) << TIC << " "
                                        << setw(8) << setprecision(1) << TIOBS << " "
                                        << setw(8) << setprecision(1) << AFCAL << " "
                                        << setw(8) << setprecision(1) << BFCAL << " "
                                        << setw(9) << setprecision(3) << AFOBS << " "
                                        << setw(9) << setprecision(3) << BFOBS << " "
                                        << setw(9) << setprecision(0) << (-K*diffractogram.XMAXINT/10) << endl;
                                }
                            }
                            else
                            {
                                file6 << setw(4) << IXX
                                    << setw(4) << IRC	<< "   "
                                    << setw(4) << IRH
                                    << setw(4) << IRK
                                    << setw(4) << IRL
                                    << setw(8) << setprecision(3) << HW
                                    << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                                    << setw(10) << setprecision(0) << TIC
                                    << setw(10) << setprecision(0) << TIOBS << endl;
                            }
                        }
                        file9 << setw(8) << setprecision(3) << refs[IX].REFS[2]+GLB_[1-1]+SHIFT
                            << " K" << setw(1) << IRC << " "
                            << setw(4) << IRH
                            << setw(4) << IRK
                            << setw(4) << IRL
                            << setw(10) << setprecision(0) << refs[IX].FMGNTD << endl;
                    }




        L481:;
                }
            }
            // version II - format for compatibility with Sakthevil's PLOT program
            for (I=1; I <= diffractogram.NPTS; ++I) file10 << diffractogram.KR[I];
            TDOBS=100.0*TDOBS/T2OBS;
            file6 << "DERIVED BRAGG R-FACTOR = " << setw(8) << setprecision(2) << TDOBS << endl;
            if (IPC != 1) file6 << "DERIVED R-F            = " << setw(8) << setprecision(2) << 100.0*RFNMR/RFDNR  << endl;
            FINAL_[ILOC_][2-(K % 2)] = TDOBS;
            if (NSIZESTRAIN  ==  9) size(K);
        }
    }




    file9
        << " NPTSZ" << setw(5) << diffractogram.NPTS << endl
        << " THMINZ" << setw(8) << setprecision(4) << diffractogram.THMIN << endl
        << " STEPZ" << setw(8) << setprecision(4) << diffractogram.STEP << endl
        << " YOBS    YCALCZ" << endl;
    for (I=1; I <= diffractogram.NPTS; ++I)
    {
        diffractogram.BK[I]=diffractogram.THMIN+static_cast<double>(I-1)*diffractogram.STEP;
        file9 << " "
            << setw(8) << setprecision(0) << diffractogram.Y[I]
        << setw(8) << setprecision(0) << diffractogram.YC[I] << endl;
        if(IPL2 != 0)
        {
            file690 << " " << setw(9) << setprecision(4) << diffractogram.THMIN+static_cast<double>(I-1)*diffractogram.STEP
                << " " << setw(9) << setprecision(0) << diffractogram.Y[I]
            << " " << setw(9) << setprecision(0) << diffractogram.YC[I]
            << " " << setw(9) << setprecision(0) << diffractogram.Y[I]-diffractogram.YC[I]
            << " " << setw(9) << setprecision(0) << ( pow(diffractogram.Y[I]-diffractogram.YC[I],2) )/diffractogram.Y[I]
            << endl;
        }
    }
    file9.close();
    file10.close();
    if(IPL2 != 0)
    {
        file69.close();
        file690.close();
    }
    //     IT BUILDS THE OBSERVED DATA FILE CORRECTED FOR ABSORPTION
    if(IPLOSS == 1 && IPBIG == 0)
    {
        file31.open("PLOTOSS.COR");
        file31 << " \"" << setw(70) << title << setw(8) << TITLE1 << endl;
        file31 << setw(6) << diffractogram.NPTS
            << setw(15) << setprecision(5) << diffractogram.STEP
            << setw(15) << setprecision(5) << diffractogram.THMIN
            << setw(15) << setprecision(5) << diffractogram.THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            file31 << setw(15) << setprecision(5) << diffractogram.Y[I] << endl;
        }
        file31.close();
    }
    //     IT BUILDS THE CALCULATED DATA FILE
    //     (BRAGG+COMPTON+DISORDINE+AMORPHOUS)
    if(IPLCAL == 1 && IPBIG == 0)
    {
        file32.open("PLOTCAL.TOT");
        file32 << " \"" << setw(70) << title << setw(8) << TITLE2 << endl;
        file32 << setw(6) << diffractogram.NPTS
            << setw(15) << setprecision(5) << diffractogram.STEP
            << setw(15) << setprecision(5) << diffractogram.THMIN
            << setw(15) << setprecision(5) << diffractogram.THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            file32 << setw(15) << setprecision(5) << diffractogram.YC[I] << endl;
        }
        file32.close();
    }
    //     IT BUILDS THE TOTAL COMPTON SCATTERING FILE
    //     FOR ALL PHASES
    if(IPLCOM == 1 && IPBIG == 0)
    {
        file33.open("PLOTCOM.TOT");
        file33 << " \"" << setw(70) << title << setw(8) << TITLE3 << endl;
        file33 << setw(6) << diffractogram.NPTS
            << setw(15) << setprecision(5) << diffractogram.STEP
            << setw(15) << setprecision(5) << diffractogram.THMIN
            << setw(15) << setprecision(5) << diffractogram.THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            file33 << setw(15) << setprecision(5) << diffractogram.BKCOM[I] << endl;
        }
        file33.close();
    }
    //     IT BUILDS THE TOTAL DISORDER SCATTERING FILE
    //     FOR ALL PHASES
    if(IPLDIS == 1 && IPBIG == 0)
    {
        file34.open("PLOTDIS.TOT");
        file34 << " \"" << setw(70) << title << setw(8) << TITLE4 << endl;
        file34 << setw(6) << diffractogram.NPTS
            << setw(15) << setprecision(5) << diffractogram.STEP
            << setw(15) << setprecision(5) << diffractogram.THMIN
            << setw(15) << setprecision(5) << diffractogram.THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            file34 << setw(15) << setprecision(5) << diffractogram.BKDIS[I] << endl;
        }
        file34.close();
    }
    //     IT BUILDS THE POLYNOMIAL BACKGROUND FILE
    if(IPLPOL == 1 && IPBIG == 0)
    {
        file37.open("PLOTPOL.TOT");
        file37 << " \"" << setw(70) << title << setw(8) << TITLE7 << endl;
        file37 << setw(6) << diffractogram.NPTS
            << setw(15) << setprecision(5) << diffractogram.STEP
            << setw(15) << setprecision(5) << diffractogram.THMIN
            << setw(15) << setprecision(5) << diffractogram.THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            file37 << setw(15) << setprecision(5) << diffractogram.BKPOL[I] << endl;
        }
        file37.close();
    }
    //     IT BUILDS THE AMORPHOUS FILE
    if(IPLAM == 1 && IPBIG == 0)
    {
        file36.open("PLOTAM.TOT");
        file36 << " \"" << setw(70) << title << setw(8) << TITLE6 << endl;
        file36 << setw(6) << diffractogram.NPTS
            << setw(15) << setprecision(5) << diffractogram.STEP
            << setw(15) << setprecision(5) << diffractogram.THMIN
            << setw(15) << setprecision(5) << diffractogram.THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            file36 << setw(15) << setprecision(5) << diffractogram.BKAM[I] << endl;
        }
        file36.close();
    }
    //     IT BUILDS THE TOTAL PLOT FILE
    if(IPBIG == 1)
    {
        file38.open("PLOTBIG.DAT");
        file38 << " \"" << setw(70) << title << setw(8) << TITLE8 << endl;

        file38 << "       ANG       OSS       CAL       +AM      +POL      +DIS      +COM       RES" << endl;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            XXXX=diffractogram.THMIN+static_cast<double>(I-1)*diffractogram.STEP;
            file38 << setw(10) << setprecision(3) << XXXX << ","
                << setw(10) << setprecision(3) << diffractogram.Y[I] << ","
                << setw(10) << setprecision(3) << diffractogram.YC[I] << ","
                << setw(10) << setprecision(3) << diffractogram.BKAM[I] << ","
                << setw(10) << setprecision(3) << diffractogram.BKAM[I]+diffractogram.BKPOL[I] << ","
                << setw(10) << setprecision(3) << diffractogram.BKAM[I]+diffractogram.BKPOL[I]+diffractogram.BKDIS[I] << ","
                << setw(10) << setprecision(3) << diffractogram.BKAM[I]+diffractogram.BKPOL[I]+diffractogram.BKDIS[I]+diffractogram.BKCOM[I] << ","
                << setw(10) << setprecision(3) << pow(diffractogram.Y[I]-diffractogram.YC[I],2)/diffractogram.Y[I]
            << endl;
        }
        file38.close();
    }
    if (IOT != 0)
    {
        NP=diffractogram.NPTS/240;
        for (N=1; N <= NP; ++N)
        {
            for (I=1; I <= 4; ++I) file6 << "2THETA   YOBS    YCALC VARIANCE ";
            file6 << endl;
            for (J=1; J <= 60; ++J)
            {
                for (I=1; I <= 4; ++I)
                {
                    file6 << setw(7) << setprecision(3) << diffractogram.BK[240*N-300+60*I+J]
                    << setw(8) << setprecision(0) << diffractogram.Y[240*N-300+60*I+J]
                    << setw(8) << setprecision(0) << diffractogram.YC[240*N-300+60*I+J]
                    << setw(8) << setprecision(0) << diffractogram.VAR[240*N-300+60*I+J];
                }
                file6 << endl;
            }
        }
        NPTS2=diffractogram.NPTS-NP*240;
        if (NPTS2 != 0)
        {
            NCOL=NPTS2/60;
            for (I=1; I <= 4; ++I) file6 << "2THETA   YOBS    YCALC VARIANCE ";
            NLINES=NPTS2-NCOL*60;
            if (NLINES != 0)
            {
                NCOL1=NCOL+1;
                NP=240*NP-60;
                for (J=1; J <= NLINES; ++J)
                {
                    for (I=1; I <= NCOL1; ++I)
                    {
                        file6 << setw(7) << setprecision(3) << diffractogram.BK[NP+60*I+J]
                        << setw(8) << setprecision(0) << diffractogram.Y[NP+60*I+J]
                        << setw(8) << setprecision(0) << diffractogram.YC[NP+60*I+J]
                        << setw(8) << setprecision(0) << diffractogram.VAR[NP+60*I+J];
                    }
                    file6 << endl;
                }
                NP=NP+NLINES;
                NLINES=60-NLINES;
            }
            for (J=1; J <= NLINES; ++J)
            {
                for (I=1; I <= NCOL1; ++I)
                {
                    file6 << setw(7) << setprecision(3) << diffractogram.BK[NP+60*I+J]
                    << setw(8) << setprecision(0) << diffractogram.Y[NP+60*I+J]
                    << setw(8) << setprecision(0) << diffractogram.YC[NP+60*I+J]
                    << setw(8) << setprecision(0) << diffractogram.VAR[NP+60*I+J];
                }
                file6 << endl;
            }
        }
    }

    if (IPL != 0)
    {
        for (I=1; I <= 12; ++I) LABEL[I]=10*I*ISCALE;
        file6 << "                              " << title << endl << endl;
        for (I=1; I <= 12; ++I) file6 << setw(10) << LABEL[I];
        file6 << endl;
        for (I=1; I <= 12; ++I) LABEL[I]=-60*IDIF+10*I*IDIF;
        for (I=1; I <= 12; ++I) file6 << setw(10) << LABEL[I];
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            IOUT = " ";
            for (J=1; J <= 120; ++J) IOUT=IOUT+ISPACE;
            IOUT[1]=IDOT;
            IOUT[61]=IDOT;
            IOUT[120]=IDOT;
            IY=static_cast<int>(diffractogram.Y[I])/ISCALE+2;
            IY=max(min(IY,120),2);
            IOUT[IY-1]=IBAR;
            IOUT[IY]=IPLUS;
            if(IY <= 119) IOUT[IY+1]=IBAR;
            IY=static_cast<int>(diffractogram.YC[I])/ISCALE+2;
            IY=max(min(IY,120),2);
            IOUT[IY]=IMINUS;
            IY=static_cast<int>(diffractogram.Y[I]-diffractogram.YC[I])/IDIF+61;
            IY=max(min(IY,120),2);
            IOUT[IY]=ISTAR;
            if((I-1 % 10) != 0)
            {
                for (J=1; J <= 120; ++J) file6 << IOUT[J];
                file6 << endl;
            }
            else
            {
                TLABEL=diffractogram.THMIN+diffractogram.STEP*static_cast<double>(I-1);
                for (J=1; J <= 113; ++J) file6 << IOUT[J];
                file6 << setw(6) << setprecision(2) << TLABEL;
                file6 << IOUT[120] << endl;
            }
        }
    }

    if (JOBTYP >= 3)
    {
        diffractogram.file4.close();
        diffractogram.file4b.open(diffractogram.file4name.data(),ios::trunc);			//  REWIND 4;
        diffractogram.file4b << setw(8) << setprecision(4) << diffractogram.THMIN
            << setw(8) << setprecision(4) << diffractogram.STEP
            << setw(8) << setprecision(4) << diffractogram.THMAX
            << endl;
        J = 0;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            ++J;
            diffractogram.file4b << setw(7) << setprecision(0) << " " << diffractogram.YC[I];
            if (J == 8)
            {
                diffractogram.file4b << endl;
                J = 0;
            }
        }
        if (diffractogram.file4b.is_open()) diffractogram.file4b.close();
    }

    if (IPLST != 0 && MAXS != 0)
    {
        INUMB = ICYRUN + 1;
        NPAGES = MAXS/12;
        if ((MAXS % 12) != 0) NPAGES = NPAGES+1;
        file6 << "PARAMETERS IN EACH CYCLE" << endl;

        if (file8o.is_open()) file8o.close();

        //     LIST PARAMETERS IN EACH CYCLE
        for (J=1; J <= NPAGES; ++J)
        {
            file8i.open(file8name.data());			// REWIND 8;
            ISTART = 1 + (J-1)*12;
            IFINIS  = min(MAXS,12 + (J-1)*12);
            file6 << "CYCLE";
            for (L=ISTART; L <= IFINIS; ++L) file6 << "   " << setw(2) << L << "     ";
            file6 << endl;
            for (K=1; K <= INUMB; ++K)
            {
                for (I=1; I <= 2*MAXS+4; ++I) file8i >> DUMMY[I];
                file6 << setw(2) << K-1 << ")  ";
                file6 << scientific;
                for (I=ISTART; I <= IFINIS; ++I) file6 <<  setw(10) << setprecision(4) << DUMMY[I] << " ";
                file6 << fixed << endl;
            }
            file8i.close();
        }
        //     LIST R-VALUES  IN EACH CYCLE
        file8i.open(file8name.data());			// REWIND 8;
        file6 << "R-VALUE VARIATION WITH CYCLE" << endl
            << "CYCLE    R-P      R-WP      S      D-W D" << endl;
        for (K=1; K <= INUMB; ++K)
        {
            for (I=1; I <= 2*MAXS+4; ++I) file8i >> DUMMY[I];
            file6 << setw(2) << K-1 << ")  ";
            for (I=2*MAXS+1; I <= 2*MAXS+4; ++I) file6 << setw(8) << setprecision(2) << DUMMY[I] << " ";
            file6 << endl;
        }
        file8i.close();

        //     LIST PARAMETER SHIFTS IN EACH CYCLE
        if (INUMB != 1)
        {
            file6 << "APPLIED PARAMETER SHIFT IN EACH CYCLE" << endl;
            for (J=1; J <= NPAGES; ++J)
            {
                file8i.open(file8name.data());			// REWIND 8;
                ISTART = 1 + MAXS+ (J-1)*12;
                IFINIS = min(2*MAXS,12 + MAXS+ (J-1)*12);
                file6 << "CYCLE";
                for (L=ISTART-MAXS; L <= IFINIS-MAXS; ++L) file6 << "   " << setw(2) << "     ";
                file6 << endl;
                for (K=1; K <= INUMB-1; ++K)
                {
                    for (I=1; I <= 2*MAXS+4; ++I) file8i >> DUMMY[I];
                    file6 << setw(2) << K << ")  ";
                    file6 << scientific;
                    for (I=ISTART; I <= IFINIS; ++I) file6 << setw(10) << setprecision(4) << DUMMY[I] << " ";
                    file6 << fixed << endl;
                }
            }
            file8i.close();
            file8o.open(file8name.data());
        }
    }

    //     CODE FOR PRINTING PARAMETERS AND STD. DEV. IN THE FINAL CYCLE
    if (IPLST == 2 && MAXS != 0)
    {
        file6 << "PARAMETERS AND STANDARD DEVIATIONS IN THE FINAL CYCLE FOR DATA BASE" << endl
            << " \" " << title << " \" " << endl;
        file6 << scientific << setw(10) << setprecision(4);
        for (I=1; I <= ILOC_; ++I)
        {
            file6 << FINAL_[I][1] << FINAL_[I][2];
        }
        file6 << fixed << endl;
    }
    return;
}

void DBWS::REFGEN(int IPHASE, double ZERO, double DIS, double TRANS, double PREFOR, int NCTR_[])
{
    const string LAU[14 +1] = {"",  "1BAR","2/M","MMM","4/M","4/MMM","3BAR   R","3BAR M R","3BAR","3BAR M 1","3BAR 1 M","6/M","6/MMM","M3","M3M"};
    int I, H1, H2, H3, I1, L1, I2, I3, IC,I1D, I2D, I3D,LX1, I12D, I13D, I23D,
        LXN,I123D,KXIS,I2DEL, I3DEL, I1MAX, I2MAX, I3MAX,IORDR1, IORDR2,IORDR3,IIPHAS;
    double ANGTTETA,SQ,RAD, SQH, POS,TLR,TAN2,XABS, TMIN, SMAX,TANX, PLOR,SMAX1,
        SHIFT,THMAX1,THMAXX;
    bool ORH1, ORH2, ORH3;

    POS = 0.0;

    RAD = M_PI/360.0;
    if ( phases[IPHASE].SYMB.NAXIS > 3 ) DBWSException("5001");
    KXIS = phases[IPHASE].SYMB.NAXIS;
    phases[IPHASE].SYMB.NAXIS = 1;
    I2D = 0;
    I3D = 0;
    I1D = 1;
    I12D = 1;
    I13D = 1;
    I123D = 1;
    I23D = 1;
    I2DEL = 1;
    I3DEL = 1;
    for (I=1; I <= 8; ++I)
    {
        if ( phases[IPHASE].SYMB.NCONT[I] == 0 ) break;
        if ( phases[IPHASE].SYMB.NCONT[I]-8192 < 0 )
        {
            if ( phases[IPHASE].SYMB.NCONT[I]-576 == 0 )
            {
                I3DEL = 2;
                I23D = 2;
            }
            else
            {
                if ( phases[IPHASE].SYMB.NCONT[I]-516 == 0 )
                {
                    I3DEL = 2;
                    I13D = 2;
                }
                else
                {
                    if ( phases[IPHASE].SYMB.NCONT[I]-68 == 0 )
                    {
                        I2DEL = 2;
                        I12D = 2;
                    }
                    else
                    {
                        if ( phases[IPHASE].SYMB.NCONT[I]-580  == 0)
                        {
                            I123D = 2;
                            I3DEL = 2;
                        }
                    }
                }
            }
        }
        else if ( phases[IPHASE].SYMB.NCONT[I]-8192 == 0 )
        {
            I3DEL = 3;
            I1D = 2;
            I123D = 3;
            break;
        }
    }
    TMIN = pow(( sin( (diffractogram.THMIN/720.0)*6.28318531 )/diffractogram.LAMDA[1]),2);
    IC = 0;
    if ( IPHASE >= 2 )
    {
        for (IIPHAS=2; IIPHAS <= IPHASE; ++IIPHAS)
        {
            IC = IC+phases[IIPHAS-1].ICR;
        }
    }
    SMAX = pow(( sin((diffractogram.THMAX/720.0)*6.28318531) /diffractogram.LAMDA[1]),2);

    //     **************************************************************
    //     SMAX IS CHANGED TO ACCOUNT FOR THE LEFT TAILS OF THE REFLECTIONS
    //     THAT ARE PRESENT AT ANGLES GREATER THAN THMAX
    //     **************************************************************
    //    PRINT *,'THMAX=',THMAX
    //      THMAX1=(U*(TAN(THMAX*RAD))**2+V*TAN(THMAX*RAD)+W+ZZZ*(1+
    //     *(TAN(THMAX*RAD))**2))
    // also incorporating the cotg^2 term  !cp may 01 97

    THMAX1 = (U_ * pow( (tan(diffractogram.THMAX*RAD)),2) +
              V_ * tan(diffractogram.THMAX*RAD)+
              W_+
              ZZZ_*
              (1+ pow( (tan(diffractogram.THMAX*RAD)),2) ) +
              UC_/(tan(diffractogram.THMAX*RAD)));
    if ( THMAX1 > 0.0 )
    {
        THMAX1= WDT_*sqrt(THMAX1)/2;
    }
    else
    {
        file6 << "  SQUARE OF FWHM NEGATIVE AT TWO-THETA" << setw(8) << setprecision(3) << POS << " PHASE NO. " << setw(4) << IPHASE << endl;
        DBWSException("SQUARE OF FWHM IS NEGATIVE");
    }
    ANGTTETA=diffractogram.THMAX+THMAX1;
    if ( (diffractogram.THMAX+THMAX1) >= 180.0 ) ANGTTETA=180.0;

    // TODO: mudar para ANGTTETA/2 * 2*pi/360
    SMAX1 = pow(( sin((ANGTTETA/720.0)*6.28318531)/diffractogram.LAMDA[1]),2);
    SMAX = SMAX1;
    THMAXX=ANGTTETA;
    if ( diffractogram.LAMDA[2] > diffractogram.LAMDA[1] ) TMIN=TMIN* pow( (diffractogram.LAMDA[1]/diffractogram.LAMDA[2]) , 2);
    if ( diffractogram.LAMDA[2] > 0.0  &&  diffractogram.LAMDA[2] < diffractogram.LAMDA[1] )SMAX=SMAX* pow((diffractogram.LAMDA[1]/diffractogram.LAMDA[2]) ,2);
    OP1(&IPHASE,NCTR_);
    file6
        << "LAUE SYMMETRY " << LAU[phases[IPHASE].SYMB.NSPGRP] << " WILL BE USED TO GENERATE INDICES" << endl
        << "       -------------------------------------------" << endl;
    I1MAX=static_cast<int>(cell_A*2.0*sqrt(SMAX));
    if ( KXIS < 3 )
    {
        I2MAX=static_cast<int>(cell_B*2.0*sqrt(SMAX));
        I3MAX=static_cast<int>(cell_C*2.0*sqrt(SMAX));
    }
    else
    {
        I2MAX = static_cast<int>(2.0*cell_C*sqrt(SMAX));
        I3MAX = static_cast<int>(2.0*cell_B*sqrt(SMAX));
    }
    L1 = phases[IPHASE].SYMB.NSPGRP;
    if ( L1 >= 13 && phases[IPHASE].SYMB.NCONT[1] == 580 ) L1=L1+2;
    if ( I12D+I13D+I23D == 5 )
    {
        I12D = 2;
        I13D = 1;
        I23D = 2;
        I2DEL = 2;
        I3DEL = 2;
    }
    if ( L1 == 6 ) I2D=1;
    if ( L1 == 7 ) I2D=2;
    I3D = I2D;
    I1 = -1;
L2200:
    I1 = 1+I1;
    if ( I1-I1MAX <= 0 )
    {
        H1 = I1;
        if ( I2D > 0 ) I2MAX=I1;
        switch (L1) {
        case 1:
            I2= -min(I2MAX,I1*I2MAX);
            break;
        case 2:
        case 3:
        case 9:
        case 11:
            I2 = 0;
            break;
        case 4:
            I2= min(I1,1);
            break;
        case 5:
        case 10:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
            I2 = I1;
            break;
        case 6:
            I2 = -2*I1;
            break;
        case 7:
            I2 = -I1/2;
            break;
        case 8:
            I2= min(-I1+1,0);
            break;
        default:
            GOTOER();
            break;
        }
    }
    else
    {
        goto L2303;
    }
    I2= I2DEL*(I2/I2DEL)+(I1 % I12D);
    goto L2221;
L2220:
    I2 = I2DEL+I2;
L2221:
    H2 = I2;
    if ( I2-I2MAX <= 0 )
    {
        switch (L1) {
        case 1:
        case 2:
        case 10:
            I3= -min(I3MAX,(I1+abs(I2))*I3MAX); //xxxxxxxxxxxxxxxxxxxxxxxxx
            I3= I3DEL*(I3/I3DEL)+  ((((I1+I1D*I2) % I123D)+I123D) % I123D) + (I1 % I13D)+ (I2 % I23D);
            if ( I3D-1 < 0 )
            {
            }
            else if ( I3D-1 == 0 )
            {
                I3MAX = I1;
                if ( I2-I1 != 0 )
                {
                    I3MAX = I3MAX-1;
                }
            }
            else
            {
                I3MAX = I2;
            }
            break;
        case 3:
        case 4:
        case 5:
        case 8:
        case 9:
        case 11:
        case 12:
            I3 = 0;
            I3= I3DEL*(I3/I3DEL)+  ((((I1+I1D*I2) % I123D)+I123D) % I123D) + (I1 % I13D)+ (I2 % I23D);
            if ( I3D-1 < 0 )
            {
            }
            else if ( I3D-1 == 0 )
            {
                I3MAX = I1;
                if ( I2-I1 != 0 )
                {
                    I3MAX = I3MAX-1;
                }
            }
            else
            {
                I3MAX = I2;
            }
            break;
        case 6:
        case 7:
            I3 = -I1-I2;
            I3= I3DEL*(I3/I3DEL)+  ((((I1+I1D*I2) % I123D)+I123D) % I123D) + (I1 % I13D)+ (I2 % I23D);
            if ( I3D-1 < 0 )
            {
            }
            else if ( I3D-1 == 0 )
            {
                I3MAX = I1;
                if ( I2-I1 != 0 )
                {
                    I3MAX = I3MAX-1;
                }
            }
            else
            {
                I3MAX = I2;
            }
            break;
        case 13:
            I3= min(I2,I1+I3DEL);
            I3= I3DEL*(I3/I3DEL)+  ((((I1+I1D*I2) % I123D)+I123D) % I123D) + (I1 % I13D)+ (I2 % I23D);
            if ( I3D-1 < 0 )
            {
            }
            else if ( I3D-1 == 0 )
            {
                I3MAX = I1;
                if ( I2-I1 != 0 )
                {
                    I3MAX = I3MAX-1;
                }
            }
            else
            {
                I3MAX = I2;
            }
            break;
        case 14:
            I3 = I2;
            I3= I3DEL*(I3/I3DEL)+  ((((I1+I1D*I2) % I123D)+I123D) % I123D) + (I1 % I13D)+ (I2 % I23D);
            if ( I3D-1 < 0 )
            {
            }
            else if ( I3D-1 == 0 )
            {
                I3MAX = I1;
                if ( I2-I1 != 0 )
                {
                    I3MAX = I3MAX-1;
                }
            }
            else
            {
                I3MAX = I2;
            }
            break;
        case 15:
            I3 = I1+2-(I2 % 2);
            if ( I1 == I2 && (I1 % 2) == 0 ) I3=I1;
            break;
        case 16:
            I3= I2+(I1 % 2);
            break;
        default:
            GOTOER();
            break;
        }
        H3 = I3;
        if ( KXIS == 3 )
        {
            H3 = I2;
            H2 = I3;
        }
        if ( I3-I3MAX <= 0 )
        {
            goto L113;
        }
        else
        {
            goto L2220;
        }
    }
    else
    {
        goto L2200;
    }



//L2240:;

L113:
    SQ = static_cast<double>(H1*H1)*cell_AL[1][1]+
         static_cast<double>(H2*H2)*cell_AL[2][2]+
         static_cast<double>(H3*H3)*cell_AL[3][3]+
         static_cast<double>(H1*H2)*cell_AL[1][2]+
         static_cast<double>(H1*H3)*cell_AL[1][3]+
         static_cast<double>(H2*H3)*cell_AL[2][3];
    SQ=SQ/4.0;
    //L2234:
    if ( SQ-SMAX <= 0 )
    {
        goto L3000;

    }
    else
    {
        I3 = I3DEL+I3;
        H3 = I3;
        if ( KXIS == 3 )
        {
            H3 = I2;
            H2 = I3;
        }
        if ( I3-I3MAX <= 0 )
        {
            goto L113;
        }
        else
        {
            goto L2220;
        }
    }
L2303:
    phases[IPHASE].ICR=IC;
    if (IPHASE >= 2)
    {
        for (IIPHAS=2; IIPHAS <= IPHASE; ++IIPHAS) phases[IPHASE].ICR = phases[IPHASE].ICR-phases[IIPHAS-1].ICR;
    }
    SORT(IPHASE);
    return;
L3000:
    if ( SQ - TMIN < 0 )
    {
        I3 = I3DEL+I3;
        H3 = I3;
        if ( KXIS == 3 )
        {
            H3 = I2;
            H2 = I3;
        }
        if ( I3-I3MAX <= 0 )
        {
            goto L113;
        }
        else
        {
            goto L2220;
        }
    }
    else
    {
        goto L3117;
    }
L3117:
    if ( IC > IRS-2 ) goto L6001;
    //     NEXT if BLOCK FOR PHASE WITH DISTINCT ORIENTATION ALONG PREF(NAXIS,I)
    if (PREFOR > 99.0)
    {
        ORH1 = static_cast<int>(phases[IPHASE].PREF[1]) == H1;
        if ( !(ORH1) && H1 != 0 && static_cast<int>(phases[IPHASE].PREF[1]) != 0) ORH1 = (H1 % static_cast<int>(phases[IPHASE].PREF[1])) == 0;
        if ( !ORH1 )
        {
            I3 = I3DEL+I3;
            H3 = I3;
            if ( KXIS == 3 )
            {
                H3 = I2;
                H2 = I3;
            }
            if ( I3-I3MAX <= 0 )
            {
                goto L113;
            }
            else
            {
                goto L2220;
            }
        }

        ORH2 = static_cast<int>(phases[IPHASE].PREF[2]) == H2;
        if ( !(ORH2) && H2 != 0 && static_cast<int>(phases[IPHASE].PREF[2]) != 0) ORH2 = (H2 % static_cast<int>(phases[IPHASE].PREF[2])) == 0;
        if ( !ORH2 )
        {
            I3 = I3DEL+I3;
            H3 = I3;
            if ( KXIS == 3 )
            {
                H3 = I2;
                H2 = I3;
            }
            if ( I3-I3MAX <= 0 )
            {
                goto L113;
            }
            else
            {
                goto L2220;
            }
        }

        ORH3 = static_cast<int>(phases[IPHASE].PREF[3]) == H3;
        if ( !(ORH3) && H3 != 0 && static_cast<int>(phases[IPHASE].PREF[3]) != 0)  ORH3 = (H3 % static_cast<int>(phases[IPHASE].PREF[3])) == 0;
        if ( !ORH3 )
        {
            I3 = I3DEL+I3;
            H3 = I3;
            if ( KXIS == 3 )
            {
                H3 = I2;
                H2 = I3;
            }
            if ( I3-I3MAX <= 0 )
            {
                goto L113;
            }
            else
            {
                goto L2220;
            }
        }

        IORDR1=0;
        IORDR2=0;
        IORDR3=0;
        if (H1 != 0) IORDR1 = static_cast<int>(H1/static_cast<int>(phases[IPHASE].PREF[1]));
        if (H2 != 0) IORDR2 = static_cast<int>(H2/static_cast<int>(phases[IPHASE].PREF[2]));
        if (H3 != 0) IORDR3 = static_cast<int>(H3/static_cast<int>(phases[IPHASE].PREF[3]));
        if (IORDR1 == IORDR2 && IORDR2 == IORDR3) goto L9257;
        if (IORDR1 == 0 && (IORDR2 == IORDR3)) goto L9257;
        if (IORDR2 == 0 && (IORDR3 == IORDR1)) goto L9257;
        if (IORDR3 == 0 && (IORDR1 == IORDR2)) goto L9257;
        if (IORDR1 == 0 && IORDR2 == 0) goto L9257;
        if (IORDR2 == 0 && IORDR3 == 0) goto L9257;
        if (IORDR3 == 0 && IORDR1 == 0) goto L9257;
    }
L9257:
    //     SEPARATE PHASE FOR ORIENTATION COMPLETE
    IHKL[1][1]=H1;
    IHKL[2][1]=H2;
    IHKL[3][1]=H3;
    AZ[1]=0.0;
    IER=0;
    SMTRY2(&IPHASE);
    if(IER != 0)
    {
        I3 = I3DEL+I3;
        H3 = I3;
        if ( KXIS == 3 )
        {
            H3 = I2;
            H2 = I3;
        }
        if ( I3-I3MAX <= 0 )
        {
            goto L113;
        }
        else
        {
            goto L2220;
        }
    }
    LXN=2;
    if(diffractogram.LAMDA[2] == 0.0 || diffractogram.LAMDA[2] == diffractogram.LAMDA[1])LXN=1;
    for (LX1=1; LX1 <= LXN; ++LX1)
    {
        SQH=SQ*diffractogram.LAMDA[LX1]*diffractogram.LAMDA[LX1];
        TAN2=SQH/(1.0-SQH);
        if (SQH >= 1.0) goto L3118;
        TANX=sqrt(TAN2);
        //     SHIFT DUE TO SAMPLE DISPLACEMENT AND TRANSPARENCY
        SHIFT =  DIS*sqrt(1-SQH)+TRANS*sqrt(1.-(1.-2.*SQH)*(1.-2.*SQH));
        POS=atan(TANX)/RAD+ ZERO + SHIFT;
        if ( POS > THMAXX  ||  POS < diffractogram.THMIN ) goto L3118;
        IC=IC+1;
        XABS=1.0;
        if(TMV_ <= 0.000001)XABS=1.0;
        PLOR=1.0/(2.0*SQH*sqrt(1.0-SQH))*XABS;
        if (INSTRM == 2)
        {
            PLOR = PLOR * (0.95+0.05*(1.-2.*SQH)*(1.-2.*SQH));
            goto L4000;
        }
        if(JOBTYP == 3)PLOR=PLOR*(1.+(1.-2.*SQH)*(1.-2.*SQH)*CTHM_);
        if(JOBTYP == 1)PLOR=PLOR*(1.+(1.-2.*SQH)*(1.-2.*SQH)*CTHM_);
L4000:
        refs[IC].IREFS=256*(256*(256*(8*IPHASE+LX1)+128+H1)+128+H2)+128+H3;
        refs[IC].FMGNTD=MULT(IPHASE,H1,H2,H3,KXIS);

        //------CALCULATE FWHM FOR PSEUDOVOIGT WITH GAUSS AND LORENTZ
        if (NPROF == _TCHZ)
        {
            refs[IC].HALFG = (U_*TAN2+V_*TANX+W_+ZZZ_*(1+TAN2));
            if (refs[IC].HALFG > 0.)
            {
                refs[IC].HALFG = sqrt(refs[IC].HALFG);
            }
            else
            {
                cout << "3" << endl;
                file6 << "  SQUARE OF FWHM NEGATIVE AT TWO-THETA" << setw(8) << setprecision(3) << POS << "FOR PHASE NO. " << setw(IPHASE) << endl;
                DBWSException("SQUARE OF FWHM IS NEGATIVE");
            }
            refs[IC].HALFL = ULOR_*TANX+VLOR_/sqrt(1.0-SQH);
            refs[IC].REFS[1] = pow(( pow(refs[IC].HALFG,5.0) +2.69269*pow(refs[IC].HALFG,4.0)*refs[IC].HALFL+2.42843*pow(refs[IC].HALFG,3.0)*pow(refs[IC].HALFL,2.0)+4.47163*pow(refs[IC].HALFG,2.0)*pow(refs[IC].HALFL,3.0)+0.07842*refs[IC].HALFG*pow(refs[IC].HALFL,4.0)+pow(refs[IC].HALFL,5.0)),0.2);
            TLR = refs[IC].HALFL/refs[IC].REFS[1];
            refs[IC].GAM = 1.36603*TLR-0.47719*TLR*TLR+0.11116* pow(TLR,3.0);
        }
        else if (NPROF == _SplitPearsonVII)
        {
            refs[IC].REFS[1] = (U_*TAN2+V_*TANX+W_);
        }
        else
        {
            // incorporating cotg^2  !cp may 01 97
            refs[IC].REFS[1]=(U_*TAN2+V_*TANX+W_+ZZZ_*(1+TAN2)+UC_/TAN2);
        }
        if (refs[IC].REFS[1] > 0.0)
        {
            refs[IC].REFS[1]= sqrt(refs[IC].REFS[1]);
        }
        else
        {
            cout << "4" << endl;
            file6 << "  SQUARE OF FWHM NEGATIVE AT TWO-THETA" << setw(8) << setprecision(3) << POS << "FOR PHASE NO. " << setw(IPHASE) << endl;
            DBWSException("SQUARE OF FWHM IS NEGATIVE");
        }
        //L7000:
        refs[IC].REFS[2]=POS;
        refs[IC].REFS[3]=PLOR;
L3118:;
    }
    I3 = I3DEL+I3;
    H3 = I3;
    if ( KXIS == 3 )
    {
        H3 = I2;
        H2 = I3;
    }
    if ( I3-I3MAX <= 0 )
    {
        goto L113;
    }
    else
    {
        goto L2220;
    }
L6001:
    file6 << "TOO MANY REFLECTIONS (" << setw(6) << IRS << "). INCREASE *IRS* IN THE PARAMETER STATEMENT IN THE SOURCE CODE " << endl;
    cout  << "TOO MANY REFLECTIONS (" << setw(6) << IRS << "). INCREASE *IRS* IN THE PARAMETER STATEMENT IN THE SOURCE CODE " << endl;
    DBWSException("");
}

//  SUBROUTINE FINDC and SUBROUTINE COMPTON subroutine DISORDER: by Canton et all.
// Added by CPS between March-May 1997
void DBWS::FINDC(int K, int NSCAT)
{
    // ----THIS SUBROUTINE ASSIGNS AUTOMATICALLY TO ATOMS
    //     THE RIGHT CONSTANTS TO CALCULATE THE COMPTON SCATTERING

    const char PIU = '+';
    const char MENO = '-';

    int I, J, L, NA, NN, NS, IOF/*, NSAVE*/;
    string NOM,NOME;

    IOF = 0;
    NS  = 0;
    //-----K = NUMBER OF PHASE UNDER CONSIDERATION
    if(K > 1)
    {
        //-----CALCULATE IOF = ALL ATOMS OF THE K-1 PHASES
        for (I = 2; I <= K; ++I) IOF = IOF + phases[I-1].AtomCount;
        //NS = NSAVE;
        NS = NSAVE;
    }


    for (I = 1; I <= phases[K].AtomCount; ++I)
    {
        for (J = 1; J <= NSCAT; ++J) if(NTYP[I+IOF] == NAM_[J]) goto L30;
        //-----212 = ALL THE POSSIBLE NAMES OF ATOMS AND IONS
        //L25:
        for (J = 1; J <= 212; ++J) if(NTYP[I+IOF] == TBXC[J]) goto L50;
        file6 << "COMPTON SCATTERING COEFFICIENT NOT FOUND FOR " << NTYP[I+IOF] << endl;
        DBWSException("COMPTON SCATTERING DATA MISSING");
L30:
        PTC_[I+IOF] = J;
        goto L9999;
L50:
        NOME = TBXC[J];
        //-----FIND NA = THE ATOMIC NUMBER OF I-TH ATOM
        if (J == 1 || J == 2 || J == 3)
        {
            NA = static_cast<int>(TBX[J][10]);
        }
        else
        {
            NA = static_cast<int>(TBX[J][10] + 1.0);
        }
        NS = NS + 1;
        PTC_[I+IOF] = NS;
        //-----PUT IN CC THE 4 COMPTON COEFFICIENTS
        for (L = 1; L <= 4; ++L) CC_[L][NS] = TCS[NA][L];
        NAM_[NS] = NOME;

        NOM = string(NOME,0,4);
        //     FIND IN WHAT COLUMN THERE IS + OR -
        for (L=1; L <= 4; ++L)
        {
            if(NOM[L] == PIU)  goto L90;
            if(NOM[L] == MENO) goto L95;
        }
        // CASE WITH NOR + NOR -
        ZEFF_[NS] = NA;
        goto L9999;
        //-----CASE WITH PLUS
L90:
        L = L + 1;
        if (L > 4) DBWSException("SOMETHING IS WRONG IN ATOMIC NAME");
        NN = NOM[L];
        ZEFF_[NS] = NA - NN;
        goto L9999;
        //-----CASE WITH MINUS
L95:
        L = L + 1;
        if (L > 4) DBWSException("SOMETHING IS WRONG IN ATOMIC NAME");
        NN = NOM[L];
        ZEFF_[NS] = NA + NN;
L9999:;
    }
    //NSAVE = NS;
    NSAVE = NS;
}

//  SUBROUTINE ABSORP and SUBROUTINE ARIA: by Canton et all. Added by cps between
//  march-may 1997
void DBWS::ABSORP(double MU, double SW, double TH, double* ABC)
{
    //-----THIS SUBROUTINE CORRECTS THE EXPERIMENTAL INTENSITIES FOR THE
    //                ABSORPTION EFFECTS AS REPORTED BY:
    //     1) H. P. KLUG & L. E. ALEXANDER, X-RAY DifFRACTION PROCEDURES,
    //        1970, PAG.487.
    //     2) A. IMMIRZI, ACTA CRYST. , 1980, B36, 2378-2385.
    //
    //     IT IS WRITTEN TAKING INTO ACCOUNT THE SYMMETRIC REFLECTION ARRANGEMENT
    //     ( KLUG & ALEXANDER, 1970, FIG. 5-52 PAG. 390) AND  THE  FINITE
    //     THICKNESS OR WIDTH OF THE SLAB, IN THIS SITUATION THE ABSORPTION
    //     IS GENERALLY LOW AND INCREASES SLIGHTLY WITH 2 THETA.
    //     MU = LINEAR ABSORPTION COEFFICIENT IN CM-1.
    //     SW = SAMPLE THICKNESS IN CM.

    double EX;
    EX  = ( 2.0 * MU * SW ) / sin(TH * 0.008726646);
    *ABC = 1.0 - exp(-EX);
}







// suboutine (QPAINIT) to compute mass fractions before starting the refinement
void DBWS::qpainit(void)
{
    int I,N,IP,IOF,ICOCO,IIPHAS,IINNMOL;
    double V0,FT, XFAC, SFIS, ARGCOS, WTOTAL;
    double VOL[99+1];
    double W[99+1];
    double XMASS[99+1];
    double FR[99+1];
    double FRP[99+1];

    file6 << "       >>> QPA before starting the refinement <<<" << endl
          << "       >>                                    <<<" << endl;
    for (IP=1; IP <= NPHASE; ++IP)
    {
        DIRECT(DCSM,DCV,&IP);
        TMASSA[IP]=0.0;
        IOF=0;
        if(IP > 1)
        {
            for (IIPHAS=2; IIPHAS <= IP; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
        }
        N=phases[IP].AtomCount;
        for (I=1; I <= N; ++I)
        {
            ICOCO=PTR_[I+IOF];
            TMASSA[IP] = TMASSA[IP] + XL_[I+IOF][5]*coeff.XMAS[ICOCO]*phases[IP].XMLTP;
        }
        XFAC = M_PI / 180.000000;
        DCV[4] = XFAC * DCV[4];
        DCSM[4][4] = DCSM[4][4] * XFAC;
        DCV[5] = XFAC * DCV[5];
        DCSM[5][5] = DCSM[5][5] * XFAC;
        DCV[6] = DCV[6] * XFAC;
        DCSM[6][6] = DCSM[6][6] * XFAC;

        //-----Calculations of VOLUME and SVZM (=W) for each phase
        ARGCOS= 1-pow((cos(DCV[4])),2)-pow((cos(DCV[5])),2)-pow((cos(DCV[6])),2) + 2 * (cos(DCV[4])) * (cos(DCV[5])) * (cos(DCV[6]));
        V0 = DCV[1] * DCV[2] * DCV[3];
        VOL[IP] = V0 * sqrt(ARGCOS);
        //		VOSQ = 0.5*VOL[IP]/ARGCOS;
        //		ARG1 = VOSQ*(2 * cos(DCV[4]) * sin(DCV[4]) - 2*sin(DCV[4]) *cos(DCV[5]) *cos(DCV[6])) * DCSM[4][4];
        //		ARG2 = VOSQ*(2 * cos(DCV[5]) * sin(DCV[5]) - 2*sin(DCV[5]) *cos(DCV[4]) *cos(DCV[6])) * DCSM[5][5];
        //		ARG3 = VOSQ*(2 * cos(DCV[6]) * sin(DCV[6]) - 2*sin(DCV[6]) *cos(DCV[4]) *cos(DCV[5])) * DCSM[6][6];
        W[IP] = phases[IP].PAR[0] * TMASSA[IP] * VOL[IP]/phases[IP].SAQF;
        file6 << "Volume("
            << setw(2) << IP << ")= "
            << setw(9) << setprecision(3) << VOL[IP]
        << " UCW= "
            << setw(7) << setprecision(2) << TMASSA[IP]
        << " U.C.Density = "
            << setw(7) << setprecision(3) << 1.66113*TMASSA[IP]/VOL[IP]
        << " gr/cm^3" << endl;
    }
    // ****** QUANTITATIVE ANALYSIS ***************
    file6 << endl;
    WTOTAL = 0.000000;
    for (I = 1; I <= NPHASE; ++I) WTOTAL = WTOTAL + W[I];
    for (I = 1; I <= NPHASE; ++I) XMASS[I] = 100.0 * W[I] / WTOTAL;
    IINNMOL = 0;
    for (I = 1; I <= NPHASE; ++I)
    {
        if (IINNMOL != 1)
        {
            if (phases[I].NMOL == 0) IINNMOL=1;
        }
    }
    if (IINNMOL == 1)
    {
        for (I = 1; I <= NPHASE; ++I)
        {
            // ** printing results
            file6 << "PHASE = "
                << setw(2) << I
                << " => %MASS = "
                << setw(6) << setprecision(2) << XMASS[I]
            << "%MOLAR = NOT COMPUTED" << endl;
        }
    }
    else
    {
        // ****    CALCULATION OF MOLAR FRACTION  ****
        FT = 0.0000000;
        for (I = 1; I <= NPHASE; ++I)
        {
            FRP[I] = XMASS[I] * phases[I].NMOL / TMASSA[I];
            FT = FT + FRP[I];
        }
        for (I = 1; I <= NPHASE; ++I)
        {
            FR[I] = 100.0 * FRP[I] / FT;
            // ** printing results
            file6 << "PHASE = "
                << setw(2) << I
                << " => %MASS = "
                << setw(6) << setprecision(2) << XMASS[I]
            << "  %MOLAR = "
                << setw(6) << setprecision(2) << FR[I] << endl;
        }
    }
    if(ISPHASE != 0)
    {
        file6 << endl << "Considering Amorphous Content:" << endl;
        SFIS=phases[ISPHASE].WTIS/XMASS[ISPHASE];
        if(SFIS > 1.0)
        {
            file6 << "PROBLEM:Amount of Internal Standard (Phase #"
                << setw(2) << ISPHASE
                << ") is less than the specified "
                << setw(6) << setprecision(2) << phases[ISPHASE].WTIS
            << "%." << endl
                << "Amorphous content not computed. Check ISWT in line 11.2 for this phase" << endl;
        }
        else
        {
            for (I=1; I <= NPHASE; ++I)
            {
                file6 << "PHASE = " << setw(2) << I << " => %MASS = " << setw(6) << setprecision(2) << XMASS[I]*SFIS << endl;
            }
            file6 << "AMORPHOUS  => %MASS = " << setw(6) << setprecision(2) <<100*(1.0-SFIS)  << endl;
        }
    }
    file6 << endl;
    return;
}









void DBWS::INPTR(void)
{
    int MULTX_;
    double Y_[192+1][3+1];
    int NCTR_[192+1];

    const string UPPER = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const string LOWER = " abcdefghijklmnopqrstuvwxyz";


    int I, J, K, N,IK, KK,IX, IOF, IRH, ICY, ICZ, ITN, KKS,
        IRL, IXX, IYY, IRK, IRC,IPTS,MLTT,IXDEL,IXRAY, ITIPO,
        ISTOP,IIPHAS,NATOMS, LCOUNT, ISTEST;
    //IBDG
    double X,TH,ABC,THX;
    string SPG,s;

    int ISIMOP_;










    IPTS = 0;


    // line 1
    getline(file5,s);
    title = s.substr(0,70);

    // line 2
    getline(file5,s);

    stringstream(s.substr( 0*4,4)) >> JOBTYP;
    stringstream(s.substr( 1*4,4)) >> NPROF;
    stringstream(s.substr( 2*4,4)) >> NPHASE;
    stringstream(s.substr( 3*4,4)) >> diffractogram.NBCKGD;
    stringstream(s.substr( 4*4,4)) >> diffractogram.NEXCRG;
    stringstream(s.substr( 5*4,4)) >> diffractogram.NSCAT;
    stringstream(s.substr( 6*4,4)) >> INSTRM;
    stringstream(s.substr( 7*4,4)) >> IPREF;
    stringstream(s.substr( 8*4,4)) >> IASYM;
    stringstream(s.substr( 9*4,4)) >> IABSR;
    stringstream(s.substr(10*4,4)) >> IDATA;
    stringstream(s.substr(11*4,4)) >> ISPHASE;

    if (ISPHASE > NPHASE)
    {
        file6 << "Internal Standard Phase does not exist. Check its number in line 2, column 12." << endl
              << "No Internal Standard will be used in the QPA" << endl;
        ISPHASE = 0;
    }

    if (NPROF == 9)
    {
        NPROF = _TCHZ;
        NSIZESTRAIN = 9;
    }
    if (diffractogram.NBCKGD == -1)
    {
        IBCKCODE=diffractogram.NBCKGD;
        diffractogram.NBCKGD = 0;
        //		IBDG   = 0;
    }
    else
    {
        FONDO  = 0;
        IBGD   = 1;
        IBCKCODE=diffractogram.NBCKGD;
    }
    if(NPHASE == 0) NPHASE=1;
    INSTRM=INSTRM+1;
    JOBTYP=JOBTYP+1;
    //NPROF=NPROF+1;

    file6 << "RIETVELD ANALYSIS PROGRAM OPENDBWS," << endl
          << "COPYRIGHT 2013 BY VEGNER UTUNI." << endl << endl;
    file6 << "PROGRAM PARAMETERS:" << endl
        << "IDSZ=" << setw(5) << IDSZ
        << "    IRS=" << setw(5) << IRS
        << "    NATS=" << setw(4) << NATS
        << "     MSZ=" << setw(3) << MSZ
        << "     NOV=" << setw(5) << NOV << endl;
    file6 << title << endl;
    cout << "RIETVELD ANALYSIS PROGRAM OPENDBWS," << endl
         << "COPYRIGHT 2013 BY VEGNER UTUNI. v3" << endl << endl;
    cout  << "PROGRAM PARAMETERS:" << endl
        << "IDSZ=" << setw(5) << IDSZ
        << "    IRS=" << setw(5) << IRS
        << "    NATS=" << setw(4) << NATS
        << "     MSZ=" << setw(3) << MSZ
        << "     NOV=" << setw(5) << NOV << endl;
    cout  << title << endl;

    switch (JOBTYP) {
    case 1:
        if (INSTRM == 2)
        {
            file6 << "COLLECTED IN SYNCHROTRON AT NSLS OR SRS" << endl;
        }
        else
        {
            file6 << "FOR X-RAY DATA" << endl;
        }
        break;
    case 2:
        if (INSTRM == 2)
        {
            file6 << "VARYING NO. OF COUNTERS AT EACH STEP" << endl;
        }
        else
        {
            file6 << "FOR NEUTRON DATA, NUCLEAR INTENSITIES ONLY" << endl;
        }
        break;
    case 3:
        file6 << "PATTERN CALCULATION,XRAY" << endl;
        break;
    case 4:
        file6 << "PATTERN CALCULATION,NEUTRON" << endl;
        break;
    default:
        DBWSException("7777");
        break;
    }

    switch (IDATA) {
    case 0:
        file6 << "READ DATA IN TRADITIONAL DBWS FORMAT" << endl;
        break;
    case 1:
        file6 << "READ DATA IN FREE FORMAT" << endl;
        break;
    case 2:
        file6 << "READ DATA IN GSAS STD FORMAT" << endl;
        break;
    case 3:
        file6 << "READ DATA IN PHILIPS UDF FORMAT" << endl;
        break;
    case 4:
        file6 << "READ DATA IN SCINTAG TXT FORMAT" << endl;
        break;
    case 5:
        file6 << "READ DATA IN SIEMENS UXD FORMAT" << endl;
        break;
    case 6:
        file6 << "READ DATA IN RIGAKU ASC FORMAT" << endl;
        break;
    default:
        break;
    }

    file6 << "NUMBER OF PHASES= " << setw(4) << NPHASE << endl
          << "NUMBER OF EXCLUDED REGIONS= " << setw(4) << diffractogram.NEXCRG << endl
          << "NUMBER OF SCATTERING SETS= " << setw(4) << diffractogram.NSCAT << endl;

    if (diffractogram.NBCKGD-1 < 0)
    {
        file6 << "BACKGROUND TO BE REFINED" << endl;
    }
    else if(diffractogram.NBCKGD-1 == 0)
    {
        file6 << "BACKGROUND DATA TO BE READ FROM FILE" << endl;
    }
    else
    {
        file6 << "BACKGROUND CORRECTION BY INTERPOLATION BETWEEN THE " << setw(4) << diffractogram.NBCKGD << " POINTS GIVEN" << endl;
    }


    switch (NPROF) {
    case _Gaussian:
        file6 << "GAUSSIAN PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    case _Lorentzian:
        file6 << "LORENTZIAN (CAUCHY) PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    case _Mod1:
        file6 << "MOD 1 LORENTZIAN PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    case _Mod2:
        file6 << "MOD 2 LORENTZIAN PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    case _SplitPearsonVII:
        file6 << "SPLIT PEARSON VII PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    case _pseudoVoigt:
        file6 << "PSEUDO-VOIGT (PV) PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    case _PearsonVII:
        file6 << "PEARSON VII PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    case _TCHZ:
        file6 << "THOMPSON-COX-HASTINGS (PV) PROFILE, NPROF = " << setw(4) << NPROF << endl;
        break;
    }

    if (IPREF == 0)
    {
        file6 << "IPREF=0, UDA-RIETVELD PREFERRED ORIENTATION FUNCTION" << endl;
    }
    else
    {
        file6 << "IPREF=1, MARCH-DOLLASE PREFERRED ORIENTATION FUNCTION" << endl;
    }

    //  SURFACE ROUGHNESS
    switch (IABSR) {
    case 1:
        file6 << "IABSR=1, CORRECTION OF SURFACE ROUGHNESS BY YOUNG" << endl;
        break;
    case 2:
        file6 << "IABSR=2, CORRECTION OF SURFACE ROUGHNESS BY, SPARKS, ET AL." << endl;
        break;
    case 3:
        file6 << "IABSR=3, CORRECTION OF SURFACE ROUGHNESS BY SUORTTI" << endl;
        break;
    case 4:
        file6 << "IABSR=4, CORRECTION OF SURFACE ROUGHNESS BY PITSCHKE, HERMANN, AND MATTERN" << endl;
        break;
    }

    //  asymmetry correction  Test for asymmetry model included !cp ap 97
    if (IASYM == 0) file6 << "IASYM=0, Usual Rietveld Asymmetry" << endl;
    if (IASYM == 1) file6 << "IASYM=1, Asymmetry by Riello et al.:',' Powder Diffraction,10,204-206,1995" << endl;
    //-----FONDO=0 :BKG EVALUATED USING STANDARD METHODS
    //-----FONDO=1 :BKG EVALUATED USING THE ISOTROPIC THERMAL FACTORS
    //-----FONDO=2 :BKG EVALUATED USING THE OVERAL THERMAL FACTORS
    //-----AIR SCATTERING,IAXX= 1 AIR SCATTERING ADDED TO THE BACKGROUND if SCAIR<>0
    //                    IAXX= -1 AIR SCATT. EVALUATED AND SUBTRACTED FROM DATA
    //                    IAXX= 0 AIR SCATT NOT EVALUATED and ibgd = 1
    //-----LINEAR ABSORP.CORR., IAS= 1  DATA CORRECTED FOR LINEAR ABSORPTION
    //
    //-----FI = AIR FRACTION INSIDE THE SAMPLE
    if(IBCKCODE == -1)
    {
        // line 2.2
        getline(file5,s);
        stringstream(s.substr( 0,4)) >> IAS;
        stringstream(s.substr( 4,4)) >> FONDO;
    }
    if (IBGD == 1)
    {
        file6 << "NO AMOPHOUS AND COMPTON CORRECTION TO THE BGD" << endl;
    }
    else
    {
        if(IAS == 0) file6 << "IAS = 0, NO LINEAR ABSORPTION CORRECTION" << endl;
        if(IAS == 1) file6 << "IAS = 1, LINEAR ABSORPTION CORRECTION IS APPLIED" << endl;
        if(FONDO == 1) file6 << "FONDO = 1,ISOTROPIC B FACTOR USED FOR BKG EVALUATION" << endl;
        if(FONDO == 2) file6 << "FONDO = 2,OVERALL Q USED FOR BKG EVALUATION" << endl;
    }


    // line 3
    getline (file5,s);
    stringstream(s.substr( 0,1)) >> IOT;
    stringstream(s.substr( 1,1)) >> IPL;
    stringstream(s.substr( 2,1)) >> IPC;
    stringstream(s.substr( 3,1)) >> MAT;
    stringstream(s.substr( 4,1)) >> NXT;

    stringstream(s.substr( 6,1)) >> LST1;
    stringstream(s.substr( 7,1)) >> LST2;
    stringstream(s.substr( 8,1)) >> LST3;
    stringstream(s.substr( 9,1)) >> IPL1;
    stringstream(s.substr(10,1)) >> IPL2;

    stringstream(s.substr(12,1)) >> IPLST;
    stringstream(s.substr(13,1)) >> IPLOSS;
    stringstream(s.substr(14,1)) >> IPLCAL;
    stringstream(s.substr(15,1)) >> IPLPOL;
    stringstream(s.substr(16,1)) >> IPLCOM;

    stringstream(s.substr(18,1)) >> IPLDIS;
    stringstream(s.substr(19,1)) >> IPLAM;
    stringstream(s.substr(20,1)) >> IPBIG;
    if(IBCKCODE != -1)
    {
        IPLOSS=0;
        IPLCAL=0;
        IPLPOL=0;
        IPLCOM=0;
        IPLDIS=0;
        IPLAM =0;
        IPBIG =0;
    }
    ISIMOP_=0;
    if (IOT != 0) file6 << "OUTPUT OBSERVED AND CALCULATED INTENSITIES ON LAST CYCLE" << endl;
    if (IPL != 0) file6 << "GENERATE LINE PRINTER PLOT" << endl;
    if (IPL2 != 0) file6 << "Generate file with obs, calc, dif intensiies, weighted diference, and Bragg peak positions" << endl;
    if (IPLST == 1) file6 << "GENERATE PARAMETER LIST" << endl;
    if (IPLST == 2) file6 << "GENERATE PARAMETER LIST" << endl << "GENERATE PARAMETERS AND STD.DEV. IN FINAL CYCLE FOR DATA BASE" << endl;
    if (IPC == 1) file6 << "OUTPUT INTENSITIES" << endl;
    if (IPC == 2) file6 << "OUTPUT ABSOLUTE VALUES OF STRUCTURE FACTORS + PHASE ANGLE" << endl;
    if (IPC == 3) file6 << "OUTPUT A AND B (CALC + OBS) STRUCTURE FACTORS" << endl;
    if (MAT != 0 && JOBTYP < 3) file6 << "OUTPUT CORRELATION MATRIX" << endl;
    if (NXT != 0 && JOBTYP < 3) file6 << "GENERATE NEW INPUT FILE" << endl;
    if (LST1 != 0) file6 << "PRINT REFLECTION LIST" << endl;
    if (NPHASE == 1)LST3=0;
    if (LST3 != 0) file6 << "Print merged reflection list" << endl;
    if (LST2 != 0) file6 << "Print corrected data" << endl;

    // line 4
    getline (file5,s);
    stringstream(s.substr( 0*8, 8)) >> diffractogram.LAMDA[1];
    stringstream(s.substr( 1*8, 8)) >> diffractogram.LAMDA[2];
    stringstream(s.substr( 2*8, 8)) >> diffractogram.RATIO[2];
    stringstream(s.substr( 3*8, 8)) >> BKPOS_;
    stringstream(s.substr( 4*8, 8)) >> WDT_;
    stringstream(s.substr( 5*8, 8)) >> CTHM_;
    stringstream(s.substr( 6*8, 8)) >> TMV_;
    stringstream(s.substr( 7*8, 8)) >> RLIM_;
    stringstream(s.substr( 8*8, 8)) >> SW;

    diffractogram.RATIO[1] = 1.0;
    diffractogram.CalcLamdaM();

    file6
        << " WAVELENGTHS= "
        << fixed
        << setw(9) << setprecision(6) << diffractogram.LAMDA[1]
        << setw(9) << setprecision(6) << diffractogram.LAMDA[2]
        << " LAMDA MEAN = "
        << setw(9) << setprecision(6) << diffractogram.LAMDAM <<endl;
    file6 << "ALPHA2:ALPHA1 RATIO = " << setw(9) << setprecision(5) << diffractogram.RATIO[2] << endl;

    IXRAY = diffractogram.SearchWavelength();

    file6 << "BASE OF PEAK = 2.0*HW*" << setw(8) << setprecision(2) << WDT_ << endl;
    file6 << "MONOCHROMATOR CORRECTION =" << setw(8) << setprecision(4) << CTHM_ << endl;
    file6 << "ABSORPTION CORRECTION COEFFICIENT = " << setw(8) << setprecision(4) << TMV_ << " CM-1" << endl
        << "SLAB-WIDTH = " << setw(8) << setprecision(4) << SW << " CM." << endl;
    if (IASYM == 0)
    {
        file6 << "RIETVELD ASYMMETRY CORRECTION FOR ANGLES LESS THAN " << setw(8) << setprecision(3) << RLIM_ << " DEGREES" << endl;
    }
    else
    {
        file6 << "ASYMMETRY CORRECTION FOR ANGLES LESS THAN "
            << setw(8) << setprecision(3) << 90.0 - RLIM_ << " DEGREES" << endl
            << "                         AND GREATER THAN "
            << setw(8) << setprecision(3) << 90.0 + RLIM_ << " DEGREES" << endl;
    }

    // line 5
    getline(file5,s);
    stringstream(s.substr(0,4)) >> MCYCLE;
    stringstream(s.substr(4,4)) >> EPS;
    stringstream(s.substr(8,4)) >> RELAX[1];
    stringstream(s.substr(12,4)) >> RELAX[2];
    stringstream(s.substr(16,4)) >> RELAX[3];
    stringstream(s.substr(20,4)) >> RELAX[4];

    if(JOBTYP > 2)
    {
        MCYCLE=1;
        stringstream(s.substr(24,8)) >> diffractogram.THMIN;
        stringstream(s.substr(32,8)) >> diffractogram.STEP;
        stringstream(s.substr(40,8)) >> diffractogram.THMAX;
    }
    ICYRUN = MCYCLE;
    file6 << "NUMBER OF CYCLES = " << setw(4) << MCYCLE << endl;
    file6 << "RELAXATION FACTORS" << endl
        << "FOR COORDINATES= " << setw(5) << setprecision(2) << RELAX[1] << endl
        << "FOR ANISOTRPIC TEMPERATURE FACTORS= " << setw(5) << setprecision(2) << RELAX[2] << endl
        << "FOR FWHM PARAMETERS= " << setw(5) << setprecision(2) << RELAX[3] << endl
        << "FOR LATTICE CONSTANTS= " << setw(5) << setprecision(2) << RELAX[4] << endl;
    file6 << "EPS-VALUE= " << setw(6) << setprecision(1) << EPS << endl;

    if(diffractogram.NBCKGD >= 2)
    {
        // line 6(*)
        for (I=1; I <= diffractogram.NBCKGD; ++I)
        {
            getline(file5,s);
            stringstream(s.substr(0,8)) >> diffractogram.POS[I];
            stringstream(s.substr(8,8)) >> diffractogram.BCK[I];
        }
        file6 << "BACKGROUND" << endl
            << "POSITION    INTENSITY" << endl;
        for (I=1; I <= diffractogram.NBCKGD; ++I)
        {
            file6 << setw(9) << setprecision(4) << diffractogram.POS[I]
            << setw(9) << setprecision(4) << diffractogram.BCK[I] << endl;
        }
    }

    if(diffractogram.NEXCRG > 0)
    {
        // line 7(*)
        for (I=1; I <= diffractogram.NEXCRG; ++I)
        {
            getline(file5,s);
            stringstream(s.substr(0,8)) >> diffractogram.ALOW[I];
            stringstream(s.substr(8,8)) >> diffractogram.AHIGH[I];
        }
        file6 << "EXCLUDED REGIONS" << endl
            << "FROM     TO" << endl;
        for (I=1; I <= diffractogram.NEXCRG; ++I)
        {
            file6 << setw(9) << setprecision(4) << diffractogram.ALOW[I]
            << setw(9) << setprecision(4) << diffractogram.AHIGH[I] << endl;
        }
    }

    if(diffractogram.NSCAT > 0)
    {
        for (I=1; I <= diffractogram.NSCAT; ++I)
        {
            if(JOBTYP == 2 || JOBTYP == 4)
            {
                // line 8.1 XRD(*)
                getline(file5,s);
                NAM_[I] = s.substr(0,4);
                stringstream(s.substr(4,8)) >> coeff.DFP[I];
                stringstream(s.substr(12,8)) >> coeff.XMAS[I];
            }
            else
            {
                // line 8.1 ND(*)
                getline(file5,s);
                NAM_[I] = s.substr(0,4);
                stringstream(s.substr(4,8)) >> coeff.DFP[I];
                stringstream(s.substr(12,8)) >> coeff.DFPP[I];
                stringstream(s.substr(20,8)) >> coeff.XMAS[I];
                K=0;
                // line 8.2 XRD(*)
                L126:
                getline(file5,s);
                for (J=1; J <= 9; ++J) stringstream(s.substr((J-1)*8,8)) >> coeff.AC[J][I];
                if (coeff.AC[1][I] == -100.0) coeff.COEF(&I,&K);
                if(coeff.AC[3][I] == 0.0)
                {
                    K=K+1;
                    coeff.POSI[K]=coeff.AC[1][I];
                    coeff.SCAT[K]=coeff.AC[2][I];
                    if (K <= 29) goto L126;
                    file6 << "TOO MANY SCATTERING TABLE ENTRIES" << endl;
                    DBWSException("7700");
                }
            }


        }
        if(JOBTYP == 1 || JOBTYP == 3)
        {
            file6 << "SCATTERING LENGTHS" << endl;
            for (I=1; I <= diffractogram.NSCAT; ++I)
            {
                file6 << "FOR " << setw(4) << NAM_[I] << "     "  << setw(10) << setprecision(6) << coeff.DFP[I] << endl;
            }
            for (I=1; I <= diffractogram.NSCAT; ++I)
            {
                for (J=1; J <= 9; ++J) coeff.AC[J][I]=0.0;
                coeff.DFPP[I]=0.0;
            }
        }
        else
        {
            file6 << "FORMFACTORS" << endl;
            for (I=1; I <= diffractogram.NSCAT; ++I)
            {
                file6 << "FOR " << setw(4) << NAM_[I]
                << " DFP=" << setw(10) << setprecision(6) << coeff.DFP[I]
                << " DFPP=" << setw(10) << setprecision(6) << coeff.DFPP[I] << endl
                    << "COEFFICIENTS= ";
                for (J=1; J <= 9; ++J) file6 << setw(10) << setprecision(6) << coeff.AC[J][I];
                file6 << endl;
            }
        }
    }


    // line 9
    getline(file5,s);
    stringstream(s.substr(0,8)) >> MAXS;
    if(MAXS > MSZ)
    {
        file6 << "* YOU HAVE DECLARED MORE CODEWORDS THAN WILL FIT INTO *" << endl
            << "* THE -MSZ- ARRAY.  EITHER DECREASE THE # OF CODEWORD *" << endl
            << "* OR  INCREASE  THE  -MSZ-  ARRAY  SIZE AND RECOMPILE *" << endl;
        cout  << "* YOU HAVE DECLARED MORE CODEWORDS THAN WILL FIT INTO *" << endl
            << "* THE -MSZ- ARRAY.  EITHER DECREASE THE # OF CODEWORD *" << endl
            << "* OR  INCREASE  THE  -MSZ-  ARRAY  SIZE AND RECOMPILE *" << endl;
        DBWSException("MAXS > MSZ");
    }
    cout << "INPUT:    CYCLES =" << setw(4) << MCYCLE << "     REFINABLE PARAMETERS =" << setw(4) <<MAXS << endl;

    //     CHECK DIMENSIONING FOR SOME EQUIVALENCED ARRAYS
    if (IDSZ < MSZ*MAXS)
    {
        file6 << "CHANGE IDSZ OR MSZ SO THAT IDSZ IS EQUAL TO OR GREATER THAN MSZ*MAXS" << endl
            << "IDSZ, MAXIMUM NO. OF DATA POINTS = " << setw(6) << IDSZ << endl
            << "MSZ,  MATRIX SIZE                = " << setw(6) << MSZ << endl
            << "MAXS, NO. OF PARAMETERS VARIED   = " << setw(6) << endl
            << "*** JUST A DIMENSIONING ERROR ***" << endl;
        DBWSException("IDSZ IS LESS THAN MSZ*MAXS");
    }
    if (JOBTYP > 2) MAXS=0;
    file6 << "NUMBER OF PARAMETERS VARIED= " << setw(5) << MAXS << endl;

    // line 10.1
    getline(file5,s);
    GLB_[1-1] = s.substr(0,8);
    GLB_[10-1] = s.substr(1*8,8);
    GLB_[11-1] = s.substr(2*8,8);
    GLB_[8-1] = s.substr(3*8,8);
    GLB_[9-1] = s.substr(4*8,8);
    GLB_[12-1] = s.substr(5*8,8);
    GLB_[13-1] = s.substr(6*8,8);

    // line 10.11
    getline(file5,s);
    stringstream(s.substr(0,8))   >> GLB_[1-1].codeword;
    stringstream(s.substr(1*8,8)) >> GLB_[10-1].codeword;
    stringstream(s.substr(2*8,8)) >> GLB_[11-1].codeword;
    stringstream(s.substr(3*8,8)) >> GLB_[8-1].codeword;
    stringstream(s.substr(4*8,8)) >> GLB_[9-1].codeword;
    stringstream(s.substr(5*8,8)) >> GLB_[12-1].codeword;
    stringstream(s.substr(6*8,8)) >> GLB_[13-1].codeword;

    file6 << "GLOBAL PARAMETERS AND CODEWORDS" << endl
        << "ZEROPOINT= "
        << setw(8) << setprecision(2) << GLB_[1-1]
    << setw(8) << setprecision(2) << GLB_[1-1].codeword << endl;

    //-----PMON1,PMON2=PARAMETER OF THE MONOCHROMATOR
    //-----IT READS PARAMETER OF THE MONOCHROMATOR AND AIR SCALE
    //-----if THE MONOCHROMATOR WORKS ON THE INCIDENT BEAM PUT :
    //-----PMON1=1;PMON2=0;
    // next 2 'READ' are for amorphous bkg codes
    // SCAIR=glb(17), FLSCAIR=aglb(17), SCAM =GLB_[20-1], FLSCAM=GLB_[20-1].A
    // PMON1=GLB_[18-1], FLMON1 =GLB_[18-1].A, PMON2=GLB_[19-1], FLMON2=GLB_[19-1].A
    //      READ(5,455,END=99999)SCAIR,FLSCAIR,SCAM,FLSCAM,
    //     *  PMON1,FLMON1,PMON2,FLMON2
    //
    //      READ(5,455,END=99999)glb(17),aglb(17),GLB_[20-1], GLB_[20-1].A,
    //     *  GLB_[18-1],GLB_[18-1].A,GLB_[19-1], GLB_[19-1].A
    //455     FORMAT(BZ,8F8.0)
    //  !cp jun 95 start
    if (IBGD != 1)
    {
        // line 10.2 and line 10.21
        getline(file5,s);
        GLB_[20-1] = s.substr(0*8,8);
        GLB_[18-1] = s.substr(1*8,8);
        GLB_[19-1] = s.substr(2*8,8);

        getline(file5,s);
        stringstream(s.substr(0*8,8)) >> GLB_[20-1].codeword;
        stringstream(s.substr(1*8,8)) >> GLB_[18-1].codeword;
        stringstream(s.substr(2*8,8)) >> GLB_[19-1].codeword;

        file6 << "AMORPHOUS SCALE and CODEWORD= "
            << setw(14) << setprecision(4) << GLB_[20-1]
        << setw(14) << setprecision(4) << GLB_[20-1].codeword << endl;
        file6 << "MONOCROMATOR BANDPASS PARAMETERS AND CODEWORDS" << endl
            << "PARAMETERS MONOC="
            << setw(8) << setprecision(4) << GLB_[18-1]
        << setw(8) << setprecision(4) << GLB_[18-1].codeword
        << "                  "
            << setw(8) << setprecision(4) << GLB_[19-1]
        << setw(8) << setprecision(4) << GLB_[19-1].codeword << endl;
    }


    if(diffractogram.NBCKGD == 0)
    {
        // line 10.4(*)
        // line 10.3 and line 10.31
        getline(file5,s);
        GLB_[2-1] = s.substr(0*9,9);
        GLB_[3-1] = s.substr(1*9,9);
        GLB_[4-1] = s.substr(2*9,9);
        GLB_[5-1] = s.substr(3*9,9);
        GLB_[6-1] = s.substr(4*9,9);
        GLB_[7-1] = s.substr(5*9,9);

        getline(file5,s);
        stringstream(s.substr(0*9,9)) >> GLB_[2-1].codeword;
        stringstream(s.substr(1*9,9)) >> GLB_[3-1].codeword;
        stringstream(s.substr(2*9,9)) >> GLB_[4-1].codeword;
        stringstream(s.substr(3*9,9)) >> GLB_[5-1].codeword;
        stringstream(s.substr(4*9,9)) >> GLB_[6-1].codeword;
        stringstream(s.substr(5*9,9)) >> GLB_[7-1].codeword;


        file6 << "BACKGROUND PARAMETERS AND CODEWORDS" << endl
            << "ORIGIN OF BACKGROUND POLYNOMIAL AT TWO-THETA = "
            << setw(8) << setprecision(3) << BKPOS_ << "DEGREES" << endl
            << setw(12) << setprecision(4) << GLB_[2-1]
            << setw(12) << setprecision(4) << GLB_[3-1]
            << setw(12) << setprecision(4) << GLB_[4-1]
            << setw(12) << setprecision(4) << GLB_[5-1]
            << setw(12) << setprecision(4) << GLB_[6-1]
            << setw(12) << setprecision(4) << GLB_[7-1] << endl
            << setw(12) << setprecision(3) << GLB_[2-1].codeword << "    "
            << setw(12) << setprecision(3) << GLB_[3-1].codeword << "    "
            << setw(12) << setprecision(3) << GLB_[4-1].codeword << "    "
            << setw(12) << setprecision(3) << GLB_[5-1].codeword << "    "
            << setw(12) << setprecision(3) << GLB_[6-1].codeword << "    "
            << setw(12) << setprecision(3) << GLB_[7-1].codeword << "    " << endl;
    }

    file6 << "DISPLACEMENT PEAKSHIFT PARAMETER AND CODEWORD"
        << setw(8) << setprecision(2) << GLB_[10-1]
    << setw(8) << setprecision(2) << GLB_[10-1].codeword << endl
        << "TRANSPARENCY PEAKSHIFT PARAMETER AND CODEWORD"
        << setw(8) << setprecision(2) << GLB_[11-1]
    << setw(8) << setprecision(2) << GLB_[11-1].codeword << endl;
    file6 << "SURFACE ROUGHNESS P PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << GLB_[8-1]
    << setw(9) << setprecision(4) << GLB_[8-1].codeword << endl
        << "SURFACE ROUGHNESS Q PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << GLB_[9-1]
    << setw(9) << setprecision(4) << GLB_[9-1].codeword << endl
        << "SURFACE ROUGHNESS R PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << GLB_[12-1]
    << setw(9) << setprecision(4) << GLB_[12-1].codeword << endl
        << "SURFACE ROUGHNESS T PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << GLB_[13-1]
    << setw(9) << setprecision(4) << GLB_[13-1].codeword << endl;


    if(JOBTYP <= 2)
    {

        if(JOBTYP == 1 && INSTRM == 2)           //-----if PATTERN CALCULATION ONLY FOR SYNCHROTRON X-RAY DATA
        {
            diffractogram.ReadSynchrotronData();
        }
        else if(JOBTYP == 2 && INSTRM == 2)           //-----if PATTERN CALCULATION ONLY FOR MULTIPLE NEUTRON DATA
        {
            diffractogram.ReadNeutronData();
        }
        else
        {
            switch (IDATA) {
            case 0:
                diffractogram.ReadDBWS();          // read DBWS formated
                break;
            case 1:
                diffractogram.ReadDBWSF();         // read DBWS formated
                break;
            case 2:
                diffractogram.GSASREAD();          // read GSAS data file (read start,stop,step and data)
                break;
            case 3:
                diffractogram.PHILIPSREAD();       // read Philips data file (read start,stop,step and data)
                break;
            case 4:
                diffractogram.scintag();           // read SCINTAG data file (read start,stop,step and data)
                break;
            case 5:
                diffractogram.SIEMENSREAD();       // read  SIEMENS UXD data file (read start,stepsize,stepcount and data)
                break;
            case 6:
                diffractogram.rigakuread();        // read  RIGAKU data file (read start,stepsize,stop,stepcount and data)
                break;
            }

            for (I=1; I <= diffractogram.NPTS; ++I)
            {
                if(diffractogram.Y[I] <= 1.0E-6) diffractogram.Y[I]=1.0;
                if(diffractogram.Y[I] > diffractogram.XMAXINT) diffractogram.XMAXINT=diffractogram.Y[I];
            }
            for (I=1; I <= diffractogram.NPTS; ++I)
            {
                //     BUILD UP AMORPHOUS-VECTOR
                //     if REQUIRED MAKE ABSORPTION CORRECTION
                //-----COMPUTE TWO-THETA
                TH = diffractogram.THMIN + static_cast<double>(I-1) * diffractogram.STEP;
                diffractogram.VAR[I]=diffractogram.Y[I];
                //-----COMPUTE SAMPLE ABSORPTION CORRECTION
                if(IAS == 1)
                {
                    ABSORP(TMV_,SW,TH,&ABC);
                    diffractogram.Y[I] = diffractogram.Y[I] / ABC;
                    diffractogram.VAR[I]=diffractogram.VAR[I]/ABC;
                }
            }
        }
    }
    else
    {
        diffractogram.NPTS = static_cast<int>((diffractogram.THMAX-diffractogram.THMIN)/diffractogram.STEP+1.5);
        if (diffractogram.NPTS > IDSZ)
        {
            file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
                << setw(5) << diffractogram.NPTS+IPTS << " POINTS WERE INPUT" << endl
                << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
            DBWSException("TOO MANY DATA POINTS");
        }
    }

    if (diffractogram.NBCKGD-1 < 0)
    {
        TH=diffractogram.THMIN-diffractogram.STEP;
        for (I=1; I <= diffractogram.NPTS; ++I)
        {
            TH=TH+diffractogram.STEP;
            THX=TH/BKPOS_-1.0;
            diffractogram.BK[I]=GLB_[2-1];
            for (J=2; J <= 6; ++J) diffractogram.BK[I]=diffractogram.BK[I] + GLB_[J+1-1] * pow(THX,J-1);
        }
    }
    else if (diffractogram.NBCKGD-1 == 0)
    {
        diffractogram.ReadFileBackground();
    }
    else
    {
        diffractogram.InitPolynomialBackground();
    }


    for (K=1; K <= NPHASE; ++K) phases[K].ICR=0;


    // start loop on phase
    for (K=1; K <= NPHASE; ++K)
    {
        // line 11.1
        getline(file5,s);
        phases[K].name = s.substr(0,50);
        file6 << "PHASE " << setw(2) << K << endl << phases[K].name << endl;

        // line 11.2
        getline(file5,s);
        phases[K].AtomCount = s.substr(0,4);
        phases[K].NMOL = s.substr(4,4);
        phases[K].SAQF = s.substr(8,7);
        stringstream(s.substr(16,4)) >> phases[K].PREF[1];
        stringstream(s.substr(20,4)) >> phases[K].PREF[2];
        stringstream(s.substr(24,4)) >> phases[K].PREF[3];
        stringstream(s.substr(28,7)) >> phases[K].WTIS;
        N=phases[K].AtomCount;
        ISTEST = ISPHASE;
        if (K == ISPHASE)
        {
            if (phases[K].WTIS == 0.0)
            {
                ISPHASE = 0;
                cout << "WARNING: Internal Standard WT% = ZERO. Check ISWT in line 11.2 for phase "
                    << setw(2) << ISTEST << "." << endl
                    << "ISPHASE (line 2) turned to ZERO. Amourphous content will not be calculated." << endl;
                file6 << "               WARNING: Wt% for internal Standard is ZERO." << endl
                    << "               Check ISWT in line 11.2 for this phase." << endl
                    << "               ISPHASE (line 2) turned to ZERO." << endl
                    << "               AMORPHOUS CONTENT WILL NOT BE CALCULATED." << endl;
            }
        }

        file6 << "NUMBER OF ATOMS= " << setw(4) << N << endl
              << "NUMBER OF FORMULA UNITS PER UNIT CELL= " << setw(4) << phases[K].NMOL << endl
              << "PARTICLE ABSORPTION FACTOR = " << setw(8) << setprecision(4) << phases[K].SAQF << endl
              << "PREFERRED ORIENTATION VECTOR= " << setw(8) << setprecision(4) << phases[K].PREF[1] << setw(8) << setprecision(4) << phases[K].PREF[2] << setw(8) << setprecision(4) << phases[K].PREF[3] << endl;
        if (K == ISPHASE)
        {
            file6 << "MASS% IN THE SAMPLE= " << setw(7) << setprecision(2) << phases[K].WTIS << endl;
        }

        // line 11.3
        getline(file5,s);
        phases[K].SYMB = s.substr(0,20);

        SPG = " " + phases[K].SYMB.get();
        for (I=1; I <= 20; ++I)
        {
            //C convert lower case to upper case
            for (IK=1; IK <= 26; ++IK) if (SPG[I] == LOWER[IK]) SPG[I]=UPPER[IK];
            // finish conversion
        }
        phases[K].SYMB.SPGP(SPG);

        // getting multiplicity of each phase !cp jun 96)
        ISIMOP_=1;
        RTMT(&MULTX_,Y_, &IPL1,NCTR_,&ISIMOP_,&K);
        phases[K].XMLTP=MLTPHASE;
        file6 << "The multiplicity of the general site is " << setw(3) << MLTPHASE << endl;
        //-----READ AND PRINT FOR EACH ATOM
        IOF=0;
        if(K > 1)
        {
            for (IIPHAS=2; IIPHAS <= K; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
        }
        file6 << "***INITIAL PARAMETERS***" << endl
            << "ATOM    M NTYP              X         Y         Z         B        So" << endl
            << "                          B11       B22       B33       B12       B13       B23" << endl;

        // line 11-4i
        // READ and WRITE, and respectives FORMAT command lines below were
        // changed to incorporate the parameter MURT(I+IOF) * !cp jun 96
        for (I=1; I <= N; ++I)
        {
            getline(file5,s);
            ATEXT[I+IOF] = s.substr(0,4);
            stringstream(s.substr(5,4)) >> MURT[I+IOF];
            NTYP[I+IOF] = s.substr(10,4);
            XL_[I+IOF][1] = s.substr(16,8);
            XL_[I+IOF][2] = s.substr(24,8);
            XL_[I+IOF][3] = s.substr(32,8);
            XL_[I+IOF][4] = s.substr(40,8);
            XL_[I+IOF][5] = s.substr(48,8);

            getline(file5,s);
            stringstream(s.substr(16,8)) >> XL_[I+IOF][1].codeword;
            stringstream(s.substr(24,8)) >> XL_[I+IOF][2].codeword;
            stringstream(s.substr(32,8)) >> XL_[I+IOF][3].codeword;
            stringstream(s.substr(40,8)) >> XL_[I+IOF][4].codeword;
            stringstream(s.substr(48,8)) >> XL_[I+IOF][5].codeword;

            getline(file5,s);
            XL_[I+IOF][6] = s.substr(0,8);
            XL_[I+IOF][7] = s.substr(8,8);
            XL_[I+IOF][8] = s.substr(16,8);
            XL_[I+IOF][9] = s.substr(24,8);
            XL_[I+IOF][10] = s.substr(32,8);
            XL_[I+IOF][11] = s.substr(40,8);

            getline(file5,s);
            stringstream(s.substr(0,8))  >> XL_[I+IOF][6].codeword;
            stringstream(s.substr(8,8))  >> XL_[I+IOF][7].codeword;
            stringstream(s.substr(16,8)) >> XL_[I+IOF][8].codeword;
            stringstream(s.substr(24,8)) >> XL_[I+IOF][9].codeword;
            stringstream(s.substr(32,8)) >> XL_[I+IOF][10].codeword;
            stringstream(s.substr(40,8)) >> XL_[I+IOF][11].codeword;
        }

        // convert lower case to upper case
        for (ITN = 1; ITN <= N; ++ITN)
        {
            for (ITIPO = 1; ITIPO <= 2; ++ITIPO)
            {
                for (IK=1; IK <= 26; ++IK) if (NTYP[ITN+IOF][ITIPO] == LOWER[IK]) NTYP[ITN+IOF][ITIPO]=UPPER[IK];
            }
        }        
        // finish conversion

        for (I=1; I <= N; ++I)
        {
            file6
                << setw(4) << ATEXT[I+IOF] << " " << setw(4) << MURT[I+IOF]
                << " " << NTYP[I+IOF] << "        "
                << setw(10) << setprecision(5) << XL_[I+IOF][1]
                << setw(10) << setprecision(5) << XL_[I+IOF][2]
                << setw(10) << setprecision(5) << XL_[I+IOF][3]
                << setw(10) << setprecision(5) << XL_[I+IOF][4]
                << setw(10) << setprecision(5) << XL_[I+IOF][5] << endl
                << "                      "
                << setw(10) << setprecision(5) << XL_[I+IOF][6]
                << setw(10) << setprecision(5) << XL_[I+IOF][7]
                << setw(10) << setprecision(5) << XL_[I+IOF][8]
                << setw(10) << setprecision(5) << XL_[I+IOF][9]
                << setw(10) << setprecision(5) << XL_[I+IOF][10]
                << setw(10) << setprecision(5) << XL_[I+IOF][11] << endl;
        }
        //  !cp jun 96 ... (CONVERT sof MULTIPLICITY)(also changed in OUTPTR)
        for (I=1; I <= N; ++I)
        {
            if(static_cast<int>(XL_[I+IOF][5].codeword/10) != 0 && XL_[I+IOF][5] == 0)
            {
                XL_[I+IOF][5]=1e-6;
            }
            XL_[I+IOF][5] = XL_[I+IOF][5] * MURT[I+IOF] / phases[K].XMLTP;
        }

        // PAR[K][21] introduced below. It is for the term cot**2 in the pv-5 FWHM !cp Aug 95
        // line 11-5, line 11-6, line 11-7, line 11-8 and line 11-9
        getline(file5,s);			// S  O_B (line 11-5)
        phases[K].PAR[0] = s.substr(0,8);
        phases[K].PAR[1] = s.substr(8,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K].PAR[0].codeword;
        stringstream(s.substr(8,8)) >> phases[K].PAR[1].codeword;

        getline(file5,s);			// FWHM (line 11-6)
        phases[K].PAR[2] = s.substr(0,8);
        phases[K].PAR[3] = s.substr(8,8);
        phases[K].PAR[4] = s.substr(16,8);
        phases[K].PAR[20] = s.substr(24,8);
        phases[K].PAR[19] = s.substr(32,8);
        phases[K].PAR[14] = s.substr(40,8);
        phases[K].PAR[15] = s.substr(48,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K].PAR[2].codeword;
        stringstream(s.substr(8,8)) >> phases[K].PAR[3].codeword;
        stringstream(s.substr(16,8)) >> phases[K].PAR[4].codeword;
        stringstream(s.substr(24,8)) >> phases[K].PAR[20].codeword;
        stringstream(s.substr(32,8)) >> phases[K].PAR[19].codeword;
        stringstream(s.substr(40,8)) >> phases[K].PAR[14].codeword;
        stringstream(s.substr(48,8)) >> phases[K].PAR[15].codeword;

        getline(file5,s);
        phases[K].PAR[5] = s.substr(0,8);
        phases[K].PAR[6] = s.substr(8,8);
        phases[K].PAR[7] = s.substr(16,8);
        phases[K].PAR[8] = s.substr(24,8);
        phases[K].PAR[9] = s.substr(32,8);
        phases[K].PAR[10] = s.substr(40,8);

        getline(file5,s);			// Unit cell (line 11-7)
        stringstream(s.substr(0,8)) >> phases[K].PAR[5].codeword;
        stringstream(s.substr(8,8)) >> phases[K].PAR[6].codeword;
        stringstream(s.substr(16,8)) >> phases[K].PAR[7].codeword;
        stringstream(s.substr(24,8)) >> phases[K].PAR[8].codeword;
        stringstream(s.substr(32,8)) >> phases[K].PAR[9].codeword;
        stringstream(s.substr(40,8)) >> phases[K].PAR[10].codeword;

        getline(file5,s);			// G1 G2 P (line 11-8)
        phases[K].PAR[11] = s.substr(0,8);
        phases[K].PAR[12] = s.substr(8,8);
        phases[K].PAR[13] = s.substr(16,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K].PAR[11].codeword;
        stringstream(s.substr(8,8)) >> phases[K].PAR[12].codeword;
        stringstream(s.substr(16,8)) >> phases[K].PAR[13].codeword;

        getline(file5,s);			// NA NB NC (line 11-91)
        phases[K].PAR[16] = s.substr(0,8);
        phases[K].PAR[17] = s.substr(8,8);
        phases[K].PAR[18] = s.substr(16,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K].PAR[16].codeword;
        stringstream(s.substr(8,8)) >> phases[K].PAR[17].codeword;
        stringstream(s.substr(16,8)) >> phases[K].PAR[18].codeword;

        getline(file5,s);			// NA NB NC HS (line 11-93)
        phases[K].PAR[23] = s.substr(0,8);
        phases[K].PAR[24] = s.substr(8,8);
        phases[K].PAR[25] = s.substr(16,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K].PAR[23].codeword;
        stringstream(s.substr(8,8)) >> phases[K].PAR[24].codeword;
        stringstream(s.substr(16,8)) >> phases[K].PAR[25].codeword;

        getline(file5,s);			// s-PVII (line 11-95);
        phases[K].PAR[26] = s.substr(0,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K].PAR[26].codeword;

        // checking for zeros if TCH-PV is being used
        if(NPROF == _TCHZ)
        {
            if(phases[K].PAR[14] == 0)phases[K].PAR[14]=1e-8;
            if(phases[K].PAR[15] == 0)phases[K].PAR[15]=1e-8;
            if(phases[K].PAR[2] == phases[K].PAR[3] && phases[K].PAR[3] == phases[K].PAR[4] && phases[K].PAR[4] == phases[K].PAR[19] && phases[K].PAR[2] == 0) phases[K].PAR[4]=1e-8;
        }
        // checking for zeros in PV #5
        if(NPROF == _pseudoVoigt)
        {
            if(phases[K].PAR[16] == 0)phases[K].PAR[16]=1e-6;
        }
        //      if(int(phases[K].PAR[19].APAR/10) != 0 && PAR[K][20] == 0)PAR[K][20]=1e-9
        // !cp Aug 95 introducing PAR[K][21]=ct
        // CHECKING FOR NON-REFINABLE PARAMETERS
        //CCC                        FOR  CT  IN TCHZ AND SPVII FUNCTIONS
        if(NPROF == _TCHZ || NPROF == _SplitPearsonVII)
        {
            phases[K].PAR[20]=0.0;
            if(phases[K].PAR[20].codeword != 0.0)
            {
                cout << "NON-REFINABLE PARAMETER TURNED ON (CT) WITH SPLIT PEARSON VII OR TCHZ PROFILE FUNCTION" << endl;
                file6 << "NON-REFINABLE PARAMETER TURNED ON (CT) WITH SPLIT PEARSON VII OR TCHZ PROFILE FUNCTION" << endl;
                DBWSException("");
            }
        }
        //CCC                             FOR  X,Y,Z IN NON-TCHZ FUNCTION
        if(NPROF != _TCHZ)
        {
            if(phases[K].PAR[19] != 0 || phases[K].PAR[15] != 0 || phases[K].PAR[14].codeword != 0)
            {
                file6 << "     NON-REFINABLE PARAMETER RESET TO ZERO (Z,X,Y) FOR" << endl
                    << "     NON TCHZ PROFILE FUNCTION" << endl;
                phases[K].PAR[19]=0.0;
                phases[K].PAR[15]=0.0;
                phases[K].PAR[14]=0.0;
            }
            if(phases[K].PAR[19].codeword != 0 || phases[K].PAR[15].codeword != 0 || phases[K].PAR[14].codeword != 0)
            {
                cout << "NON-REFINABLE PARAMETER TURNED ON (Z,X,Y) WITH NOT TCHZ PROFILE FUNCTION" << endl;
                file6 << "NON-REFINABLE PARAMETER TURNED ON (Z,X,Y) WITH NOT TCHZ PROFILE FUNCTION" << endl;
                DBWSException("");
            }
        }
        //CCCC                            FOR  RIET_ASYM,X,Y,Z,CT IN SPVII
        if (NPROF == _SplitPearsonVII)
        {
            for (KKS=1; KKS <= 3; ++KKS)
            {
                if(phases[K].PAR[13+KKS-1] != 0.0)
                {
                    phases[K].PAR[13+KKS-1]=0.0;
                    file6 << " RIET_ASYM X Y RESET TO ZERO FOR SPVII FUNCTION" << endl;
                }
                if(phases[K].PAR[13+KKS-1].codeword != 0.0)
                {
                    cout << "NON-REFINEABLE PARAMETER TURNED ON (X,Y,R/RCF_ASYM)" << endl
                        << ", WITH SPLIT PEARSON VII PROFILE FUNCTION" << endl;
                    file6 << "NON-REFINEABLE PARAMETER TURNED ON (X,Y,R/RCF_ASYM)" << endl
                        << ", WITH SPLIT PEARSON VII PROFILE FUNCTION" << endl;
                    DBWSException("");
                }
            }
            for (KKS=1; KKS <= 2; ++KKS)
            {
                if ( phases[K].PAR[19+KKS-1] != 0.0 )
                {
                    phases[K].PAR[19+KKS-1]=0.0;
                    file6 << "CT,Z  RESET TO ZERO FOR SPVII FUNCTION" << endl;
                }
                if ( phases[K].PAR[19+KKS-1].codeword != 0.0 )
                {
                    cout << "NON-REFINEABLE PARAMETER TURNED ON (CT,Z), WITH SPLIT PEARSON VII PROFILE FUNCTION" << endl;
                    file6 << "NON-REFINEABLE PARAMETER TURNED ON (CT,Z), WITH SPLIT PEARSON VII PROFILE FUNCTION" << endl;
                    DBWSException("");
                }
            }
        }
        //CCCCC             FOR HIGH SIDE PARAMETERS IN NON-PVII PROFILE FUNCTION
        if (NPROF != _SplitPearsonVII)
        {
            for (KK=1; KK <= 4; ++KK)
            {
                if ( phases[K].PAR[23+KK-1] != 0.0 )
                {
                    phases[K].PAR[23+KK-1]=0.0;
                    file6 << "NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY RESET TO ZERO FOR A, NON-SPLIT PEARSON VII FUNCTION" << endl;
                }
                if(phases[K].PAR[23+KK-1] != 0.0)
                {
                    cout << "NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY ONLY REFINABLE FOR A SPLIT PEARSON VII FUNCTION" << endl;
                    file6 << "NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY ONLY REFINABLE FOR A SPLIT PEARSON VII FUNCTION" << endl;
                }
            }
        }
        //CCCCC END OF CHECKING NON-REFINABLE PARAMETERS
        file6 << "OVERALL SCALE FACTOR=" << scientific << setw(12) << setprecision(6) << phases[K].PAR[0] << fixed << endl
            << "OVERALL TEMP. FACTOR=" << setw(12) << setprecision(5) << phases[K].PAR[1] << endl;

        file6 << "DIRECT CELL PARAMETERS:" << endl
        << "a = " << setw(9) << setprecision(4) << phases[K].PAR[5] << endl
        << "b = " << setw(9) << setprecision(4) << phases[K].PAR[6] << endl
        << "c = " << setw(9) << setprecision(4) << phases[K].PAR[7] << endl
        << "α = " << setw(9) << setprecision(4) << phases[K].PAR[8] << endl
        << "β = " << setw(9) << setprecision(4) << phases[K].PAR[9] << endl
        << "γ = " << setw(9) << setprecision(4) << phases[K].PAR[10] << endl;

        file6 << "PREFERRED ORIENTATION PARAMETERS="
            << setw(7) << setprecision(3) << phases[K].PAR[11]
        << setw(7) << setprecision(3) << phases[K].PAR[12] << endl
        << "ASYMMETRY PARAMETER="
            << setw(8) << setprecision(4) << phases[K].PAR[13] << endl;

        //-----CHECK FOR SPLIT PEARSON PROFILE
        if (NPROF == _SplitPearsonVII)
        {
            file6 << "LOW SIDE EXPONENT COEFFICIENTS="
                << setw(12) << setprecision(4) << phases[K].PAR[16]
            << setw(12) << setprecision(4) << phases[K].PAR[17]
            << setw(12) << setprecision(4) << phases[K].PAR[18]
            << "HIGH SIDE EXPONENT COEFFICIENTS="
                << setw(12) << setprecision(4) << phases[K].PAR[23]
            << setw(12) << setprecision(4) << phases[K].PAR[24]
            << setw(12) << setprecision(4) << phases[K].PAR[25] << endl
            << "SPLIT PEARSON VII ASYMMETRY PARAMETER="
                << setw(8) << setprecision(4) << phases[K].PAR[26] << endl;
        }
        else
        {
            //-----if NOT THE SPLIT PEARSON VII PROFILE
            file6 << "MIXING PARAMETERS = "
                << scientific
                << setw(10) << setprecision(3) << phases[K].PAR[16] << " "
                << setw(10) << setprecision(3) << phases[K].PAR[17] << " "
                << setw(10) << setprecision(3) << phases[K].PAR[18] << fixed << endl;
        }
        file6 << "FWHM PARAMETERS (U,V,W,CT,Z,X,Y)="
            << setw(9) << setprecision(4) << phases[K].PAR[2]
        << setw(9) << setprecision(4) << phases[K].PAR[3]
        << setw(9) << setprecision(4) << phases[K].PAR[4]
        << setw(9) << setprecision(4) << phases[K].PAR[20]
        << setw(9) << setprecision(4) << phases[K].PAR[19]
        << setw(9) << setprecision(4) << phases[K].PAR[14]
        << setw(9) << setprecision(4) << phases[K].PAR[15] << endl;
        cell_A=phases[K].PAR[5];
        cell_B=phases[K].PAR[6];
        cell_C=phases[K].PAR[7];
        cell_ALPHA=phases[K].PAR[8];
        cell_BETA=phases[K].PAR[9];
        cell_GAMMA=phases[K].PAR[10];
        CELL2(K);
        // ************************************** !cp ap 97 (from It code)
        if(FONDO == 1 || FONDO == 2)
        {
            phases[K].SCABKG = phases[K].GCOM * phases[K].PAR[0];
        }
        // ************************************************************
        file6 << "CELL VOLUME PHASE(" << setw(2) << K << " ) = " << setw(12) << setprecision(4) << phases[K].VOLI << endl;
        for (I=1; I <= 6; ++I) SAVE[K][I]=phases[K].PAR[I+5-1];

        // change to reci
        for (I=1; I <= 3; ++I) phases[K].PAR[I+5-1]=cell_AL[I][I];			// a,b and c cell parameters
        phases[K].PAR[8]=cell_AL[2][3];			// alpha cell
        phases[K].PAR[9]=cell_AL[1][3];			// beta cell
        phases[K].PAR[10]=cell_AL[1][2];			// gamma cell
        file6 << "***Coding of variables***" << endl
            << "ATOM                        X         Y         Z         B         So" << endl
            << "                          B11       B22       B33       B12       B13       B23" << endl;
        //-----PRINT CODEWORDS FOR ATOMIC PARAMETERS
        for (I=1; I <= N; ++I)
        {
            file6 << setw(4) << ATEXT[I+IOF]
            << "                 "
            << setw(10) << setprecision(2) << XL_[I+IOF][1].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][2].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][3].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][4].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][5].codeword << endl
                << "                     "
            << setw(10) << setprecision(2) << XL_[I+IOF][6].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][7].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][8].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][9].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][10].codeword
            << setw(10) << setprecision(2) << XL_[I+IOF][11].codeword << endl;
        }
        for (I=1; I <= N; ++I)
        {
            for (J=1; J <= 11; ++J)
            {
                X=XL_[I+IOF][J].codeword;
                IYY=static_cast<int>(abs(X)/10.0);
                if(IYY > MSZ)
                {
                    DBWSException("MATRIX SIZE IS TOO SMALL");
                }
                XL_[I+IOF][J].L=IYY;
                XL_[I+IOF][J].codeword=(abs(X)-10.*static_cast<double>(IYY))*sign(X);
            }
        }
        //-----PRINT CODEWORDS FOR PROFILE PARAMETERS
        file6 << "OVERALL SCALE FACTOR=" << setw(8) << setprecision(2) << phases[K].PAR[0].codeword << endl
            << "OVERALL TEMP. FACTOR=" << setw(8) << setprecision(2) << phases[K].PAR[1].codeword << endl
            << "DIRECT CELL PARAMETERS="
            << setw(8) << setprecision(2) << phases[K].PAR[5].codeword
        << setw(8) << setprecision(2) << phases[K].PAR[6].codeword
        << setw(8) << setprecision(2) << phases[K].PAR[7].codeword
        << setw(8) << setprecision(2) << phases[K].PAR[8].codeword
        << setw(8) << setprecision(2) << phases[K].PAR[9].codeword
        << setw(8) << setprecision(2) << phases[K].PAR[10].codeword << endl
            << "PREFERRED ORIENTATION PARAMETERS="
            << setw(8) << setprecision(2) << phases[K].PAR[11].codeword
        << setw(8) << setprecision(2) << phases[K].PAR[12].codeword << endl
            << "ASYMMETRY PARAMETER="
            << setw(8) << setprecision(4) << phases[K].PAR[13].codeword << endl;

        // !cp ap 97 (from It code)
        if(FONDO == 1 && (phases[K].PAR[1] != 0.0 || phases[K].PAR[1].codeword != 0.0))
        {
            DBWSException("WITH FONDO=1 YOU MUST USE ISOTROPIC THERMAL FACTORS");
        }
        if(FONDO == 2 && phases[K].PAR[1] == 0.0 && phases[K].PAR[1].codeword == 0.0)
        {
            DBWSException("if FONDO=2 YOU MUST USE THE OVERALL THERMAL FACTOR");
        }
        //-----CHECK FOR SPLIT PEARSON PROFILE
        if (NPROF == _SplitPearsonVII)
        {
            file6 << "LOW SIDE EXPONENT COEFFICIENTS="
                << setw(8) << setprecision(2) << phases[K].PAR[16].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[17].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[18].codeword << endl
                << "HIGH SIDE EXPONENT COEFFICIENTS="
                << setw(8) << setprecision(2) << phases[K].PAR[23].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[24].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[25].codeword << endl
                << "SPLIT PEARSON VII ASSYMETRY PARAMETER="
                << setw(8) << setprecision(2) << phases[K].PAR[26].codeword << endl;
        }
        else
        {
            //-----if NOT THE SPLIT PEARSON VII PROFILE
            file6
                << "MIXING PARAMETERS = "
                << setw(10) << setprecision(3) << phases[K].PAR[16].codeword
                << setw(10) << setprecision(3) << phases[K].PAR[17].codeword
                << setw(10) << setprecision(3) << phases[K].PAR[18].codeword << endl;
        }
        file6
            << "FWHM PARAMETERS (U,V,W,CT,Z,X,Y)="
            << setw(8) << setprecision(2) << phases[K].PAR[2].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[3].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[4].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[20].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[19].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[14].codeword
            << setw(8) << setprecision(2) << phases[K].PAR[15].codeword << endl;
        for (I=1; I <= 27; ++I)
        {
            X=phases[K].PAR[I-1].codeword;
            IYY=static_cast<int>(abs(X)/10.0);
            phases[K].PAR[I-1].L=IYY;
            phases[K].PAR[I-1].codeword=(abs(X)-10.*static_cast<double>(IYY))*sign(X);
        }
        //L151:
        LOOKUP(K,phases[K].AtomCount,diffractogram.NSCAT,IXRAY);
        if(FONDO == 1 || FONDO == 2) FINDC(K,diffractogram.NSCAT);
        U_=phases[K].PAR[2];
        V_=phases[K].PAR[3];
        W_=phases[K].PAR[4];
        ZZZ_ = phases[K].PAR[19];
        UC_ = phases[K].PAR[20];
        if (NPROF == _TCHZ)
        {
            ULOR_=phases[K].PAR[14];
            VLOR_=phases[K].PAR[15];
        }
        REFGEN(K,GLB_[1-1],GLB_[10-1],GLB_[11-1],phases[K].PAR[11],NCTR_);
        ISIMOP_=0;
        RTMT(&MULTX_,Y_,&IPL1,NCTR_,&ISIMOP_,&K);
        ICY=1;
        ICZ=0;
        for (IIPHAS=1; IIPHAS <= K; ++IIPHAS) ICZ = ICZ + phases[IIPHAS].ICR;
        if(K >= 2) ICY=1+ICZ-phases[K].ICR;
        if(LST1 != 1)goto L479;
        IXDEL=0;
        for (IXX=ICY; IXX <= ICZ; ++IXX)
        {
            IX=IXX-IXDEL;
            IRL=(refs[IXX].IREFS % 256)-128;
            IRK=((refs[IXX].IREFS/256) % 256)-128;
            IRH=((refs[IXX].IREFS/(256*256)) % 256)-128;
            IRC=(refs[IXX].IREFS/(256*256*256)) % 8;
            MLTT = static_cast<int>(refs[IXX].FMGNTD);
            if (NPROF == _TCHZ)
            {
                if((IX-1 % 60) == 0) file6 << "NO.  CODE    H   K   L  MULT   HW     POSN      FACTOR       HWL       HWG     ETA" << endl;
                file6 << setw(4) << IX
                    << setw(4) << IRC << "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(6) << MLTT
                    << setw(8) << setprecision(3) << refs[IXX].REFS[1]
                << setw(8) << setprecision(3) << refs[IXX].REFS[2]
                << setw(10) << setprecision(6) << refs[IXX].REFS[3]
                << setw(8) << setprecision(3) << refs[IXX].HALFL
                << setw(8) << setprecision(3) << refs[IXX].HALFG
                << setw(8) << setprecision(3) << refs[IXX].GAM << endl;
            }
            else if (NPROF == _SplitPearsonVII)
            {
                refs[IXX].FWHM[1]=2.0*(refs[IXX].REFS[1])*phases[K].PAR[26]/(1.0+phases[K].PAR[26]);
                refs[IXX].FWHM[2]=2.0*(refs[IXX].REFS[1])/(1.0+phases[K].PAR[26]);
                if((IX-1 % 60) == 0) file6 << "NO.  CODE    H   K   L  MULT      HWL    HWH     FWHM   POSN    FACTOR" << endl;
                file6 << setw(4) << IX
                    << setw(4) << IRC << "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(6) << MLTT
                    << setw(8) << setprecision(3) << refs[IXX].FWHM[1]
                << setw(8) << setprecision(3) << refs[IXX].FWHM[2]
                << setw(8) << setprecision(3) << refs[IXX].REFS[1]
                << setw(8) << setprecision(3) << refs[IXX].REFS[2]
                << setw(10) << setprecision(6) << refs[IXX].REFS[3] << endl;
            }
            else
            {
                if ( (IX-1 % 60) == 0 ) file6 << "NO.  CODE    H   K   L  MULT     HW     POSN    FACTOR" << endl;
                file6 << setw(4) << IX
                    << setw(4) << IRC << "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(6) << MLTT
                    << setw(8) << setprecision(3) << refs[IXX].REFS[1]
                << setw(8) << setprecision(3) << refs[IXX].REFS[2]
                << setw(10) << setprecision(6) << refs[IXX].REFS[3] << endl;
            }
        }
L479:
        for (IX=ICY; IX <= ICZ; ++IX)
        {
            refs[IX].REFS[3]=refs[IX].REFS[3]*refs[IX].FMGNTD; // double(MLTT(IX))
        }
    }
    // end of great loop on phases

    // FIND CODEWORDS FOR GLOBAL PARAMETERS: LGLB(I) & AGLB(I)
    for (I=1; I <= 20; ++I)
    {
        X=GLB_[I-1].codeword;
        IYY=static_cast<int>(abs(X)/10.0);
        GLB_[I-1].L=IYY;
        GLB_[I-1].codeword=(abs(X)-10.*static_cast<double>(IYY))*sign(X);
    }
    qpainit();
    if (JOBTYP <= 2)
    {
        NATOMS = 0;
        for (K=1; K <= NPHASE; ++K) NATOMS = NATOMS + phases[K].AtomCount;
        //     COUNT NO. OF USES OF SAME LOCATION IN NORMAL MATRIX
        for (I=1; I <= MAXS; ++I)
        {
            LCOUNT = 0;
            for (J=1; J <= 20; ++J) if ( I == GLB_[J-1].L ) LCOUNT = LCOUNT + 1;
            for (K=1; K <= NPHASE; ++K)
            {
                for (J=1; J <= 27; ++J) if ( I == phases[K].PAR[J-1].L )  LCOUNT = LCOUNT + 1;
            }
            for (K=1; K <= NATOMS; ++K)
            {
                for (J=1; J <= 11; ++J) if(I == XL_[K][J].L) LCOUNT = LCOUNT + 1;
            }
            if (LCOUNT >= 2) file6 << "******* CODEWORD " << setw(3) << I << " is used " << setw(3) << LCOUNT << " times **" << endl;
        }
        //     END COUNT NO. OF USES OF SAME LOCATION IN NORMAL MATRIX
        //     CHECK FOR HOLES IN THE NORMAL MATRIX
        ISTOP=0;
        for (I=1; I <= MAXS; ++I)
        {
            for (J=1; J <= 20; ++J) if (I == GLB_[J-1].L) goto L7200;
            //  HOLES IN PHASE PARAMETER
            for (K=1; K <= NPHASE; ++K)
            {
                for (J=1; J <= 27; ++J) if(I == phases[K].PAR[J-1].L) goto L7200;
            }
            // HOLES IN ATOMS PARAMETERS
            for (K=1; K <= NATOMS; ++K)
            {
                for (J=1; J <= 11; ++J) if(I == XL_[K][J].L) goto L7200;
            }
            file6 << "     ***** HOLE IN THE MATRIX ******. ELEMENT " << setw(3) << I << " IN THE NORMAL MATRIX IS MISSING" << endl;
            cout << "     ***** HOLE IN THE MATRIX ******. ELEMENT " << setw(3) << I << " IN THE NORMAL MATRIX IS MISSING" << endl;
            ISTOP = 1;
    L7200:;
        }
        if (ISTOP == 1) DBWSException("");
    }
}

//ccccc subroutine inpam: To read the amorphous data file !cp may 10 97
void DBWS::inpam(void)
{
    int I;
    double TH,SW1, ABC, TMV1, STEP1, THMIN1, THMAX1;
    string s,DATAID1;


    // unit 11 = amorphous file
    file11.open(file11name.data());

    //-----CONTROL if THERE IS THE AMORPHOUS FILE .
    //-----AND if IT IS ON THE SAME POINTS OF THE DATA FILE
    getline(file11,s);
    stringstream(s.substr(0,8)) >> THMIN1;
    stringstream(s.substr(8,8)) >> STEP1;
    stringstream(s.substr(16,8)) >> THMAX1;
    stringstream(s.substr(24,8)) >> TMV1;
    stringstream(s.substr(32,8)) >> SW1;
    DATAID1 = s.substr(40,16);
    file6 << "DATA AMORPHOUS " << DATAID1 << endl;
    if (THMIN1 != diffractogram.THMIN)
    {
        file6 << "AMORPHOUS THMIN="
            << setw(8) << setprecision(2) << THMIN1
            << "IS DIFFERENT FROM THE DATA THMIN="
            << setw(8) << setprecision(2) << diffractogram.THMIN << endl;
    }
    if (STEP1 != diffractogram.STEP)
    {
        file6 << "AMORPHOUS STEP="
            << setw(8) << setprecision(2) << STEP1
            << "IS DIFFERENT FROM THE DATA STEP="
            << setw(8) << setprecision(2) << diffractogram.STEP << endl;
    }
    if (THMAX1 != diffractogram.THMAX)
    {
        file6 << "AMORPHOUS THMAX="
            << setw(8) << setprecision(2) << THMAX1
            << "IS DIFFERENT FROM THE DATA THMAX="
            << setw(8) << setprecision(2) << diffractogram.THMAX << endl;
    }
    //      read the rest of the file in free format
    for (I=1; I <= diffractogram.NPTS; ++I) file11 >> diffractogram.AMORPHOUS[I];
    for (I = 1; I <= diffractogram.NPTS; ++I)
    {
        TH=diffractogram.THMIN+(I-1)*diffractogram.STEP;
        //------AMORPHOUS CORRECTION FOR AIR SCATTERING
        //       CALL ARIA(TMV1,SW1,FI1,TH,SCA)
        //       AMORPHOUS(I)=0.0
        if (IAS == 1)
        {
            ABSORP(TMV1,SW1,TH,&ABC);
            diffractogram.AMORPHOUS[I] = diffractogram.AMORPHOUS[I] / ABC;
        }
    }
    return;
    //L99999:
    //	file6 << "END OF FILE ENCOUNTERED AT START ON TAPE = 11 OF PHASE AMORPHUS" << endl;
    //	DBWSException("*** END OF FILE ENCOUNTERED FOR AMORPHOUS FILE ***");
}



int main(int argc, char *argv[])
{
    DBWS app;

    if (argc > 0)
    {
        // unit 4=data file
        app.diffractogram.file4name = string(argv[1]);
        app.diffractogram.file4.open(app.diffractogram.file4name.data());

        if (argc > 1)
        {
            // unit 5 Input control File
            app.file5name = string(argv[2]);
            app.file5.open(app.file5name.data());

            if (argc > 2)
            {
                // unit 6 = output file
                app.file6name = string(argv[3]);
                app.file6.open(app.file6name.data());


                // unit 8 = scratch
                app.file8name = "unit8";
                app.file8o.open(app.file8name.data(),ios::trunc);			//open (8,file='unit8',form='unformatted',status='unknown')

                app.run();



                if ( app.diffractogram.file4.is_open() ) app.diffractogram.file4.close();
                if ( app.file5.is_open() ) app.file5.close();
                if ( app.file6.is_open() ) app.file6.close();
                if ( app.file8o.is_open() ) app.file8o.close();    //close(8,status='DELETE');
                if ( app.file8i.is_open() ) app.file8i.close();
                if ( app.file11.is_open() ) app.file11.close();
                if ( app.file20.is_open() ) app.file20.close();

            }
            else
            {
                cout << "falta o nome do arquivo OUT" << endl;
            }
        }
        else
        {
            cout << "falta o nome do arquivo ICF" << endl;
        }
    }
    else
    {
        cout << "Falta o nome do arquivo RIT" << endl;
    }








    return 0;
}
