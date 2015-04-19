/*


RTMT   COMMON/BLNK1/N,NMAX,MULTX,Y(192,3),XLT,IXB,FILL(192)               blnk1_1*
OPERTR COMMON/BLNK1/N(2),L,X(192),Y(192),Z(192),XLT,IXB,FILL2(192)        blnk1_2*
SYMOPR COMMON/BLNK1/N(2),L,X(192),Y(192),Z(192),XLT,KXB,NCTR(192)         blnk1_3*
CEL000 COMMON/BLNK1/N(2),L,X(192),Y(192),Z(192),XLT,IXB,FILL(192)         blnk1_4*
OP1    COMMON/BLNK1/N(581),NCTR_(192)                                     blnk1_5*
INPTR  COMMON/BLNK1/NJNK(22),FILL2(751)                                   blnk1_6


RTMT   COMMON/BLNK1/N    ,NMAX ,MULTX,Y(192,1),Y(192,2),Y(192,3),XLT,IXB,FILL(192)
OPERTR COMMON/BLNK1/N1   ,N2   ,L    ,X(192)  ,Y(192)  ,Z(192)  ,XLT,IXB,FILL2(192)
SYMOPR COMMON/BLNK1/N1   ,N2   ,L    ,X(192)  ,Y(192)  ,Z(192)  ,XLT,KXB,NCTR(192)
CEL000 COMMON/BLNK1/N1   ,N2   ,L    ,X(192)  ,Y(192)  ,Z(192)  ,XLT,IXB,FILL(192)
OP1    COMMON/BLNK1/N1   ,N2   ,N3   ,N(192)  ,N(192)  ,N(192)  ,N  ,  N,NCTR_(192)
INPTR  COMMON/BLNK1/NJNK1,NJNK2,NJNK3,NJNK(22-3),FILL2(751)


RTMT   COMMON/BLNK1/            MULTX,Y(192,1),Y(192,2),Y(192,3)
OPERTR COMMON/BLNK1/           ,L    ,X(192)  ,Y(192)  ,Z(192)  ,    IXB,FILL2(192)
SYMOPR COMMON/BLNK1/           ,L    ,X(192)  ,Y(192)  ,Z(192)  ,XLT,KXB,NCTR(192)
CEL000 COMMON/BLNK1/           ,L    ,X(192)  ,Y(192)  ,Z(192)  ,
OP1    COMMON/BLNK1/                                                     NCTR_(192)


*/
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





class DC
{
public:
    double SAVE[99+1][6+1];
    int NSAVE;
};

class REFLS
{
public:
    int IREFS[IRS+1];
    double REFS[IRS+1][3+1];
    double FMGNTD[IRS+1];
    int ICR[99+1];
    double HALFL[IRS+1];
    double HALFG[IRS+1];
    double GAM[IRS+1];
    double FWHM[IRS+1][2+1];
};

class G1
{
public:
    double THMIN;
    double STEP;
    double THMAX;
    double U;
    double V;
    double W;
    double LAMDA[2+1];
    double TMV;
    double CTHM;
    double RLIM;
    double BKPOS; // sbx
    double WDT;
    double ULOR;
    double VLOR;
    double ZZZ;
    double UC;
};

class CNTRLS
{
public:
    int JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,IPLCOM,IPLDIS,IPLAM,IPBIG,MCYCLE;
    double EPS;
    int MAXS, INSTRM,MAXSX,MCYCLX,ICYRUN,IPREF,IABSR;
    int FONDO;
    int IAS,IASYM;
    double SW;
    int IBGD,IDATA,I2D94;

    int ISPHASE; // Internal Standard Phase
};

class PARAMS
{
public:
    Parameter GLB_[20];  // GLOBAL PARAMETERS

    double XL[NATS][11 +1];
    int LP[NATS][11 +1];
    double A[NATS][11 +1];    
    int PTR[NATS +1];

    double RELAX[4 +1];
    double RATIO[2 +1];
};

class PARAC
{
public:
    string ATEXT[NATS+1];
    string NTYP[NATS+1];
};

class JNK
{
public:
    double ALOW[100+1];
    double AHIGH[100+1];
    double POS[100+1];
    double BCK[100+1];
    int NBCKGD;
    int NEXCRG;
    int NSCAT; //nscatx
};

class ALLP
{
public:
    double FINAL[NFINAL+1][2+1];
    int ILOC;
};


class COEFC
{
public:
    string NAM[16+1];
};



class CNV
{
public:
    int N;
    double SINLAM[30+1];
    double F[30+1];

    void CV1(double AB[],double V[]);
    void CV2(double AB[],double* G);
    void STEEP(double X[], int N, int M);
};

class COEFF
{
public:
    double AC[10+1][16+1];
    double POSI[30+1];
    double SCAT[30+1];
    double DFP[16+1];
    double DFPP[16+1];
    double XMAS[16+1];
    CNV cnv;

    void COEF(int* J, int* K);
};

class BKGSCALE
{
public:
    double SCABKG[99+1];
};

class DATAX
{
public:
    double Y[IDSZ+1];
    double VAR[IDSZ+1];
    double YC[IDSZ+1];
    int KR[IDSZ+1];  // phs
    double BK[IDSZ+1];
    int NPTS;
    double AMORPHOUS[IDSZ+1];
};

class RTMTX
{
public:
    //int IVEC[99+1][192+1];
    //int MLTPHS[99+1];
    //int ICNTPHS[99+1];
};

class SIZESTRAIN
{
public:
    double SIZEG[15+1];
    double STRAING[15+1];
    double SIZEL[15+1];
    double STRAINL[15+1];
    double SIZ[15+1];
    double STRAIN[15+1];
    int NSIZESTRAIN;
};



class VOLUME
{
public:
    double VOLI[99+1];
    double GCOM[99+1];
};

class G2
{
public:
    double S1,S2,SS2,S3,SS4,D1,D2,D4,R1,R2,R3,R2NOBK,R3NOBK;
};

class F1
{
public:
    double RJAC[MSZ+1][MSZ+1]; // SMM
    double VX[MSZ+1]; // V1, V
};

class G3
{
public:
    double COND;
    int IORD1;
    int IORD2;
    double TH;
    int NUM;
};

class FONDI
{
public:
    double BKCOM[IDSZ+1];
    double BKDIS[IDSZ+1];
    double BKAM[IDSZ+1];
    double BKPOL[IDSZ+1];
};

class MULTIP
{
public:
    double TMASSA[99+1];
    int MLTPHASE;
    double XMLTP[99+1];
    int MURT[NATS+1];
    //double SAQF[99+1];
    //double WTIS[99+1];
};

class SIMOPER
{
public:
    int ISIMOP;
};

class HKLCTL
{
public:
    int IHKL[3+1][48+1];
    double AZ[48+1];
    int NC1[10+1][7+1][99+1];
    int ICHKL[99+1];
    int N1HKL[99+1];
    int IER;
};

class CELLX
{
public:
    double A; // AA
    double B;
    double C;
    double ALPHA;
    double BETA;
    double GAMMA;
    double AL[3+1][3+1];
};

class PRFX
{
public:
    double DELT;
    double TL;
    double GAM1;
    double GAM2;
    double PRFDER;
    int IPH;
    double DELTA;
};

class PVII
{
public:
    double TF1,TF2,TF4,TF6,TF8,TF9,C4;
};

class SPVII
{
public:
    double RL,DA1L,DA2L,DA3L,DA4L,DA5L,DA6L,DA7L,
        RH,DA1H,DA2H,DA3H,DA4H,DA5H,DA6H,DA7H;
};

class STRUPHASE
{
public:
    double APHASE[IRS+1];
    double TAVIX[IRS+1];
    double SRIX[IRS+1];
};

class MAXINT
{
public:
    double XMAXINT;
};

class CODEBCK
{
public:
    int IBCKCODE;
};

class G4
{
public:
    double TANN[NOV+1];
    double DERSTO[NOV+1][MSZ+1];
};

class COMP
{
public:
    double CC[4+1][16+1];
    double ZEFF[16+1];
    int PTC[NATS+1];
};

class ATFAT
{
public:
    double FI2[2+1];
};

class DIRCV
{
public:
    double DCSM[6+1][6+1];
    double DCV[6+1];
};

class LABELS
{
public:
    int LB1,LB2,LB3;
    string TEST;
};

class DBWSException {
public:
    DBWSException(string m)
    {
        msg = m;
    };
    ~DBWSException() {};
    const string Show() const
    {
        return msg;
    }
private:
    string msg;
};

class DBWS
{
public:
    DBWS(void);
    virtual ~DBWS(void);
    void run(void);

public:
    string title;

    Phase phases[99+1];

    DC* dc;
    REFLS* refls;
    G1* g1;
    CNTRLS* cntrls;
    PARAMS* params;
    PARAC* parac;
    JNK* jnk;
    ALLP* allp;
    COEFC* coefc;
    COEFF* coeff;
    BKGSCALE* bkgscale;
    DATAX* datax;
    RTMTX* rtmtx;
    SIZESTRAIN* sizestrain;
    VOLUME* volume;
    G2* g2;
    F1* f1;
    G3* g3;
    FONDI* fondi;
    MULTIP* multip;
    SIMOPER* simoper;
    HKLCTL* hklctl;
    CELLX* cellx;
    PRFX* prfx;
    PVII* pvii;
    SPVII* spvii;
    STRUPHASE* struphase;
    MAXINT* maxint;
    CODEBCK* codebck;
    G4* g4;
    COMP* comp;
    ATFAT* atfat;
    DIRCV* dircv;
    LABELS* labels;

    ifstream file3;
    string file3name;

    ifstream file4;
    ofstream file4b;

    string file4name;

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

    ofstream file53; // at INPTR

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
    void SUMMAT(int IPM, double CSK[], double DISK[], double DYCDD[], double ISODER[], double TOTCS);
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
    void RTMT(int* MULTX_, double Y_[][3+1], int* IPRT, int NCTR_[], int* IPHASE);
    void OP1(int* IPHASE, int NCTR_[]);
    void LOOKUP(int K, int N, int NSCAT, int IXRAY, int JOB);
    void CELL2(int NPHASE, double LAMDAM);
    double MULT(int IPHASE, int IH, int IK, int IL, int KXIS);
    void REWRIT(int ISCALE, int IDIF);
    void size(int K);
    void WRITE94(int ISCALE, int IDIF);
    void EXPUT(void);
    void REFGEN(int IPHASE, double ZERO, double DIS, double TRANS, double PREFOR, int NCTR_[]);
    void FINDC(int K, int NSCAT);
    void ABSORP(double MU, double SW, double TH, double* ABC);
    void GSASREAD(void);
    void PHILIPSREAD(void);
    void qpainit(void);
    void rigakuread(void);
    void readasc(void);
    void scintag(void);
    void SIEMENSREAD(void);
    void INPTR(void);
    void inpam(void);

};

void CNV::CV1(double AB[],double V[])
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

void CNV::CV2(double AB[],double* G)
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

void CNV::STEEP(double X[], int N, int M)
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
        cnv.SINLAM[I]=POSI[I]*POSI[I];
        cnv.F[I]=SCAT[I];
    }
    cnv.N=*K;
    NN = min(cnv.N,9);
    for (I=1; I <= NN; ++I)
    {
        AB[I]=I;
    }
    NN1=NN+1;
    for (I=NN1; I <= 10; ++I)
    {
        AB[I]=0.0;
    }
    cnv.STEEP(AB,3,10);
    if (AB[3] == 0.0) AB[3]=-1E-6;
    for (I=1; I <= 9; ++I)
    {
        AC[I][*J]=AB[I];
    }
}




DBWS::DBWS(void)
{
    dc = new DC();
    refls = new REFLS();
    g1 = new G1();
    cntrls = new CNTRLS();
    params = new PARAMS();
    parac = new PARAC();
    jnk = new JNK();
    allp = new ALLP();
    coefc = new COEFC();
    coeff = new COEFF();
    bkgscale = new BKGSCALE();
    datax = new DATAX();
    rtmtx = new RTMTX();
    sizestrain = new SIZESTRAIN();


    volume = new VOLUME();
    g2 = new G2();
    f1 = new F1();
    g3 = new G3();
    fondi = new FONDI();
    multip = new MULTIP();
    simoper = new SIMOPER();
    hklctl = new HKLCTL();
    cellx = new CELLX();
    prfx = new PRFX();
    pvii = new PVII();
    spvii = new SPVII();
    struphase = new STRUPHASE();
    maxint = new MAXINT();
    codebck = new CODEBCK();
    g4 = new G4();
    comp = new COMP();
    atfat = new ATFAT();
    dircv = new DIRCV();
    labels = new LABELS();
}

DBWS::~DBWS(void)
{
    delete dc;
    delete refls;
    delete g1;
    delete cntrls;
    delete params;
    delete parac;
    delete jnk;
    delete allp;
    delete coefc;
    delete coeff;
    delete bkgscale;
    delete datax;
    delete rtmtx;
    delete sizestrain;
    delete volume;
    delete g2;
    delete f1;
    delete g3;
    delete fondi;
    delete multip;
    delete simoper;
    delete hklctl;
    delete cellx;
    delete prfx;
    delete pvii;
    delete spvii;
    delete struphase;
    delete maxint;
    delete codebck;
    delete g4;
    delete comp;
    delete atfat;
    delete dircv;
    delete labels;
}

void DBWS::run(void)
{
    try
    {
        INPTR();
        // Canton et all code starts here !cp may 03 97
        //-----OPEN,if NECESSARY, FILE CONTAINING AMORPHOUS SCATTERING
        if (params->GLB_[20-1] != 0.0 || params->GLB_[20-1].codeword != 0.0 ) inpam();
        // Canton et all code stops here

        ITER();
        if ( cntrls->MAXS > 0 )
        {
            cntrls->MCYCLX =cntrls->MCYCLE;
            cntrls->MAXSX  = cntrls->MAXS;
            cntrls->MCYCLE = 1;
            cntrls->MAXS   = 0;
            //-------LAST CALL TO ITER if MCYCLE = 1 & MAXS = 0
            ITER();
            cntrls->MCYCLE = cntrls->MCYCLX;
            cntrls->MAXS   = cntrls->MAXSX;
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
        for (IIPHAS=2; IIPHAS <= IPHASE; ++IIPHAS) IOF = IOF + refls->ICR[IIPHAS-1];
    }
    if (IPHASE == 1000)
    {
        IC=0;
        for (IIPHAS=1; IIPHAS <= 99; ++IIPHAS) IC  = IC  + refls->ICR[IIPHAS];
    }
    else
    {
        IC=refls->ICR[IPHASE];
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
                R=refls->REFS[INP1][2]-refls->REFS[INP][2];
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
    for (I=1; I <= IC; ++I) ITEMP[I]=refls->IREFS[I+IOF];
    for (I=1; I <= IC; ++I)
    {
        LX=L[I];
        refls->IREFS[I+IOF]=ITEMP[LX];
    }
    for (J=1; J <= 3; ++J)
    {
        for (I=1; I <= IC; ++I) TEMP[I]=refls->REFS[I+IOF][J];
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refls->REFS[I+IOF][J]=TEMP[LX];
        }
    }
    for (I=1; I <= IC; ++I) TEMP[I]=refls->FMGNTD[I+IOF];
    for (I=1; I <= IC; ++I)
    {
        LX=L[I];
        refls->FMGNTD[I+IOF]=TEMP[LX];
    }
    if (cntrls->NPROF == _TCHZ)
    {
        for (I=1; I <= IC; ++I) TEMP[I]=refls->HALFL[I+IOF];
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refls->HALFL[I+IOF]=TEMP[LX];
        }
        for (I=1; I <= IC; ++I) TEMP[I]=refls->HALFG[I+IOF];
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refls->HALFG[I+IOF]=TEMP[LX];
        }
        for (I=1; I <= IC; ++I) TEMP[I]=refls->GAM[I+IOF];
        for (I=1; I <= IC; ++I)
        {
            LX=L[I];
            refls->GAM[I+IOF]=TEMP[LX];
        }
    }
}

void DBWS::ASSIGN_(void)
{
    int I,J,K,IN1,IN2,IRC,ICX,IRK,IRL,IRH,MIN,MAX,IPHAS;
    double PX, WT,YUN, RMIN, RMAX;
    int KRR[2+1][IDSZ+1];

    ICX=0;
    for (J=1; J <= cntrls->NPHASE; ++J)
    {
        ICX = ICX + refls->ICR[J];
    }
    if (cntrls->NPHASE > 1) SORT(1000);
    for (J=1; J <= 2; ++J)
    {
        for (K=1; K <= IDSZ; ++K) KRR[J][K]=0;
    }
    if(cntrls->LST3 == 0) goto L7;
    for (I=1; I <= ICX; ++I)
    {
        IRL=( refls->IREFS[I] % 256)-128;
        IRK=((refls->IREFS[I]/256) % 256)-128;
        IRH=((refls->IREFS[I]/(256*256)) % 256)-128;
        IRC=(refls->IREFS[I]/(256*256*256)) % 8;
        IPHAS=refls->IREFS[I]/(256*256*256*8);
        if ( (I-1 % 60) == 0) file6 << "NO.  CODE    H   K   L  PHASE  HW       POSN" << endl;
        file6 << setw(4) << I
            << setw(4) << IRC << "   "
            << setw(4) << IRH
            << setw(4) << IRK
            << setw(4) << IRL
            << setw(6) << IPHAS
            << setw(8) << setprecision(3) << refls->REFS[I][1]
            << setw(8) << setprecision(3) << refls->REFS[I][2] << endl;
    }
    L7:
    for(I=1; I <= ICX; ++I)
    {
        if (cntrls->NPROF == _SplitPearsonVII)
        {
            RMIN=refls->REFS[I][2]-g1->WDT*refls->FWHM[I][1];
            RMAX=refls->REFS[I][2]+g1->WDT*refls->FWHM[I][2];
        }
        else
        {
            RMIN=refls->REFS[I][2]-g1->WDT*refls->REFS[I][1];
            RMAX=refls->REFS[I][2]+g1->WDT*refls->REFS[I][1];
        }
        MIN= static_cast<int>( (RMIN-g1->THMIN)/g1->STEP+1.5 );
        MAX= static_cast<int>( (RMAX-g1->THMIN)/g1->STEP+1.5 );
        MIN=max(MIN,1);
        MIN=min(MIN,datax->NPTS);
        MAX=min(MAX,datax->NPTS);
        for (K=MIN; K <= MAX; ++K)
        {
            KRR[2][K]=I;
            if(KRR[1][K] == 0)KRR[1][K]=I;
        }
    }
    for (J=1; J <= jnk->NEXCRG; ++J)
    {
        if (jnk->AHIGH[J] <= g1->THMIN) goto L482;
        IN1= static_cast<int>( (jnk->ALOW[J]-g1->THMIN)/g1->STEP+1.5 );
        IN1=max(IN1,1);
        IN2= static_cast<int>( (jnk->AHIGH[J]-g1->THMIN)/g1->STEP+1.5 );
        IN2=min(IN2,datax->NPTS);
        for(I=IN1; I <= IN2; ++I)
        {
            KRR[2][I]=1;
            KRR[1][I]=0;
        }
        if(IN2 == datax->NPTS)goto L484;
        L482:;
    }
    L484:
    if(cntrls->LST2 == 0)goto L530;
    file6 << "PATTERN FROM"
          << setw(8) << setprecision(4) << g1->THMIN
          << " TO" << setw(8) << setprecision(4) << g1->THMAX
          << " IN STEPS OF"  << setw(8) << setprecision(4) << g1->STEP << "DEGREES" << endl;
    file6 << "POSN      I+B     B     I       100*W   K11  K21" << endl;
    for(J=1; J <= datax->NPTS; ++J)
    {
        PX=g1->THMIN + static_cast<double>(J-1) * g1->STEP;
        YUN=datax->Y[J]-datax->BK[J];
        if(cntrls->JOBTYP < 3) WT=1.0/datax->VAR[J];
        file6 << setw(8) << setprecision(4) << PX
              << setw(7) << setprecision(0) << datax->Y[J]
              << setw(7) << setprecision(0) << datax->BK[J]
              << setw(7) << setprecision(0) << YUN
              << setw(9) << setprecision(4) << WT
              << setw(5) << KRR[1][J]
              << setw(5) << KRR[1][J] << endl;
    }
    L530:;
    for (I=1; I <= datax->NPTS; ++I)
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
    for (I=1; I <= datax->NPTS; ++I) datax->KR[I]=KRR[1][I]+IRS *KRR[2][I];
    for (I=1; I <= datax->NPTS; ++I) if(datax->KR[I] != 0 && datax->KR[I] != IRS )goto L603;
    DBWSException("NO REFLECTIONS FOUND");
    // test for detecting the asymmetry model required    !cp ap 16 97
    // new code included in line 2 of ICF
    // iasym = 0 (usual Rietveld asymmetry)
    // iasym = 1 (new asymmetry. Riello, Canton & Fagherazzi.PD 10,3,204-206,1997)
    L603:
    if (cntrls->IASYM == 0)
    {
        // THE FOLLOWING TEST IS REALLY ONLY VALID FOR THE SINGLE PHASE CASE
        if(refls->REFS[ KRR[1][I] ][2] >= g1->RLIM && phases[1-1].PAR[13].L != 0)
        {
            file6 << "ASYMMETRY PARAMETER USAGE INVALID" << endl;
            DBWSException("ASYMMETRY PARAMETER USAGE INVALID");
        }
    }
    else
    {
        if(refls->REFS[ KRR[1][I] ][2] >= (90.0-g1->RLIM) && phases[1-1].PAR[13].L != 0)
        {
            file6 << "ASYMMETRY PARAMETER USAGE INVALID" << endl;
            DBWSException("ASYMMETRY PARAMETER USAGE INVALID");
        }
    }
    for (I=1; I <= cntrls->NPHASE; ++I)
    {
        if(phases[I-1].PAR[11] == 0.0 && phases[I-1].PAR[12] == 0.0 && phases[I-1].PAR[12].L != 0)
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

    g2->S1=0.0;
    S1NOBK=0.0;
    g2->S2=0.0;
    S2NOBK = 0.0;
    g2->SS2=0.0;
    g2->S3=0.0;
    //S3NOBK = 0.0;
    g2->SS4=0.0;
    g2->D1=0.0;
    g2->D2=0.0;
    g2->D4=0.0;
    SDELW = 1.0E+25;
    SDELWP= 1.0E+25;
    D1NOBK = 0.0;
    D2NOBK = 0.0;
    for (I=1; I <= datax->NPTS; ++I)
    {
        if((jnk->NBCKGD != 0 && datax->KR[I] == 0) || datax->KR[I] == IRS )goto L10;
        if(jnk->NEXCRG > 0)
        {
            for (IEXC=1; IEXC <= jnk->NEXCRG; ++IEXC)
            {
                NPTLOW = static_cast<int>( (jnk->ALOW[IEXC]-g1->THMIN)/g1->STEP + 1.5 );
                NPTHI  = static_cast<int>( (jnk->AHIGH[IEXC]-g1->THMIN)/g1->STEP + 1.5 );
                if (I >= NPTLOW && I <= NPTHI) goto L10;
            }
        }
        if (SDELW < 1.0e+24) SDELWP = SDELW;
        DEL=datax->Y[I]-datax->BK[I]-datax->YC[I];
        //DELW = DEL/sqrt(datax->VAR[I]);
        SDELW = DEL;
        g2->S1 = g2->S1 + abs(DEL);
        S1NOBK = S1NOBK + DEL*(datax->Y[I]-datax->BK[I])/datax->Y[I];                              //numerador para calcular Rp(-bcg)
        g2->S2=g2->S2+DEL*DEL/datax->VAR[I];
        S2NOBK=S2NOBK+DEL*DEL*(datax->Y[I]-datax->BK[I])*(datax->Y[I]-datax->BK[I])/(datax->Y[I]*datax->Y[I])/datax->VAR[I];   //numerador para calcular Rwp(-bcg)
        g2->SS2 = g2->SS2 + DEL*DEL;
        if (SDELWP < 1.0e+24) g2->SS4=g2->SS4+(SDELW-SDELWP)*(SDELW-SDELWP);
        L10:;
    }
    for (I=1; I <= datax->NPTS; ++I)
    {
        if((jnk->NBCKGD != 0 && datax->KR[I] == 0) || datax->KR[I] == IRS )goto L15;
        if(jnk->NEXCRG > 0)
        {
            for (IEXC=1; IEXC <= jnk->NEXCRG; ++IEXC)
            {
                NPTLOW = static_cast<int>( (jnk->ALOW[IEXC]-g1->THMIN)/g1->STEP + 1.5 );
                NPTHI  = static_cast<int>( (jnk->AHIGH[IEXC]-g1->THMIN)/g1->STEP + 1.5 );
                if (I >= NPTLOW && I <= NPTHI) goto L15;
            }
        }
        g2->D1=g2->D1+datax->Y[I];
        D1NOBK = D1NOBK+(datax->Y[I]-datax->BK[I]);                                //denominador para calcular Rp(-bcg)
        g2->D2=g2->D2+datax->Y[I]*datax->Y[I]/datax->VAR[I];
        D2NOBK = D2NOBK + (datax->Y[I]-datax->BK[I])*(datax->Y[I]-datax->BK[I])/datax->VAR[I];          //denominador para calcular Rwp(-bcg)
        g2->S3=g2->S3+datax->BK[I]+datax->YC[I];
        L15:;
    }
    g2->R1=0.0;
    g2->R2=100.0*g2->S1/g2->D1;
    g2->R2NOBK = 100.0*S1NOBK/D1NOBK;                                //Rp without background 03nov00 (no writable)
    g2->R3=100.0*sqrt(g2->S2/g2->D2);
    g2->R3NOBK = 100.0*sqrt(S2NOBK/D2NOBK);                           //Rwp without background 03nov00 (no writable)
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
        S  =  STH / g1->LAMDA[ICX];
        S2 =  S * S;
        S4 = S2 * S2;
        CSP[ICX] = 0.0;
        for (I = 1; I <= NK; ++I)
        {
            J = comp->PTC[I+IOF];
            FNUM    =   1.0 + comp->CC[1][J] * S2 + comp->CC[2][J] * S4;
            FDEN    = pow( ( 1.0 + comp->CC[3][J] * S2 + comp->CC[4][J] * S4 ) ,2);
            CSP[ICX]=CSP[ICX] + comp->ZEFF[J] * params->XL[I+IOF][5] * (1.0-(FNUM/FDEN));
        }
        CSP[ICX] = CSP[ICX] * static_cast<double>(IRL);
        //-----UPDATE S
        S  = 2.0 * S;
        S2 =   S * S;
        //-----COMPUTE ASS1
        ASS1 =   1.0 + 0.75 * HMC * g1->LAMDA[ICX] * S2;
        //-----COMPUTE BDF
        BDF  = pow(( 1.0 + 0.50 * HMC * g1->LAMDA[ICX] * S2 ) , 2);
        //-----COMPUTE ASS2 = EXP (-DELTA(MU)*D)
        //      WHERE:
        //      MU       = ABSORPTION COEFFICIENT OF AIR
        //      DELTA(MU)= VARIATION OF MU WITH LAMDA, WAS CALCULATED BY RIELLO
        //                 USING REGRESSION ANALYSIS ON ESTIMATED MU BY ASSUMING
        //                 AIR COMPOSITION 20% O2 AND 80% N2 AT 300�K AND 1 ATM.:
        //                 MU = (E ** 3.40089)*(1.0E-04)*( LAMDA ** 2.79287)
        //      D = 17.3 CM. DISTANCE SPECIMENT-COUNTER OR RADIUS OF THE CAMERA
        ASS2 = exp(-1.5*HMC*17.3*(29.99E-04)* pow(g1->LAMDA[ICX],3.79287)*S2);
        //-----COMPUTE ASS3
        //     FORMULA MUST BE CHANGED FOR OTHER RADIATION AND MONOCHROMATOR
        //-----  ASS3 IS A LORENTZIAN FUNCTION
        ASS3=1/(1+params->GLB_[18-1]* pow(S,params->GLB_[19-1]));
        //----- ASS3 IS A GAUSSIAN FUNCTION
        //       ASS3=EXP(-params->GLB_[18-1]*S2)
        //     ASS3=0.38*EXP(-4.0*S2)+0.62*EXP(-0.15*S2)
        CSP[ICX] = CSP[ICX] * ASS2 * ASS3  / ( ASS1 * BDF );
    }
    //     COMPUTE CISK = COMPTON INTENSITY SCATTERED BY THE K-TH PHASE
    *CISK = params->RATIO[1] * CSP[1] + params->RATIO[2] * CSP[2];
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
        STHL2[ICX] = pow(( STH / g1->LAMDA[ICX] ) , 2);
        if(FONDO == 1)
        {
            //-----COMPUTE FI2 = SUM OF SQUARE SCATTERING FACTORS DUE TO ALL ATOMS
            //                   IN THE CELL AT LAMDA(ICX) CORRECTED FOR THE ISOTROPIC
            //                    THERMAL FACTORS.
            atfat->FI2[ICX] = 0.0;
            for (I = 1; I <= NK; ++I)
            {
                NI = params->PTR[I+IOF];
                FI = 0.0;
                coeff->AC[10][NI] = 0.0;
                for (II = 1; II <= 9; II = II + 2) FI = FI + coeff->AC[II][NI]*exp(-coeff->AC[II+1][NI]*STHL2[ICX]);
                FI = FI + coeff->DFP[NI];
                atfat->FI2[ICX] += params->XL[I+IOF][5] *
                    static_cast<double>(IRL)*( 1.0 - exp(-params->XL[I+IOF][4]*2.0*STHL2[ICX]) )* (pow(FI,2)+pow(coeff->DFPP[NI],2));
                //-----NEXT LINES EVALUATE THE DERIVATES
                if(params->LP[I+IOF][4] != 0 && IDERIV == 2)
                {
                    DERDIS[I+IOF][ICX] = params->XL[I+IOF][5] * 2.0 * STHL2[ICX] *
                        exp(-params->XL[I+IOF][4]*2.0*STHL2[ICX])*static_cast<double>(IRL)*(pow(FI,2)+pow(coeff->DFPP[NI],2));
                }
                else
                {
                    DERDIS[I+IOF][ICX]=0.0;
                }
            }
            DIS[ICX] = atfat->FI2[ICX];
        }
        if(FONDO == 2)
        {
            //-----COMPUTE FI2 = SUM OF SQUARE SCATTERING FACTORS DUE TO ALL ATOMS
            //                   IN THE CELL AT LAMDA(ICX) CORRECTED FOR THE OVERALL
            //                    THERMAL FACTORS.
            atfat->FI2[ICX] = 0.0;
            for (I = 1; I <= NK; ++I)
            {
                NI = params->PTR[I+IOF];
                FI = 0.0;
                coeff->AC[10][NI] = 0.0;
                for (II = 1; II <= 9; II = II + 2) FI += coeff->AC[II][NI]*exp(-coeff->AC[II+1][NI]*STHL2[ICX]);
                FI += coeff->DFP[NI];
                atfat->FI2[ICX] += params->XL[I+IOF][5] * (pow(FI,2)+pow(coeff->DFPP[NI],2));
            }
            atfat->FI2[ICX] += static_cast<double>(IRL);
            DIS[ICX] = ( 1.0 - exp(-phases[K-1].PAR[1]*2.0*STHL2[ICX]) ) * atfat->FI2[ICX];
            LK =phases[K-1].PAR[1].L;
            if(LK != 0 && IDERIV == 2)
            {
                DER[ICX] = 2.0 * STHL2[ICX] * atfat->FI2[ICX] * exp(-phases[K-1].PAR[1]*2.0*STHL2[ICX]);
            }
            else
            {
                DER[ICX] = 0.0;
            }
        }
    }
    //-----COMPUTE SDK = SCATTERING DISORDER DUE TO THE K-TH PHASE
    *SDK = params->RATIO[1] * DIS[1] + params->RATIO[2] * DIS[2];
    //-----COMPUTE  DERIVATIVE OF YC RESPECT TO ISOTROPIC THERMAL FACTORS
    //                   XL(I+IOF,4) IN THE K-TH PHASE
    if (FONDO == 1)
    {
        for (I=1; I <= NK; ++I)
        {
            if(params->LP[I+IOF][4] != 0 && IDERIV == 2)
            {
                DERISO[I+IOF]=params->RATIO[1] * DERDIS[I+IOF][1] + params->RATIO[2] * DERDIS[I+IOF][2];
            }
            else
            {
                DERISO[I+IOF]=0.0;
            }
        }
    }
    //-----COMPUTE  DERIVATIVE OF YC RESPECT TO OVERALL THERMAL FACTOR
    //                   params->PAR[K][2] IN THE K-TH PHASE
    if(FONDO == 2)
    {
        LK =phases[K-1].PAR[1].L;				// TODO: Coloquei talez esteja incorreto, verificar depois de retirar os gotos
        if(LK != 0 && IDERIV == 2)
        {
            *DYC = params->RATIO[1] * DER[1] + params->RATIO[2] * DER[2];
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
    if ( hklctl->IHKL[1][1] < 0 )
    {
        goto L11;
    }
    else if ( hklctl->IHKL[1][1] == 0 )
    {
        goto L9;
    }
    else
    {
        goto L14;
    }
L9:
    if ( hklctl->IHKL[2][1] < 0 )
    {
        goto L12;
    }
    else if ( hklctl->IHKL[2][1] == 0 )
    {
        goto L10;
    }
    else
    {
        goto L14;
    }
L10:
    if ( hklctl->IHKL[3][1] < 0 )
    {
        goto L13;
    }
    else
    {
        goto L14;
    }
L11:
    hklctl->IHKL[1][1] = -hklctl->IHKL[1][1];
L12:
    hklctl->IHKL[2][1] = -hklctl->IHKL[2][1];
L13:
    hklctl->IHKL[3][1] = -hklctl->IHKL[3][1];
L14:
    IPH[1] = hklctl->IHKL[3][1]+512*(hklctl->IHKL[2][1]+512*hklctl->IHKL[1][1]);
    CI = 1.0;
    hklctl->ICHKL[*IPHASE] = 1;
    hklctl->AZ[1] = 0.0;
    LD[1][1] = 0;
    LD[2][1] = 0;
    LD[3][1] = 0;
    NCH[1] = 0;
    NCK[1] = 0;
    NCL[1] = 1;
    if ( hklctl->N1HKL[*IPHASE] == 0 ) goto L1000;
    L = 2;
    for (I=1; I <= hklctl->N1HKL[*IPHASE]; ++I)
    {
        CI = CI*2.0;
        if ( hklctl->NC1[I][1][*IPHASE] > 0 ) CI = CI*1.5;
        for (J=1; J <= hklctl->ICHKL[*IPHASE]; ++J)
        {
            NCO = hklctl->NC1[I][1][*IPHASE]+1;
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
            if( ((hklctl->IHKL[1][1]-hklctl->IHKL[2][1]-hklctl->IHKL[3][1]) % 3) != 0 ) goto L1002;
            goto L900;
L20:
            NCH[L] =  NCH[J] ^ hklctl->NC1[I][2][*IPHASE];
            if ( (hklctl->NC1[I][2][*IPHASE] % 2) != 0 ) NCH[L] = NCK[J] ^ hklctl->NC1[I][2][*IPHASE];
            M = 1+(NCH[L] % 2);
            MS = 1-2*((NCH[L]/2) % 2);
            N = abs(NCL[J]);
            hklctl->IHKL[M][L] = hklctl->IHKL[1][N]*MS;
            M = 1+(hklctl->NC1[I][2][*IPHASE] % 2);
            MS = 1-2*((hklctl->NC1[I][2][*IPHASE]/2) % 2);
            NCX = hklctl->NC1[I][3][*IPHASE]+1;
            LD[1][L] = KD[NCX]+LD[M][J]*MS;
            NCK[L] =  NCK[J] ^ hklctl->NC1[I][4][*IPHASE];
            if ( (hklctl->NC1[I][4][*IPHASE] % 2) != 0 ) NCK[L]= NCH[J] ^ hklctl->NC1[I][4][*IPHASE];
            M = 2-(NCK[L] % 2);
            MS = 1-2*((NCK[L]/2) % 2);
            hklctl->IHKL[M][L] = hklctl->IHKL[2][N]*MS;
            M = 2-(hklctl->NC1[I][4][*IPHASE] % 2);
            MS = 1-2*((hklctl->NC1[I][4][*IPHASE]/2) % 2);
            NCX = hklctl->NC1[I][5][*IPHASE]+1;
            LD[2][L] = KD[NCX]+LD[M][J]*MS;
            MS = (1-2*hklctl->NC1[I][6][*IPHASE]) *  sign(NCL[J]);
            NCL[L] = MS*NCL[N];
            hklctl->IHKL[3][L] = hklctl->IHKL[3][N]*MS;
            //L70:
            NCO = hklctl->NC1[I][7][*IPHASE]+1;
            if ( NCO > 5 ) NCO=6;
            MS = 1-2*hklctl->NC1[I][6][*IPHASE];
            if ( N > 1 && NCK[J] >= 4 && M == 1 ) MS=-MS;
            LD[3][L] = KD[NCO]+LD[3][J]*MS;
            //L80:
            IXIT = 0;
L81:
            IPH[L] = hklctl->IHKL[3][L]+512*(hklctl->IHKL[2][L]+512*hklctl->IHKL[1][L]);
            //L87:
            LA = LD[1][L]*hklctl->IHKL[1][1]+LD[2][L]*hklctl->IHKL[2][1]+LD[3][L]*hklctl->IHKL[3][1];
            hklctl->AZ[L] = static_cast<double>(LA)/12.0;
            if ( phases[*IPHASE].SYMB.NC != 0 ) goto L89;
            if ( phases[*IPHASE].SYMB.NSPGRP < 12 ) goto L89;
            if ( hklctl->NC1[I][1][*IPHASE] > 0 ) goto L89;
            if ( (hklctl->NC1[I][2][*IPHASE] % 2) == 1 ) goto L89;
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
            if ( fmod(abs(hklctl->AZ[L]-hklctl->AZ[M-1]),1.0) != 0.0 ) goto L1002;
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
            K1JD = hklctl->IHKL[1][JD];
            hklctl->IHKL[1][L] = hklctl->IHKL[2][JD];
            hklctl->IHKL[2][L] = hklctl->IHKL[3][JD];
            hklctl->IHKL[3][L] = K1JD;
            JD = L;
            goto L81;
L800:
            IXIT = 4;
            JD = J;
            NCO = hklctl->NC1[I][6][*IPHASE];
L850:
            IXIT = IXIT-2;
            NCH[L] = 0;
            NCK[L] = 4;
            NCL[L] = L;
            LD[1][L] = 0;
            LD[2][L] = 0;
            hklctl->IHKL[1][L] = -hklctl->IHKL[1][JD]-hklctl->IHKL[2][JD];
            hklctl->IHKL[2][L] = hklctl->IHKL[1][JD];
            hklctl->IHKL[3][L] = hklctl->IHKL[3][JD];
            LD[3][L] = KF[NCO]+LD[3][JD];
            if ( (NCH[J] % 2) != 0 ) LD[3][L]=KE[NCO]+LD[3][JD];
            JD = L;
            goto L81;
L890:;
        }
L900:
        hklctl->ICHKL[*IPHASE] = L-1;
    }
L1000:
    hklctl->IER = 0;
    if ( phases[*IPHASE].SYMB.NC != 0 ) goto L1001;
    for (M=2; M <= hklctl->ICHKL[*IPHASE]; ++M)
    {
        if ( IPH[M] > 0 ) goto L1003;
        hklctl->IHKL[1][M] = -hklctl->IHKL[1][M];
        hklctl->IHKL[2][M] = -hklctl->IHKL[2][M];
        hklctl->IHKL[3][M] = -hklctl->IHKL[3][M];
L1003:;
    }
L1001:
    return;
L1002:
    hklctl->IER = 1;
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
    //double SINTL[NOV+1];
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
    prfx->IPH = refls->IREFS[NN]/(256*256*256*8);  // 256/256/256/8
    IOF = 0;
    if ( prfx->IPH > 1 )
    {
        for (IIPHAS=2; IIPHAS <= prfx->IPH; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
    }
    IRL = phases[prfx->IPH].MLTPHS;
    N = phases[prfx->IPH].AtomCount;
    ICENT = phases[prfx->IPH].ICNTPHS;
    //-----ZEROIZE THE DERIVATIVES OF THIS REFLECTION W.R.T. TO PARAMETERS
    NX=0;
    for (IIPHAS=1; IIPHAS <= cntrls->NPHASE; ++IIPHAS) NX = NX+phases[IIPHAS].AtomCount;
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
    PAC = phases[prfx->IPH-1].PAR[11] != 0.0  ||  phases[prfx->IPH-1].PAR[11].L != 0  ||
          phases[prfx->IPH-1].PAR[12] != 0.0  ||  phases[prfx->IPH-1].PAR[12].L != 0;
    for (I=1; I <= 3; ++I) AL[I][I]=phases[prfx->IPH-1].PAR[I+5-1];
    AL[3][2]=phases[prfx->IPH-1].PAR[8];
    AL[2][3]=AL[3][2];
    AL[3][1]=phases[prfx->IPH-1].PAR[9];
    AL[1][3]=AL[3][1];
    AL[2][1]=phases[prfx->IPH-1].PAR[10];
    AL[1][2]=AL[2][1];
    AH=phases[prfx->IPH-1].PAR[2];
    BH=phases[prfx->IPH-1].PAR[3];
    CH=phases[prfx->IPH-1].PAR[4];
    DH=phases[prfx->IPH-1].PAR[19];
    if (cntrls->NPROF == _pseudoVoigt) EH=phases[prfx->IPH-1].PAR[20];
    if (cntrls->NPROF == _TCHZ)
    {
        AH2=phases[prfx->IPH-1].PAR[14];
        BH2=phases[prfx->IPH-1].PAR[15];
    }
    B1=phases[prfx->IPH-1].PAR[11];
    B2=phases[prfx->IPH-1].PAR[12];
    NM=(NN % NOV)+1;
    ICX=(refls->IREFS[NN]/(256*256*256)) % 8;
    SLABDA=g1->LAMDA[ICX]*g1->LAMDA[ICX]/4.0;
    N=phases[prfx->IPH].AtomCount;
    refls->FMGNTD[NN]=0.0;
    PAKNN=0.0;
    PRECOR=1.0;
    TR=0.0;
    HNN[3]=(refls->IREFS[NN] % 256)-128;
    HNN[2]=((refls->IREFS[NN]/256) % 256)-128;
    HNN[1]=((refls->IREFS[NN]/(256*256)) % 256)-128;

    for (I=1; I <= 3; ++I) H[I]=HNN[I];
    //-----CALCULATION OF TEMP.FACTOR,POSITION AND FWHM
    SS=0.0;
    for (I=1; I <= 3; ++I)
    {
        hklctl->IHKL[I][1] = HNN[I];
        for (J=I; J <= 3; ++J) SS = HNN[I]*AL[I][J]*HNN[J]+SS;                         // GSAS DH2
    }
    if ( PAC )
    {
        TT = 0.0;
        for (I=1; I <= 3; ++I)
        {
            for (J=I; J <= 3; ++J) TT = TT + phases[prfx->IPH].PREF[I]*AL[I][J] * phases[prfx->IPH].PREF[J];             // GSAS DP2
        }
        SMTRY2(&prfx->IPH);
        PRECOR = 0.0;
        PREXP = 0.0;
        DPRECOR = 0.0;
        for (IJ=1; IJ <= hklctl->ICHKL[prfx->IPH]; ++IJ)
        {
            PAK = 0.0;
            for (I=1; I <= 3; ++I)
            {
                for (J=I; J <= 3; ++J) PAK = phases[prfx->IPH].PREF[I]*AL[I][J]*static_cast<double>(hklctl->IHKL[J][IJ])+PAK;     // GSAS CA
            }
            PAK = PAK*PAK/(TT*SS);
            PAKNN = pow((PI/2.0),2);
            if (PAK != 0) PAKNN=pow(atan(sqrt(abs((1.0-PAK)/PAK))),2);
            if (cntrls->IPREF == 0)
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
        PREXP = PREXP/static_cast<double>(hklctl->ICHKL[prfx->IPH]);
        PRECOR = PRECOR/static_cast<double>(hklctl->ICHKL[prfx->IPH]);
        DPRECOR = DPRECOR/static_cast<double>(hklctl->ICHKL[prfx->IPH]);
    }
    //L13:
    //SINTL[NM] = sqrt(SS);
    SSNN = 0.25*SS;
    TAV = exp(-2.0*phases[prfx->IPH-1].PAR[1]*SSNN)*PRECOR;
    SINTH = SLABDA*SS;
    COSTH = 1.0-SINTH;
    TANTH = sqrt(SINTH/COSTH);
    g4->TANN[NM] = TANTH;
    //   Correction of microabsorption
    if (cntrls->IABSR == 1)
    {
        ISITH=(1.0)/(sqrt(SINTH));
        SR = params->GLB_[12-1]*(1.0-params->GLB_[8-1]*exp(-params->GLB_[9-1])+params->GLB_[8-1]*exp(-params->GLB_[9-1]/sqrt(SINTH)))+(1.0-params->GLB_[12-1])*(1+params->GLB_[13-1]*(asin(sqrt(SINTH)))-1.5707963268);
    }
    else if (cntrls->IABSR == 2)
    {
        ISITH = sqrt(SINTH);
        SR = 1.0-params->GLB_[13-1]*(asin(ISITH)-1.5707963268);
    }
    else if (cntrls->IABSR == 3)
    {
        ISITH=(1.0)/(sqrt(SINTH));
        SR = 1.0-params->GLB_[8-1]*exp(-params->GLB_[9-1])+params->GLB_[8-1]*exp(-params->GLB_[9-1]*ISITH);
    }
    else if (cntrls->IABSR == 4)
    {
        ISITH=(1.0)/(sqrt(SINTH));
        SR = 1.0-params->GLB_[8-1]*params->GLB_[9-1]*(1.0-params->GLB_[9-1])-ISITH*params->GLB_[8-1]*params->GLB_[9-1]*(1.0-params->GLB_[9-1]*ISITH);
    }

    refls->REFS[NN][2]=atan(TANTH)*PI;
    refls->HALFG[NN]=(AH*TANTH*TANTH+BH*TANTH+CH+DH/COSTH+EH/(TANTH*TANTH));
    if (refls->HALFG[NN] > 0.0)
    {
        refls->HALFG[NN] = sqrt(refls->HALFG[NN]);
    }
    else
    {
        cout << "1 = " << HNN[1] << " " << HNN[2] << " " << HNN[3]  << endl;

        file6 << "   SQUARE OF FWHM NEGATIVE AT TWO-THETA=" << setw(8) << setprecision(3) << refls->REFS[NN][2] << " FOR PHASE NO. " << setw(4) << prfx->IPH << endl;
        cout << "   SQUARE OF FWHM NEGATIVE AT TWO-THETA=" << setw(8) << setprecision(3) << refls->REFS[NN][2] << " FOR PHASE NO. " << setw(4) << prfx->IPH << endl;
        DBWSException("SQUARE OF FWHM IS NEGATIVE");
    }
    if (cntrls->NPROF == _TCHZ)
    {
        refls->HALFL[NN] = AH2*TANTH+BH2*sqrt(1.0+TANTH*TANTH);
        BB = pow((pow(refls->HALFG[NN],5.0)+2.69269*pow(refls->HALFG[NN],4.0)*refls->HALFL[NN]+2.42843*pow(refls->HALFG[NN],3.0)*pow(refls->HALFL[NN],2.0)+4.47163*pow(refls->HALFG[NN],2.0)*pow(refls->HALFL[NN],3.0) +0.07842*refls->HALFG[NN]*pow(refls->HALFL[NN],4.0) + pow(refls->HALFL[NN],5.0)),0.2);
        TLR = refls->HALFL[NN]/BB;
        refls->GAM[NN] = 1.36603*TLR-0.47719*TLR*TLR+0.11116*pow(TLR,3.0);
    }
    else
    {
        BB = refls->HALFG[NN];
    }
    prfx->TL=BB;
    refls->REFS[NN][1]=BB;
    BB=BB*BB;
    //-----VERT=.TRUE. if ASYMMETRY CORRECTION IS TO BE CALCULATED
    VERT=false;
    if(cntrls->IASYM == 0)
    {
        if(refls->REFS[NN][2] <= g1->RLIM && cntrls->NPROF != _SplitPearsonVII) VERT=true;
    }
    else
    {
        if (abs(refls->REFS[NN][2]-90.0) >= g1->RLIM) VERT=true;
    }
    //-----CALCULATION OF COS(H.X),SIN(H.X) AND TEMP. FACTOR FOR EACH ATOM
    //L3:
    for (I=1; I <= N; ++I)
    {
        SNXI=0.0;
        SAI=0.0;
        SBI=0.0;
        for (J=1; J <= 11; ++J) XI[J]=params->XL[I+IOF][J];
        for (IR=1; IR <= IRL; ++IR)
        {
            for (J=1; J <= 3; ++J)
            {
                IV=phases[prfx->IPH].IVEC[IR] /32768 / static_cast<int>(pow(32,3-J));
                IV=(IV % 32);
                SM[J][1]=IV/9-1;
                SM[J][2]=((IV/3) % 3)-1;
                SM[J][3]=(IV % 3)-1;
                T[J] = static_cast<double>( ((phases[prfx->IPH].IVEC[IR]/
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
            for (J=1; J <= 3; ++J) ARG=H[J]*XI[J]+ARG;
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
                if(ICENT == 1)SUMBX[I][JJ]=SUMBX[I][JJ]+H[JJ]*COSA;
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
            //L10:;
        }
        TEMP[I]=exp(-params->XL[I+IOF][4]*SSNN);
        SA[I]=SAI;
        SB[I]=SBI;
        NI=params->PTR[I+IOF];
        BNI=coeff->DFPP[NI];
        FFX=coeff->DFP[NI];
        coeff->AC[10][NI]=0.0;
        for (II=1; II <= 9; II = II + 2) FFX=FFX+coeff->AC[II][NI]*exp(-coeff->AC[II+1][NI]*SSNN);
        SNEXI=FFX*params->XL[I+IOF][5]*TEMP[I];
        SNXI=BNI*params->XL[I+IOF][5]*TEMP[I];
        //-----CALCULATE A AND B OF F
        AV=SNEXI*SAI+AV;
        BV=SNEXI*SBI+BV;
        SNEX[I]=2.0*SNEXI*TAV*params->RATIO[ICX];
        SNX[I]=2.0*SNXI*TAV*params->RATIO[ICX];
        CV=CV+SNXI*SAI;
        DV=DV+SNXI*SBI;
    }
    FNN=params->RATIO[ICX]*(CV*CV+AV*AV+DV*DV+BV*BV)*TAV*SR;

    // PREPARING PHASE and struc fact TO BE PRINTED  !cp june 98
    struphase->TAVIX[NN]=TAV;
    struphase->SRIX[NN]=SR;
    if(AV == 0)AV=1E-6;
    struphase->APHASE[NN]=atan(BV/AV);
    if(AV < 0.0) struphase->APHASE[NN] = struphase->APHASE[NN]+3.1415927359;
    if(BV < 0 && AV == 0) struphase->APHASE[NN] =1.5*3.1415927359;
    if(BV > 0 && AV == 0) struphase->APHASE[NN] =0.5*3.1415927359;
    //L120:
    refls->FMGNTD[NN]=FNN;
    if(cntrls->MAXS == 0) return;

    //-----CALCULATE DERIVATIVES
    for (I=1; I <= N; ++I)
    {
        SNEXI=SNEX[I];
        SNXI=SNX[I];
        SAI=SA[I];
        SBI=SB[I];
        for (J=1; J <= 11; ++J)
        {
            K=params->LP[I+IOF][J];
            if(K == 0) goto L22;
            if(J > 5) goto L221;
            if(J > 3)goto L29;
            SUMA=SUMAX[I][J];
            SUMB=SUMBX[I][J];
            DER=-((AV*SUMA-BV*SUMB)*SNEXI+(CV*SUMA-DV*SUMB)*SNXI)*6.2831853071;
            goto L26;
L29:
            if(J > 4)goto L31;
            DER=-((SAI*AV+SBI*BV)*SNEXI+(SAI*CV+SBI*DV)*SNXI)*SSNN;
            goto L26;
L31:
            DER=((SAI*AV+SBI*BV)*SNEXI+(SAI*CV+SBI*DV)*SNXI)/params->XL[I+IOF][5];
            goto L26;
L221:
            SUMA=SUMAX[I][J-2];
            SUMB=SUMBX[I][J-2];
            DER=-((AV*SUMA+BV*SUMB)*SNEXI+(CV*SUMA+DV*SUMB)*SNXI);
            if(J >= 9) DER = 2.0*DER;
L26:
            DERIV[K] = sign(params->A[I+IOF][J])*DER+DERIV[K];
L22:;
        }
    }
    //-----CALCULATE DERIVATIVES
    //-----Preferred Orientation Derivatives
    K = phases[prfx->IPH-1].PAR[11].L;
    if ( K != 0 )
    {
        DERIV[K] = DERIV[K]+FNN*DPRECOR;
        //         print '(3x,i3,3f10.5)',k,fnn,dprecor,deriv(k)         ! ********
    }

    //////////////////////////////////////////////////////////
    // TODO: remover este isso. Aideia é dividir a rotina INPTR em duas e fazer estes teste por lá

    K = phases[prfx->IPH-1].PAR[12].L;
    if (cntrls->IPREF == 0)
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
    if (cntrls->IABSR == 1)
    {
        K = params->GLB_[13-1].L;
        SRD = (1.0-params->GLB_[12-1])*(asin(sqrt(SINTH))-1.5707963268);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = params->GLB_[12-1].L;
        SRD = 1.0 -params->GLB_[8-1]*exp(-params->GLB_[9-1])+params->GLB_[8-1]*exp(-params->GLB_[9-1]/sqrt(SINTH))-1.0-params->GLB_[13-1]*(asin(sqrt(SINTH))-1.5707963268);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = params->GLB_[9-1].L;
        SRD = params->GLB_[8-1]*params->GLB_[12-1]*(exp(-params->GLB_[9-1])-exp(-params->GLB_[9-1]/sqrt(SINTH))/sqrt(SINTH));
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = params->GLB_[8-1].L;
        SRD = -params->GLB_[12-1]*(exp(-params->GLB_[9-1])+exp(-params->GLB_[9-1]/sqrt(SINTH)));
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
    }
    else if (cntrls->IABSR == 2)
    {
        ///////////////////////////////
        // TODO: Eliminar este if
        KKL = params->GLB_[12-1].L;
        K  = params->GLB_[9-1].L;
        KL = params->GLB_[8-1].L;
        if (K != 0 || KL != 0 || KKL != 0)
        {
            file6 << "P AND/OR Q AND/OR R ARE NOT REFINABLE PARAMETERS FOR IABSR=2" << endl;
            cout << "P AND/OR Q AND/OR R ARE NOT REFINABLE PARAMETERS FOR IABSR=2" << endl;
            DBWSException("");
        }
        ///////////////////////////

        K = params->GLB_[13-1].L;
        SRD = 1.5707963268 - asin(sqrt(SINTH));
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
    }
    else if (cntrls->IABSR == 3)
    {
        /////////////////////////////////////////////////////////
        // TODO: Eliminar este if!!!
        KKL = params->GLB_[13-1].L;
        K = params->GLB_[12-1].L;
        if (K != 0 || KKL != 0)
        {
            file6 << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            cout << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            DBWSException("");
        }
        ////////////////////////////////////////////////////////////

        K = params->GLB_[9-1].L;
        SRD = params->GLB_[8-1]*exp(-params->GLB_[9-1])-params->GLB_[8-1]*ISITH*exp(-params->GLB_[9-1]*ISITH);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN;
        K = params->GLB_[8-1].L;
        SRD = -exp(-params->GLB_[9-1])+exp(-params->GLB_[9-1]*ISITH);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN;
    }
    else if (cntrls->IABSR == 4)
    {
        //////////////////////////////////////////
        // TODO: Eliminar este if!!!!
        KKL = params->GLB_[13-1].L;
        K = params->GLB_[12-1].L;
        if (K != 0 || KKL != 0)
        {
            file6 << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            cout << "R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR CHOICE" << endl;
            DBWSException("");
        }
        /////////////////////////////

        K = params->GLB_[9-1].L;
        SRD = params->GLB_[8-1]*(2*params->GLB_[9-1]-1)+params->GLB_[8-1]*ISITH*(2*params->GLB_[9-1]*ISITH-1);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
        K = params->GLB_[8-1].L;
        SRD = params->GLB_[9-1]*(params->GLB_[9-1]-1)-params->GLB_[9-1]*ISITH*(1-params->GLB_[9-1]*ISITH);
        if (K != 0) DERIV[K] = DERIV[K] + SRD*FNN/SR;
    }
    //----Overall Temperature and Scale Factor
    K=phases[prfx->IPH-1].PAR[1].L;
    if(K != 0) DERIV[K]=DERIV[K]-2.0*SSNN*FNN;
    K=phases[prfx->IPH-1].PAR[0].L;
    if(K != 0) DERIV[K]=DERIV[K]+FNN/phases[prfx->IPH-1].PAR[0];
    SINTH=FNN*PI*SLABDA/(sqrt(SINTH*COSTH)*BB);
    SS=FNN/prfx->TL;
    X=TANTH*TANTH;
    //-----Broadening Derivatives
    if (cntrls->NPROF == _SplitPearsonVII || cntrls->NPROF == _TCHZ) goto L9212;
    for (J=3; J <= 5; ++J)
    {
        K=phases[prfx->IPH-1].PAR[J-1].L;
        if(K == 0) goto L78;
        DERIV[K]=X*SS+DERIV[K];
L78:
        X=X/TANTH;
    }
    K=phases[prfx->IPH-1].PAR[20].L;
    if(K == 0) goto L9212;
    // for cot^2 case                                    !cp nov 29 96
    DERIV[K]=SS/(TANTH*TANTH)+DERIV[K];
    //-----Split Pearson VII Broadening Derivatives
L9212:
    if (cntrls->NPROF == _SplitPearsonVII)
    {
        if (prfx->DELTA < 0.0)
        {
            DA3 = spvii->DA3L;
            DA1 = spvii->DA1L;
        }
        else
        {
            DA3 = spvii->DA3H;
            DA1 = spvii->DA1H;
        }
        for (J=3; J <= 5; ++J)
        {
            K=phases[prfx->IPH-1].PAR[J-1].L;
            if(K == 0) goto L780;
            DERIV[K]=DERIV[K]+X*SS*((DA3*prfx->DELT/(1.0+DA1*prfx->DELT))-(1.0/prfx->TL));
L780:
            X=X/TANTH;
        }
    }
    //L9211:

    //-----TCHZ Broadening Derivatives
    if (cntrls->NPROF == _TCHZ)
    {
        prfx->TL = refls->REFS[NN][1];
        TLG = refls->HALFG[NN];
        TLL = refls->HALFL[NN];
        DHDHG = 0.2/pow(prfx->TL,4.0)*(5.*pow(TLG,4.0)+10.77076*pow(TLG,3.0)*TLL+ 7.28529*TLG*TLG*TLL*TLL+8.94326*TLG*pow(TLL,3.0) + 0.07842*pow(TLL,4.0));
        DHDHL = 0.2/pow(prfx->TL,4.0)*(2.69269*pow(TLG,4.0)+ 4.85686*pow(TLG,3.0)*TLL +13.41489*TLG*TLG*TLL*TLL + 0.31368*TLG*pow(TLL,3.0)+5.*pow(TLL,4.0));
        for (J=3; J <= 5; ++J)
        {
            K=phases[prfx->IPH-1].PAR[J-1].L;
            if(K == 0)goto L9078;
            DERIV[K]=DHDHG*X*SS+DERIV[K];
L9078:
            X=X/TANTH;
        }
        K=phases[prfx->IPH-1].PAR[19].L;
        if  (K != 0) DERIV[K]=DHDHG*SS/COSTH+DERIV[K];
        K = phases[prfx->IPH-1].PAR[14].L;
        if (K == 0) goto L9213;
        DERIV[K] = 2.0*FNN*DHDHL*TANTH+DERIV[K];
L9213:
        K = phases[prfx->IPH-1].PAR[15].L;
        if (K == 0) goto L9214;
        DERIV[K] = 2.0*FNN*DHDHL/sqrt(COSTH) + DERIV[K];
L9214:;
    }
    //-----Profile Shape Derivatives
    K = phases[prfx->IPH-1].PAR[16].L;
    if(K != 0 && (cntrls->NPROF == _pseudoVoigt || cntrls->NPROF == _PearsonVII)) DERIV[K]=DERIV[K]+ FNN;
    K = phases[prfx->IPH-1].PAR[17].L;
    if(K != 0 && cntrls->NPROF == _pseudoVoigt) DERIV[K]=DERIV[K]+ FNN * refls->REFS[NN][2];
    if(K != 0 && cntrls->NPROF == _PearsonVII) DERIV[K]=DERIV[K]+ FNN / refls->REFS[NN][2];
    K = phases[prfx->IPH-1].PAR[18].L;
    if(K != 0 && cntrls->NPROF == _PearsonVII)DERIV[K]=DERIV[K]+FNN/refls->REFS[NN][2] /refls->REFS[NN][2];
    //-----Split Pearson VII Shape Derivative
    if (cntrls->NPROF == _SplitPearsonVII)
    {
        K = phases[prfx->IPH-1].PAR[16].L;
        if (K != 0.0)
        {
            if (prfx->DELTA < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+spvii->DA1L*prfx->DELT/BB)+ spvii->DA7L*prfx->DELT/BB/(1.0+spvii->DA1L*prfx->DELT/BB));
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*spvii->DA6L;
            }
        }
        K = phases[prfx->IPH-1].PAR[17].L;
        if (K != 0.0)
        {
            if (prfx->DELTA < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+spvii->DA1L*prfx->DELT/BB)+ spvii->DA7L*prfx->DELT/BB/(1.0+spvii->DA1L*prfx->DELT/BB))/refls->REFS[NN][2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*spvii->DA6L/refls->REFS[NN][2];
            }
        }
        K = phases[prfx->IPH-1].PAR[18].L;
        if (K != 0.0)
        {
            if (prfx->DELTA < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+spvii->DA1L*prfx->DELT/BB)+spvii->DA7L*prfx->DELT/BB/(1.0+spvii->DA1L*prfx->DELT/BB))/refls->REFS[NN][2]/refls->REFS[NN][2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*spvii->DA6L/refls->REFS[NN][2]/refls->REFS[NN][2];
            }
        }
        K = phases[prfx->IPH-1].PAR[23].L;
        if (K != 0.0)
        {
            if (prfx->DELTA < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*spvii->DA6H;
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+spvii->DA1H*prfx->DELT/BB)+spvii->DA7H*prfx->DELT/BB/(1.0+spvii->DA1H*prfx->DELT/BB));
            }
        }
        K = phases[prfx->IPH-1].PAR[24].L;
        if (K != 0.0)
        {
            if (prfx->DELTA < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*spvii->DA6H/refls->REFS[NN][2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+spvii->DA1H*prfx->DELT/BB)+spvii->DA7H*prfx->DELT/BB/(1.0+spvii->DA1H*prfx->DELT/BB))/refls->REFS[NN][2];
            }
        }
        K = phases[prfx->IPH-1].PAR[25].L;
        if (K != 0.0)
        {
            if (prfx->DELTA < 0.0)
            {
                DERIV[K] = DERIV[K]+FNN*spvii->DA6H/refls->REFS[NN][2]/refls->REFS[NN][2];
            }
            else
            {
                DERIV[K] = DERIV[K]+FNN*(-log(1.0+spvii->DA1H*prfx->DELT/BB)+spvii->DA7H*prfx->DELT/BB/(1.0+spvii->DA1H*prfx->DELT/BB))/refls->REFS[NN][2]/refls->REFS[NN][2];
            }
        }
    }
    //-----Split Pearson VII Asymmetry Derivative
    K = phases[prfx->IPH-1].PAR[26].L;
    if (prfx->DELTA < 0.0)
    {
        DA1 = spvii->DA1L;
        DA4 = spvii->DA4L;
        DA5 = spvii->DA5L;
    }
    else
    {
        DA1 = spvii->DA1H;
        DA4 = spvii->DA4H;
        DA5 = spvii->DA5H;
    }
    if (K != 0 && cntrls->NPROF == _SplitPearsonVII)
    {
        DERIV[K]=DERIV[K]+FNN*(DA4+DA5*prfx->DELT/BB/(1.0+DA1*prfx->DELT/BB));
    }
    //-----Zero, Displacement, and Transparancy Derivatives
    K=params->GLB_[1-1].L;
    if(K != 0)DERIV[K]=DERIV[K]+2.0*FNN/BB;
    K=params->GLB_[10-1].L;
    if (K != 0) DERIV[K]=DERIV[K]+2.0*FNN/BB*sqrt(COSTH);
    K=params->GLB_[11-1].L;
    if (K != 0) DERIV[K]=DERIV[K]+2.0*FNN/BB*sin(refls->REFS[NN][2]/57.2958);
    //c-----Lattice Parameter Derivatives
    for (J=1; J <= 6; ++J)
    {
        K=phases[prfx->IPH-1].PAR[J+5-1].L;
        if(K == 0) goto L79;
        if(J < 4) X=HNN[J]*HNN[J];
        if(J == 4) X=HNN[2]*HNN[3];
        if(J == 5) X=HNN[1]*HNN[3];
        if(J == 6) X=HNN[1]*HNN[2];
        DERIV[K]= X*SINTH + DERIV[K];
L79:;
    }
    //-Asymmetry Derivative.  Test for asymmetry model included !cp may 97
    K=phases[prfx->IPH-1].PAR[13].L;
    if((K != 0) && VERT)
    {
        if (cntrls->IASYM == 0)
        {
            DERIV[K]=-FNN/TANTH+DERIV[K];
        }
        else
        {
            TANTHE=TANTH;
            if (TANTHE >= 1.0) TANTHE=tan(atan(TANTHE-3.14159265359/2));
            DERIV[K]=-FNN/TANTH+DERIV[K];
        }
    }
    //-----STORE DERIVATIVES FOR LIMO REFLECTIONS AT A TIME
    for (I=1; I <= cntrls->MAXS; ++I) g4->DERSTO[NM][I]=DERIV[I];
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
    const double PIG = 3.1415926;
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
    pvii->TF1=sqrt( pow(2.0,(1.0/T))-1.0);
    pvii->C4=2.0*RG*pvii->TF1/sqrt(PIG);
    pvii->TF2=AL2* pow(2.0,(1.0/T))/(2.0*pow(pvii->TF1,2));
    pvii->TF4=4.0*AL2* pow(2.0,(1.0/T))/T;
    pvii->TF6=4.0*T*pow(pvii->TF1,2);
    pvii->TF8=DG/RG-pvii->TF2/pow(T,2);
    pvii->TF9=4.0*pow(pvii->TF1,2);
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

    DL = gamma(spvii->RL-0.5)/(sqrt( pow(2.0,(1.0/spvii->RL))-1.0)*gamma(spvii->RL));
    DH = gamma(spvii->RH-0.5)/(sqrt( pow(2.0,(1.0/spvii->RH))-1.0)*gamma(spvii->RH));
    DLD = gamma(spvii->RL-0.5+0.001)/(sqrt(pow(2.0,(1.0/(spvii->RL+0.001)))-1.0)*gamma(spvii->RL+0.001));
    DHD = gamma(spvii->RH-0.5+0.001)/(sqrt( pow(2.0,(1.0/(spvii->RH+0.001))) -1.0)*gamma(spvii->RH+0.001));
    spvii->DA1L = ( pow(((1.0+A)/A),2.0))*( pow(2.0,(1.0/spvii->RL))-1.0);
    spvii->DA1H = (pow((1.0+A),2.0))*(pow(2.0,(1/spvii->RH))-1.0);
    spvii->DA2L = 1.128379167*(1.0+A)*(1.0/(A*DL+DH))/W;
    spvii->DA2H = spvii->DA2L;
    spvii->DA3L = 2.0*spvii->RL*spvii->DA1L/W;
    spvii->DA3H = 2.0*spvii->RH*spvii->DA1H/W;
    spvii->DA4L = (1.0/(1.0+A))-(DL/(A*DL+DH));
    spvii->DA4H = spvii->DA4L;
    spvii->DA5L = 2.0*spvii->RL*pow(((1+A)/A),2.0)*(1.0+A)/(pow(A,3));
    spvii->DA5H = -2.0*spvii->RH*pow((1+A),2.0)*(1+A);
    spvii->DA6L = -A*1000.0*(DLD-DL)/(A*DL+DH);
    spvii->DA6H = -1000.0*(DHD-DH)/(A*DL+DH);
    spvii->DA7L = log(2.0)*pow(((1.0+A)/A),2.0)*(pow(2.0,(1.0/spvii->RL))-1.0+1.0)/spvii->RL;
    spvii->DA7H = log(2.0)*pow((1.0+A),2.0)*(pow(2.0,(1.0/spvii->RH))-1.0+1.0)/spvii->RH;
}

double DBWS::PROFIL(int N, double X)
{
    double r;

    switch (N) {
    case _Gaussian:
        r=0.939437279*exp(-2.772588722*X);
        prfx->PRFDER=2.772588722;
        break;
    case _Lorentzian:
        r=0.636619772/(1.0+4.0*X);
        prfx->PRFDER=4.0/(1.+4.0*X);
        break;
    case _Mod1:
        r=0.819449653/ pow((1.0+1.656854248*X),2.0);
        prfx->PRFDER=3.313708496/(1.0+1.656854248*X);
        break;
    case _Mod2:
        r=0.766420937/pow((1.0+2.349604208*X),1.5);
        prfx->PRFDER=3.5244063/(1.0+2.349604208*X);
        break;
    case _SplitPearsonVII:
        if (prfx->DELTA < 0.0)
        {
            r = pow((1.0+spvii->DA1L*X),(-spvii->RL)) *spvii->DA2L;
            prfx->PRFDER = spvii->DA1L*spvii->RL/(1.0+spvii->DA1L*X);
        }
        else
        {
            r = pow((1.0+spvii->DA1H*X),(-spvii->RH)) *spvii->DA2H;
            prfx->PRFDER = spvii->DA1H*spvii->RH/(1.0+spvii->DA1H*X);
        }
        break;
    case _pseudoVoigt:
        r=prfx->GAM1*0.636619772/(1.0+4.0*X)+(1-prfx->GAM1)*0.939437279*exp(-2.772588722*X);
        prfx->PRFDER=(prfx->GAM1*2.546479088/ pow((1.0+4.0*X),2) +(1-prfx->GAM1)*2.6046732048*exp(-2.772588722*X))/r;
        break;
    case _PearsonVII:
        r=pvii->C4/pow((1.0+4.0*( pow(2,(1.0/prfx->GAM1)) -1)*X) , prfx->GAM1);
        prfx->PRFDER=pvii->TF6*prfx->GAM1/(1.0+4.0*( pow(2,(1.0/prfx->GAM1))-1)*X);
        break;
    case _TCHZ:
        r=prfx->GAM1*0.636619772/(1.0+4.0*X)+(1-prfx->GAM1)*0.939437279*exp(-2.772588722*X);
        prfx->PRFDER=(prfx->GAM1*2.546479088/ pow((1.0+4.0*X),2) +(1-prfx->GAM1)*2.6046732048*exp(-2.772588722*X))/r;
        break;
    default:
        DBWSException("ILLEGAL PROFILE FUNCTION REQUEST");
    }
    return r;
}

void DBWS::SUMMAT(int IPM, double CSK[], double DISK[], double DYCDD[], double ISODER[], double TOTCS)
{
    int I, J, K, M,II, IL, KK, LK, KM,LK1, LK2,IOF,KRP1, IISO;
    double X, Z, X1, BB,YX,DER, TLL, TLR,ASS5,ESSE,OMEGA, YCALC,
        TANNJ, SHIFT, OMEGA8, LAMDAM,DERMON,PRTEMP;
    bool VERT;
    double DERIV[MSZ+1];

    SHIFT = 0.0;
    for (J=1; J <= MSZ; ++J) DERIV[J]=0.0;
    YCALC=0.0;
    IL=0;
    //-----CALCULATE THE CONTRIBUTION OF THE REFLECTIONS ORD1 TO ORD2 TO THE
    //-----DERIVATIVES W.R.T. THE PROFILE INTENSITY YOBS
    if( g3->IORD1 == 0) goto L12;
    for (I=g3->IORD1; I <= g3->IORD2; ++I)
    {
        prfx->IPH=refls->IREFS[I]/(256*256*256*8);
        IL=IL+1;
        J=(I % NOV)+1;
        // test for asymmetry function !cp ap 12 97
        if (cntrls->IASYM == 0)
        {
            VERT = refls->REFS[I][2] <= g1->RLIM;
        }
        else
        {
            VERT = abs(refls->REFS[I][2]-90.0) >= g1->RLIM;
        }
        SHIFT = params->GLB_[10-1] * cos(refls->REFS[I][2]/2.0/57.2958) + params->GLB_[11-1] * sin(refls->REFS[I][2]/57.2958);
        prfx->DELTA=g3->TH-refls->REFS[I][2]-params->GLB_[1-1]-SHIFT;
        //TANTH=tan(g3->TH*3.14159265359/360.0);
        prfx->DELT=prfx->DELTA*prfx->DELTA;
        prfx->TL=refls->REFS[I][1];
        if (cntrls->NPROF == _pseudoVoigt) prfx->GAM1 = phases[prfx->IPH-1].PAR[16] + phases[prfx->IPH-1].PAR[17] * refls->REFS[I][2];
        if (cntrls->NPROF == _PearsonVII) prfx->GAM1 = phases[prfx->IPH-1].PAR[16] + phases[prfx->IPH-1].PAR[17] / refls->REFS[I][2] + phases[prfx->IPH-1].PAR[18]/refls->REFS[I][2]/refls->REFS[I][2];
        if (cntrls->NPROF == _PearsonVII) PRSVII(prfx->GAM1);
        if (cntrls->NPROF == _SplitPearsonVII)
        {
            spvii->RL=phases[prfx->IPH-1].PAR[16]+(phases[prfx->IPH-1].PAR[17]+phases[prfx->IPH-1].PAR[18]/refls->REFS[I][2])/refls->REFS[I][2];
            spvii->RH=phases[prfx->IPH-1].PAR[23]+(phases[prfx->IPH-1].PAR[24]+phases[prfx->IPH-1].PAR[25]/refls->REFS[I][2])/refls->REFS[I][2];
            mspvii(phases[prfx->IPH-1].PAR[26],prfx->TL);
        }
        if(cntrls->NPROF == _TCHZ)
        {
            TLL = refls->HALFL[I];
            prfx->GAM1 = refls->GAM[I];
            TLR = TLL/prfx->TL;
        }
        BB=prfx->TL*prfx->TL;
        //-----NEXT LINE IS NECESSEARY FOR 2 PHASES WITH VERY DifFERENT FWHM.
        if (prfx->DELT/BB > g1->WDT*g1->WDT) goto L33;
        if (VERT)
        {
            //       test for asymmetry model               !cp may 01 97
            if (cntrls->IASYM == 0)
            {
                YX=prfx->DELT*sign(prfx->DELTA);
                Z=1.0-phases[prfx->IPH-1].PAR[13]*YX/g4->TANN[J];
                if ( Z <= 0.0 ) Z=0.0001;
            }
            else
            {
                YX=sign(prfx->DELTA)*prfx->DELTA/(2*prfx->TL);
                TANNJ=g4->TANN[J];
                if (TANNJ >= 1.0) TANNJ=tan(atan(TANNJ)-3.14159265359/2);
                Z=(phases[prfx->IPH-1].PAR[13]/TANNJ) * (2.0*(prfx->DELTA/(2*prfx->TL))*exp(-YX));
                Z=1+Z;
                if ( Z <= 0.0 ) Z=0.0001;
            }
        }
        else
        {
            Z=1.0;
        }
        //L5:
        PRTEMP = PROFIL(cntrls->NPROF,prfx->DELT/BB);
        if (cntrls->NPROF == _SplitPearsonVII)
        {
            OMEGA = refls->REFS[I][3]*Z*PRTEMP*phases[prfx->IPH-1].PAR[0];
        }
        else
        {
            OMEGA = refls->REFS[I][3]*Z*PRTEMP*phases[prfx->IPH-1].PAR[0]/prfx->TL;
        }
        YCALC = YCALC+OMEGA*refls->FMGNTD[I];
        if ( cntrls->JOBTYP > 2 ) goto L33;
        X = prfx->PRFDER*2.0* prfx->DELT/BB-1.0;
        for (K=1; K <= cntrls->MAXS; ++K)
        {
            DER=1.0;
            //-----Broadening Coeficients Derivatives
            if(cntrls->NPROF != _SplitPearsonVII)
            {
                for (M=3; M <= 5; ++M) if(phases[prfx->IPH-1].PAR[M-1].L == K) DER=X/prfx->TL/2.0;
            }
            if (phases[prfx->IPH-1].PAR[19].L == K) DER=X/prfx->TL/2.0;
            X1=0.0;
            //-----Asymmetry Derivative
            if (VERT)
            {
                if (cntrls->IASYM == 0)
                {
                    X1=phases[prfx->IPH-1].PAR[13]*sign(prfx->DELTA)*BB/g4->TANN[J]/Z;
                }
                else
                {
                    X1=-phases[prfx->IPH-1].PAR[13]*exp(-YX)*(prfx->TL/(2*prfx->DELTA)-sign(prfx->DELTA)*1.0/4)/TANNJ/Z;
                }
            }
            //-----Zero, Displacement, and Transparancy Derivative
            if ( params->GLB_[1-1].L == K ) DER=prfx->DELTA*(prfx->PRFDER+X1);
            if ( params->GLB_[10-1].L == K ) DER=prfx->DELTA*(prfx->PRFDER+X1);
            if ( params->GLB_[11-1].L == K ) DER=prfx->DELTA*(prfx->PRFDER+X1);
            if ( (phases[prfx->IPH-1].PAR[13].L == K)  &&  VERT )
            {
                if ( cntrls->IASYM == 0 )
                {
                    DER = YX/Z;
                }
                else
                {
                    DER = -2.0*(prfx->DELTA/(2.0*prfx->TL))*exp(-YX)/Z;
                }
            }
            if ( cntrls->NPROF == _TCHZ ) goto L8;
            //-----Pseudo-Voigt Shape Derivatives
            if(cntrls->NPROF == _pseudoVoigt)
            {
                KRP1=phases[prfx->IPH-1].PAR[16].L;
                if (K == KRP1) DER=(0.636619772/(1.0+4.0*prfx->DELT/BB)-0.939437279*exp(-2.772588722*prfx->DELT/BB))/PRTEMP;
                KRP1=phases[prfx->IPH-1].PAR[17].L;
                if (K == KRP1) DER=(0.636619772/(1.0+4.0*prfx->DELT/BB)-0.939437279*exp(-2.772588722*prfx->DELT/BB))/PRTEMP;
            }
            //-----Pearson VII Shape Derivatives
            if (cntrls->NPROF == _PearsonVII)
            {
                KRP1=phases[prfx->IPH-1].PAR[16].L;
                if(K == KRP1) DER=-log(1.0+pvii->TF9*prfx->DELT/BB)+pvii->TF4*(prfx->DELT/BB)/(1.0+pvii->TF9*prfx->DELT/BB)+pvii->TF8;
                KRP1=phases[prfx->IPH-1].PAR[17].L;
                if(K == KRP1)   DER=-log(1.0+pvii->TF9*prfx->DELT/BB)+pvii->TF4*(prfx->DELT/BB)/(1.0+pvii->TF9*prfx->DELT/BB)+pvii->TF8;
                KRP1=phases[prfx->IPH-1].PAR[18].L;
                if(K == KRP1) DER=-log(1.0+pvii->TF9*prfx->DELT/BB)+pvii->TF4*(prfx->DELT/BB)/(1.0+pvii->TF9*prfx->DELT/BB)+pvii->TF8;
            }
            //-----Lattice Parameter Derivatives
L8:
            for (M=6; M <= 11; ++M) if(phases[prfx->IPH-1].PAR[M-1].L == K) DER=(prfx->PRFDER+X1)*prfx->DELTA;
            DERIV[K]=g4->DERSTO[J][K]*DER*OMEGA+DERIV[K];
            //L3:;
        }
        //----TCHZ Profile Derivatives
        if(cntrls->NPROF == _TCHZ)
        {
            OMEGA8 = Z*phases[prfx->IPH-1].PAR[0]*refls->REFS[I][3]/prfx->TL;
            for (K = 1; K <= cntrls->MAXS; ++K)
            {
                for (M=3; M <= 5; ++M)
                {
                    if (phases[prfx->IPH-1].PAR[M-1].L != K) goto L1001;
                    DERIV[K] = DERIV[K]+ OMEGA8*g4->DERSTO[J][K]/2.0*(0.939437279*exp(-2.772588722*prfx->DELT/BB) - 0.636619772/(1.0+4.0*prfx->DELT/BB)) * (1.36603*TLR/prfx->TL-0.95438*TLR*TLR/prfx->TL+0.33348 * pow(TLR,3.0)/prfx->TL);
L1001:;
                }
                if (phases[prfx->IPH-1].PAR[19].L == K) DERIV[K] = DERIV[K]+ OMEGA8*g4->DERSTO[J][K]/2.0*(0.939437279*exp(-2.772588722*prfx->DELT/BB) - 0.636619772/(1.0+4.0*prfx->DELT/BB)) * (1.36603*TLR/prfx->TL-0.95438*TLR*TLR/prfx->TL+0.33348*pow(TLR,3.0)/prfx->TL);
                if (phases[prfx->IPH-1].PAR[14].L != K) goto L1002;
                DERIV[K] = DERIV[K]+ OMEGA8*(0.939437279*exp(-2.772588722*prfx->DELT/BB) -0.636619772/(1.0+4.0*prfx->DELT/BB)) *((1.36603*TLR/prfx->TL-0.95438*TLR*TLR/prfx->TL+0.33348*pow(TLR,3.0)/prfx->TL)* g4->DERSTO[J][K]/2.0 - refls->FMGNTD[I]*g4->TANN[J]*(1.36603/prfx->TL-0.95438*TLR/prfx->TL+0.33348*TLR*TLR/prfx->TL));
L1002:
                if (phases[prfx->IPH-1].PAR[15].L != K) goto L1003;
                DERIV[K] = DERIV[K]+ OMEGA8*(0.939437279*exp(-2.772588722*prfx->DELT/BB) -0.636619772/(1.0+4.0*prfx->DELT/BB)) *((1.36603*TLR/prfx->TL-0.95438*TLR*TLR/prfx->TL+0.33348*pow(TLR,3.0)/prfx->TL)* g4->DERSTO[J][K]/2.0 - refls->FMGNTD[I]*sqrt(1+g4->TANN[J]*g4->TANN[J])*(1.36603/prfx->TL-0.95438*TLR/prfx->TL+0.33348*TLR*TLR/prfx->TL));
L1003:;
            }
        }
L33:;
    }
    //-----FORM SUMS
L12:
    if(jnk->NBCKGD != 0)goto L11;
    for (II=2; II <= 7; ++II)
    {
        if(params->GLB_[II-1].L == 0)goto L10;
        KM=params->GLB_[II-1].L;
        if(II == 2)DERIV[KM]=DERIV[KM]+1.0;
        if(II == 2)goto L10;
        DERIV[KM]=DERIV[KM]+ pow( ((g1->THMIN+static_cast<double>(IPM-1)*g1->STEP)/g1->BKPOS-1.0) , (II-2));
L10:;
    }
L11:
    datax->YC[IPM]=YCALC;
    if (cntrls->JOBTYP > 2) goto L20;
    for (K = 1; K <= cntrls->NPHASE; ++K)
    {
        LK1  = phases[K-1].PAR[0].L;
        LK2  = phases[K-1].PAR[1].L;
        //
        //-----UPDATING GLOBAL SCALE DERIVATE FOR BKG CONTRIBUTE
        //
        if(cntrls->FONDO == 1 || cntrls->FONDO == 2)
        {
            if(LK1 != 0) DERIV[LK1]=DERIV[LK1]+volume->GCOM[K]*CSK[K]+volume->GCOM[K]*DISK[K];
        }
        //
        //-----UPDATING DERIVATE OF Q OVERALL FOR BKG CONTRIBUTE
        //
        if(cntrls->FONDO == 2)
        {
            if(LK2 != 0) DERIV[LK2] = DERIV[LK2] +  DYCDD[K];
        }
        //
        //-----UPDATING DERIVATE OF ISOTROPIC THERMAL PARAMETERS FOR BKG CONTRIBUTE
        //
        if(cntrls->FONDO == 1)
        {
            IOF = 0;
            if(K > 1)
            {
                for (I = 2; I <= K; ++I) IOF = IOF + phases[I-1].AtomCount;
            }
            for (I = 1; I <= phases[K].AtomCount; ++I)
            {
                IISO=params->LP[I+IOF][4];
                if(IISO != 0) DERIV[IISO] = DERIV[IISO] + ISODER[I+IOF];
                ISODER[I+IOF]=0.0;
            }
        }
    }
    //-----DYC RESPECT TO AMORPHOUS SCALE FACTOR
    LK = params->GLB_[20-1].L;
    if(LK != 0) DERIV[LK] = DERIV[LK] + datax->AMORPHOUS[IPM];
    //-----MONOCHROMATOR PARAMETERS DERIVATIVES
    LAMDAM=(g1->LAMDA[1]*params->RATIO[1]+g1->LAMDA[2]*params->RATIO[2])/(params->RATIO[1]+params->RATIO[2]);
    ESSE=2*sin((g3->TH-params->GLB_[1-1]-SHIFT)*0.00872665)/LAMDAM;
    LK=params->GLB_[18-1].L;
    if ( LK != 0 )
    {
        //-------NEXT LINE FOR A LORENTZIAN MONOCHROMATOR BASS-BAND  FUNCTION
        ASS5=1/(1+params->GLB_[18-1]*  pow(ESSE,params->GLB_[19-1]) );
        DERMON=-( pow(ESSE,params->GLB_[19-1]) )/ ( pow((1+params->GLB_[18-1]*  pow(ESSE,params->GLB_[19-1]) ),2) );
        //--------------------------------------------------------------------
        DERIV[LK]=DERIV[LK]+TOTCS/ASS5*DERMON;
    }
    //   !cp ap 20 97
    LK=params->GLB_[19-1].L;
    ASS5=1/(1+params->GLB_[18-1]*  pow(ESSE,params->GLB_[19-1]) );
    if ( LK != 0 )
    {
        DERMON=-params->GLB_[18-1]*log(ESSE)*( pow(ESSE,params->GLB_[19-1]) )/ pow( (1+params->GLB_[18-1]*(  pow(ESSE,params->GLB_[19-1]) )) , 2);
        DERIV[LK]=DERIV[LK]+TOTCS/ASS5*DERMON;
    }
    //-----FORM THE UPPER TRIANGULAR OF "RJAC" = MATRIX OF NORMAL EQUATIONS
    prfx->DELTA = datax->Y[IPM]-datax->BK[IPM]-YCALC;
    for (J=1; J <= cntrls->MAXS; ++J)
    {
        X = DERIV[J]/(datax->VAR[IPM]);
        f1->VX[J] = f1->VX[J]+prfx->DELTA*X;
        for (KK=J; KK <= cntrls->MAXS; ++KK) f1->RJAC[J][KK] = f1->RJAC[J][KK]+X*DERIV[KK];
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
        E[I+3]=E[I+3]/(2.0*sqrt(S[I+3]*(1.0-S[I+3])))*180.0/3.14159265359;
        S[I+3]=atan(sqrt((1.0-S[I+3])/S[I+3]))*180.0/3.14159265359;
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
        V[J]=phases[*IPH-1].PAR[J+5-1];
        K=phases[*IPH-1].PAR[J+5-1].L;
        X=phases[*IPH-1].PAR[J+5-1].codeword;
        if(K == 0)X=0.0;
        for (L=J; L <= 6; ++L)
        {
            M=phases[*IPH-1].PAR[L+5-1].L;
            SM[L][J]=0.0;
            if((M == 0) || (K == 0))goto L1;
            SM[L][J]=f1->RJAC[M][K]*X*phases[*IPH-1].PAR[L+5-1].codeword;
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
        << "R-P        = " << setw(8) << setprecision(2) << g2->R2 << "%" << endl
        << "R-WP       = " << setw(8) << setprecision(2) << g2->R3 << "%     "
        << "R-WP(Background Removed) = " << setw(8) << setprecision(2) << g2->R3NOBK << "%" << endl;
    if ( cntrls->MAXS == 0  &&  cntrls->MAXSX != 0 )
    {
        X=100.0*sqrt((static_cast<double>(g3->NUM)-static_cast<double>(cntrls->MAXSX))*1.0/g2->D2);
    }
    else
    {
        X=100.0*sqrt((static_cast<double>(g3->NUM)-static_cast<double>(cntrls->MAXS))*1.0/g2->D2);
    }
    if ( cntrls->MAXS == 0  &&  cntrls->MCYCLE == 1 )
    {
        OUTSCR(cntrls->MCYCLX,g2->R2,g2->R3,X);
    }
    else
    {
        OUTSCR(ICYCLE-1,g2->R2,g2->R3,X);
    }
    file6
        << "R-EXPECTED = " << setw(8) << setprecision(2) << X << "%" << endl
        << "S          = " << setw(8) << setprecision(2) << sqrt(g2->S2/static_cast<double>(g3->NUM-cntrls->MAXS))
        << "     SQRT(RESIDUAL/N-P)GOODNESS OF FIT" << endl
        << "D - W D    = " << setw(8) << setprecision(2) << g2->SS4/g2->SS2 << "     UNWEIGHTED DURBIN-WATSON STATISTIC D" << endl;
    I=g3->NUM-cntrls->MAXS;
    file6
        << "N-P        = " << setw(8) << I << endl

        << endl
        << scientific
        << "SUMYDIF     = " << setw(14) << setprecision(8) << g2->S1 << endl
        << "SUMYOBS     = " << setw(14) << setprecision(8) << g2->D1 << endl
        << "SUMYCALC    = " << setw(14) << setprecision(8) << g2->S3 << endl
        << "SUMWYOBSSQ  = " << setw(14) << setprecision(8) << g2->D2 << endl
        << "RESIDUAL    = " << setw(14) << setprecision(8) << g2->S2 << endl
        << "CONDITION   = " << setw(14) << setprecision(8) << g3->COND << endl
        << fixed << endl;

    DUMMY[2*MSZ+1] = g2->R2;
    DUMMY[2*MSZ+2] = g2->R3;
    DUMMY[2*MSZ+3] = sqrt(g2->S2/static_cast<double>(g3->NUM));
    DUMMY[2*MSZ+4] = g2->SS4/g2->SS2;

    //     FINAL PARAMETERS AND R-FACTORS
    if (cntrls->IPLST != 0 && cntrls->MAXS == 0)
    {
        for (IP=1; IP <= cntrls->NPHASE; ++IP)
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
                    KM=params->LP[I+IOF][J];
                    if(KM != 0)
                    {
                        //  !cp jun 96 start
                        if (J ==  5)
                        {
                            DUMMY[KM] = params->XL[I+IOF][J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                        }
                        else
                        {
                            // !cp jun 96 stop
                            DUMMY[KM] = params->XL[I+IOF][J];
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
                    KM=params->LP[I+IOF][J];
                    if(KM != 0)
                    {
                        //  !cp jun 96 start
                        if (J ==  5)
                        {
                            DUMMY[KM] = params->XL[I+IOF][J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                        }
                        // !cp jun 96 stop
                        DUMMY[KM]  =params->XL[I+IOF][J];
                    }
                }
            }
            for (J=1; J <= 27; ++J)
            {
                KM=phases[IP-1].PAR[J-1].L;
                if (KM != 0) DUMMY[KM] = phases[IP-1].PAR[J-1];
            }
            DIRECT(dircv->DCSM,dircv->DCV,&IP);
            for (I=1; I <= 6; ++I)
            {
                KM = phases[IP-1].PAR[I+5-1].L;
                MATCH = 0;
                if (KM != 0)
                {
                    if (I >= 4)
                    {
                        for (MMM=1; MMM <= 3; ++MMM) if (KM == phases[IP-1].PAR[MMM+5-1].L) MATCH=1;
                    }
                    if (MATCH == 0) DUMMY[KM] = dc->SAVE[IP][I];
                }
            }
        }
        for (J=1; J <= 20; ++J)
        {
            if(cntrls->NPROF == _TCHZ && (J == 14 || J == 15 || J == 16)) goto L6019;
            if(cntrls->NPROF == _PearsonVII && (J == 14 || J == 15 || J == 16)) goto L6019;
            if(cntrls->NPROF == _pseudoVoigt && (J == 14 || J == 15 || J == 16)) goto L6019;
            if(cntrls->NPROF == _SplitPearsonVII && (J == 14 || J == 15 || J == 16)) goto L6019;
            KM=params->GLB_[J-1].L;
            if(KM != 0) DUMMY[KM] = params->GLB_[J-1];
L6019:;
        }
        for (I=1; I <= cntrls->MAXSX; ++I) file8o << DUMMY[I];
        for (I=MSZ+1; I <= MSZ+cntrls->MAXSX; ++I) file8o << DUMMY[I];
        for (I=2*MSZ+1; I <= 2*MSZ+4; ++I) file8o << DUMMY[I];
        for (I=1; I <= 4; ++I) if ((I % 2) == 1) allp->ILOC = allp->ILOC + 1;
        //L9487:
        allp->FINAL[allp->ILOC][2-(I % 2)] = DUMMY[2*MSZ+I];
    }
    if (cntrls->MAXS == 0 && cntrls->MCYCLE == 1) return;
    allp->ILOC = 0;
    for (I=1; I <= NFINAL; ++I)
    {
        for (J=1; J <= 2; ++J) allp->FINAL[I][J] = 0.0;
    }
    file6 << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CYCLE NUMBER=" << setw(4) << ICYCLE << endl;
    for (IP=1; IP <= cntrls->NPHASE; ++IP)
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
        multip->TMASSA[IP]=0.0;
        STMASSA[IP]=0.0;
        for (I=1; I <= N; ++I)
        {
            for (J=1; J <= 5; ++J)
            {
                allp->ILOC = allp->ILOC + 1;
                KM=params->LP[I+IOF][J];
                if(KM != 0)goto L6;
                SY[J]=0.0;
                SZ[J]=0.0;
                //  !cp jun 97 start
                if(J == 5)
                {
                    allp->FINAL[allp->ILOC][1] = params->XL[I+IOF][J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                }
                else
                {
                    allp->FINAL[allp->ILOC][1] = params->XL[I+IOF][J];
                }
                //  !cp jun 97 stop
                goto L5;
L6:
                SZ[J]=sqrt(abs(f1->RJAC[KM][KM]));
                SY[J]=f1->VX[KM]*params->A[I+IOF][J]*params->RELAX[1];
                //  !cp jun 96 start
                if(J == 5)
                {
                    DUMMY[KM] = params->XL[I+IOF][J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                    params->XL[I+IOF][J]=params->XL[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ]  = SY[J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                    allp->FINAL[allp->ILOC][1] = params->XL[I+IOF][J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                    allp->FINAL[allp->ILOC][2] = SZ[J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                }
                else
                {
                    // !cp jun 96 stop
                    DUMMY[KM]  =params->XL[I+IOF][J];
                    params->XL[I+IOF][J]=params->XL[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ] = SY[J];
                    allp->FINAL[allp->ILOC][1] = params->XL[I+IOF][J];
                    allp->FINAL[allp->ILOC][2] = SZ[J];
                    SMASS[IP]=SZ[5];
                }
L5:;
            }
            // !cp jun 96 start
            //L4:
            file6
                << setw(4) << parac->ATEXT[I+IOF] << "  "
                << setw(9) << setprecision(6) << params->XL[I+IOF][1]
                << setw(9) << setprecision(6) << SY[1]
                << setw(9) << setprecision(6) << SZ[1]
                << setw(9) << setprecision(6) << params->XL[I+IOF][2]
                << setw(9) << setprecision(6) << SY[2]
                << setw(9) << setprecision(6) << SZ[2]
                << setw(9) << setprecision(6) << params->XL[I+IOF][3]
                << setw(9) << setprecision(6) << SY[3]
                << setw(9) << setprecision(6) << SZ[3]
                << setw(9) << setprecision(6) << params->XL[I+IOF][4]
                << setw(9) << setprecision(6) << SY[4]
                << setw(9) << setprecision(6) << SZ[4]
                << setw(9) << setprecision(6) << params->XL[I+IOF][5]*multip->XMLTP[IP]/multip->MURT[I+IOF]
                << setw(9) << setprecision(6) << SY[5]*multip->XMLTP[IP]/multip->MURT[I+IOF]
                << setw(9) << setprecision(6) << SZ[5]*multip->XMLTP[IP]/multip->MURT[I+IOF] << endl;
        }
        file6
            << "ATOM       B11      DB11      SB11       B22      DB22      SB22       B33      DB33      SB33" << endl
            << "           B12      DB12      SB12       B13      DB13      SB13       B23      DB23      SB23" << endl;
        for (I=1; I <= N; ++I)
        {
            for (J=6; J <= 11; ++J)
            {
                allp->ILOC = allp->ILOC + 1;
                KM=params->LP[I+IOF][J];
                if(KM != 0)goto L9;
                SY[J]=0.0;
                SZ[J]=0.0;
                allp->FINAL[allp->ILOC][1] = params->XL[I+IOF][J];
                goto L8;
L9:
                SZ[J]=sqrt(abs(f1->RJAC[KM][KM]));
                SY[J]=f1->VX[KM]*params->A[I+IOF][J]*params->RELAX[2];
                //  !cp jun 96 start
                if(J ==  5)
                {
                    DUMMY[KM] = params->XL[I+IOF][J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                    params->XL[I+IOF][J]=params->XL[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ] = SY[J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                    allp->FINAL[allp->ILOC][1] = params->XL[I+IOF][J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                    allp->FINAL[allp->ILOC][2] = SZ[J]*multip->XMLTP[IP]/multip->MURT[I+IOF];
                }
                else
                {
                    DUMMY[KM]  =params->XL[I+IOF][J];
                    params->XL[I+IOF][J]=params->XL[I+IOF][J]+SY[J];
                    DUMMY[KM+MSZ] = SY[J];
                    allp->FINAL[allp->ILOC][1] = params->XL[I+IOF][J];
                    allp->FINAL[allp->ILOC][2] = SZ[J];
                }
                // !cp jun 96 stop
L8:;
            }
            //L7:
            file6 << setw(4) << parac->ATEXT[I+IOF]
                << setw(10) << setprecision(6) << params->XL[I+IOF][6]
                << setw(10) << setprecision(6) << SY[6]
                << setw(10) << setprecision(6) << SZ[6]
                << setw(10) << setprecision(6) << params->XL[I+IOF][7]
                << setw(10) << setprecision(6) << SY[7]
                << setw(10) << setprecision(6) << SZ[7]
                << setw(10) << setprecision(6) << params->XL[I+IOF][8]
                << setw(10) << setprecision(6) << SY[8]
                << setw(10) << setprecision(6) << SZ[8] << endl
                << "    "
                << setw(10) << setprecision(6) << params->XL[I+IOF][9]
                << setw(10) << setprecision(6) << SY[9]
                << setw(10) << setprecision(6) << SZ[9]
                << setw(10) << setprecision(6) << params->XL[I+IOF][10]
                << setw(10) << setprecision(6) << SY[10]
                << setw(10) << setprecision(6) << SZ[10]
                << setw(10) << setprecision(6) << params->XL[I+IOF][11]
                << setw(10) << setprecision(6) << SY[11]
                << setw(10) << setprecision(6) << SZ[11] << endl;
        }
        for (J=1; J <= 27; ++J)
        {
            if ( J == 12 )
            {
                allp->ILOC = allp->ILOC + 1;
                allp->FINAL[allp->ILOC][1] = phases[IP].PREF[1];
                allp->FINAL[allp->ILOC][2] = phases[IP].PREF[2];
                allp->ILOC = allp->ILOC + 1;
                allp->FINAL[allp->ILOC][1] = phases[IP].PREF[3];
            }
            allp->ILOC = allp->ILOC + 1;
            if ( J == 6 ) ILOC1 = allp->ILOC;
            KM = phases[IP-1].PAR[J-1].L;
            if ( KM == 0 )
            {
                SY[J] = 0.0;
                SZ[J] = 0.0;
                allp->FINAL[allp->ILOC][1] = phases[IP-1].PAR[J-1];
            }
            else
            {
                SZ[J] = sqrt(abs(f1->RJAC[KM][KM]));
                SY[J] = f1->VX[KM]*phases[IP-1].PAR[J-1].codeword*params->RELAX[3];
                DUMMY[KM] = phases[IP-1].PAR[J-1];
                phases[IP-1].PAR[J-1] = phases[IP-1].PAR[J-1]+SY[J];
                DUMMY[KM+MSZ] = SY[J];
                allp->FINAL[allp->ILOC][1] = phases[IP-1].PAR[J-1];
                allp->FINAL[allp->ILOC][2] = SZ[J];
            }
        }
        DIRECT(dircv->DCSM,dircv->DCV,&IP);
        ILOC2 = allp->ILOC;
        allp->ILOC = ILOC1-1;
        for (I=1; I <= 6; ++I)
        {
            allp->ILOC = allp->ILOC + 1;
            KM = phases[IP-1].PAR[I+5-1].L;
            MATCH=0;
            SY[I+5]=dircv->DCV[I]-dc->SAVE[IP][I];
            if (KM != 0)
            {
                if (I >= 4)
                {
                    for (MMM=1; MMM <= 3; ++MMM) if (KM == phases[IP-1].PAR[MMM+5-1].L) MATCH=1;
                }
                if (MATCH == 0)
                {
                    DUMMY[KM] = dc->SAVE[IP][I];
                    DUMMY[KM+MSZ] = SY[I+5];
                }
            }
            allp->FINAL[allp->ILOC][1] = dircv->DCV[I];
            allp->FINAL[allp->ILOC][2] = dircv->DCSM[I][I];
            dc->SAVE[IP][I]=dircv->DCV[I];
        }
        allp->ILOC = ILOC2;
        if(dircv->DCV[1] == dircv->DCV[2] && dircv->DCV[2] != dircv->DCV[3])
        {
            if(dircv->DCV[4] == dircv->DCV[5] && dircv->DCV[5] != 90.0)
            {
                if(dircv->DCV[6] == 90.0 || dircv->DCV[6] == 120.0)
                {
                    RATIODCV=dircv->DCV[3]/dircv->DCV[1];
                    SRATIODCV=(dircv->DCV[1]*dircv->DCSM[3][3]+dircv->DCV[3]*dircv->DCSM[1][1])/(dircv->DCV[3] * dircv->DCV[3]);
                    file6
                        << "OVERALL SCALE FACTOR=" << scientific
                        << setw(9) << setprecision(3) << phases[IP-1].PAR[0]
                        << setw(9) << setprecision(3) << SY[1]
                        << setw(9) << setprecision(3) << SZ[1] << endl << fixed
                        << "OVERALL TEMP. FACTOR="
                        << setw(9) << setprecision(4) << phases[IP-1].PAR[1]
                        << setw(9) << setprecision(4) << SY[2]
                        << setw(9) << setprecision(4) << SZ[2] << endl
                        << "CELL PARAMETERS="
                        << setw(11) << setprecision(6) << dircv->DCV[1]
                        << setw(11) << setprecision(6) << SY[1+5]
                        << setw(11) << setprecision(6) << dircv->DCSM[1][1] << endl
                        << setw(11) << setprecision(6) << dircv->DCV[2]
                        << setw(11) << setprecision(6) << SY[2+5]
                        << setw(11) << setprecision(6) << dircv->DCSM[2][2] << endl
                        << setw(11) << setprecision(6) << dircv->DCV[3]
                        << setw(11) << setprecision(6) << SY[3+5]
                        << setw(11) << setprecision(6) << dircv->DCSM[3][3] << endl
                        << setw(11) << setprecision(4) << dircv->DCV[4]
                        << setw(11) << setprecision(4) << SY[4+5]
                        << setw(11) << setprecision(4) << dircv->DCSM[4][4] << endl
                        << setw(11) << setprecision(4) << dircv->DCV[5]
                        << setw(11) << setprecision(4) << SY[5+5]
                        << setw(11) << setprecision(4) << dircv->DCSM[5][5] << endl
                        << setw(11) << setprecision(4) << dircv->DCV[6]
                        << setw(11) << setprecision(4) << SY[6+5]
                        << setw(11) << setprecision(4) << dircv->DCSM[6][6] << endl << endl
                        << "            c/a= "  << setw(11) << setprecision(6) << RATIODCV << "  +/-" << setw(11) << setprecision(6) << SRATIODCV << endl << endl
                        << "PREFERRED ORIENTATION PARAMETERS="
                        << " "
                        << setw(8) << setprecision(5) << phases[IP-1].PAR[11]
                        << setw(8) << setprecision(5) << SY[12]
                        << setw(8) << setprecision(5) << SZ[12]
                        << " "
                        << setw(8) << setprecision(5) << phases[IP-1].PAR[12]
                        << setw(8) << setprecision(5) << SY[13]
                        << setw(8) << setprecision(5) << SZ[13] << endl
                        << "ASYMMETRY PARAMETER="
                        << " "
                        << setw(8) << setprecision(4) << phases[IP-1].PAR[13]
                        << setw(8) << setprecision(4) << SY[14]
                        << setw(8) << setprecision(4) << SZ[14] << endl
                        << "LORENTZIAN HALF WIDTH PARAMS (X AND Y) "
                        << " "
                        << setw(8) << setprecision(5) << phases[IP-1].PAR[14]
                        << setw(8) << setprecision(5) << SY[15]
                        << setw(8) << setprecision(5) << SZ[15] << endl
                        << " "
                        << setw(8) << setprecision(5) << phases[IP-1].PAR[15]
                        << setw(8) << setprecision(5) << SY[16]
                        << setw(8) << setprecision(5) << SZ[16] << endl;
                    goto L1032;
                }
            }
        }
        file6 << "OVERALL SCALE FACTOR=" << scientific
            << setw(9) << setprecision(3) << phases[IP-1].PAR[0]
        << setw(9) << setprecision(3) << SY[1]
        << setw(9) << setprecision(3) << SZ[1] << fixed << endl
            << "OVERALL TEMP. FACTOR="
            << setw(9) << setprecision(4) << phases[IP-1].PAR[1]
        << setw(9) << setprecision(4) << SY[2]
        << setw(9) << setprecision(4) << SZ[2] << endl
            << "CELL PARAMETERS=" << endl
            << setw(11) << setprecision(6) << dircv->DCV[1]
        << setw(11) << setprecision(6) << SY[1+5]
        << setw(11) << setprecision(6) << dircv->DCSM[1][1] << endl
            << setw(11) << setprecision(6) << dircv->DCV[2]
        << setw(11) << setprecision(6) << SY[2+5]
        << setw(11) << setprecision(6) << dircv->DCSM[2][2] << endl
            << setw(11) << setprecision(6) << dircv->DCV[3]
        << setw(11) << setprecision(6) << SY[3+5]
        << setw(11) << setprecision(6) << dircv->DCSM[3][3] << endl
            << setw(11) << setprecision(4) << dircv->DCV[4]
        << setw(11) << setprecision(4) << SY[4+5]
        << setw(11) << setprecision(4) << dircv->DCSM[4][4] << endl
            << setw(11) << setprecision(4) << dircv->DCV[5]
        << setw(11) << setprecision(4) << SY[5+5]
        << setw(11) << setprecision(4) << dircv->DCSM[5][5] << endl
            << setw(11) << setprecision(4) << dircv->DCV[6]
        << setw(11) << setprecision(4) << SY[6+5]
        << setw(11) << setprecision(4) << dircv->DCSM[6][6] << endl
            << "PREFERRED ORIENTATION PARAMETERS=" << endl
            << " "
            << setw(8) << setprecision(5) << dircv->DCV[1]
        << setw(8) << setprecision(5) << SY[1+5]
        << setw(8) << setprecision(5) << dircv->DCSM[1][1] << endl
            << " "
            << setw(8) << setprecision(5) << dircv->DCV[2]
        << setw(8) << setprecision(5) << SY[2+5]
        << setw(8) << setprecision(5) << dircv->DCSM[2][2] << endl
            << " "
            << setw(8) << setprecision(5) << dircv->DCV[3]
        << setw(8) << setprecision(5) << SY[3+5]
        << setw(8) << setprecision(5) << dircv->DCSM[3][3] << endl
            << "ASYMMETRY PARAMETER=" << endl
            << " "
            << setw(8) << setprecision(4) << dircv->DCV[4]
        << setw(8) << setprecision(4) << SY[4+5]
        << setw(8) << setprecision(4) << dircv->DCSM[4][4] << endl
            << "LORENTZIAN HALF WIDTH PARAMS (X AND Y)" << endl
            << " "
            << setw(8) << setprecision(5) << dircv->DCV[5]
        << setw(8) << setprecision(5) << SY[5+5]
        << setw(8) << setprecision(5) << dircv->DCSM[5][5] << endl
            << " "
            << setw(8) << setprecision(5) << dircv->DCV[6]
        << setw(8) << setprecision(5) << SY[6+5]
        << setw(8) << setprecision(5) << dircv->DCSM[6][6] << endl
            << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP-1].PAR[11]
        << setw(8) << setprecision(5) << SY[12]
        << setw(8) << setprecision(5) << SZ[12] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP-1].PAR[12]
        << setw(8) << setprecision(5) << SY[13]
        << setw(8) << setprecision(5) << SZ[13] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP-1].PAR[13]
        << setw(8) << setprecision(5) << SY[14]
        << setw(8) << setprecision(5) << SZ[14] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP-1].PAR[14]
        << setw(8) << setprecision(5) << SY[15]
        << setw(8) << setprecision(5) << SZ[15] << endl
            << " "
            << setw(8) << setprecision(5) << phases[IP-1].PAR[15]
        << setw(8) << setprecision(5) << SY[16]
        << setw(8) << setprecision(5) << SZ[16] << endl;
L1032:
        //-----CHECK FOR THE SPLIT PEARSON VII PROFILE
        if(cntrls->NPROF == _SplitPearsonVII)
        {
            file6 << "LOW SIDE EXPONENT PARAMETERS (NA, NB, NC)=" << endl
                //<< scientific
                << " "
                << setw(10) << setprecision(4) << phases[IP-1].PAR[16]
                << setw(10) << setprecision(4) << SY[17]
                << setw(10) << setprecision(4) << SZ[17] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP-1].PAR[17]
                << setw(10) << setprecision(4) << SY[18]
                << setw(10) << setprecision(4) << SZ[18] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP-1].PAR[18]
                << setw(10) << setprecision(4) << SY[19]
                << setw(10) << setprecision(4) << SZ[19] << endl
                << "HIGH SIDE EXPONENT PARAMETERS (NA, NB, NC)=" << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP-1].PAR[23]
                << setw(10) << setprecision(4) << SY[24]
                << setw(10) << setprecision(4) << SZ[24] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP-1].PAR[24]
                << setw(10) << setprecision(4) << SY[25]
                << setw(10) << setprecision(4) << SZ[25] << endl
                << " "
                << setw(10) << setprecision(4) << phases[IP-1].PAR[25]
                << setw(10) << setprecision(4) << SY[26]
                << setw(10) << setprecision(4) << SZ[26] << endl
                << fixed
                << "SPLIT PEARSON VII ASSYMETRY PARAMETER="
                << " "
                << setw(10) << setprecision(6) << phases[IP-1].PAR[26]
                << setw(10) << setprecision(6) << SY[27]
                << setw(10) << setprecision(6) << SZ[27] << endl;
        }
        else
        {
            //-----if NOT A SPLIT PEARSON VII PROFILE
            file6 << "MIXING PARAMETERS" << endl
                << "NA= "
                //<< scientific
                << setw(10) << setprecision(3) << phases[IP-1].PAR[16]
                << setw(10) << setprecision(3) << SY[17]
                << setw(10) << setprecision(3) << SZ[17] << endl
                << "NB= "
                << setw(10) << setprecision(3) << phases[IP-1].PAR[17]
                << setw(10) << setprecision(3) << SY[18]
                << setw(10) << setprecision(3) << SZ[18] << endl
                << "NC= "
                << setw(10) << setprecision(3) << phases[IP-1].PAR[18]
                << setw(10) << setprecision(3) << SY[19]
                << setw(10) << setprecision(3) << SZ[19] << fixed << endl;
        }
        file6 << "FWHM PARAMETERS=" << endl
            << "U = "
            << setw(10) << setprecision(6) << phases[IP-1].PAR[2]
            << setw(10) << setprecision(6) << SY[3]
            << setw(10) << setprecision(6) << SZ[3] << endl
            << "V = "
            << setw(10) << setprecision(6) << phases[IP-1].PAR[3]
            << setw(10) << setprecision(6) << SY[4]
            << setw(10) << setprecision(6) << SZ[4] << endl
            << "W = "
            << setw(10) << setprecision(6) << phases[IP-1].PAR[4]
            << setw(10) << setprecision(6) << SY[5]
            << setw(10) << setprecision(6) << SZ[5] << endl
            << "CT= "
            << setw(10) << setprecision(6) << phases[IP-1].PAR[20]
            << setw(10) << setprecision(6) << SY[21]
            << setw(10) << setprecision(6) << SZ[21] << endl
            << "Z = "
            << setw(10) << setprecision(6) << phases[IP-1].PAR[19]
            << setw(10) << setprecision(6) << SY[20]
            << setw(10) << setprecision(6) << SZ[20] << endl;

        //-----Modification introduced by Carlos O. Paiva-Santos to perform
        //-----Quantitative phase analysis, 03/94. Added 05/94, T.S. Moss
        //-----CHANGES TO INCORPORATE THE REFINED OCCUPANCY. Paiva-Santos (Feb-Mar/95)
        for (I=1; I <= N; ++I)
        {
            ICOCO=params->PTR[I+IOF];
            multip->TMASSA[IP] = multip->TMASSA[IP] + params->XL[I+IOF][5]*coeff->XMAS[ICOCO]*multip->XMLTP[IP];
            STMASSA[IP] = STMASSA[IP] + SMASS[IP]*coeff->XMAS[ICOCO]*multip->XMLTP[IP];
        }
        XFAC = 3.141592654 / 180.000000;
        dircv->DCV[4] = XFAC * dircv->DCV[4];
        dircv->DCSM[4][4] = dircv->DCSM[4][4] * XFAC;
        dircv->DCV[5] = XFAC * dircv->DCV[5];
        dircv->DCSM[5][5] = dircv->DCSM[5][5] * XFAC;
        dircv->DCV[6] = dircv->DCV[6] * XFAC;
        dircv->DCSM[6][6] = dircv->DCSM[6][6] * XFAC;
        //-----Calculations of VOLUME and SVZM (=W) for each phase
        //-----and respectives standard deviations
        //-----New standard deviation code introduced in nov 96 !cp
        ARGCOS= 1-pow((cos(dircv->DCV[4])),2)-pow((cos(dircv->DCV[5])),2)-pow((cos(dircv->DCV[6])),2) + 2 * (cos(dircv->DCV[4])) * (cos(dircv->DCV[5])) * (cos(dircv->DCV[6]));
        V0 = dircv->DCV[1] * dircv->DCV[2] * dircv->DCV[3];
        VOL[IP] = V0 * sqrt(ARGCOS);
        VOSQ = 0.5*VOL[IP]/ARGCOS;
        ARG1 = VOSQ*(2 * cos(dircv->DCV[4]) * sin(dircv->DCV[4]) - 2*sin(dircv->DCV[4]) *cos(dircv->DCV[5]) *cos(dircv->DCV[6])) * dircv->DCSM[4][4];
        ARG2 = VOSQ*(2 * cos(dircv->DCV[5]) * sin(dircv->DCV[5]) - 2*sin(dircv->DCV[5]) *cos(dircv->DCV[4]) *cos(dircv->DCV[6])) * dircv->DCSM[5][5];
        ARG3 = VOSQ*(2 * cos(dircv->DCV[6]) * sin(dircv->DCV[6]) - 2*sin(dircv->DCV[6]) *cos(dircv->DCV[4]) *cos(dircv->DCV[5])) * dircv->DCSM[6][6];
        DVOL[IP] = sqrt(pow((VOL[IP] * dircv->DCSM[1][1] / dircv->DCV[1]),2) + pow((VOL[IP] * dircv->DCSM[2][2] / dircv->DCV[2]),2) + pow((VOL[IP] * dircv->DCSM[3][3] / dircv->DCV[3]),2) + pow(ARG1,2) + pow(ARG2,2) + pow(ARG3,2));
        // standard deviations are calculed below                      !cp nov 96
        W[IP] = phases[IP-1].PAR[0] * multip->TMASSA[IP] * VOL[IP]/phases[IP].NMOL;
        DW[IP] = (SZ[1]/phases[IP-1].PAR[0]) + (DVOL[IP]/VOL[IP]) + (STMASSA[IP]/multip->TMASSA[IP])/phases[IP].SAQF;
        //   end of std
        file6 << "Volume= "
            << setw(9) << setprecision(3) << VOL[IP]
        << "(+/-)"
            << setw(7) << setprecision(3) << DVOL[IP]
        << " UCW= "
            << setw(7) << setprecision(2) << multip->TMASSA[IP]
        << " U.C.Density = "
            << setw(7) << setprecision(3) << 1.66113*multip->TMASSA[IP]/VOL[IP]
        << " gr/cm^3" << endl
            << "                     _____________________________" << endl;
        //L10:;
    }
    // ****** QUANTITATIVE ANALYSIS ***************
    WTOTAL = 0.000000;
    //DWTOTAL = 0.000000;
    for (I = 1; I <= cntrls->NPHASE; ++I) WTOTAL = WTOTAL + W[I];
    for (I = 1; I <= cntrls->NPHASE; ++I)
    {
        //      !cp nov 10 96
        if (cntrls->NPHASE == 1)
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
    for (I = 1; I <= cntrls->NPHASE; ++I)
    {
        if (IINNMOL == 1) goto L2713;
        if (phases[I].NMOL == 0) IINNMOL=1;
L2713:;
    }
    if (IINNMOL == 1)
    {
        for (I = 1; I <= cntrls->NPHASE; ++I)
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
        for (I = 1; I <= cntrls->NPHASE; ++I)
        {
            FRP[I] = XMASS[I] * phases[I].NMOL / multip->TMASSA[I];
            FT = FT + FRP[I];
        }
        for (I = 1; I <= cntrls->NPHASE; ++I)
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
    if (cntrls->ISPHASE != 0)
    {
        file6 << "Considering Amorphous Content:" << endl;
        SFIS=phases[cntrls->ISPHASE].WTIS/XMASS[cntrls->ISPHASE];
        if(SFIS > 1.0)
        {
            file6 << "PROBLEM:Amount of Internal Standard (Phase #"
                << setw(2) << cntrls->ISPHASE
                << ") is less than the specified "
                << setw(6) << setprecision(2) << phases[cntrls->ISPHASE].WTIS
            << "%." << endl
                << "Amorphous content not computed. Check ISWT in line 11.2 for this phase" << endl;
            goto L2720;
        }
        for (I=1; I <= cntrls->NPHASE; ++I)
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
        allp->ILOC = allp->ILOC +1;
        KM=params->GLB_[J-1].L;
        if(KM != 0)goto L20;
        SY[J]=0.0;
        SZ[J]=0.0;
        allp->FINAL[allp->ILOC][1] = params->GLB_[J-1];
        goto L19;
L20:
        SZ[J]=sqrt(abs(f1->RJAC[KM][KM]));
        SY[J]=f1->VX[KM]*params->GLB_[J-1].codeword*params->RELAX[4];
        DUMMY[KM] = params->GLB_[J-1];
        params->GLB_[J-1]=params->GLB_[J-1]+SY[J];
        DUMMY[KM+MSZ]  = SY[J];
        allp->FINAL[allp->ILOC][1] = params->GLB_[J-1];
        allp->FINAL[allp->ILOC][2] = SZ[J];
L19:;
    }
    if (allp->ILOC  >  NFINAL) DBWSException("Parameter NFINAL in PARAM.INC file too small");
    file6 << "GLOBAL PARAMETERS" << endl
        << "ZEROPOINT (ZER)             :"
        << setw(8) << setprecision(4) << params->GLB_[1-1]
    << setw(8) << setprecision(4) << SY[1]
    << setw(8) << setprecision(4) << SZ[1] << endl;
    file6 << "SAMPLE DISPLACEMENT (DISP)  :"
        << setw(8) << setprecision(4) << params->GLB_[10-1]
    << setw(8) << setprecision(4) << SY[10]
    << setw(8) << setprecision(4) << SZ[10] << endl
        << "SAMPLE TRANSPARENCY (TRANSP):"
        << setw(8) << setprecision(4) << params->GLB_[11-1]
    << setw(8) << setprecision(4) << SY[11]
    << setw(8) << setprecision(4) << SZ[11] << endl
        << "ROUGHNESS PARAMETERS        :"
        << "             P              :"
        << setw(8) << setprecision(4) << params->GLB_[8-1]
    << setw(8) << setprecision(4) << SY[8]
    << setw(8) << setprecision(4) << SZ[8] << endl
        << "             Q              :"
        << setw(8) << setprecision(4) << params->GLB_[9-1]
    << setw(8) << setprecision(4) << SY[9]
    << setw(8) << setprecision(4) << SZ[9] << endl
        << "             R              :"
        << setw(8) << setprecision(4) << params->GLB_[12-1]
    << setw(8) << setprecision(4) << SY[12]
    << setw(8) << setprecision(4) << SZ[12] << endl
        << "             T              :"
        << setw(8) << setprecision(4) << params->GLB_[13-1]
    << setw(8) << setprecision(4) << SY[13]
    << setw(8) << setprecision(4) << SZ[13] << endl;
    // !cp ap 20 97  !from It. codes
    file6 << "AMORPHOUS SCALE (SCAM):"
        << setw(11) << setprecision(4) << params->GLB_[20-1]
    << setw(11) << setprecision(4) << SY[20]
    << setw(11) << setprecision(4) << SZ[20] << endl;
    // !cp ap 20 97   !from It. codes
    file6 << "MONOCROMATOR BANDPASS PARAMETERS (PMONI)" << endl
        << setw(8) << setprecision(4) << params->GLB_[18-1]
    << setw(8) << setprecision(4) << SY[18]
    << setw(8) << setprecision(4) << SZ[18]
    << setw(8) << setprecision(4) << params->GLB_[19-1]
    << setw(8) << setprecision(4) << SY[19]
    << setw(8) << setprecision(4) << SZ[19] << endl;
    if (jnk->NBCKGD != 0) goto L90;
    file6 << "BACKGROUND PARAMETERS" << endl //<< scientific
        << setw(12) << setprecision(6) << params->GLB_[2-1]
    << setw(12) << setprecision(6) << SY[2]
    << setw(12) << setprecision(6) << SZ[2] << endl
        << setw(12) << setprecision(6) << params->GLB_[3-1]
    << setw(12) << setprecision(6) << SY[3]
    << setw(12) << setprecision(6) << SZ[3] << endl
        << setw(12) << setprecision(6) << params->GLB_[4-1]
    << setw(12) << setprecision(6) << SY[4]
    << setw(12) << setprecision(6) << SZ[4] << endl
        << setw(12) << setprecision(6) << params->GLB_[5-1]
    << setw(12) << setprecision(6) << SY[5]
    << setw(12) << setprecision(6) << SZ[5] << endl
        << setw(12) << setprecision(6) << params->GLB_[6-1]
    << setw(12) << setprecision(6) << SY[6]
    << setw(12) << setprecision(6) << SZ[6] << endl
        << setw(12) << setprecision(6) << params->GLB_[7-1]
    << setw(12) << setprecision(6) << SY[7]
    << setw(12) << setprecision(6) << SZ[7] << endl;
L90:
    for (I=1; I <= cntrls->MAXS; ++I) file8o << DUMMY[I];
    for (I=MSZ+1; I <= MSZ+cntrls->MAXS; ++I) file8o << DUMMY[I];
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
    double CSK[99+1];
    double DISK[99+1];
    double DYCDD[99+1];
    double VLAST[MSZ+1];
    double DERISO[NATS+1];
    double ISODER[NATS+1];


    for (I=1; I <= MSZ; ++I) VLAST[I]=0.0;

    //     **** start cycles ****
    for (IX=1; IX <= cntrls->MCYCLE; ++ IX)
    {
        ASSIGN_();
        for (I=1; I <= MSZ; ++I)
        {
            f1->VX[I]=0.0;
            for (J=1; J <= MSZ; ++J) f1->RJAC[I][J]=0.0;
        }
        for (I=1; I <= NATS; ++I)
        {
            if (cntrls->IBGD == 1)
            {
                ISODER[I] = 1;
            }
            else
            {
                ISODER[I]=0;
            }
        }
        g3->NUM=0;
        IPREV=1;
        if(jnk->NBCKGD != 0)goto L30;
        g3->TH=g1->THMIN-g1->STEP;
        for (I=1; I <= datax->NPTS; ++I)
        {
            g3->TH=g3->TH+g1->STEP;
            THX=g3->TH/g1->BKPOS-1.0;
            datax->BK[I]=params->GLB_[2-1];
            for (J=2; J <= 6; ++J) datax->BK[I]=datax->BK[I]+params->GLB_[J+1-1] * pow(THX,J-1);
            if(cntrls->MCYCLE == 1 && cntrls->MAXS == 0 && cntrls->IPLPOL == 1) fondi->BKPOL[I]=datax->BK[I];
        }
L30:
        g3->TH=g1->THMIN-g1->STEP;

        //           ****** START GREAT LOOP ON NPTS STEP-POINTS ******
        for (IPM=1; IPM <= datax->NPTS; ++IPM)
        {
            if ( (IPM % (datax->NPTS/80+1)) == 0) cout << ".";
            g3->TH=g3->TH+g1->STEP;
            // !cp jun 97 test for bck option
            if (cntrls->IBGD == 1)
            {
                CSK[1]    = 1;
                DISK[1]   = 1;
                DYCDD[1]  = 1;
                TOTCS  = 1;
                goto L4334;
            }
            //-----COMPUTE FPOL = POLARIZATION FACTOR FOR THIS IPM POINT
            //     TH  = 2.0 * THETA
            //     STH = SIN(THETA)
            STH = sin(g3->TH * 0.008726646);
            FPOL = (1.0 + pow(1.0 - 2.0 * STH * STH,2) * g1->CTHM ) / 2.0;
            if(datax->KR[IPM] == IRS)
            {
                // ----IDERIV = 1 MEANS ONLY BACKGROUND (FOR THE STEP-POINTS THAT
                //                                       ARE IN EXCLUDED REGIONS)
                IDERIV = 1;
            }
            else
            {
                //-----IDERIV = 2 MEANS BACKGROUND + DERIVATIVES (FOR STEP-POINTS THAT ARE
                //                                                IN INCLUDED REGIONS)
                IDERIV = 2;
            }
            //-----if NECESSARY COMPUTE TOTCS :
            //     TOTCS = TOTAL COMPTON INTENSITY SCATTERED BY ALL CRYSTALLINE
            //             PHASES AT THE IPM-TH POINT OF THE X-RAY PATTERN.
            TOTCS = 0.0;
            if ( cntrls->FONDO == 1  ||  cntrls->FONDO == 2 )
            {
                for (K = 1; K <= cntrls->NPHASE; ++K)
                {
                    COMPTON(K,STH,&CISK);
                    CSK[K] = CISK * FPOL;
                    TOTCS = TOTCS + bkgscale->SCABKG[K] * CSK[K];
                }
                //-----COMPTON UPDATE
                if ( cntrls->MCYCLE == 1  &&  cntrls->MAXS == 0  &&  cntrls->IPLCOM == 1 )   fondi->BKCOM[IPM] = TOTCS;
                datax->BK[IPM] = datax->BK[IPM] + TOTCS;
            }
            //-----if NECESSARY CALL DIS = SUBROUTINE TO COMPUTE THERMAL AND LATTICE
            //                             DISORDER SCATTERING SDK  AND  DERIVATIVES
            //                             DYC IN THE K-TH PHASE
            //     TOTDS = TOTAL DISORDER SCATTERING
            TOTDS = 0.0;
            if (cntrls->FONDO == 1 || cntrls->FONDO == 2)
            {
                for (K = 1; K <= cntrls->NPHASE; ++K)
                {
                    IOF = 0;
                    if(K > 1)
                    {
                        for (I = 2; I <= K; ++I) IOF = IOF + phases[I-1].AtomCount;
                    }
                    for (I=1; I <= NATS; ++I) DERISO[I]=0;
                    DISORDER(K,STH,IDERIV,&SDK,&DYC,cntrls->FONDO,DERISO);
                    //int K, double STH, int IDERIV, double* SDK, double* DYC, int FONDO, double DERISO[])
                    DISK[K]  =   SDK * FPOL;
                    TOTDS    = TOTDS + bkgscale->SCABKG[K] * DISK[K];
                    //     UPDATING DERIVATE OF ISOTROPIC THERMAL FACTORS
                    if(cntrls->FONDO == 1)
                    {
                        for (I=1; I <= phases[K].AtomCount; ++I) ISODER[IOF+I]= ISODER[IOF+I]+bkgscale->SCABKG[K]*FPOL*DERISO[IOF+I];
                    }
                    //     UPDATING DERIVATE OF OVERALL THERMAL FACTOR
                    if(cntrls->FONDO == 2) DYCDD[K]=DYC*bkgscale->SCABKG[K]*FPOL;
                }
                //-----DISORDER UPDATE
                if(cntrls->MCYCLE == 1  &&  cntrls->MAXS == 0  &&  cntrls->IPLDIS == 1 ) fondi->BKDIS[IPM] = TOTDS;
                datax->BK[IPM] = datax->BK[IPM] + TOTDS;
            }
            //------AMORPHOUS EVALUATIONS
            TOTAS =0.0;
            TOTAS = params->GLB_[20-1] * datax->AMORPHOUS[IPM];
            //-----AMORPHOUS UPDATE
            if ( cntrls->MCYCLE == 1  &&  cntrls->MAXS == 0  &&  cntrls->IPLAM == 1 ) fondi->BKAM[IPM] = TOTAS;
            datax->BK[IPM] = datax->BK[IPM] + TOTAS;
L4334:
            if(datax->KR[IPM] == IRS )goto L9;
            g3->IORD1=datax->KR[IPM] % IRS;
            g3->IORD2=(datax->KR[IPM]/IRS) % IRS;
            if ( (g3->IORD2 == 0  ||  g3->IORD1 == 0)  &&  jnk->NBCKGD != 0 ) goto L9;
            //L15:
            g3->NUM=g3->NUM+1;
            IORDLIM2=g3->IORD2;
            if(IPREV > IORDLIM2) goto L40;
            for (J=IPREV; J <= IORDLIM2; ++J) CALCUL(J);
L40:
            IPREV=max(IPREV,IORDLIM2+1);
            //L108:
            SUMMAT(IPM,CSK,DISK,DYCDD,ISODER,TOTCS);
L9:;
        }
        cout << ":" << endl;
        if ( cntrls->JOBTYP > 2 ) return;
        for (I=1; I <= cntrls->MAXS; ++I)
        {
            for (J=I; J <= cntrls->MAXS; ++J) f1->RJAC[J][I] = f1->RJAC[I][J];
        }
        g3->COND = DPINV(f1->RJAC,f1->VX,&cntrls->MAXS);
        CHISQ();
        OOM=g2->S2 / static_cast<double>(g3->NUM-cntrls->MAXS);
        for (I=1; I <= cntrls->MAXS; ++I)
        {
            for (J=1; J <= cntrls->MAXS; ++J) f1->RJAC[I][J] = f1->RJAC[I][J]*OOM;
        }
        //              CODE TO ATTEMPT TO STABILIZE OSCILLATIONS
        for (I=1; I <= cntrls->MAXS; ++I)
        {
            if( sign(f1->VX[I]) == sign(VLAST[I])) goto L60;
            if(abs(f1->VX[I]) > 1.2*abs(VLAST[I])) goto L60;
            if(abs(f1->VX[I]) < 0.8*abs(VLAST[I])) goto L60;
            f1->VX[I]=f1->VX[I]/2.0;
L60:
            VLAST[I]=f1->VX[I];
        }
        OUTPTR(IX);
        for (I=1; I <= cntrls->MAXS; ++I)
        {
            X1=sqrt(abs(f1->RJAC[I][I]))*cntrls->EPS;
            if(abs(f1->VX[I]) > X1)goto L10;
        }
        if (cntrls->MAXS > 0)
        {
            file6 << endl
                << "          ***** EPSED OUT *****" << endl;
        }
        if (cntrls->MAXS > 0) cntrls->ICYRUN = IX;
        goto L20;
L10:;
    }
    //                 ***** END CYCLES *****
L20:
    //     CODE FOR PRINTING CORRELATION MATRIX MOVED FROM EXPUT (65-85)
    if(cntrls->MAT == 0 || cntrls->MAXS == 0) return;
    file6 << endl << endl << "CORRELATION MATRIX=" << endl;
    IA=1;
    LIM=19;
L38:
    LIM=min(cntrls->MAXS,LIM);
    for (I=IA; I <= LIM; ++I) file6 << setw(6) << I;
    file6 << endl;
    for (I=1; I <= cntrls->MAXS; ++I)
    {
        L=0;
        X=f1->RJAC[I][I];
        for (J=IA; J <= LIM; ++J)
        {
            L=L+1;
            YY=f1->RJAC[J][J]*X;
            KS[L]= static_cast<int>( 100.0*f1->RJAC[I][J]/(sqrt(abs(YY)) * sign(YY) )+0.5 );
        }
        file6 << setw(5) << I;
        for (J=1; J <= L; ++J) file6 << setw(6) << KS[J];
        file6 << endl;
    }
    if(LIM >= cntrls->MAXS) goto L37;
    IA=LIM+1;
    LIM=LIM+19;
    goto L38;
L37:;
}

void DBWS::CEL000(int* MULTX_, double Y_[][3+1])
{
    //
    //                             THIS SR IS USED TO FORCE A POINT XYZ TO
    //                           LIE IN THE UNIT CELL AT 000.
    //

    Y_[*MULTX_-1][1] = fmod(Y_[*MULTX_-1][1]+8.0,1.0);
    if ( Y_[*MULTX_-1][1]-0.99999 < 0 ) goto L2; else goto L1;
L1:
    Y_[*MULTX_-1][1] = 0.0;
L2:
    Y_[*MULTX_-1][2] = fmod(Y_[*MULTX_-1][2]+8.0,1.0);
    if ( Y_[*MULTX_-1][2]-0.99999 < 0 ) goto L4; else goto L3;
L3:
    Y_[*MULTX_-1][2] = 0.0;
L4:
    Y_[*MULTX_-1][3] = fmod(Y_[*MULTX_-1][3]+8.0,1.0);
    if ( Y_[*MULTX_-1][3]-0.99999 < 0 ) goto L6; else goto L5;
L5:
    Y_[*MULTX_-1][3] = 0.0;
L6:;
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
        goto L100;
        break;
    case 2:
        goto L200;
        break;
    case 3:
        goto L300;
        break;
    }
    GOTOER();
L100:
    switch (I1) {
    case 1:
        goto L101;
        break;
    case 2:
        goto L102;
        break;
    case 3:
        goto L103;
        break;
    case 4:
        goto L104;
        break;
    }
    GOTOER();
L101:
    Y_[*MULTX_][1] = Y_[*L2][1];
    goto L105;
L102:
    Y_[*MULTX_][1] = Y_[*L2][2];
    goto L105;
L103:
    Y_[*MULTX_][1] = -Y_[*L2][1];
    goto L105;
L104:
    Y_[*MULTX_][1] = -Y_[*L2][2];
L105:
    switch (I2) {
    case 1:
        goto L110;
        break;
    case 2:
        goto L107;
        break;
    case 3:
        goto L108;
        break;
    case 4:
        goto L106;
        break;
    }
    GOTOER();
L106:
    Y_[*MULTX_][1] = Y_[*MULTX_][1]+0.25;
L107:
    Y_[*MULTX_][1] = Y_[*MULTX_][1]+0.25;
L108:
    Y_[*MULTX_][1] = Y_[*MULTX_][1]+0.25;
L110:
    switch (I3) {
    case 1:
        goto L111;
        break;
    case 2:
        goto L112;
        break;
    case 3:
        goto L113;
        break;
    case 4:
        goto L114;
        break;
    }
    GOTOER();
L111:
    Y_[*MULTX_][2] = Y_[*L2][2];
    goto L115;
L112:
    Y_[*MULTX_][2] = Y_[*L2][1];
    goto L115;
L113:
    Y_[*MULTX_][2] = -Y_[*L2][2];
    goto L115;
L114:
    Y_[*MULTX_][2] = -Y_[*L2][1];
L115:
    switch (I4) {
    case 1:
        goto L120;
        break;
    case 2:
        goto L117;
        break;
    case 3:
        goto L118;
        break;
    case 4:
        goto L116;
        break;
    }
    GOTOER();
L116:
    Y_[*MULTX_][2] = Y_[*MULTX_][2]+0.25;
L117:
    Y_[*MULTX_][2] = Y_[*MULTX_][2]+0.25;
L118:
    Y_[*MULTX_][2] = Y_[*MULTX_][2]+0.25;
L120:
    Y_[*MULTX_][3] = Y_[*L2][3] * static_cast<double>(1-2*I5);
    switch (I6) {
    case 1:
        goto L125;
        break;
    case 2:
        goto L122;
        break;
    case 3:
        goto L123;
        break;
    case 4:
        goto L121;
        break;
    }
    GOTOER();
L121:
    Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.25;
L122:
    Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.25;
L123:
    Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.25;
L125:
    switch (I7) {
    case 1:
        goto L130;
        break;
    case 2:
        goto L127;
        break;
    case 3:
        goto L126;
        break;
    }
    GOTOER();
L126:
    Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.333333333;
L127:
    Y_[*MULTX_][3] = Y_[*MULTX_][3]+0.333333333;
L130:
    *MULTX_ = *MULTX_+1;
    CEL000(MULTX_,Y_);
    return;
L200:
    Y_[*MULTX_][1] = -Y_[*L2][2];
    Y_[*MULTX_][2] = Y_[*L2][1]-Y_[*L2][2];
    Y_[*MULTX_][3] = Y_[*L2][3]+0.333333333 * static_cast<double>(I2-2);
    *MULTX_ = *MULTX_+1;
    CEL000(MULTX_,Y_);
    Y_[*MULTX_][1] = -Y_[*MULTX_-1][2];
    Y_[*MULTX_][2] = -Y_[*L2][1];
    Y_[*MULTX_][3] = Y_[*L2][3]-0.333333333 * static_cast<double>(I2-2);
    goto L130;
L300:
    I1 = *L2;
L301:
    Y_[*MULTX_][1] = Y_[I1][3];
    Y_[*MULTX_][2] = Y_[I1][1];
    Y_[*MULTX_][3] = Y_[I1][2];
    if ( I1-*L2 <= 0 ) goto L302; else goto L130;
L302:
    I1 = *MULTX_;
    *MULTX_ = *MULTX_+1;
    CEL000(MULTX_,Y_);
    goto L301;
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
    //L130:
L200:
    *IXB_ = max( (*IXB_ & 18874367), (*IXB_ & 16777215));
    if (  (*IXB_ & 384) == 384 ) *IXB_ |= 512;
}

void DBWS::RTMT(int* MULTX_, double Y_[][3+1], int* IPRT, int NCTR_[],  int* IPHASE)
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
    multip->MLTPHASE=phases[*IPHASE].SYMB.MULTP;
    for (I=1; I <= phases[*IPHASE].SYMB.MULTP; ++I)
    {
        phases[*IPHASE].IVEC[I] = 0;
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
            phases[*IPHASE].IVEC[I] = phases[*IPHASE].IVEC[I]+(
              9*(static_cast<int>(VCTR[K][1]+1.5))+
              3*(static_cast<int>(VCTR[K][2]+1.5))+
              static_cast<int>(VCTR[K][3]+1.5))*32768 * static_cast<int>(pow(32,3-K))+
              (static_cast<int>(VCTR[K][4]*12.0+.5)+16) * static_cast<int>(pow(32,3-K));
        }
        if ( *IPRT > 0 )
        {
            if(simoper->ISIMOP != 1)
            {
                if ( (I-1 % 15) == 0) file6 << "The operations of the space group are" << endl;
                file6
                    << " ("	<< setw(3) << setprecision(0) << VCTR[1][1]	<< setw(3) << setprecision(0) << VCTR[1][2]	<< setw(3) << setprecision(0) << VCTR[1][3]	<< " )   ( X )   ("	<< setw(6) << setprecision(3) << VCTR[1][4] << " )   ( X2 )" << endl
                    << " ("	<< setw(3) << setprecision(0) << VCTR[2][1]	<< setw(3) << setprecision(0) << VCTR[2][2]	<< setw(3) << setprecision(0) << VCTR[2][3]	<< " ) * ( Y ) + (" << setw(6) << setprecision(3) << VCTR[2][4] << " ) = ( Y2 )" << endl
                    << " ("	<< setw(3) << setprecision(0) << VCTR[3][1]	<< setw(3) << setprecision(0) << VCTR[3][2]	<< setw(3) << setprecision(0) << VCTR[3][3]	<< " )   ( Z )   (" << setw(6) << setprecision(3) << VCTR[3][4] << " )   ( Z2 )"
                    << "          VEC(" << setw(3) << I << ") =" << setw(15) << phases[*IPHASE].IVEC[I] << endl << endl;
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
    for (I=1; I <= 7; ++I) hklctl->NC1[1][I][*IPHASE] = 0;
    NCTR_[1] = 0;
    for (I=1; I <= 8; ++I)
    {
        if ( phases[*IPHASE].SYMB.NCONT[I] <= 0 ) goto L101;
        if ( phases[*IPHASE].SYMB.NCONT[I] < 8192 ) goto L80;
        if ( phases[*IPHASE].SYMB.NCONT[I] <= 8192 ) goto L60;		// TODO: Possivel erro aqui estava 08192
        if ( phases[*IPHASE].SYMB.NCONT[I] == 16384 ) goto L70;
        phases[*IPHASE].SYMB.NPOL = 4;
        hklctl->NC1[I][1][*IPHASE] = 2;
L50:
        for (J=2; J <= 7; ++J) hklctl->NC1[I][J][*IPHASE] = 0;
        hklctl->NC1[I][6][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/4) % 32;
        goto L90;
L60:
        hklctl->NC1[I][1][*IPHASE] = 3;
        phases[*IPHASE].SYMB.NPOL = 1;
        goto L50;
L70:
        hklctl->NC1[I][1][*IPHASE] = 1;
        goto L50;
L80:
        hklctl->NC1[I][1][*IPHASE] = 0;
        hklctl->NC1[I][2][*IPHASE] = phases[*IPHASE].SYMB.NCONT[I] % 4;
        hklctl->NC1[I][3][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/4) % 4;
        hklctl->NC1[I][4][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/16) % 4;
        hklctl->NC1[I][5][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/64) % 4;
        hklctl->NC1[I][6][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/256) % 2;
        hklctl->NC1[I][7][*IPHASE] = (phases[*IPHASE].SYMB.NCONT[I]/512) % 16;
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
    hklctl->N1HKL[*IPHASE] = I-1;
    if ( hklctl->N1HKL[*IPHASE] < 2 ) goto L105;
    for (N2=2; N2 <= hklctl->N1HKL[*IPHASE]; ++N2) if ( hklctl->NC1[N2-1][1][*IPHASE] == 1 ) goto L103;
    goto L105;
L103:
    for (I=N2; I <= hklctl->N1HKL[*IPHASE]; ++I)
    {
        for (J=1; J <= 7; ++J) hklctl->NC1[I-1][J][*IPHASE] = hklctl->NC1[I][J][*IPHASE];
    }
    hklctl->NC1[hklctl->N1HKL[*IPHASE]][1][*IPHASE] = 1;
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

void DBWS::LOOKUP(int K, int N, int NSCAT, int IXRAY, int JOB)
{
    int I,J,L,NS,IOF,NSL,IIPHAS,TBXPTR;

    IOF=0;
    NS=0;
    if (K > 1)
    {
        for (IIPHAS=2; IIPHAS <= K; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
    }
    if (K > 1) NS=dc->NSAVE;
    if (JOB == 2 || JOB == 4) goto L70;
    for (I=1; I <= N; ++I)
    {
        NS=max(NSCAT,NS);
        for (J=1; J <= NS; ++J) if (parac->NTYP[I+IOF] == coefc->NAM[J]) goto L30;
        for (J=1; J <= 212; ++J) if (parac->NTYP[I+IOF] == TBXC[J]) goto L50;
        file6 << " SCATTERING COEFFICIENTS NOT FOUND FOR " << parac->NTYP[I+IOF] << endl;
        //41      FORMAT(' SCATTERING COEFFICIENTS NOT FOUND FOR ',A4)
        DBWSException("SCATTERING DATA MISSING");
L30:
        params->PTR[I+IOF]=J;
        goto L10;
L50:
        NS=NS+1;
        params->PTR[I+IOF]=NS;
        for (L=1; L <= 9; ++ L) coeff->AC[L][NS]=TBX[J][L];
        coefc->NAM[NS]=TBXC[J];
        TBXPTR=static_cast<int>(TBX[J][10]);			//TBXPTR=static_cast<int>(TBX[J][10]+0.5);
        coeff->DFP[NS]=TBD[TBXPTR][IXRAY];
        coeff->DFPP[NS]=TBD[TBXPTR][IXRAY+10];
        coeff->XMAS[NS]=TBM[TBXPTR];
L10:;
    }
    dc->NSAVE=NS;
    return;
L70:
    for (I=1; I <= N; ++I)
    {
        NS=max(NSCAT,NS);
        for (J=1; J <= NS; ++J) if(parac->NTYP[I+IOF] == coefc->NAM[J])goto L100;
        for (J=1; J <= 85; ++J) if(parac->NTYP[I+IOF] == TABNC[J])goto L120;
        file6 << " SCATTERING LENGTHS NOT FOUND FOR " << parac->NTYP[I+IOF] << endl;
        //111     FORMAT(34H SCATTERING LENGTHS NOT FOUND FOR ,A4)
        DBWSException("SCATTERING DATA MISSING");
L100:
        params->PTR[I+IOF]=J;
        goto L80;
L120:
        NS=NS+1;
        params->PTR[I+IOF]=NS;
        coeff->DFP[NS]=TABN[J];
        coefc->NAM[NS]=TABNC[J];
        if (J  >=  61  &&  J  <=  81) goto L81;
        if (J  ==  82) goto L82;
        if (J  >=  83  &&  J  <=  85) goto L83;
        goto L84;
L81:
        NSL = J+2;
        goto L180;
L82:
        NSL = 90;
        goto L180;
L83:
        NSL = J+9;
        goto L180;
L84:
        NSL = J;
L180:
        coeff->XMAS[NS]=TBM[NSL];
L80:;
    }
    dc->NSAVE=NS;
}

void DBWS::CELL2(int NPHASE, double LAMDAM)
{
    const string LAU[14+1] = {"",  "1BAR","2/M","MMM","4/M","4/MMM","3BAR   R","3BAR M R","3BAR","3BAR M 1","3BAR 1 M","6/M","6/MMM","M3","M3M"};

    double COSA,COSB,COSC,SINA,SINB,SINC,SINASR,SINBSR,COSASR,COSBSR,COSCSR;


    file6 << "THE LAUE SYMMETRY IS " << LAU[phases[NPHASE].SYMB.NSPGRP] << endl;
    if ( phases[NPHASE].SYMB.NAXIS > 3 ) DBWSException("5001");
    switch (phases[NPHASE].SYMB.NSPGRP) {
    case 1:
        goto L1002;
        break;
    case 2:
        goto L1071;
        break;
    case 3:
        goto L1072;
        break;
    case 4:
        goto L1073;
        break;
    case 5:
        goto L1073;
        break;
    case 6:
        goto L1074;
        break;
    case 7:
        goto L1074;
        break;
    case 8:
        goto L1075;
        break;
    case 9:
        goto L1075;
        break;
    case 10:
        goto L1075;
        break;
    case 11:
        goto L1075;
        break;
    case 12:
        goto L1075;
        break;
    case 13:
        goto L1076;
        break;
    case 14:
        goto L1076;
        break;
    }
    GOTOER();
L1071:
    if ( phases[NPHASE].SYMB.NAXIS != 1 ) cellx->ALPHA=90.0;
    if ( phases[NPHASE].SYMB.NAXIS != 2 ) cellx->BETA=90.0;
    if ( phases[NPHASE].SYMB.NAXIS != 3) cellx->GAMMA=90.0;
    goto L1002;
L1073:
    cellx->B = cellx->A;
L1072:
    cellx->ALPHA = 90.0;
    cellx->BETA = 90.0;
    cellx->GAMMA = 90.0;
    goto L1002;
L1074:
    cellx->B = cellx->A;
    cellx->C = cellx->A;
    cellx->BETA = cellx->ALPHA;
    cellx->GAMMA = cellx->ALPHA;
    goto L1002;
L1075:
    cellx->B = cellx->A;
    cellx->ALPHA = 90.;
    cellx->BETA = 90.0;
    cellx->GAMMA = 120.0;
    goto L1002;
L1076:
    cellx->B = cellx->A;
    cellx->C = cellx->A;
    goto L1072;
L1002:
    COSA = cos(6.28318531*cellx->ALPHA/360.0);
    //L105:
    COSB = cos(6.28318531*cellx->BETA/360.0);
    //L1052:
    COSC = cos(6.28318531*cellx->GAMMA/360.0);
    //L1053:
    SINA = sin(6.28318531*cellx->ALPHA/360.0);		// SINA=sqrt(1.0-pow(COSA,2));
    SINB = sin(6.28318531*cellx->BETA/360.0);		// SINB=sqrt(1.0-pow(COSB,2));
    SINC = sin(6.28318531*cellx->GAMMA/360.0);		// SINC=sqrt(1.0-pow(COSC,2));




    //
    //-----THE VOLUME OF THE CELL IN THE DIRECT SPACE OF THE REAL CELL OF
    //-----THE K-TH PHASE=NPHASE IS CALCULATED
    //-----AND THE RAY GMAX OF THE SPHERE WITH THE VOLUME EQUIVALENT TO THE
    //-----BRILLOUIN CELL 4/3*PI*GMAX**3=1/VOL

    // Volume of cell in direct space
    volume->VOLI[NPHASE] = cellx->A * cellx->B * cellx->C * sqrt(1.0 - pow(COSA,2) - pow(COSB,2) - pow(COSC,2) + 2.0 * COSA * COSB * COSC );

    volume->GCOM[NPHASE] = 0.877298169 * volume->VOLI[NPHASE] / pow(LAMDAM,3);

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
    cellx->AL[1][1]=1.0/(cellx->A * SINC * SINBSR);

    // b*, a in space reciprocal. Giacovazzo, p67 table 2.1
    cellx->AL[2][2]=1.0/(cellx->B * SINC * SINASR);

    // c*, a in space reciprocal. Giacovazzo, p67 table 2.1
    cellx->AL[3][3]=1.0/(cellx->C * SINA * SINBSR);

    cellx->AL[2][3]=cellx->AL[2][2]*cellx->AL[3][3]*COSASR*2.0;
    cellx->AL[1][3]=cellx->AL[3][3]*cellx->AL[1][1]*COSBSR*2.0;
    cellx->AL[1][2]=cellx->AL[2][2]*cellx->AL[1][1]*COSCSR*2.0;

    cellx->AL[1][1]=cellx->AL[1][1]*cellx->AL[1][1];
    cellx->AL[2][2]=cellx->AL[2][2]*cellx->AL[2][2];
    cellx->AL[3][3]=cellx->AL[3][3]*cellx->AL[3][3];
}

double DBWS::MULT(int IPHASE, int IH, int IK, int IL, int KXIS)
{
    double P;

    P = 2.0;
    switch (phases[IPHASE].SYMB.NSPGRP)
    {
    case 1:
        goto L220;
        break;
    case 2:
        goto L100;
        break;
    case 3:
        goto L110;
        break;
    case 4:
        goto L120;
        break;
    case 5:
        goto L130;
        break;
    case 6:
        goto L150;
        break;
    case 7:
        goto L180;
        break;
    case 8:
        goto L140;
        break;
    case 9:
        goto L170;
        break;
    case 10:
        goto L160;
        break;
    case 11:
        goto L170;
        break;
    case 12:
        goto L190;
        break;
    case 13:
        goto L200;
        break;
    case 14:
        goto L210;
        break;
    }
    GOTOER();
L100:
    if ( IH <= 0 ) goto L220; else goto L101;
L101:
    if ( KXIS == 3 )
    {
        if ( IL > 0  &&  IH+abs(IK) > 0 ) P=2.0*P;
    }
    else
    {
        if ( IK > 0  &&  IH+abs(IL) > 0 ) P=2.0*P;
    }
    goto L220;
L110:
    P = 1.0;
    if ( IH <= 0 ) goto L112; else goto L111;
L111:
    P = 2.0*P;
L112:
    if ( IK <= 0) goto L114; else goto L113;
L113:
    P = 2.0*P;
L114:
    if ( IL <= 0 ) goto L220; else goto L115;
L115:
    P = 2.0*P;
    goto L220;
L120:
    if ( IK <= 0 ) goto L220; else goto L121;
L121:
    P = 2.0*P;
    if ( IL <= 0 ) goto L220; else goto L122;
L122:
    P = 2.0*P;
    goto L220;
L130:
    if ( IK <= 0) goto L220; else goto L131;
L131:
    P = 2.0*P;
    if ( IH <= 0 ) goto L134; else goto L132;
L132:
    if ( IH-IK < 0 ) goto L133; else goto L134;
L133:
    P = 2.0*P;
L134:
    if ( IL <= 0 ) goto L220; else goto L135;
L135:
    P = 2.0*P;
    goto L220;
L140:
    if ( IH+abs(IK) <= 0 ) goto L142; else goto L141;
L141:
    P = 3.0*P;
L142:
    if ( IL <= 0 ) goto L220; else goto L143;
L143:
    P = 2.0*P;
    goto L220;
L150:
    if ( IH+IK+IL == 0 ) goto L152; else goto L151;
L151:
    P = 2.0*P;
L152:
    if ( IH-IK == 0 ) goto L154; else goto L153;
L153:
    P = 3.0*P;
    goto L220;
L154:
    if ( IH-IL == 0) goto L220; else goto L153;
L160:
    if ( IK <= 0 ) goto L220; else goto L161;
L161:
    P = 3.0*P;
    if ( IK-IH <= 0 ) goto L220; else goto L162;
L162:
    if ( IH <= 0 ) goto L220; else goto L135;
L170:
    if ( IL <= 0 ) goto L172; else goto L171;
L171:
    P = 2.0*P;
L172:
    if ( IH <= 0 ) goto L174; else goto L173;
L173:
    P = 3.0*P;
    if ( IK <= 0 ) goto L220; else goto L135;
L174:
    if ( IK <= 0 ) goto L220; else goto L153;
L180:
    if ( IH-IK == 0 ) goto L181; else goto L182;
L181:
    if ( IH-IL == 0 ) goto L220; else goto L183;
L182:
    if ( IK-IL == 0 ) goto L183; else goto L184;
L183:
    P = P/2.0;
L184:
    P = 3.0*P;
    if ( IH+IK+IL == 0 ) goto L220; else goto L135;
L190:
    if ( IK <= 0 ) goto L194; else goto L191;
L191:
    P = 6.0*P;
    if ( IH == 0 ) goto L194; else goto L192;
L192:
    if ( IK-IH == 0 ) goto L194; else goto L193;
L193:
    P = 2.0*P;
L194:
    if ( IL == 0 ) goto L220; else goto L135;
L200:
    if ( IH-IK == 0 ) goto L201; else goto L203;
L201:
    if ( IK-IL == 0 ) goto L202; else goto L203;
L202:
    P = 4.0*P;
    goto L220;
L203:
    P = 3.0*P;
    if ( IK <= 0 ) goto L205; else goto L204;
L204:
    P = 2.0*P;
L205:
    if ( IH <= 0 ) goto L220; else goto L135;
L210:
    if ( IH == 0 ) goto L214; else goto L211;
L211:
    P = 4.0*P;
    if ( IH-IL == 0 ) goto L220; else goto L212;
L212:
    P = 3.0*P;
    if ( IH-IK == 0 ) goto L220; else goto L213;
L213:
    if ( IK-IL == 0 ) goto L220; else goto L135;
L214:
    P = 3.0*P;
    if ( IK == 0 ) goto L220; else goto L215;
L215:
    P = 2.0*P;
    if ( IK-IL == 0 ) goto L220; else goto L135;
    //-----THIS LINE ADDED TO CORRECT THE PROBLEMS WITH THE -3, -3M1, 6/M,
    //-----AND 6/MMM SPACE GROUPS THAT WERE OFF BY A FACTOR OF TWO
L220:
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
    cntrls->JOBTYP=cntrls->JOBTYP-1;
    //cntrls->NPROF=cntrls->NPROF-1;
    NPLOF = cntrls->NPROF;
    if (sizestrain->NSIZESTRAIN == 9) NPLOF = 9;
    cntrls->INSTRM = cntrls->INSTRM-1;

    // line 2.1
    // line 2.1 changed due to size-strain calculation (NsizeStrain)
    file5b << setw(4) << cntrls->JOBTYP
        << setw(4) << NPLOF
        << setw(4) << cntrls->NPHASE
        << setw(4) << codebck->IBCKCODE
        << setw(4) << jnk->NEXCRG
        << setw(4) << jnk->NSCAT
        << setw(4) << cntrls->INSTRM
        << setw(4) << cntrls->IPREF
        << setw(4) << cntrls->IASYM
        << setw(4) << cntrls->IABSR
        << setw(4) << cntrls->IDATA
        << setw(4) << cntrls->ISPHASE
        << setw(4) << cntrls->I2D94
        << "     LINE 2.1" << endl;

    // line 2.2
    if(codebck->IBCKCODE == -1)
    {
        file5b << setw(4) << cntrls->IAS
            << setw(4) << cntrls->FONDO
            << "                                                 LINE 2.2" << endl;
    }

    // line 3
    file5b
        << setw(1) << cntrls->IOT
        << setw(1) << cntrls->IPL
        << setw(1) << cntrls->IPC
        << setw(1) << cntrls->MAT
        << setw(1) << cntrls->NXT
        << " "
        << setw(1) << cntrls->LST1
        << setw(1) << cntrls->LST2
        << setw(1) << cntrls->LST3
        << setw(1) << cntrls->IPL1
        << setw(1) << cntrls->IPL2
        << " "
        << setw(1) << cntrls->IPLST
        << setw(1) << cntrls->IPLOSS
        << setw(1) << cntrls->IPLCAL
        << setw(1) << cntrls->IPLPOL
        << setw(1) << cntrls->IPLCOM
        << " "
        << setw(1) << cntrls->IPLDIS
        << setw(1) << cntrls->IPLAM
        << setw(1) << cntrls->IPBIG
        << "                                    LINE 3" << endl;

    // line 4
    file5b
        << setw(8) << setprecision(5) << g1->LAMDA[1]
        << setw(8) << setprecision(5) << g1->LAMDA[2]
        << setw(8) << setprecision(5) << params->RATIO[2]
        << setw(8) << setprecision(4) << g1->BKPOS
        << setw(8) << setprecision(4) << g1->WDT
        << setw(8) << setprecision(4) << g1->CTHM
        << setw(8) << setprecision(4) << g1->TMV
        << setw(8) << setprecision(4) << g1->RLIM
        << setw(8) << setprecision(4) << cntrls->SW
        << endl;

    // line 5
    file5b
        << setw(4) << cntrls->MCYCLE
        << setw(4) << setprecision(2) << cntrls->EPS
        << setw(4) << setprecision(2) << params->RELAX[1]
        << setw(4) << setprecision(2) << params->RELAX[2]
        << setw(4) << setprecision(2) << params->RELAX[3]
        << setw(4) << setprecision(2) << params->RELAX[4]
        << "                                 CYCLS EPS RELAX P_CALC" << endl;
    if(jnk->NBCKGD < 2)goto L120;

    // line 6(*)
    for (I=1; I <= jnk->NBCKGD; ++I) file5b << setw(8) << setprecision(2) << jnk->POS[I] << jnk->BCK[I] << endl;
L120:
    if(jnk->NEXCRG <= 0)goto L122;

    // line 7(*)
    for (I=1; I <= jnk->NEXCRG; ++I)
        file5b
            << setw(8) << setprecision(2) << jnk->ALOW[I]
            << setw(8) << setprecision(2) << jnk->AHIGH[I]
            << "                                         EXCLUDED REGION" << endl;
L122:
    if (jnk->NSCAT <= 0) goto L124;
    for (I=1; I <= jnk->NSCAT; ++I)
    {
        if (cntrls->JOBTYP  ==  1 || cntrls->JOBTYP == 3) goto L1228;
        //C line 8.1 XRD (*)
        file5b
            << setw(4) << coefc->NAM[I]
            << setw(8) << setprecision(4) << coeff->DFP[I]
            << setw(8) << setprecision(4) << coeff->DFPP[I]
            << setw(8) << setprecision(4) << coeff->XMAS[I]
            << "                             SCATTERING SET" << setw(2) << jnk->NSCAT << endl;
        goto L126;
        // line 8.1 ND(*)
L1228:
        file5b
            << setw(4) << coefc->NAM[I]
            << setw(8) << setprecision(4) << coeff->DFP[I]
            << setw(8) << setprecision(4) << coeff->XMAS[I]
            << "                                     SCATTERING SET " << setw(2) << jnk->NSCAT << endl;
        // line 8.2 XRD(*)
L126:
        if (cntrls->JOBTYP == 0 || cntrls->JOBTYP == 2)
        {
            for (J=1; J <= 9; ++J) file5b << setw(8) << setprecision(5) << coeff->AC[J][I];
            file5b << endl;
        }
    }
L124:;

    // line 9
    file5b << setw(8) << cntrls->MAXS << "                                                 PARAMS REFINED" << endl;
    N=0;
    for (IIPHAS=1; IIPHAS <= cntrls->NPHASE; ++IIPHAS) N=N+phases[IIPHAS].AtomCount;
    for (I=1; I <= N; ++I)
    {
        for (J=1; J <= 11; ++J) params->A[I][J]=sign(params->A[I][J])*(static_cast<double>(10*params->LP[I][J])+abs(params->A[I][J]));
    }
    for (I=1; I <= cntrls->NPHASE; ++I)
    {
        for (J=1; J <= 6; ++J) phases[I-1].PAR[J+5-1]=dc->SAVE[I][J];
        for (J=1; J <= 27; ++J) phases[I-1].PAR[J-1].codeword=sign(phases[I-1].PAR[J-1].codeword)*(static_cast<double>(10*phases[I-1].PAR[J-1].L)+abs(phases[I-1].PAR[J-1].codeword));
    }
    for (J=1; J <= 20; ++J) params->GLB_[J-1].codeword=sign(params->GLB_[J-1].codeword)*(static_cast<double>(10*params->GLB_[J-1].L)+abs(params->GLB_[J-1].codeword));

    // line 10.1
    file5b
        << setw(8) << setprecision(4) << params->GLB_[1-1]
        << setw(8) << setprecision(4) << params->GLB_[10-1]
        << setw(8) << setprecision(4) << params->GLB_[11-1]
        << setw(8) << setprecision(4) << params->GLB_[8-1]
        << setw(8) << setprecision(4) << params->GLB_[9-1]
        << setw(8) << setprecision(4) << params->GLB_[12-1]
        << setw(8) << setprecision(4) << params->GLB_[13-1] << " ZER DISP TRANS p q r t" << endl;

    // line 10.11
    file5b
        << setw(8) << setprecision(4) << params->GLB_[1-1]
        << setw(8) << setprecision(4) << params->GLB_[10-1].codeword
        << setw(8) << setprecision(4) << params->GLB_[11-1].codeword
        << setw(8) << setprecision(4) << params->GLB_[8-1].codeword
        << setw(8) << setprecision(4) << params->GLB_[9-1].codeword
        << setw(8) << setprecision(4) << params->GLB_[12-1].codeword
        << setw(8) << setprecision(4) << params->GLB_[13-1].codeword << " CODEWORDS" << endl;
    if (cntrls->IBGD == 1) goto L4600;
    // line 10.2
    file5b
        << setw(8) << setprecision(4) << params->GLB_[20-1]
        << setw(8) << setprecision(4) << params->GLB_[18-1]
        << setw(8) << setprecision(4) << params->GLB_[19-1]
        << "                                 AM MON1 MON2" << endl;

    // line 10.21
    file5b
        << setw(8) << setprecision(4) << params->GLB_[20-1].codeword
        << setw(8) << setprecision(4) << params->GLB_[18-1].codeword
        << setw(8) << setprecision(4) << params->GLB_[19-1].codeword
        << "                                 CODEWORS" << endl;
L4600:
    if(jnk->NBCKGD == 0)
    {
        file5b
            << setw(9) << setprecision(2) << params->GLB_[2-1]
            << setw(9) << setprecision(2) << params->GLB_[3-1]
            << setw(9) << setprecision(2) << params->GLB_[4-1]
            << setw(9) << setprecision(2) << params->GLB_[5-1]
            << setw(9) << setprecision(2) << params->GLB_[6-1]
            << setw(9) << setprecision(2) << params->GLB_[7-1]
            << "   BACKGROUND" << endl;
        file5b
            << setw(9) << setprecision(2) << params->GLB_[2-1].codeword
            << setw(9) << setprecision(2) << params->GLB_[3-1].codeword
            << setw(9) << setprecision(2) << params->GLB_[4-1].codeword
            << setw(9) << setprecision(2) << params->GLB_[5-1].codeword
            << setw(9) << setprecision(2) << params->GLB_[6-1].codeword
            << setw(9) << setprecision(2) << params->GLB_[7-1].codeword
            << "   CODEWORDS" << endl;
    }
    //L477:
    for (K=1; K <= cntrls->NPHASE; ++K)
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
        for (ISOF=1; ISOF <= N; ++ISOF) params->XL[ISOF+IOF][5]=params->XL[ISOF+IOF][5]*multip->XMLTP[K]/multip->MURT[ISOF+IOF];
        // !cp oct 96 #6 murt parametrs included below. FORMAT modified...
        for (I=1; I <= N; ++I)
        {
            file5b
                << setw(4) << parac->ATEXT[I+IOF] << " "
                << setw(4) << multip->MURT[I+IOF] << " "
                << setw(4) << parac->NTYP[I+IOF] << "  "
                << setw(8) << setprecision(5) << params->XL[I+IOF][1]
                << setw(8) << setprecision(5) << params->XL[I+IOF][2]
                << setw(8) << setprecision(5) << params->XL[I+IOF][3]
                << setw(8) << setprecision(5) << params->XL[I+IOF][4]
                << setw(8) << setprecision(5) << params->XL[I+IOF][5]
                << "  LBL M NTYP x y z B So" << endl
                << "                "
                << setw(8) << setprecision(5) << params->A[I+IOF][1]
                << setw(8) << setprecision(5) << params->A[I+IOF][2]
                << setw(8) << setprecision(5) << params->A[I+IOF][3]
                << setw(8) << setprecision(5) << params->A[I+IOF][4]
                << setw(8) << setprecision(5) << params->A[I+IOF][5]
                << "  CODEWORDS" << endl
                << setw(8) << setprecision(5) << params->XL[I+IOF][6]
                << setw(8) << setprecision(5) << params->XL[I+IOF][7]
                << setw(8) << setprecision(5) << params->XL[I+IOF][8]
                << setw(8) << setprecision(5) << params->XL[I+IOF][9]
                << setw(8) << setprecision(5) << params->XL[I+IOF][10]
                << setw(8) << setprecision(5) << params->XL[I+IOF][11]
                << "          BETAS" << endl
                << setw(8) << setprecision(2) << params->A[I+IOF][6]
                << setw(8) << setprecision(2) << params->A[I+IOF][7]
                << setw(8) << setprecision(2) << params->A[I+IOF][8]
                << setw(8) << setprecision(2) << params->A[I+IOF][9]
                << setw(8) << setprecision(2) << params->A[I+IOF][10]
                << setw(8) << setprecision(2) << params->A[I+IOF][11]
            << "          CODEWORDS" << endl;
        }
        file5b
            << scientific << setw(8) << setprecision(2) << phases[K-1].PAR[0]
            << fixed
            << setw(8) << setprecision(4) << phases[K-1].PAR[1]
            << "                                         SCALE Bo(OVERALL)" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[0].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[1].codeword << endl
            << setw(8) << setprecision(5) << phases[K-1].PAR[2]
            << setw(8) << setprecision(5) << phases[K-1].PAR[3]
            << setw(8) << setprecision(5) << phases[K-1].PAR[4]
            << setw(8) << setprecision(5) << phases[K-1].PAR[20]
            << setw(8) << setprecision(5) << phases[K-1].PAR[19]
            << setw(8) << setprecision(5) << phases[K-1].PAR[14]
            << setw(8) << setprecision(5) << phases[K-1].PAR[15]
            << " U V W CT Z X Y" << endl

            << setw(8) << setprecision(2) << phases[K-1].PAR[2].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[3].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[4].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[20].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[19].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[14].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[15].codeword << endl

            << setw(8) << setprecision(4) << phases[K-1].PAR[5]
            << setw(8) << setprecision(4) << phases[K-1].PAR[6]
            << setw(8) << setprecision(4) << phases[K-1].PAR[7]
            << setw(8) << setprecision(4) << phases[K-1].PAR[8]
            << setw(8) << setprecision(4) << phases[K-1].PAR[9]
            << setw(8) << setprecision(4) << phases[K-1].PAR[10]
            << "         CELL PARAMETERS" << endl

            << setw(8) << setprecision(2) << phases[K-1].PAR[5].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[6].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[7].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[8].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[9].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[10].codeword << endl

            << setw(8) << setprecision(5) << phases[K-1].PAR[11]
            << setw(8) << setprecision(5) << phases[K-1].PAR[12]
            << setw(8) << setprecision(5) << phases[K-1].PAR[13]
            << "                                 PREF1 PREF2 R/RCF_ASYM" << endl

            << setw(8) << setprecision(2) << phases[K-1].PAR[11].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[12].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[13].codeword << endl

            << setw(8) << setprecision(4) << phases[K-1].PAR[16]
            << setw(8) << setprecision(4) << phases[K-1].PAR[17]
            << setw(8) << setprecision(4) << phases[K-1].PAR[18]
            << "                                 NA NB NC (MIX_PARAMS)" << endl

            << setw(8) << setprecision(2) << phases[K-1].PAR[16].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[17].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[18].codeword << endl

            << setw(8) << setprecision(4) << phases[K-1].PAR[23]
            << setw(8) << setprecision(4) << phases[K-1].PAR[24]
            << setw(8) << setprecision(4) << phases[K-1].PAR[25]
            << "                                 NA NB NC (HIGH SIDE)" << endl

            << setw(8) << setprecision(2) << phases[K-1].PAR[23].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[24].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[25].codeword << endl

            << setw(8) << setprecision(4) << phases[K-1].PAR[26] << "" << endl

            << setw(8) << setprecision(2) << phases[K-1].PAR[26].codeword << endl;
    }
    if (cntrls->IPL != 0) file5b << setw(8) << ISCALE << IDIF << "                                         LINE PRINTER INFO" << endl;
    if (file5b.is_open()) file5b.close();
    //L151:
    cntrls->JOBTYP=cntrls->JOBTYP+1;
}

// subroutine to compute size&strain (NsizeStrain)
void DBWS::size(int K)
{
    double A7, B7, C7, D7,HD, UA, XA, YA, ZA, HS, UI, XI, YI, ZI, HDG, HDL, HSG, HSL,RADO;
    int ISIZETEST,ISTRAINTEST;
    ifstream file17;
    string s;

    RADO=3.1415927/180.0;
    A7 = 2.69269;
    B7 = 2.42843;
    C7 = 4.47163;
    D7 = 0.07842;
    file17.open("instr.dat");
    getline(file17,s);				// na primeira linha de padrao.dat sera apenas para comentario
    // a segunda linha ira ler U,Z,X,Y da FWHM do padrao
    file17 >> UI >> ZI >> XI >> YI;
    file17.close();
    UA = phases[K-1].PAR[2];
    ZA = phases[K-1].PAR[19];
    XA = phases[K-1].PAR[14];
    YA = phases[K-1].PAR[15];
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
        sizestrain->SIZEG[K] =  g1->LAMDA[1] / HSG;
    }
    else
    {
        HSG = 0.0;
        sizestrain->SIZEG[K] = 99999.0;
    }

    //
    //        COMPUTE SIZE BASED ON Loren-COS Y parameter [sizeL(k)]
    //  params->PAR[K][16] = Y (Lorentz)
    //
    if(YA  >  YI)
    {
        HSL = (YA - YI) * RADO;
        sizestrain->SIZEL[K] =  g1->LAMDA[1] / HSL;
    }
    else
    {
        HSL = 0.0;
        sizestrain->SIZEL[K] = 99999.0;
    }
    //
    //        COMPUTE WEIGHTED SIZE BASED on Gauss & Lorentz [siz(k)]
    //
    if(ISIZETEST == 0)
    {
        HS  = pow((pow(HSG,5) + A7 * pow(HSG,4) * HSL + B7 * pow(HSG,3) * pow(HSL,2) + C7 * pow(HSG,2) * pow(HSL,3) + D7 * HSG * pow(HSL,4) + pow(HSL,5)),0.2);
        sizestrain->SIZ[K] = g1->LAMDA[1] / HS;
    }
    else
    {
        sizestrain->SIZ[K] = 0.0;
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
        sizestrain->STRAING[K] =  HDG / 2;
    }
    else
    {
        HDG = 0.0;
        sizestrain->STRAING[K] = 0.0;
    }
    //
    // compute strain based on Lorentz-X * [strainL(k)]
    //
    if(XA  >  XI)
    {
        HDL = (XA - XI) * RADO;
        sizestrain->STRAINL[K] =  HDL / 2;
    }
    else
    {
        HDL = 0.0;
        sizestrain->STRAINL[K] = 0.0;
    }
    // compute weighted strain: based on Gauss-U and Lorentz-X * [strain(k)]
    //
    if (ISTRAINTEST == 0)
    {
        HD  = pow(( pow(HDG,5) + A7 * pow(HDG,4) * HDL + B7 * pow(HDG,3) * pow(HDL,2) + C7 * pow(HDG,2) * pow(HDL,3) + D7 * HDG * pow(HDL,4) + pow(HDL,5)),0.2);
        sizestrain->STRAIN[K] = HD / 2;
    }
    else
    {
        sizestrain->STRAIN[K] = 0.0;
    }
    file6 << "Gauss(phase " << setw(2) << K << "): Size =" //<< scientific
        << setw(9) << setprecision(4) << sizestrain->SIZEG[K]
    << " Strain ="
        << setw(9) << setprecision(3) << sizestrain->STRAING[K] << endl
        << "Lorentz(phase " << setw(2) << K << "): Size ="
        << setw(9) << setprecision(4) << sizestrain->SIZEL[K]
    << " Strain ="
        << setw(9) << setprecision(3) << sizestrain->STRAINL[K] << endl
        << "Weighted(phase " << setw(K) << "): Size ="
        << setw(9) << setprecision(4) << sizestrain->SIZ[K]
    << " Strain ="
        << setw(9) << setprecision(3) << sizestrain->STRAIN[K] << endl;
L20:
    return;
}

//      SUBROUTINE WRIT94 Gera uma saida para o 9411, pra poder Usar no ATOMS
void DBWS::WRITE94(int ISCALE, int IDIF)
{
    int I, J, K, N, IOF, IIPHAS;

    //REWIND 53;

    // line 1
    file53 << setw(70) << title << endl;
    cntrls->JOBTYP=cntrls->JOBTYP-1;
    cntrls->INSTRM = cntrls->INSTRM-1;

    // line 2.1
    file53 << setw(9) << setprecision(4) << cntrls->JOBTYP
        << setw(9) << setprecision(4) << cntrls->NPROF
        << setw(9) << setprecision(4) << cntrls->NPHASE
        << setw(9) << setprecision(4) << codebck->IBCKCODE
        << setw(9) << setprecision(4) << jnk->NEXCRG
        << setw(9) << setprecision(4) << jnk->NSCAT
        << setw(9) << setprecision(4) << cntrls->INSTRM
        << setw(9) << setprecision(4) << cntrls->IPREF
        << setw(9) << setprecision(4) << cntrls->IABSR
        << "                     LINE 2.1" << endl;

    // line 3
    file53 << setw(1) << cntrls->IOT
        << setw(1) << cntrls->IPL
        << setw(1) << cntrls->IPC
        << setw(1) << cntrls->MAT
        << setw(1) << cntrls->NXT
        << setw(1) << cntrls->LST1
        << setw(1) << cntrls->LST2
        << setw(1) << cntrls->LST3
        << setw(1) << cntrls->IPL1
        << setw(1) << cntrls->IPL2
        << setw(1) << cntrls->IPLST
        << "                                              LINE 3" << endl;

    // line 4
    file53 << setw(8) << setprecision(5) << g1->LAMDA[1]
    << setw(8) << setprecision(5) << g1->LAMDA[2]
    << setw(8) << setprecision(5) << params->RATIO[2]
    << setw(8) << setprecision(4) << g1->BKPOS
        << setw(8) << setprecision(4) << g1->WDT
        << setw(8) << setprecision(4) << g1->CTHM
        << setw(8) << setprecision(4) << g1->TMV
        << setw(8) << setprecision(4) << g1->RLIM << endl;

    // line 5
    file53 << setw(4) << cntrls->MCYCLE
        << setw(4) << setprecision(2) << cntrls->EPS
        << setw(4) << setprecision(2) << params->RELAX[1]
    << setw(4) << setprecision(2) << params->RELAX[2]
    << setw(4) << setprecision(2) << params->RELAX[3]
    << setw(4) << setprecision(2) << params->RELAX[4]
    << "                                 CYCLS EPS RELAX P_CALC" << endl;
    if (jnk->NBCKGD < 2) goto L120;

    // line 6(*)
    for (I=1; I <= jnk->NBCKGD; ++I)
    {
        file53 << setw(8) << setprecision(2) << jnk->POS[I]
        << setw(8) << setprecision(2) << jnk->BCK[I] << endl;
    }
L120:
    if(jnk->NEXCRG <= 0)goto L122;
    // line 7(*)
    for (I=1; I <= jnk->NEXCRG; ++I)
    {
        file53 << setw(8) << setprecision(2) << jnk->ALOW[I]
        << setw(8) << setprecision(2) << jnk->AHIGH[I] << endl;
    }
L122:
    if (jnk->NSCAT <= 0) goto L124;
    for (I=1; I <= jnk->NSCAT; ++I)
    {
        if (cntrls->JOBTYP  ==  1 || cntrls->JOBTYP == 3) goto L1228;
        // line 8.1 XRD (*)
        file53 << setw(4) << coefc->NAM[I]
        << setw(8) << setprecision(4) << coeff->DFP[I]
        << setw(8) << setprecision(4) << coeff->DFPP[I]
        << setw(8) << setprecision(4) << coeff->XMAS[I]
        << "                             SCATTERING SET " << setw(2) << jnk->NSCAT << endl;
        goto L126;
L1228:
        // line 8.1 ND(*)
        file53 << setw(4) << coefc->NAM[I]
        << setw(8) << setprecision(4) << coeff->DFP[I]
        << setw(8) << setprecision(4) << coeff->XMAS[I]
        << "                                     SCATTERING SET " << setw(2) << jnk->NSCAT << endl;
L126:
        // line 8.2 XRD(*)
        if(cntrls->JOBTYP == 0 || cntrls->JOBTYP == 2)
        {
            for (J=1; J <= 9; ++J) file53 << setw(8) << setprecision(5) << coeff->AC[J][I];
            file53 << endl;
        }
    }
L124:

    // line 9
    file53 << setw(8) << cntrls->MAXS << "                                                 PARAMS REFINED" << endl;
    N=0;
    for (IIPHAS=1; IIPHAS <= cntrls->NPHASE; ++IIPHAS) N=N+phases[IIPHAS].AtomCount;
    for (I=1; I <= N; ++I)
    {
        for (J=1; J <= 11; ++J) params->A[I][J]=sign(params->A[I][J])*(static_cast<double>(10*params->LP[I][J])+abs(params->A[I][J]));
    }
    for (I=1; I <= cntrls->NPHASE; ++I)
    {
        for (J=1; J <= 6; ++J) phases[I-1].PAR[J+5-1]=dc->SAVE[I][J];
        for (J=1; J <= 27; ++J) phases[I-1].PAR[J-1].codeword=sign(phases[I-1].PAR[J-1].codeword)*(static_cast<double>(10*phases[I-1].PAR[J-1].L)+abs(phases[I-1].PAR[J-1].codeword));
    }
    for (J=1; J <= 20; ++J) params->GLB_[J-1].codeword=sign(params->GLB_[J-1].codeword)*(static_cast<double>(10*params->GLB_[J-1].L)+abs(params->GLB_[J-1].codeword));

    // line 10.1
    file53 << setw(8) << setprecision(4) << params->GLB_[1-1]
    << setw(8) << setprecision(4) << params->GLB_[1-1].codeword
    << setw(8) << setprecision(4) << params->GLB_[10-1]
    << setw(8) << setprecision(4) << params->GLB_[10-1].codeword
    << setw(8) << setprecision(4) << params->GLB_[11-1]
    << setw(8) << setprecision(4) << params->GLB_[11-1].codeword
    << "         ZER DISP TRANS + CODEWORS" << endl
        << setw(9) << setprecision(4) << 0.0
        << setw(9) << setprecision(4) << 0.0
        << setw(9) << setprecision(4) << 0.0
        << setw(9) << setprecision(4) << 0.0
        << setw(9) << setprecision(4) << 0.0
        << setw(9) << setprecision(4) << 0.0
        << "   CODEWORDS" << endl;
    //L477:
    for (K=1; K <= cntrls->NPHASE; ++K)
    {
        IOF=0;
        if(K > 1)
        {
            for (IIPHAS=2; IIPHAS <= K; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
        }
        // line 11.1
        file53 << setw(50) << phases[K].name << "       PHASE NUMBER " << setw(2) << K << endl;

        // line 11.2
        file53 << setw(4) << phases[K].AtomCount
        << setw(4) << phases[K].NMOL
        << "        "
            << setw(4) << setprecision(1) << phases[K].PREF[1]
        << setw(4) << setprecision(1) << phases[K].PREF[2]
        << setw(4) << setprecision(1) << phases[K].PREF[3]
        << "                             #ATMS #FU PREFDIR" << endl;
        N=phases[K].AtomCount;

        // line 11.3
        file53 << setw(20) << phases[K].SYMB
        << "                                     SPACE GROUP" << endl;
        for (I=1; I <= N; ++I)
        {
            file53 << setw(4) << parac->ATEXT[I+IOF]
            << setw(4) << parac->NTYP[I+IOF]
            << "        "
                << setw(8) << setprecision(5) << params->XL[I+IOF][1]
            << setw(8) << setprecision(5) << params->XL[I+IOF][2]
            << setw(8) << setprecision(5) << params->XL[I+IOF][3]
            << setw(8) << setprecision(5) << params->XL[I+IOF][4]
            << setw(8) << setprecision(5) << params->XL[I+IOF][5]
            << "  LBL NTYP x y z B So" << endl
                << "                "
                << setw(8) << setprecision(2) << params->A[I+IOF][1]
            << setw(8) << setprecision(2) << params->A[I+IOF][2]
            << setw(8) << setprecision(2) << params->A[I+IOF][3]
            << setw(8) << setprecision(2) << params->A[I+IOF][4]
            << setw(8) << setprecision(2) << params->A[I+IOF][5]
            << "  CODEWORDS" << endl
                << setw(8) << setprecision(5) << params->XL[I+IOF][6]
            << setw(8) << setprecision(5) << params->XL[I+IOF][7]
            << setw(8) << setprecision(5) << params->XL[I+IOF][8]
            << setw(8) << setprecision(5) << params->XL[I+IOF][9]
            << setw(8) << setprecision(5) << params->XL[I+IOF][10]
            << setw(8) << setprecision(5) << params->XL[I+IOF][11]
            << "          BETAS" << endl
                << setw(8) << setprecision(2) << params->A[I+IOF][6]
            << setw(8) << setprecision(2) << params->A[I+IOF][7]
            << setw(8) << setprecision(2) << params->A[I+IOF][8]
            << setw(8) << setprecision(2) << params->A[I+IOF][9]
            << setw(8) << setprecision(2) << params->A[I+IOF][10]
            << setw(8) << setprecision(2) << params->A[I+IOF][11]
            << "          CODEWORDS" << endl;
        }
        file53 << scientific
            << setw(8) << setprecision(2) << phases[K-1].PAR[0]
        << fixed
            << setw(8) << setprecision(4) << phases[K-1].PAR[1]
        << "                                    SCALE Bo(OVERALL)" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[0].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[1].codeword << endl
            << setw(8) << setprecision(5) << phases[K-1].PAR[2]
        << setw(8) << setprecision(5) << phases[K-1].PAR[3]
        << setw(8) << setprecision(5) << phases[K-1].PAR[4]
        << setw(8) << setprecision(5) << phases[K-1].PAR[20]
        << setw(8) << setprecision(5) << phases[K-1].PAR[19]
        << setw(8) << setprecision(5) << phases[K-1].PAR[14]
        << setw(8) << setprecision(5) << phases[K-1].PAR[15]
        << " U V W CT Z X Y" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[2].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[3].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[4].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[20].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[19].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[14].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[15].codeword << endl
            << setw(8) << setprecision(4) << phases[K-1].PAR[5]
        << setw(8) << setprecision(4) << phases[K-1].PAR[6]
        << setw(8) << setprecision(4) << phases[K-1].PAR[7]
        << setw(8) << setprecision(4) << phases[K-1].PAR[8]
        << setw(8) << setprecision(4) << phases[K-1].PAR[9]
        << setw(8) << setprecision(4) << phases[K-1].PAR[10]
        << "         CELL PARAMETERS" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[5].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[6].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[7].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[8].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[9].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[10].codeword << endl
            << setw(8) << setprecision(5) << phases[K-1].PAR[11]
        << setw(8) << setprecision(5) << phases[K-1].PAR[12]
        << setw(8) << setprecision(5) << phases[K-1].PAR[13]
        << "                                 PREF1 PREF2 R/RCF_ASYM" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[11].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[12].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[13].codeword << endl
            << setw(8) << setprecision(4) << phases[K-1].PAR[16]
        << setw(8) << setprecision(4) << phases[K-1].PAR[17]
        << setw(8) << setprecision(4) << phases[K-1].PAR[18]
        << "                                 NA NB NC (MIX_PARAMS)" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[16].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[17].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[18].codeword << endl
            << setw(8) << setprecision(4) << phases[K-1].PAR[23]
        << setw(8) << setprecision(4) << phases[K-1].PAR[24]
        << setw(8) << setprecision(4) << phases[K-1].PAR[25]
        << "                                 NA NB NC (HIGH SIDE)" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[23].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[24].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[25].codeword << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[26]
        << "                                                 PEARSON ASYM.FACTOR" << endl
            << setw(8) << setprecision(2) << phases[K-1].PAR[26].codeword << endl;
    }

    if (file5.is_open()) file5.close();
    file5b.open(file5name.data());
    if (cntrls->IPL != 0) file5b << setw(8) << ISCALE << setw(8) << IDIF << "                                         LINE PRINTER INFO" << endl;
    file5.close();

    //L151:
    //                       JOBTYP  must be reproduced !!!
    cntrls->JOBTYP=cntrls->JOBTYP+1;
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
    //if (cntrls->NXT != 0 || cntrls->JOBTYP >= 3) cntrls->NPROF=cntrls->NPROF-1;
    RAD = 45.0/atan(1.0);
    //	ZER = params->GLB_[1-1];
    if (cntrls->IPL != 0)
    {
        getline(file5,s);
        stringstream(s.substr(0,8)) >> ISCALE;
        stringstream(s.substr(8,8)) >> IDIF;
    }
    //L10000:
    if (cntrls->NXT != 0 && cntrls->JOBTYP < 3) REWRIT(ISCALE,IDIF);
    if (cntrls->IPL2 != 0)
    {
        file69.open("bragg.dat");
        file690.open("xy-int.dat");
        file690 << "  2theta_i      y_o       y_c        yo-yc  w(yo-yc)^2" << endl;
    }
    if (cntrls->JOBTYP > 2) goto L37;
    file6 << endl
        << "AVERAGE INTENSITY DIFFERENCE FOR PATTERN" << endl
        << "GIVEN FOR BLOCKS OF 20 OBSERVATIONS." << endl;
    for (I=1; I <= datax->NPTS; I = I + 200)
    {
        for (J=1; J <= 10; ++J)
        {
            J20=20*(J-1);
            DEL[J]=0;
            for (K=1; K <= 20; ++K)
            {
                J20IK=J20+I+K-1;
                DEL[J]=DEL[J]+datax->Y[J20IK]-datax->BK[J20IK]-datax->YC[J20IK];
            }
            KI[J]=J+I/20;
            DEL[J]=DEL[J]/20.0;
            if (KI[J]*20 >= datax->NPTS) goto L15;
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
L37:
    for (I=1; I <= datax->NPTS; ++I)
    {
        datax->YC[I]=datax->YC[I]+datax->BK[I];
        if (cntrls->JOBTYP > 2) datax->Y[I]=datax->YC[I];
    }
    if (cntrls->IPC == 0) goto L36;
    //INQUIRE(UNIT=6,NAME=OUTFILE);
    //LBKSL = 0;
    //for (I=1; I <= 80; ++I) if ( OUTFILE[I] == "\\" ) LBKSL=I;
    //OUTFILELBKSL+1:) = 'plotinfo.bin';

    file9.open("plotinfo");
    file10.open("plotinfo.bin");
    for (I = 1; I <= datax->NPTS; ++I) file10 << datax->BK[I];
    file9 << title << endl;
    file9 << "NO. OF PHASESZ  " << cntrls->NPHASE << endl
        << "NO. OF REFLECTIONS IN EACH PHASEQ  ";
    for (IIPHAS=1; IIPHAS <= cntrls->NPHASE; ++IIPHAS) file9 << setw(4) << refls->ICR[IIPHAS];
    file9 << endl;
    file9 << "BRAGG POSITIONSZ" << endl;
    for (K=1; K <= cntrls->NPHASE; ++K)
    {
        for (I=1; I <= datax->NPTS; ++I) datax->KR[I] = 0;
        if ((K % 2) == 1) allp->ILOC = allp->ILOC + 1;
        file6 << "PHASE NO. = " << K << "     PHASE NAME " << phases[K].name << endl;
        file9 << phases[K].name;
        T2OBS=0.0;
        TDOBS=0.0;
        RFNMR=0.0;
        RFDNR=0.0;
        ICZ=0;
        for (IIPHAS=1; IIPHAS <= cntrls->NPHASE; ++IIPHAS) ICZ = ICZ + refls->ICR[IIPHAS];
        IXX=0;
        IXXX=0;
        for (IX=1; IX <= ICZ; ++IX)
        {
            if( K != refls->IREFS[IX]/(256*256*256*8) )goto L481;
            SHIFT = params->GLB_[10-1] * cos(refls->REFS[IX][2]/2.0/57.2958) + params->GLB_[11-1] * sin(refls->REFS[IX][2]/57.2958);
            IXX=IXX+1;
            for (J=1; J <= jnk->NEXCRG; ++J)
            {
                //-----CHECK FOR SPLIT PEARSON VII PROFILE
                if (cntrls->NPROF == _SplitPearsonVII)
                {
                    if ((refls->REFS[IX][2]+SHIFT+params->GLB_[1-1]) >= (jnk->ALOW[J]-g1->WDT*refls->FWHM[IX][1]) && (refls->REFS[IX][2]+SHIFT+params->GLB_[1-1]) <= (jnk->AHIGH[J]+g1->WDT*refls->FWHM[IX][2]))  goto L481;
                    //-----FOR ALL OTHER PROFILES
                }
                else
                {
                    if ((refls->REFS[IX][2]+SHIFT+params->GLB_[1-1]) >= (jnk->ALOW[J]-g1->WDT*refls->REFS[IX][1]) && (refls->REFS[IX][2]+SHIFT+params->GLB_[1-1]) <= (jnk->AHIGH[J]+g1->WDT*refls->REFS[IX][1])) goto L481;
                }
            }
            IXXX=IXXX+1;                           //test !cp 29 jun 98
            IRL=(refls->IREFS[IX] % 256)-128;
            IRK=((refls->IREFS[IX]/256) % 256)-128;
            IRH=((refls->IREFS[IX]/(256*256)) % 256)-128;
            IRC=(refls->IREFS[IX]/(256*256*256)) % 8;
            //-----CHECK FOR THE SPLIT PEARSON VII PROFILE
            //-----if SO CHANGE THE PROFILE LIMITS
            if (cntrls->NPROF == _SplitPearsonVII)
            {
                RMIN=refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT-g1->WDT*refls->FWHM[IX][1];
                RMAX=refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT+g1->WDT*refls->FWHM[IX][2];
            }
            else
            {
                //-----FOR ALL OTHER PROFILES
                RMIN=refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT-g1->WDT*refls->REFS[IX][1];
                RMAX=refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT+g1->WDT*refls->REFS[IX][1];
            }
            MIN=static_cast<int>((RMIN-g1->THMIN)/g1->STEP+1.5);
            MAX=static_cast<int>((RMAX-g1->THMIN)/g1->STEP+1.5);
            MIN=max(MIN,1);
            MIN=min(MIN,datax->NPTS);
            MAX=min(MAX,datax->NPTS);
            //-----PATCH TO CALCULATE R-BRAGG
            prfx->TL=refls->REFS[IX][1];
            VERT=refls->REFS[IX][2] <= g1->RLIM;
            if (cntrls->NPROF == _pseudoVoigt)
            {
                prfx->GAM1=phases[K-1].PAR[16] + phases[K-1].PAR[17] * refls->REFS[IX][2];
            }
            else if (cntrls->NPROF == _PearsonVII)
            {
                prfx->GAM1=phases[K-1].PAR[16]+(phases[K-1].PAR[17]+phases[K-1].PAR[18]/refls->REFS[IX][2])/refls->REFS[IX][2];
                PRSVII(prfx->GAM1);
            }
            else if (cntrls->NPROF == _TCHZ)
            {
                prfx->GAM1=refls->GAM[IX];
                TLG=refls->HALFG[IX];
                TLL=refls->HALFL[IX];
            }
            else if (cntrls->NPROF == _SplitPearsonVII)
            {
                spvii->RL=phases[K-1].PAR[16]+(phases[K-1].PAR[17]+phases[K-1].PAR[18]/refls->REFS[IX][2])/refls->REFS[IX][2];
                spvii->RH=phases[K-1].PAR[23]+(phases[K-1].PAR[24]+phases[K-1].PAR[25]/refls->REFS[IX][2])/refls->REFS[IX][2];
                mspvii(phases[K-1].PAR[26],prfx->TL);
            }
            BB=prfx->TL*prfx->TL;
            TIOBS=0.0;
            TIC=0.0;
            refls->FMGNTD[IX] = refls->FMGNTD[IX]*refls->REFS[IX][3]*phases[K-1].PAR[0];
            for (IS=MIN; IS <= MAX; ++IS)
            {
                TH=g1->THMIN+static_cast<double>(IS-1)*g1->STEP;
                if(jnk->NEXCRG > 0)
                {
                    for (IEXC=1; IEXC <= jnk->NEXCRG; ++IEXC)
                    {
                        if (TH >= jnk->ALOW[IEXC] && TH <= jnk->AHIGH[IEXC]) goto L410;
                    }
                }
                TH2 = 0.5*TH/RAD;
                prfx->DELTA=TH-refls->REFS[IX][2]-params->GLB_[1-1]-SHIFT;
                prfx->DELT=prfx->DELTA*prfx->DELTA;
                //     NEXT LINE IS NECESSEARY FOR 2 PHASES WITH VERY DifFERENT FWHM.
                if (prfx->DELT/BB > g1->WDT*g1->WDT) goto L410;
                if( !VERT )goto L4;
                if (cntrls->IASYM == 0)
                {
                    YX=prfx->DELT * sign(prfx->DELTA);
                    Z=1.0-phases[K-1].PAR[13]*YX/tan(TH2);
                }
                else
                {
                    YX=sign(prfx->DELTA)*prfx->DELTA/(2*prfx->TL);
                    if (TH2 > (45.0/RAD)) TH2 = TH2-(90.0/RAD);
                    Z=(phases[K-1].PAR[13]/tan(TH2)) * (2.0*(prfx->DELTA/(2*prfx->TL))*exp(-YX));
                    Z=1+Z;
                }
                if(Z <= 0.0)Z=0.0001;
                goto L5;
L4:
                Z=1.0;
L5:
                OMEGA= Z*PROFIL(cntrls->NPROF,prfx->DELT/BB);
                if (cntrls->NPROF == _SplitPearsonVII)
                {
                    datax->KR[IS] = datax->KR[IS] + static_cast<int>(OMEGA*refls->FMGNTD[IX]);
                }
                else
                {
                    datax->KR[IS] = datax->KR[IS] + static_cast<int>(OMEGA*refls->FMGNTD[IX]/prfx->TL);
                }
                TIC   = TIC + OMEGA;
                TIOBS=TIOBS+OMEGA*(datax->Y[IS]-datax->BK[IS])/ max((datax->YC[IS]-datax->BK[IS]),1.0);
L410:;
            }
            TIC   = TIC   * refls->FMGNTD[IX]/prfx->TL;
            TIOBS = TIOBS * refls->FMGNTD[IX]/prfx->TL;
            // ! FROM ITALIAN CODE  !cp ap 97
            //        ***************************************************************
            //        NEXT LINE IS FOR NOT EVALUATING R-BRAGG FOR REFLECTIONS
            //        WHICH ARE OUTSIDE THE MEASUREMENTS RANGE BUT WHOSE TAILS ARE
            //        IN THE PATTERN
            //        ***************************************************************
            if((refls->REFS[IX][2]+SHIFT+params->GLB_[1-1]) <= g1->THMAX)
            {
                T2OBS = T2OBS+TIOBS;
                TDOBS = TDOBS + abs(TIOBS-TIC);
                if (cntrls->IPC == 2 || cntrls->IPC == 3)
                {
                    TFCAL= sqrt(abs(TIC/refls->REFS[IX][3]/phases[K-1].PAR[0]/struphase->TAVIX[IX]/struphase->SRIX[IX]));
                    TFOBS= sqrt(abs(TIOBS/refls->REFS[IX][3]/phases[K-1].PAR[0]/struphase->TAVIX[IX]/struphase->SRIX[IX]));
                    AFCAL=TFCAL*cos(struphase->APHASE[IX]);
                    BFCAL=TFCAL*sin(struphase->APHASE[IX]);
                    AFOBS=TFOBS*cos(struphase->APHASE[IX]);
                    BFOBS=TFOBS*sin(struphase->APHASE[IX]);
                    //        Line below also from italian code   !cp ap 97
                    if(refls->REFS[IX][2] <= g1->THMAX)
                    {
                        RFNMR = RFNMR + abs(TFOBS-TFCAL);
                        RFDNR = RFDNR + abs(TFOBS);
                    }
                }
            }
            if (cntrls->NPROF == _TCHZ) goto L9222;
            if (cntrls->NPROF == _SplitPearsonVII) goto L9221;
            if( (IXXX-1 % 60) == 0)
            {
                // alterei aqui para criar arquivo com dados para determinacao de estruturas e graficos ! cp 22abril01
                if(cntrls->IPC == 1)
                {
                    file6 << "NO.  CODE    H   K   L   HW     POSN      ICALC     IOBS" << endl;
                    if(cntrls->IPL2 != 0) file69 << "   H   K   L   HW     POSN      ICALC       IOBS   Bragg_Pos" << endl;
                }
                if(cntrls->IPC == 2)
                {
                    file6 << "NO.  CODE     H   K   L    FWHM    POSN    ICALC      IOBS      FCALC     FOBS   PHASE_C" << endl;
                    if(cntrls->IPL2 != 0) file69 << "   H   K   L     FWHM    POSN    ICALC    IOBS      FCALC      FOBS   PHASE_C  Bragg_Pos" << endl;
                }
                if(cntrls->IPC == 3)
                {
                    file6 << "NO.  CODE     H   K   L     FWHM    POSN    ICALC      IOBS     A_CALC    B_CALC     A_OBS     B_OBS" << endl;
                    if(cntrls->IPL2 != 0) file69 << "   H   K   L     FWHM    POSN    ICALC       IOBS     A_CAL  B_CALC      A_OBS     B_OBS  Bragg_Pos" << endl;
                }
            }
            HW=refls->REFS[IX][1];
            if (cntrls->IPC == 1)
            {
                file6 << setw(4) << IXX
                    << setw(4) << IRC	<< "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(8) << setprecision(3) << HW
                    << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                    << setw(10) << setprecision(0) << TIC
                    << setw(10) << setprecision(0) << TIOBS << endl;
                if(IRC == 1)
                {
                    if (cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << HW << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(0) << TIC << " "
                            << setw(8) << setprecision(0) << TIOBS << " "
                            << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
                    }
                }
            }
            else if (cntrls->IPC == 2)
            {
                if (IRC == 1)
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC	<< "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << HW
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(10) << setprecision(3) << TFCAL
                        << setw(10) << setprecision(3) << TFOBS
                        << setw(8) << setprecision(2) << RAD*struphase->APHASE[IX];
                    if (cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << HW << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(1) << TIC << " "
                            << setw(8) << setprecision(1) << TIOBS << " "
                            << setw(8) << setprecision(1) << TFCAL << " "
                            << setw(8) << setprecision(1) << TFOBS << " "
                            << setw(9) << setprecision(4) << RAD*struphase->APHASE[IX] << " "
                            << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
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
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS << endl;
                }
            }
            else if (cntrls->IPC == 3)
            {
                if (IRC == 1)
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC	<< "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << HW
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(10) << setprecision(3) << AFCAL
                        << setw(10) << setprecision(3) << BFCAL
                        << setw(10) << setprecision(3) << AFOBS
                        << setw(10) << setprecision(3) << BFOBS << endl;
                    if (cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << HW << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(1) << TIC << " "
                            << setw(8) << setprecision(1) << TIOBS << " "
                            << setw(8) << setprecision(1) << AFCAL << " "
                            << setw(8) << setprecision(1) << BFCAL << " "
                            << setw(9) << setprecision(3) << AFOBS << " "
                            << setw(9) << setprecision(3) << BFOBS << " "
                            << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
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
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS << endl;
                }
            }
            file9 << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                << " K" << setw(1) << IRC << " "
                << setw(4) << IRH
                << setw(4) << IRK
                << setw(4) << IRL
                << setw(10) << setprecision(0) << refls->FMGNTD[IX] << endl;
            goto L481;
L9221:
            if( (IXXX-1 % 60) == 0)
            {
                if(cntrls->IPC == 1)
                {
                    file6 << "NO.  CODE    H   K   L      HWL     HWH  POSN       ICALC     IOBS  PHASE" << endl;
                    if(cntrls->IPL2 != 0) file69 << "  H   K   L      HWL       HWH      POSN     ICALC      IOBS    PHASE  Bragg_Pos" << endl;
                }
                if(cntrls->IPC == 2)
                {
                    file6 << "NO.  CODE    H   K   L      HWL     HWH  POSN      ICALC     IOBS        FCALC    FOBS       PHASE_C" << endl;
                    if(cntrls->IPL2 != 0) file69 << "  H   K   L      HWL       HWH      POSN     ICALC      IOBS    FCALC     FOBS     PHASE_C  Bragg_Pos" << endl;
                }
                if(cntrls->IPC == 3)
                {
                    file6 << "NO.  CODE    H   K   L      HWL     HWH  POSN      ICALC     IOBS        A_CALC     B_CALC     A_OBS     B_OBS" << endl;
                    if(cntrls->IPL2 != 0) file69 << "  H   K   L      HWL      HWH      POSN     ICALC      IOBS   A_CALC   B_CALC   A_OBS    B_OBS  Bragg_Pos" << endl;
                }
            }
            if(cntrls->IPC == 1)
            {
                file6 << setw(4) << IXX
                    << setw(4) << IRC << "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(8) << setprecision(3) << refls->FWHM[IX][1]
                    << setw(8) << setprecision(3) << refls->FWHM[IX][2]
                    << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                    << setw(10) << setprecision(3) << TIC
                    << setw(10) << setprecision(3) << TIOBS << endl;
                if(IRC == 1)
                {
                    if (cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << refls->FWHM[IX][1] << " "
                            << setw(8) << setprecision(4) << refls->FWHM[IX][2] << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(0) << TIC << " "
                            << setw(8) << setprecision(0) << TIOBS << " "
                            << setw(8) << setprecision(4) << (-K*maxint->XMAXINT/10) << endl;
                    }
                }
            }
            else if(cntrls->IPC == 2)
            {
                if (IRC == 1)
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC << "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << refls->FWHM[IX][1]
                        << setw(8) << setprecision(3) << refls->FWHM[IX][2]
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(10) << setprecision(3) << TFCAL
                        << setw(10) << setprecision(3) << TFOBS << " "
                        << setw(8) << setprecision(2) << RAD*struphase->APHASE[IX] << endl;
                    if (cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << refls->FWHM[IX][1] << " "
                            << setw(8) << setprecision(4) << refls->FWHM[IX][2] << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(0) << TIC << " "
                            << setw(8) << setprecision(0) << TIOBS << " "
                            << setw(8) << setprecision(4) << TFCAL << " "
                            << setw(8) << setprecision(4) << TFOBS << " "
                            << setw(9) << setprecision(4) << RAD*struphase->APHASE[IX]
                            << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
                    }
                }
                else
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC << "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << refls->FWHM[IX][1]
                        << setw(8) << setprecision(3) << refls->FWHM[IX][2]
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC << endl;
                }
            }
            else if(cntrls->IPC == 3)
            {
                if (IRC == 1)
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC << "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << refls->FWHM[IX][1]
                        << setw(8) << setprecision(3) << refls->FWHM[IX][2]
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(10) << setprecision(3) << AFCAL
                        << setw(10) << setprecision(3) << BFCAL
                        << setw(10) << setprecision(3) << AFOBS
                        << setw(10) << setprecision(3) << BFOBS << endl;
                    if (cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << refls->FWHM[IX][1] << " "
                            << setw(8) << setprecision(4) << refls->FWHM[IX][2] << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(0) << TIC << " "
                            << setw(8) << setprecision(0) << TIOBS << " "
                            << setw(8) << setprecision(2) << AFCAL << " "
                            << setw(8) << setprecision(2) << BFCAL << " "
                            << setw(8) << setprecision(2) << AFOBS
                            << setw(8) << setprecision(2) << BFOBS
                            << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
                    }
                }
                else
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC << "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << refls->FWHM[IX][1]
                        << setw(8) << setprecision(3) << refls->FWHM[IX][2]
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS << endl;
                }
            }
            file9 << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " K"
                << setw(1) << IRC << " "
                << setw(4) << IRH
                << setw(4) << IRK
                << setw(4) << IRL
                << setw(10) << setprecision(4) << refls->FMGNTD[IX] << endl;
            goto L481;
L9222:
            if((IXXX-1 % 60) == 0)
            {
                if (cntrls->IPC == 1)
                {
                    file6 << "NO. CODE     H   K   L     HW     POSN     ICALC    IOBS       HG      HL      ETA   PHASE_C" << endl;
                    if(cntrls->IPL2 != 0) file69 << "  H   K   L     HW      POSN     ICALC    IOBS    HGHL      ETA   PHASE_C" << endl;
                }
                if (cntrls->IPC == 2)
                {
                    file6 << "NO. CODE     H   K   L     HW     POSN    ICALC      IOBS      HG      HL      ETA     FCALC     FOBS   PHASE_C" << endl;
                    if(cntrls->IPL2 != 0) file69 << "  H   K   L     HW       POSN    ICALC    IOBS    HG     HL      ETA   FCALC     FOBS   PHASE_C" << endl;
                }
                if (cntrls->IPC == 3)
                {
                    file6 << "NO. CODE     H   K   L     HW     POSN    ICALC      IOBS      HG      HL      ETA    A_CALC    B_CALC     A_OBS     B_OBS" << endl;
                    if(cntrls->IPL2 != 0) file69 << "  H   K   L     HW       POSN    ICALC    IOBS    HG     HL      ETA   A_CALC    B_CALC     A_OBS     B_OBS" << endl;
                }
            }
            HW=refls->REFS[IX][1];
            if (cntrls->IPC == 1)
            {
                file6 << setw(4) << IXX
                    << setw(4) << IRC << "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(8) << setprecision(3) << HW
                    << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                    << setw(10) << setprecision(0) << TIC
                    << setw(10) << setprecision(0) << TIOBS
                    << setw(8) << setprecision(3) << TLG
                    << setw(8) << setprecision(3) << TLL
                    << setw(8) << setprecision(3) << prfx->GAM1
                    << setw(10) << setprecision(3) << RAD*struphase->APHASE[IX] << endl;
                if(cntrls->IPL2 != 0)
                {
                    file69 << setw(3) << IRH << " "
                        << setw(3) << IRK << " "
                        << setw(3) << IRL << " "
                        << setw(8) << setprecision(4) << HW << " "
                        << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                        << setw(8) << setprecision(0) << TIC << " "
                        << setw(8) << setprecision(4) << TIOBS << " "
                        << setw(8) << setprecision(4) << TLG << " "
                        << setw(8) << setprecision(4) << TLL << " "
                        << setw(8) << setprecision(4) << prfx->GAM1 << " "
                        << setw(9) << setprecision(4) << RAD*struphase->APHASE[IX]
                        << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
                }
            }
            else if(cntrls->IPC == 2)
            {
                if (IRC == 1)
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC << "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << HW
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(8) << setprecision(3) << TLG
                        << setw(8) << setprecision(3) << TLL
                        << setw(8) << setprecision(3) << prfx->GAM1
                        << setw(10) << setprecision(3) << TFCAL
                        << setw(10) << setprecision(3) << TFOBS << " "
                        << setw(8) << setprecision(2) << struphase->APHASE[IX] << endl;
                    if(cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << HW << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(0) << TIC << " "
                            << setw(8) << setprecision(0) << TIOBS << " "
                            << setw(8) << setprecision(4) << TLG << " "
                            << setw(8) << setprecision(4) << TLL << " "
                            << setw(8) << setprecision(4) << prfx->GAM1 << " "
                            << setw(8) << setprecision(4) << TFCAL << " "
                            << setw(8) << setprecision(4) << TFOBS << " "
                            << setw(9) << setprecision(4) << struphase->APHASE[IX] << " "
                            << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
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
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(8) << setprecision(3) << TLG
                        << setw(8) << setprecision(3) << TLL
                        << setw(8) << setprecision(3) << prfx->GAM1 << endl;
                }
            }
            else if(cntrls->IPC == 3)
            {
                if (IRC == 1)
                {
                    file6 << setw(4) << IXX
                        << setw(4) << IRC << "   "
                        << setw(4) << IRH
                        << setw(4) << IRK
                        << setw(4) << IRL
                        << setw(8) << setprecision(3) << HW
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(8) << setprecision(3) << TLG
                        << setw(8) << setprecision(3) << TLL
                        << setw(8) << setprecision(3) << prfx->GAM1
                        << setw(10) << setprecision(4) << AFCAL
                        << setw(10) << setprecision(4) << BFCAL
                        << setw(10) << setprecision(4) << AFOBS
                        << setw(10) << setprecision(4) << BFOBS << endl;
                    if(cntrls->IPL2 != 0)
                    {
                        file69 << setw(3) << IRH << " "
                            << setw(3) << IRK << " "
                            << setw(3) << IRL << " "
                            << setw(8) << setprecision(4) << HW << " "
                            << setw(8) << setprecision(4) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT << " "
                            << setw(8) << setprecision(0) << TIC << " "
                            << setw(8) << setprecision(0) << TIOBS << " "
                            << setw(8) << setprecision(4) << TLG << " "
                            << setw(8) << setprecision(4) << TLL << " "
                            << setw(8) << setprecision(4) << prfx->GAM1 << " "
                            << setw(8) << setprecision(4) << AFCAL << " "
                            << setw(8) << setprecision(4) << BFCAL << " "
                            << setw(8) << setprecision(4) << AFOBS << " "
                            << setw(8) << setprecision(4) << BFOBS << " "
                            << setw(9) << setprecision(0) << (-K*maxint->XMAXINT/10) << endl;
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
                        << setw(8) << setprecision(3) << refls->REFS[IX][2]+params->GLB_[1-1]+SHIFT
                        << setw(10) << setprecision(0) << TIC
                        << setw(10) << setprecision(0) << TIOBS
                        << setw(8) << setprecision(3) << TLG
                        << setw(8) << setprecision(3) << TLL
                        << setw(8) << setprecision(3) << prfx->GAM1 << endl;
                }
            }
L481:;
        }
        // version II - format for compatibility with Sakthevil's PLOT program
        for (I=1; I <= datax->NPTS; ++I) file10 << datax->KR[I];
        TDOBS=100.0*TDOBS/T2OBS;
        file6 << "DERIVED BRAGG R-FACTOR = " << setw(8) << setprecision(2) << TDOBS << endl;
        if (cntrls->IPC != 1) file6 << "DERIVED R-F            = " << setw(8) << setprecision(2) << 100.0*RFNMR/RFDNR  << endl;
        allp->FINAL[allp->ILOC][2-(K % 2)] = TDOBS;
        if (sizestrain->NSIZESTRAIN  ==  9) size(K);
        //L400:
    }
L36:
    file9
        << " NPTSZ" << setw(5) << datax->NPTS << endl
        << " THMINZ" << setw(8) << setprecision(4) << g1->THMIN << endl
        << " STEPZ" << setw(8) << setprecision(4) << g1->STEP << endl
        << " YOBS    YCALCZ" << endl;
    for (I=1; I <= datax->NPTS; ++I)
    {
        datax->BK[I]=g1->THMIN+static_cast<double>(I-1)*g1->STEP;
        file9 << " "
            << setw(8) << setprecision(0) << datax->Y[I]
        << setw(8) << setprecision(0) << datax->YC[I] << endl;
        if(cntrls->IPL2 != 0)
        {
            file690 << " " << setw(9) << setprecision(4) << g1->THMIN+static_cast<double>(I-1)*g1->STEP
                << " " << setw(9) << setprecision(0) << datax->Y[I]
            << " " << setw(9) << setprecision(0) << datax->YC[I]
            << " " << setw(9) << setprecision(0) << datax->Y[I]-datax->YC[I]
            << " " << setw(9) << setprecision(0) << ( pow(datax->Y[I]-datax->YC[I],2) )/datax->Y[I]
            << endl;
        }
    }
    file9.close();
    file10.close();
    if(cntrls->IPL2 != 0)
    {
        file69.close();
        file690.close();
    }
    //     IT BUILDS THE OBSERVED DATA FILE CORRECTED FOR ABSORPTION
    if(cntrls->IPLOSS == 1 && cntrls->IPBIG == 0)
    {
        file31.open("PLOTOSS.COR");
        file31 << " \"" << setw(70) << title << setw(8) << TITLE1 << endl;
        file31 << setw(6) << datax->NPTS
            << setw(15) << setprecision(5) << g1->STEP
            << setw(15) << setprecision(5) << g1->THMIN
            << setw(15) << setprecision(5) << g1->THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= datax->NPTS; ++I)
        {
            file31 << setw(15) << setprecision(5) << datax->Y[I] << endl;
        }
        file31.close();
    }
    //     IT BUILDS THE CALCULATED DATA FILE
    //     (BRAGG+COMPTON+DISORDINE+AMORPHOUS)
    if(cntrls->IPLCAL == 1 && cntrls->IPBIG == 0)
    {
        file32.open("PLOTCAL.TOT");
        file32 << " \"" << setw(70) << title << setw(8) << TITLE2 << endl;
        file32 << setw(6) << datax->NPTS
            << setw(15) << setprecision(5) << g1->STEP
            << setw(15) << setprecision(5) << g1->THMIN
            << setw(15) << setprecision(5) << g1->THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= datax->NPTS; ++I)
        {
            file32 << setw(15) << setprecision(5) << datax->YC[I] << endl;
        }
        file32.close();
    }
    //     IT BUILDS THE TOTAL COMPTON SCATTERING FILE
    //     FOR ALL PHASES
    if(cntrls->IPLCOM == 1 && cntrls->IPBIG == 0)
    {
        file33.open("PLOTCOM.TOT");
        file33 << " \"" << setw(70) << title << setw(8) << TITLE3 << endl;
        file33 << setw(6) << datax->NPTS
            << setw(15) << setprecision(5) << g1->STEP
            << setw(15) << setprecision(5) << g1->THMIN
            << setw(15) << setprecision(5) << g1->THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= datax->NPTS; ++I)
        {
            file33 << setw(15) << setprecision(5) << fondi->BKCOM[I] << endl;
        }
        file33.close();
    }
    //     IT BUILDS THE TOTAL DISORDER SCATTERING FILE
    //     FOR ALL PHASES
    if(cntrls->IPLDIS == 1 && cntrls->IPBIG == 0)
    {
        file34.open("PLOTDIS.TOT");
        file34 << " \"" << setw(70) << title << setw(8) << TITLE4 << endl;
        file34 << setw(6) << datax->NPTS
            << setw(15) << setprecision(5) << g1->STEP
            << setw(15) << setprecision(5) << g1->THMIN
            << setw(15) << setprecision(5) << g1->THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= datax->NPTS; ++I)
        {
            file34 << setw(15) << setprecision(5) << fondi->BKDIS[I] << endl;
        }
        file34.close();
    }
    //     IT BUILDS THE POLYNOMIAL BACKGROUND FILE
    if(cntrls->IPLPOL == 1 && cntrls->IPBIG == 0)
    {
        file37.open("PLOTPOL.TOT");
        file37 << " \"" << setw(70) << title << setw(8) << TITLE7 << endl;
        file37 << setw(6) << datax->NPTS
            << setw(15) << setprecision(5) << g1->STEP
            << setw(15) << setprecision(5) << g1->THMIN
            << setw(15) << setprecision(5) << g1->THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= datax->NPTS; ++I)
        {
            file37 << setw(15) << setprecision(5) << fondi->BKPOL[I] << endl;
        }
        file37.close();
    }
    //     IT BUILDS THE AMORPHOUS FILE
    if(cntrls->IPLAM == 1 && cntrls->IPBIG == 0)
    {
        file36.open("PLOTAM.TOT");
        file36 << " \"" << setw(70) << title << setw(8) << TITLE6 << endl;
        file36 << setw(6) << datax->NPTS
            << setw(15) << setprecision(5) << g1->STEP
            << setw(15) << setprecision(5) << g1->THMIN
            << setw(15) << setprecision(5) << g1->THMAX
            << setw(5) << 1 << endl;
        for (I=1; I <= datax->NPTS; ++I)
        {
            file36 << setw(15) << setprecision(5) << fondi->BKAM[I] << endl;
        }
        file36.close();
    }
    //     IT BUILDS THE TOTAL PLOT FILE
    if(cntrls->IPBIG == 1)
    {
        file38.open("PLOTBIG.DAT");
        file38 << " \"" << setw(70) << title << setw(8) << TITLE8 << endl;

        file38 << "       ANG       OSS       CAL       +AM      +POL      +DIS      +COM       RES" << endl;
        for (I=1; I <= datax->NPTS; ++I)
        {
            XXXX=g1->THMIN+static_cast<double>(I-1)*g1->STEP;
            file38 << setw(10) << setprecision(3) << XXXX << ","
                << setw(10) << setprecision(3) << datax->Y[I] << ","
                << setw(10) << setprecision(3) << datax->YC[I] << ","
                << setw(10) << setprecision(3) << fondi->BKAM[I] << ","
                << setw(10) << setprecision(3) << fondi->BKAM[I]+fondi->BKPOL[I] << ","
                << setw(10) << setprecision(3) << fondi->BKAM[I]+fondi->BKPOL[I]+fondi->BKDIS[I] << ","
                << setw(10) << setprecision(3) << fondi->BKAM[I]+fondi->BKPOL[I]+fondi->BKDIS[I]+fondi->BKCOM[I] << ","
                << setw(10) << setprecision(3) << pow(datax->Y[I]-datax->YC[I],2)/datax->Y[I]
            << endl;
        }
        file38.close();
    }
    //L12045:
    if (cntrls->IOT == 0) goto L30;
    NP=datax->NPTS/240;
    //xxxxxxxxxxxxxxxxxxxxxxx
    for (N=1; N <= NP; ++N)
    {
        for (I=1; I <= 4; ++I) file6 << "2THETA   YOBS    YCALC VARIANCE ";
        file6 << endl;
        for (J=1; J <= 60; ++J)
        {
            for (I=1; I <= 4; ++I)
            {
                file6 << setw(7) << setprecision(3) << datax->BK[240*N-300+60*I+J]
                << setw(8) << setprecision(0) << datax->Y[240*N-300+60*I+J]
                << setw(8) << setprecision(0) << datax->YC[240*N-300+60*I+J]
                << setw(8) << setprecision(0) << datax->VAR[240*N-300+60*I+J];
            }
            file6 << endl;
        }
    }
    NPTS2=datax->NPTS-NP*240;
    if (NPTS2 == 0) goto L30;
    NCOL=NPTS2/60;
    for (I=1; I <= 4; ++I) file6 << "2THETA   YOBS    YCALC VARIANCE ";
    NLINES=NPTS2-NCOL*60;
    if (NLINES == 0) goto L33;
    NCOL1=NCOL+1;
    NP=240*NP-60;
    for (J=1; J <= NLINES; ++J)
    {
        for (I=1; I <= NCOL1; ++I)
        {
            file6 << setw(7) << setprecision(3) << datax->BK[NP+60*I+J]
            << setw(8) << setprecision(0) << datax->Y[NP+60*I+J]
            << setw(8) << setprecision(0) << datax->YC[NP+60*I+J]
            << setw(8) << setprecision(0) << datax->VAR[NP+60*I+J];
        }
        file6 << endl;
    }
    NP=NP+NLINES;
    NLINES=60-NLINES;
L33:
    for (J=1; J <= NLINES; ++J)
    {
        for (I=1; I <= NCOL1; ++I)
        {
            file6 << setw(7) << setprecision(3) << datax->BK[NP+60*I+J]
            << setw(8) << setprecision(0) << datax->Y[NP+60*I+J]
            << setw(8) << setprecision(0) << datax->YC[NP+60*I+J]
            << setw(8) << setprecision(0) << datax->VAR[NP+60*I+J];
        }
        file6 << endl;
    }
L30:
    if (cntrls->IPL == 0) goto L45;
    for (I=1; I <= 12; ++I) LABEL[I]=10*I*ISCALE;
    file6 << "                              " << title << endl << endl;
    for (I=1; I <= 12; ++I) file6 << setw(10) << LABEL[I];
    file6 << endl;
    for (I=1; I <= 12; ++I) LABEL[I]=-60*IDIF+10*I*IDIF;
    for (I=1; I <= 12; ++I) file6 << setw(10) << LABEL[I];
    for (I=1; I <= datax->NPTS; ++I)
    {
        //L1:
        IOUT = " ";
        for (J=1; J <= 120; ++J) IOUT=IOUT+ISPACE;
        IOUT[1]=IDOT;
        IOUT[61]=IDOT;
        IOUT[120]=IDOT;
        //L200:
        IY=static_cast<int>(datax->Y[I])/ISCALE+2;
        IY=max(min(IY,120),2);
        IOUT[IY-1]=IBAR;
        IOUT[IY]=IPLUS;
        if(IY <= 119) IOUT[IY+1]=IBAR;
        IY=static_cast<int>(datax->YC[I])/ISCALE+2;
        IY=max(min(IY,120),2);
        IOUT[IY]=IMINUS;
        IY=static_cast<int>(datax->Y[I]-datax->YC[I])/IDIF+61;
        IY=max(min(IY,120),2);
        IOUT[IY]=ISTAR;
        if((I-1 % 10) != 0) goto L55;
        TLABEL=g1->THMIN+g1->STEP*static_cast<double>(I-1);
        for (J=1; J <= 113; ++J) file6 << IOUT[J];
        file6 << setw(6) << setprecision(2) << TLABEL;
        file6 << IOUT[120] << endl;
        goto L60;
L55:
        for (J=1; J <= 120; ++J) file6 << IOUT[J];
        file6 << endl;
L60:;
    }
L45:
    if (cntrls->JOBTYP < 3) goto L72;

    file4.close();
    file4b.open(file4name.data(),ios::trunc);			//  REWIND 4;
    file4b << setw(8) << setprecision(4) << g1->THMIN
        << setw(8) << setprecision(4) << g1->STEP
        << setw(8) << setprecision(4) << g1->THMAX
        << endl;
    J = 0;
    for (I=1; I <= datax->NPTS; ++I)
    {
        ++J;
        file4b << setw(7) << setprecision(0) << " " << datax->YC[I];
        if (J == 8)
        {
            file4b << endl;
            J = 0;
        }
    }
    if (file4b.is_open()) file4b.close();
L72:
    if (cntrls->IPLST != 0 && cntrls->MAXS != 0)
    {
        INUMB = cntrls->ICYRUN + 1;
        NPAGES = cntrls->MAXS/12;
        if ((cntrls->MAXS % 12) != 0) NPAGES = NPAGES+1;
        file6 << "PARAMETERS IN EACH CYCLE" << endl;

        if (file8o.is_open()) file8o.close();

        //     LIST PARAMETERS IN EACH CYCLE
        for (J=1; J <= NPAGES; ++J)
        {
            file8i.open(file8name.data());			// REWIND 8;
            ISTART = 1 + (J-1)*12;
            IFINIS  = min(cntrls->MAXS,12 + (J-1)*12);
            file6 << "CYCLE";
            for (L=ISTART; L <= IFINIS; ++L) file6 << "   " << setw(2) << L << "     ";
            file6 << endl;
            for (K=1; K <= INUMB; ++K)
            {
                for (I=1; I <= 2*cntrls->MAXS+4; ++I) file8i >> DUMMY[I];
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
            for (I=1; I <= 2*cntrls->MAXS+4; ++I) file8i >> DUMMY[I];
            file6 << setw(2) << K-1 << ")  ";
            for (I=2*cntrls->MAXS+1; I <= 2*cntrls->MAXS+4; ++I) file6 << setw(8) << setprecision(2) << DUMMY[I] << " ";
            file6 << endl;
        }
        file8i.close();

        //C TESTE PARA GERAR SAIDA PARA 9411
        if(cntrls->I2D94 != 0)
        {
            WRITE94(ISCALE,IDIF);
        }
        //     LIST PARAMETER SHIFTS IN EACH CYCLE
        if (INUMB == 1) goto L88880;
        file6 << "APPLIED PARAMETER SHIFT IN EACH CYCLE" << endl;
        for (J=1; J <= NPAGES; ++J)
        {
            file8i.open(file8name.data());			// REWIND 8;
            ISTART = 1 + cntrls->MAXS+ (J-1)*12;
            IFINIS = min(2*cntrls->MAXS,12 + cntrls->MAXS+ (J-1)*12);
            file6 << "CYCLE";
            for (L=ISTART-cntrls->MAXS; L <= IFINIS-cntrls->MAXS; ++L) file6 << "   " << setw(2) << "     ";
            file6 << endl;
            for (K=1; K <= INUMB-1; ++K)
            {
                for (I=1; I <= 2*cntrls->MAXS+4; ++I) file8i >> DUMMY[I];
                if (K == 1) goto L7964;
L7964:
                file6 << setw(2) << K << ")  ";
                file6 << scientific;
                for (I=ISTART; I <= IFINIS; ++I) file6 << setw(10) << setprecision(4) << DUMMY[I] << " ";
                file6 << fixed << endl;
            }
        }
        file8i.close();

        file8o.open(file8name.data());

    }
L88880:
    //     CODE FOR PRINTING PARAMETERS AND STD. DEV. IN THE FINAL CYCLE
    if (cntrls->IPLST == 2 && cntrls->MAXS != 0)
    {
        file6 << "PARAMETERS AND STANDARD DEVIATIONS IN THE FINAL CYCLE FOR DATA BASE" << endl
            << " \" " << title << " \" " << endl;
        file6 << scientific << setw(10) << setprecision(4);
        for (I=1; I <= allp->ILOC; ++I)
        {
            file6 << allp->FINAL[I][1] << allp->FINAL[I][2];
        }
        file6 << fixed << endl;
    }
    return;
    //L99990:
    //	DBWSException("ERROR IN WRITING TO FILE 6 IN SUBROUTINE EXPUT");
    //L99991:
    //	DBWSException("ERROR WRITING PARAMS,ST.DEV. IN SUBROUTINE EXPUT");
    //L99992:
    //	DBWSException("ERROR IN READING FROM UNIT 8 IN SUBROUTINE EXPUT");
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

    RAD = 3.14159265359/360.0;
    if ( phases[IPHASE].SYMB.NAXIS > 3 ) DBWSException("5001");
    KXIS = phases[IPHASE].SYMB.NAXIS;
    phases[IPHASE].SYMB.NAXIS = 1;
    //L1002:
    //J = 0;
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
        if ( phases[IPHASE].SYMB.NCONT[I] == 0 ) goto L9000;
        if ( phases[IPHASE].SYMB.NCONT[I]-8192 < 0 )
        {
            goto L8002;
        }
        else if ( phases[IPHASE].SYMB.NCONT[I]-8192 == 0 )
        {
            goto L8001;
        }
        else
        {
            goto L8000;
        }
        L8001:
        I3DEL = 3;
        I1D = 2;
        I123D = 3;
        goto L9000;
        L8002:
        if ( phases[IPHASE].SYMB.NCONT[I]-576 == 0 ) goto L8003; else goto L8004;
        L8003:
        I3DEL = 2;
        I23D = 2;
        goto L8000;
        L8004:
        if ( phases[IPHASE].SYMB.NCONT[I]-516 == 0 ) goto L8009; else goto L8008;
        L8009:
        I3DEL = 2;
        I13D = 2;
        goto L8000;
        L8008:
        if ( phases[IPHASE].SYMB.NCONT[I]-68 == 0 ) goto L8013; else goto L8012;
        L8012:
        if ( phases[IPHASE].SYMB.NCONT[I]-580  == 0) goto L8016; else goto L8015;
        L8013:
        I2DEL = 2;
        I12D = 2;
        goto L8000;
        L8016:
        I123D = 2;
        I3DEL = 2;
        L8015:;
        L8000:;
    }
    L9000:
    //N = J-1;
    //L97:
    TMIN = pow(( sin( (g1->THMIN/720.0)*6.28318531 )/g1->LAMDA[1]),2);
    IC = 0;
    if ( IPHASE >= 2 )
    {
        for (IIPHAS=2; IIPHAS <= IPHASE; ++IIPHAS)
        {
            IC = IC+refls->ICR[IIPHAS-1];
        }
    }
    SMAX = pow(( sin((g1->THMAX/720.0)*6.28318531) /g1->LAMDA[1]),2);

    //     **************************************************************
    //     SMAX IS CHANGED TO ACCOUNT FOR THE LEFT TAILS OF THE REFLECTIONS
    //     THAT ARE PRESENT AT ANGLES GREATER THAN THMAX
    //     **************************************************************
    //    PRINT *,'THMAX=',THMAX
    //      THMAX1=(U*(TAN(THMAX*RAD))**2+V*TAN(THMAX*RAD)+W+ZZZ*(1+
    //     *(TAN(THMAX*RAD))**2))
    // also incorporating the cotg^2 term  !cp may 01 97

    THMAX1 = (g1->U * pow( (tan(g1->THMAX*RAD)),2) +
              g1->V * tan(g1->THMAX*RAD)+
              g1->W+
              g1->ZZZ*
              (1+ pow( (tan(g1->THMAX*RAD)),2) ) +
              g1->UC/(tan(g1->THMAX*RAD)));
    if ( THMAX1 > 0.0 )
    {
        THMAX1= g1->WDT*sqrt(THMAX1)/2;
    }
    else
    {
        file6 << "  SQUARE OF FWHM NEGATIVE AT TWO-THETA" << setw(8) << setprecision(3) << POS << " PHASE NO. " << setw(4) << IPHASE << endl;
        DBWSException("SQUARE OF FWHM IS NEGATIVE");
    }
    ANGTTETA=g1->THMAX+THMAX1;
    if ( (g1->THMAX+THMAX1) >= 180.0 ) ANGTTETA=180.0;

    // TODO: mudar para ANGTTETA/2 * 2*pi/360
    SMAX1 = pow(( sin((ANGTTETA/720.0)*6.28318531)/g1->LAMDA[1]),2);
    SMAX = SMAX1;
    THMAXX=ANGTTETA;
    if ( g1->LAMDA[2] > g1->LAMDA[1] ) TMIN=TMIN* pow( (g1->LAMDA[1]/g1->LAMDA[2]) , 2);
    if ( g1->LAMDA[2] > 0.0  &&  g1->LAMDA[2] < g1->LAMDA[1] )SMAX=SMAX* pow((g1->LAMDA[1]/g1->LAMDA[2]) ,2);
    OP1(&IPHASE,NCTR_);
    file6
        << "LAUE SYMMETRY " << LAU[phases[IPHASE].SYMB.NSPGRP] << " WILL BE USED TO GENERATE INDICES" << endl
        << "       -------------------------------------------" << endl;
    I1MAX=static_cast<int>(cellx->A*2.0*sqrt(SMAX));
    if ( KXIS < 3 )
    {
        I2MAX=static_cast<int>(cellx->B*2.0*sqrt(SMAX));
        I3MAX=static_cast<int>(cellx->C*2.0*sqrt(SMAX));
    }
    else
    {
        I2MAX = static_cast<int>(2.0*cellx->C*sqrt(SMAX));
        I3MAX = static_cast<int>(2.0*cellx->B*sqrt(SMAX));
    }
    L1 = phases[IPHASE].SYMB.NSPGRP;
    if ( L1 >= 13 && phases[IPHASE].SYMB.NCONT[1] == 580 ) L1=L1+2;
    if ( I12D+I13D+I23D != 5 ) goto L2131;
    I12D = 2;
    I13D = 1;
    I23D = 2;
    I2DEL = 2;
    I3DEL = 2;
L2131:
    if ( L1 == 6 ) I2D=1;
    if ( L1 == 7 ) I2D=2;
    I3D = I2D;
    I1 = -1;
L2200:
    I1 = 1+I1;
    if ( I1-I1MAX <= 0 ) goto L2201; else goto L2303;
L2201:
    H1 = I1;
    if ( I2D > 0 ) I2MAX=I1;
    switch (L1) {
    case 1:
        goto L2203;
        break;
    case 2:
        goto L2204;
        break;
    case 3:
        goto L2204;
        break;
    case 4:
        goto L2205;
        break;
    case 5:
        goto L2206;
        break;
    case 6:
        goto L2208;
        break;
    case 7:
        goto L2209;
        break;
    case 8:
        goto L2207;
        break;
    case 9:
        goto L2204;
        break;
    case 10:
        goto L2206;
        break;
    case 11:
        goto L2204;
        break;
    case 12:
        goto L2206;
        break;
    case 13:
        goto L2206;
        break;
    case 14:
        goto L2206;
        break;
    case 15:
        goto L2206;
        break;
    case 16:
        goto L2206;
        break;
    }
    GOTOER();
L2203:
    I2= -min(I2MAX,I1*I2MAX);
    goto L2210;
L2204:
    I2 = 0;
    goto L2210;
L2205:
    I2= min(I1,1);
    goto L2210;
L2206:
    I2 = I1;
    goto L2210;
L2207:
    I2= min(-I1+1,0);
    goto L2210;
L2208:
    I2 = -2*I1;
    goto L2210;
L2209:
    I2 = -I1/2;
L2210:
    I2= I2DEL*(I2/I2DEL)+(I1 % I12D);
    goto L2221;
L2220:
    I2 = I2DEL+I2;
L2221:
    H2 = I2;
    if ( I2-I2MAX <= 0 ) goto L2222; else goto L2200;
L2222:
    switch (L1) {
    case 1:
        goto L2223;
        break;
    case 2:
        goto L2223;
        break;
    case 3:
        goto L2224;
        break;
    case 4:
        goto L2224;
        break;
    case 5:
        goto L2224;
        break;
    case 6:
        goto L2225;
        break;
    case 7:
        goto L2225;
        break;
    case 8:
        goto L2224;
        break;
    case 9:
        goto L2224;
        break;
    case 10:
        goto L2223;
        break;
    case 11:
        goto L2224;
        break;
    case 12:
        goto L2224;
        break;
    case 13:
        goto L2226;
        break;
    case 14:
        goto L2227;
        break;
    case 15:
        goto L2236;
        break;
    case 16:
        goto L2235;
        break;
    }
    GOTOER();
L2223:
    I3= -min(I3MAX,(I1+abs(I2))*I3MAX); //xxxxxxxxxxxxxxxxxxxxxxxxx
    goto L2228;
L2224:
    I3 = 0;
    goto L2228;
L2225:
    I3 = -I1-I2;
    goto L2228;
L2226:  I3= min(I2,I1+I3DEL);
    goto L2228;
L2227:
    I3 = I2;
L2228:  I3= I3DEL*(I3/I3DEL)+  ((((I1+I1D*I2) % I123D)+I123D) % I123D) + (I1 % I13D)+ (I2 % I23D);
    if ( I3D-1 < 0 )
    {
        goto L2232;
    }
    else if ( I3D-1 == 0 )
    {
        goto L2229;
    }
    else
    {
        goto L2231;
    }
L2229:
    I3MAX = I1;
    if ( I2-I1 == 0 ) goto L2232; else goto L2230;
L2230:
    I3MAX = I3MAX-1;
    goto L2232;
L2231:
    I3MAX = I2;
    goto L2232;
L2235:
    I3= I2+(I1 % 2);
    goto L2232;
L2236:
    I3 = I1+2-(I2 % 2);
    if ( I1 == I2 && (I1 % 2) == 0 ) I3=I1;
    goto L2232;
L2233:
    I3 = I3DEL+I3;
L2232:
    H3 = I3;
    if ( KXIS != 3 ) goto L2240;
    H3 = I2;
    H2 = I3;
L2240:;
    if ( I3-I3MAX <= 0 ) goto L113; else goto L2220;
L113:
    SQ = static_cast<double>(H1*H1)*cellx->AL[1][1]+
         static_cast<double>(H2*H2)*cellx->AL[2][2]+
         static_cast<double>(H3*H3)*cellx->AL[3][3]+
         static_cast<double>(H1*H2)*cellx->AL[1][2]+
         static_cast<double>(H1*H3)*cellx->AL[1][3]+
         static_cast<double>(H2*H3)*cellx->AL[2][3];
    SQ=SQ/4.0;
    //L2234:
    if ( SQ-SMAX <= 0 ) goto L3000; else goto L2233;
L2303:
    refls->ICR[IPHASE]=IC;
    if (IPHASE >= 2)
    {
        for (IIPHAS=2; IIPHAS <= IPHASE; ++IIPHAS) refls->ICR[IPHASE] = refls->ICR[IPHASE]-refls->ICR[IIPHAS-1];
    }
    SORT(IPHASE);
    return;
L3000:
    if ( SQ - TMIN < 0 ) goto L2233; else goto L3117;
L3117:
    if ( IC > IRS-2 ) goto L6001;
    //     NEXT if BLOCK FOR PHASE WITH DISTINCT ORIENTATION ALONG PREF(NAXIS,I)
    if (PREFOR > 99.0)
    {
        ORH1 = static_cast<int>(phases[IPHASE].PREF[1]) == H1;
        if ( !(ORH1) && H1 != 0 && static_cast<int>(phases[IPHASE].PREF[1]) != 0) ORH1 = (H1 % static_cast<int>(phases[IPHASE].PREF[1])) == 0;
        if ( !ORH1 ) goto L2233;

        ORH2 = static_cast<int>(phases[IPHASE].PREF[2]) == H2;
        if ( !(ORH2) && H2 != 0 && static_cast<int>(phases[IPHASE].PREF[2]) != 0) ORH2 = (H2 % static_cast<int>(phases[IPHASE].PREF[2])) == 0;
        if ( !ORH2 ) goto L2233;

        ORH3 = static_cast<int>(phases[IPHASE].PREF[3]) == H3;
        if ( !(ORH3) && H3 != 0 && static_cast<int>(phases[IPHASE].PREF[3]) != 0)  ORH3 = (H3 % static_cast<int>(phases[IPHASE].PREF[3])) == 0;
        if ( !ORH3 ) goto L2233;

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
    hklctl->IHKL[1][1]=H1;
    hklctl->IHKL[2][1]=H2;
    hklctl->IHKL[3][1]=H3;
    hklctl->AZ[1]=0.0;
    hklctl->IER=0;
    SMTRY2(&IPHASE);
    if(hklctl->IER != 0)goto L2233;
    LXN=2;
    if(g1->LAMDA[2] == 0.0 || g1->LAMDA[2] == g1->LAMDA[1])LXN=1;
    for (LX1=1; LX1 <= LXN; ++LX1)
    {
        SQH=SQ*g1->LAMDA[LX1]*g1->LAMDA[LX1];
        TAN2=SQH/(1.0-SQH);
        if (SQH >= 1.0) goto L3118;
        TANX=sqrt(TAN2);
        //     SHIFT DUE TO SAMPLE DISPLACEMENT AND TRANSPARENCY
        SHIFT =  DIS*sqrt(1-SQH)+TRANS*sqrt(1.-(1.-2.*SQH)*(1.-2.*SQH));
        POS=atan(TANX)/RAD+ ZERO + SHIFT;
        if ( POS > THMAXX  ||  POS < g1->THMIN ) goto L3118;
        IC=IC+1;
        XABS=1.0;
        if(g1->TMV <= 0.000001)XABS=1.0;
        PLOR=1.0/(2.0*SQH*sqrt(1.0-SQH))*XABS;
        if (cntrls->INSTRM == 2)
        {
            PLOR = PLOR * (0.95+0.05*(1.-2.*SQH)*(1.-2.*SQH));
            goto L4000;
        }
        if(cntrls->JOBTYP == 3)PLOR=PLOR*(1.+(1.-2.*SQH)*(1.-2.*SQH)*g1->CTHM);
        if(cntrls->JOBTYP == 1)PLOR=PLOR*(1.+(1.-2.*SQH)*(1.-2.*SQH)*g1->CTHM);
L4000:
        refls->IREFS[IC]=256*(256*(256*(8*IPHASE+LX1)+128+H1)+128+H2)+128+H3;
        refls->FMGNTD[IC]=MULT(IPHASE,H1,H2,H3,KXIS);

        //------CALCULATE FWHM FOR PSEUDOVOIGT WITH GAUSS AND LORENTZ
        if (cntrls->NPROF == _TCHZ)
        {
            refls->HALFG[IC] = (g1->U*TAN2+g1->V*TANX+g1->W+g1->ZZZ*(1+TAN2));
            if (refls->HALFG[IC] > 0.)
            {
                refls->HALFG[IC] = sqrt(refls->HALFG[IC]);
            }
            else
            {
                cout << "3" << endl;
                file6 << "  SQUARE OF FWHM NEGATIVE AT TWO-THETA" << setw(8) << setprecision(3) << POS << "FOR PHASE NO. " << setw(IPHASE) << endl;
                DBWSException("SQUARE OF FWHM IS NEGATIVE");
            }
            refls->HALFL[IC] = g1->ULOR*TANX+g1->VLOR/sqrt(1.0-SQH);
            refls->REFS[IC][1] = pow(( pow(refls->HALFG[IC],5.0) +2.69269*pow(refls->HALFG[IC],4.0)*refls->HALFL[IC]+2.42843*pow(refls->HALFG[IC],3.0)*pow(refls->HALFL[IC],2.0)+4.47163*pow(refls->HALFG[IC],2.0)*pow(refls->HALFL[IC],3.0)+0.07842*refls->HALFG[IC]*pow(refls->HALFL[IC],4.0)+pow(refls->HALFL[IC],5.0)),0.2);
            TLR = refls->HALFL[IC]/refls->REFS[IC][1];
            refls->GAM[IC] = 1.36603*TLR-0.47719*TLR*TLR+0.11116* pow(TLR,3.0);
        }
        else if (cntrls->NPROF == _SplitPearsonVII)
        {
            refls->REFS[IC][1] = (g1->U*TAN2+g1->V*TANX+g1->W);
        }
        else
        {
            // incorporating cotg^2  !cp may 01 97
            refls->REFS[IC][1]=(g1->U*TAN2+g1->V*TANX+g1->W+g1->ZZZ*(1+TAN2)+g1->UC/TAN2);
        }
        if (refls->REFS[IC][1] > 0.0)
        {
            refls->REFS[IC][1]= sqrt(refls->REFS[IC][1]);
        }
        else
        {
            cout << "4" << endl;
            file6 << "  SQUARE OF FWHM NEGATIVE AT TWO-THETA" << setw(8) << setprecision(3) << POS << "FOR PHASE NO. " << setw(IPHASE) << endl;
            DBWSException("SQUARE OF FWHM IS NEGATIVE");
        }
        //L7000:
        refls->REFS[IC][2]=POS;
        refls->REFS[IC][3]=PLOR;
L3118:;
    }
    goto L2233;
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
        NS = dc->NSAVE;
    }


    for (I = 1; I <= phases[K].AtomCount; ++I)
    {
        for (J = 1; J <= NSCAT; ++J) if(parac->NTYP[I+IOF] == coefc->NAM[J]) goto L30;
        //-----212 = ALL THE POSSIBLE NAMES OF ATOMS AND IONS
        //L25:
        for (J = 1; J <= 212; ++J) if(parac->NTYP[I+IOF] == TBXC[J]) goto L50;
        file6 << "COMPTON SCATTERING COEFFICIENT NOT FOUND FOR " << parac->NTYP[I+IOF] << endl;
        DBWSException("COMPTON SCATTERING DATA MISSING");
L30:
        comp->PTC[I+IOF] = J;
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
        comp->PTC[I+IOF] = NS;
        //-----PUT IN CC THE 4 COMPTON COEFFICIENTS
        for (L = 1; L <= 4; ++L) comp->CC[L][NS] = TCS[NA][L];
        coefc->NAM[NS] = NOME;

        NOM = string(NOME,0,4);
        //     FIND IN WHAT COLUMN THERE IS + OR -
        for (L=1; L <= 4; ++L)
        {
            if(NOM[L] == PIU)  goto L90;
            if(NOM[L] == MENO) goto L95;
        }
        // CASE WITH NOR + NOR -
        comp->ZEFF[NS] = NA;
        goto L9999;
        //-----CASE WITH PLUS
L90:
        L = L + 1;
        if (L > 4) DBWSException("SOMETHING IS WRONG IN ATOMIC NAME");
        NN = NOM[L];
        comp->ZEFF[NS] = NA - NN;
        goto L9999;
        //-----CASE WITH MINUS
L95:
        L = L + 1;
        if (L > 4) DBWSException("SOMETHING IS WRONG IN ATOMIC NAME");
        NN = NOM[L];
        comp->ZEFF[NS] = NA + NN;
L9999:;
    }
    //NSAVE = NS;
    dc->NSAVE = NS;
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


// SUBROUTINE GSASREAD * READ GSAS FORMATTED DATA FILE
void DBWS::GSASREAD(void)
{
    int I, J, K, L4, L5, L10, L11, L12, L13;
    int LAB[20+1];
    int NCTR[IDSZ+1];
    string s,TEST,DATAID;

    getline(file4,s);
    DATAID = s.substr(0,72);
    J=0;
    K=1;
L10:
    getline(file4,s);
    TEST = " " + s.substr(0,72);
    if( TEST.substr(1,4) != "BANK") goto L10;
L31:
    for (I=K; I <= 72; ++I)
    {
        if( TEST[I] == ' ')
        {
            J=J+1;
            LAB[J]=I-1;
            goto L21;
        }
    }

L21:
    for (K=I; K <= 66; ++K)
    {
        if(TEST[K] != ' ')
        {
            J=J+1;
            LAB[J]=K;
            goto L31;
        }
    }
    L4=LAB[4];
    L5=LAB[5];
    L10=LAB[10];
    L11=LAB[11];
    L12=LAB[12];
    L13=LAB[13];
    stringstream(TEST.substr(L4,L5-L4+1)) >> datax->NPTS;			// npts
    stringstream(TEST.substr(L10,L11-L10+1)) >> g1->THMIN;			// start*100
    stringstream(TEST.substr(L12,L13-L12+1)) >> g1->STEP;			// step*100
    g1->THMIN=g1->THMIN/100;
    g1->STEP=g1->STEP/100;
    g1->THMAX =g1->THMIN+datax->NPTS*g1->STEP;
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    file6 << "    DATA ID " << setw(56) << DATAID << endl;
    file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;

    if (datax->NPTS > IDSZ)
    {
        file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= datax->NPTS; ++I)
    {
        file4 >> NCTR[I] >> datax->Y[I];			//format(10(i2,f6.0))
    }
    for (I=1; I <= datax->NPTS; ++I)
    {
        if (NCTR[I] == 0) NCTR[I]=1;
        if (datax->Y[I] == 0.0) datax->Y[I]=0.0;		// TODO: ???????
    }
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
}

// SUBROUTINE TO READ PHILIPS UDF DATA FILE
void DBWS::PHILIPSREAD(void)
{
    int I,LB1,LB2,LB3;
    string s,TEST,DATAID,ATHMIN,ATHMAX;

    getline(file4,s);
    DATAID = s.substr(0,56);
L10:
    getline(file4,s);
    TEST = " " + s.substr(0,56);
    for (I=1; I <= 56; ++I)
    {
        if(TEST[I] == ',')
        {
            LB1=I;
            goto L31;
        }
    }
L31:
    if( TEST.substr(1,LB1-1) != "DataAngleRange") goto L10;
    //                            DataAngleRange
    for (I=LB1+1; I <= 56; ++I)
    {
        if (TEST[I] == ',')
        {
            LB2=I;
            goto L41;
        }
    }

L41:
    for (I=LB2+1; I <= 56; ++I)
    {
        if(TEST[I] == ',')
        {
            LB3=I;
            goto L51;
        }
    }

L51:
    ATHMIN=TEST.substr(LB1+1,LB2-1 - (LB1+1)+1);
    ATHMAX=TEST.substr(LB2+1,LB3-1 - (LB2+1)+1);
    getline(file4,s);
    TEST = s.substr(0,13);
    stringstream(s.substr(13,8)) >> g1->STEP;
L30:
    getline(file4,s);
    TEST = " " + s.substr(0,56);
    for (I=1; I <= 56; ++I)
    {
        if( TEST[I] == ',') LB1=I;
    }

    if( TEST.substr(1,LB1-1-1+1) != "RawScan") goto L30;
    //                            RawScan
    stringstream(ATHMIN) >> g1->THMIN;
    stringstream(ATHMAX) >> g1->THMAX;

    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    file6 << "    DATA ID " << setw(56) << DATAID << endl;
    file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    datax->NPTS=static_cast<int>((g1->THMAX-g1->THMIN)/g1->STEP+1.5);
    if (datax->NPTS > IDSZ)
    {
        file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= datax->NPTS; ++I) file4 >> datax->Y[I];
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
}

// suboutine (QPAINIT) to compute mass fractions before starting the refinement
void DBWS::qpainit(void)
{
    int I,N,IP,IOF,ICOCO,IIPHAS,IINNMOL;
    double V0,FT, XFAC, SFIS, ARGCOS, WTOTAL;
    //VOSQ
    //ARG1, ARG2, ARG3,
    double VOL[99+1];
    double W[99+1];
    double XMASS[99+1];
    double FR[99+1];
    double FRP[99+1];

    file6 << "       >>> QPA before starting the refinement <<<" << endl
        << "       >>                                    <<<" << endl;
    for (IP=1; IP <= cntrls->NPHASE; ++IP)
    {
        DIRECT(dircv->DCSM,dircv->DCV,&IP);
        multip->TMASSA[IP]=0.0;
        IOF=0;
        if(IP > 1)
        {
            for (IIPHAS=2; IIPHAS <= IP; ++IIPHAS) IOF = IOF + phases[IIPHAS-1].AtomCount;
        }
        N=phases[IP].AtomCount;
        for (I=1; I <= N; ++I)
        {
            ICOCO=params->PTR[I+IOF];
            multip->TMASSA[IP] = multip->TMASSA[IP] + params->XL[I+IOF][5]*coeff->XMAS[ICOCO]*multip->XMLTP[IP];
        }
        XFAC = 3.141592654 / 180.000000;
        dircv->DCV[4] = XFAC * dircv->DCV[4];
        dircv->DCSM[4][4] = dircv->DCSM[4][4] * XFAC;
        dircv->DCV[5] = XFAC * dircv->DCV[5];
        dircv->DCSM[5][5] = dircv->DCSM[5][5] * XFAC;
        dircv->DCV[6] = dircv->DCV[6] * XFAC;
        dircv->DCSM[6][6] = dircv->DCSM[6][6] * XFAC;

        //-----Calculations of VOLUME and SVZM (=W) for each phase
        ARGCOS= 1-pow((cos(dircv->DCV[4])),2)-pow((cos(dircv->DCV[5])),2)-pow((cos(dircv->DCV[6])),2) + 2 * (cos(dircv->DCV[4])) * (cos(dircv->DCV[5])) * (cos(dircv->DCV[6]));
        V0 = dircv->DCV[1] * dircv->DCV[2] * dircv->DCV[3];
        VOL[IP] = V0 * sqrt(ARGCOS);
        //		VOSQ = 0.5*VOL[IP]/ARGCOS;
        //		ARG1 = VOSQ*(2 * cos(dircv->DCV[4]) * sin(dircv->DCV[4]) - 2*sin(dircv->DCV[4]) *cos(dircv->DCV[5]) *cos(dircv->DCV[6])) * dircv->DCSM[4][4];
        //		ARG2 = VOSQ*(2 * cos(dircv->DCV[5]) * sin(dircv->DCV[5]) - 2*sin(dircv->DCV[5]) *cos(dircv->DCV[4]) *cos(dircv->DCV[6])) * dircv->DCSM[5][5];
        //		ARG3 = VOSQ*(2 * cos(dircv->DCV[6]) * sin(dircv->DCV[6]) - 2*sin(dircv->DCV[6]) *cos(dircv->DCV[4]) *cos(dircv->DCV[5])) * dircv->DCSM[6][6];
        W[IP] = phases[IP-1].PAR[0] * multip->TMASSA[IP] * VOL[IP]/phases[IP].SAQF;
        file6 << "Volume("
            << setw(2) << IP << ")= "
            << setw(9) << setprecision(3) << VOL[IP]
        << " UCW= "
            << setw(7) << setprecision(2) << multip->TMASSA[IP]
        << " U.C.Density = "
            << setw(7) << setprecision(3) << 1.66113*multip->TMASSA[IP]/VOL[IP]
        << " gr/cm^3" << endl;
    }
    // ****** QUANTITATIVE ANALYSIS ***************
    file6 << endl;
    WTOTAL = 0.000000;
    for (I = 1; I <= cntrls->NPHASE; ++I) WTOTAL = WTOTAL + W[I];
    for (I = 1; I <= cntrls->NPHASE; ++I) XMASS[I] = 100.0 * W[I] / WTOTAL;
    IINNMOL = 0;
    for (I = 1; I <= cntrls->NPHASE; ++I)
    {
        if (IINNMOL == 1) goto L2713;
        if (phases[I].NMOL == 0) IINNMOL=1;
L2713:;
    }
    if (IINNMOL == 1)
    {
        for (I = 1; I <= cntrls->NPHASE; ++I)
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
        for (I = 1; I <= cntrls->NPHASE; ++I)
        {
            FRP[I] = XMASS[I] * phases[I].NMOL / multip->TMASSA[I];
            FT = FT + FRP[I];
        }
        for (I = 1; I <= cntrls->NPHASE; ++I)
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
    if(cntrls->ISPHASE != 0)
    {
        file6 << endl << "Considering Amorphous Content:" << endl;
        SFIS=phases[cntrls->ISPHASE].WTIS/XMASS[cntrls->ISPHASE];
        if(SFIS > 1.0)
        {
            file6 << "PROBLEM:Amount of Internal Standard (Phase #"
                << setw(2) << cntrls->ISPHASE
                << ") is less than the specified "
                << setw(6) << setprecision(2) << phases[cntrls->ISPHASE].WTIS
            << "%." << endl
                << "Amorphous content not computed. Check ISWT in line 11.2 for this phase" << endl;
            goto L2720;
        }
        for (I=1; I <= cntrls->NPHASE; ++I)
        {
            file6 << "PHASE = " << setw(2) << I << " => %MASS = " << setw(6) << setprecision(2) << XMASS[I]*SFIS << endl;
        }
        file6 << "AMORPHOUS  => %MASS = " << setw(6) << setprecision(2) <<100*(1.0-SFIS)  << endl;
    }
L2720:
    file6 << endl;
    return;
}

// SUBROUTINE TO READ RIGAKU DATA FILE
void  DBWS::rigakuread(void)
{
    int I;
    string s,TEST,ASTEP,ANPTS,ATHMIN,ATHMAX;

L10:
    getline(file4,s);
    TEST = s.substr(0,72);
    if(TEST.substr(1,6) == "*START")
    {
        for (I=1; I <= 56; ++I)
        {
            if(TEST[I] == '=')
            {
                ATHMIN=TEST.substr(I+1,I+11-(I+1)+1);
                goto L10;
            }
        }
    }
    if(TEST.substr(1,5) == "*STOP")
    {
        for (I=1; I <= 56; ++I)
        {
            if(TEST[I] == '=')
            {
                ATHMAX=TEST.substr(I+1,I+11-(I+1)+1);
                goto L10;
            }
        }
    }
    if(TEST.substr(1,5) == "*STEP")
    {
        for (I=1; I <= 56; ++I)
        {
            if(TEST[I] == '=')
            {
                ASTEP=TEST.substr(I+1,I+8-(I+1)+1);
                goto L10;
            }
        }
    }
    if(TEST.substr(1,6) == "*COUNT" && TEST.substr(1,7) != "*COUNTE")
    {
        for (I=1; I <= 56; ++I)
        {
            if(TEST[I] == '=')
            {
                ANPTS=TEST.substr(I+1,I+8-(I+1)+1);
                goto L20;
            }
        }
    }
    goto L10;
L20:
    stringstream(ATHMIN) >> g1->THMIN;
    stringstream(ASTEP) >> g1->STEP;
    stringstream(ATHMAX) >> g1->THMAX;
    stringstream(ANPTS) >> datax->NPTS;
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    //file6 << "    DATA ID " << setw(56) << DATAID << endl;
    file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    if (datax->NPTS > IDSZ)
    {
        file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    // read the data
    for (I=1; I <= datax->NPTS; ++I) file4 >> datax->Y[I];
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
    //L99999:
    //	DBWSException("IS THE FILE NOT RIGAKU FORMATED?");
}

void DBWS::readasc(void)
{
    int I;

    for (I=1; I <= 56; ++I)
    {
        if( labels->TEST[I] == ':')
        {
            labels->LB1=I;
            goto L31;
        }
    }
L31:
    for (I=labels->LB1+1; I <= 56; ++I)
    {
        if( labels->TEST[I] != ' ')
        {
            labels->LB2=I;
            goto L41;
        }
    }
L41:
    for (I=labels->LB2; I <= 56; ++I)
    {
        if( labels->TEST[I] == ' ')
        {
            labels->LB3=I-1;
            goto L51;
        }
    }
L51:;
}

// subroutine to read SCINTAG TXT file
void DBWS::scintag(void)
{
    int I;
    double X1,X2,X3;
    string s,ASTEP, DATAID, ATHMIN, ATHMAX;

    getline(file4,s);
    DATAID = " "+s.substr(0,56);
L10:
    getline(file4,s);
    labels->TEST = " "+s.substr(0,56);
    if( labels->TEST.substr(1,12) != "Start Angle:") goto L10;
    //                         Start Angle:
    readasc();
    ATHMIN=labels->TEST.substr(labels->LB2,labels->LB3-labels->LB2+1);
L20:
    getline(file4,s);
    labels->TEST = " "+s.substr(0,56);
    if(labels->TEST.substr(1,11) != "Stop Angle:")goto L20;
    //                         Stop Angle:
    readasc();
    ATHMAX=labels->TEST.substr(labels->LB2,labels->LB3-labels->LB2+1);
L30:
    getline(file4,s);
    labels->TEST = " "+s.substr(0,56);
    if(labels->TEST.substr(1,10) != "Step Size:")goto L30;
    //                         Step Size:
    readasc();
    ASTEP=labels->TEST.substr(labels->LB2,labels->LB3-labels->LB2+1);

L70:
    getline(file4,s);
    if( labels->TEST.substr(1,5) != "Range")goto L70;
    stringstream(ATHMIN) >> g1->THMIN;
    stringstream(ATHMAX) >> g1->THMAX;
    stringstream(ASTEP) >> g1->STEP;
    datax->NPTS=static_cast<int>((g1->THMAX-g1->THMIN)/g1->STEP+1.5);
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    file6 << "    DATA ID " << setw(56) << DATAID << endl;
    file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;

    if (datax->NPTS > IDSZ)
    {
        file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= datax->NPTS; ++I) file4 >> X1 >> datax->Y[I] >> X2 >> X3;
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
}


// SUBROUTINE TO READ SIEMENS UXD DATA FILE
void DBWS::SIEMENSREAD(void)
{
    int I;
    string s,TEST, ASTEP, ANPTS, DATAID, ATHMIN;

    //character*72 test,dataid
    //character*15 athmin,astep,anpts
    //COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ),BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
    //COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT,ULOR,VLOR,ZZZ,UC

    getline(file4,s);
    DATAID = s.substr(0,72);
L10:
    getline(file4,s);
    TEST = s.substr(0,72);
    if(TEST.substr(2,9-2+1) == "STEPSIZE")
    {
        ASTEP=TEST.substr(12,20-12+1);
        goto L10;
    }
    else if(TEST.substr(2,7-2+1) == "2THETA")
    {
        ATHMIN=TEST.substr(10,16-10+1);
        goto L10;
    }else if(TEST.substr(2,10-2+1) == "STEPCOUNT")
    {
        ANPTS=TEST.substr(13,19-13+1);
        goto L10;
    }else if(TEST.substr(2,7-2+1) == "COUNTS")
    {
        goto L20;
    }
    else
    {
        goto L10;
    }
L20:
    stringstream(ATHMIN) >> g1->THMIN;
    stringstream(ASTEP) >> g1->STEP;
    stringstream(ANPTS) >> datax->NPTS;
    g1->THMAX=g1->THMIN + g1->STEP*datax->NPTS;
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    file6 << "    DATA ID " << setw(56) << DATAID << endl;
    file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << g1->THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << g1->THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << g1->STEP << endl;
    if (datax->NPTS > IDSZ)
    {
        file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= datax->NPTS; ++I) file4 >> datax->Y[I];
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
    //L99999:
    //	DBWSException("IS THE FILE NOT UXD SIEMENS FORMAT?");
}

void DBWS::INPTR(void)
{
    int MULTX_;
    double Y_[192+1][3+1];
    int NCTR_[192+1];

    const string UPPER = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const string LOWER = " abcdefghijklmnopqrstuvwxyz";
    // Adding the code from Ian Madsen (10 feb 98)  (Au,Cr,Fe,Co,Cu,
    //                                               Mo,Ag,Ta,W ,Au)
    const double XRYZ[10+1] = { 0.0,
        2.748510,2.289620, 1.935970, 1.788965, 1.540520 ,
        0.709260,0.559360, 0.215947, 0.209010, 0.180195 };

    int I, J, K, N,IK, KK,IX, N2X,NBC, IOF, IRH, ICY, ICZ,NBX, ITN, KKS,
        IRL, IXX, IYY, IRK, NXX, IRC,NINC,IPTS,MLTT,IXDEL,IXRAY, ITIPO,
        ISTOP,NRANGE,IIPHAS,NATOMS, LCOUNT, ISTEST;
    //IBDG
    double X,TH,ABC,THX,DIFB,TAUK,BSTEP,FLEET1, FLEET2, OFSTI0, OFSTI1,
        LAMDAM, CHMBRI,ANGMIN,ANGMAX,STPTIM;
    string SPG,DATE,DATAID,s,s1;

    IPTS = 0;


    // line 1
    getline(file5,s);
    title = s.substr(0,70);

    // line 2
    getline(file5,s);
    stringstream(s.substr( 0*4,4)) >> cntrls->JOBTYP;
    stringstream(s.substr( 1*4,4)) >> cntrls->NPROF;
    stringstream(s.substr( 2*4,4)) >> cntrls->NPHASE;
    stringstream(s.substr( 3*4,4)) >> jnk->NBCKGD;
    stringstream(s.substr( 4*4,4)) >> jnk->NEXCRG;
    stringstream(s.substr( 5*4,4)) >> jnk->NSCAT;
    stringstream(s.substr( 6*4,4)) >> cntrls->INSTRM;
    stringstream(s.substr( 7*4,4)) >> cntrls->IPREF;
    stringstream(s.substr( 8*4,4)) >> cntrls->IASYM;
    stringstream(s.substr( 9*4,4)) >> cntrls->IABSR;
    stringstream(s.substr(10*4,4)) >> cntrls->IDATA;
    stringstream(s.substr(11*4,4)) >> cntrls->ISPHASE;
    stringstream(s.substr(12*4,4)) >> cntrls->I2D94;

    if (cntrls->ISPHASE > cntrls->NPHASE)
    {
        file6 << "Internal Standard Phase does not exist. Check its number in line 2, column 12." << endl
            << "No Internal Standard will be used in the QPA" << endl;
        cntrls->ISPHASE = 0;
    }

    //     open file to write +/- dbws9006 & dbws9411 format (readable by ATOMS and ZORTEP)
    if(cntrls->I2D94 == 1) file53.open("ICF94.ICF");
    if (cntrls->NPROF == 9)
    {
        cntrls->NPROF = _TCHZ;
        sizestrain->NSIZESTRAIN = 9;
    }
    if (jnk->NBCKGD == -1)
    {
        codebck->IBCKCODE=jnk->NBCKGD;
        jnk->NBCKGD = 0;
        //		IBDG   = 0;
    }
    else
    {
        cntrls->FONDO  = 0;
        cntrls->IBGD   = 1;
        codebck->IBCKCODE=jnk->NBCKGD;
    }
    if(cntrls->NPHASE == 0) cntrls->NPHASE=1;
    cntrls->INSTRM=cntrls->INSTRM+1;
    cntrls->JOBTYP=cntrls->JOBTYP+1;
    //cntrls->NPROF=cntrls->NPROF+1;
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
    if(cntrls->JOBTYP == 1) file6 << "FOR X-RAY DATA" << endl;
    if(cntrls->JOBTYP == 1 && cntrls->INSTRM == 2) file6 << "COLLECTED IN SYNCHROTRON AT NSLS OR SRS" << endl;
    if(cntrls->JOBTYP == 2) file6 << "FOR NEUTRON DATA, NUCLEAR INTENSITIES ONLY" << endl;
    if(cntrls->JOBTYP == 2 && cntrls->INSTRM == 2) file6 << "VARYING NO. OF COUNTERS AT EACH STEP" << endl;
    if(cntrls->JOBTYP == 3) file6 << "PATTERN CALCULATION,XRAY" << endl;
    if(cntrls->JOBTYP == 4) file6 << "PATTERN CALCULATION,NEUTRON" << endl;
    if(cntrls->JOBTYP < 1 || cntrls->JOBTYP > 4) DBWSException("7777");
    if(cntrls->IDATA == 0) file6 << "READ DATA IN TRADITIONAL DBWS FORMAT" << endl;
    if(cntrls->IDATA == 1) file6 << "READ DATA IN FREE FORMAT" << endl;
    if(cntrls->IDATA == 2) file6 << "READ DATA IN GSAS STD FORMAT" << endl;
    if(cntrls->IDATA == 3) file6 << "READ DATA IN PHILIPS UDF FORMAT" << endl;
    if(cntrls->IDATA == 4) file6 << "READ DATA IN SCINTAG TXT FORMAT" << endl;
    if(cntrls->IDATA == 5) file6 << "READ DATA IN SIEMENS UXD FORMAT" << endl;
    if(cntrls->IDATA == 6) file6 << "READ DATA IN RIGAKU ASC FORMAT" << endl;
    file6 << "NUMBER OF PHASES= " << setw(4) << cntrls->NPHASE << endl
        << "NUMBER OF EXCLUDED REGIONS= " << setw(4) << jnk->NEXCRG << endl
        << "NUMBER OF SCATTERING SETS= " << setw(4) << jnk->NSCAT << endl;

    if (jnk->NBCKGD-1 < 0)
    {
        file6 << "BACKGROUND TO BE REFINED" << endl;
    }
    else if(jnk->NBCKGD-1 == 0)
    {
        file6 << "BACKGROUND DATA TO BE READ FROM FILE" << endl;
    }
    else
    {
        file6 << "BACKGROUND CORRECTION BY INTERPOLATION BETWEEN THE " << setw(4) << jnk->NBCKGD << " POINTS GIVEN" << endl;
    }


    switch (cntrls->NPROF) {
    case _Gaussian:
        file6 << "GAUSSIAN PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    case _Lorentzian:
        file6 << "LORENTZIAN (CAUCHY) PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    case _Mod1:
        file6 << "MOD 1 LORENTZIAN PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    case _Mod2:
        file6 << "MOD 2 LORENTZIAN PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    case _SplitPearsonVII:
        file6 << "SPLIT PEARSON VII PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    case _pseudoVoigt:
        file6 << "PSEUDO-VOIGT (PV) PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    case _PearsonVII:
        file6 << "PEARSON VII PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    case _TCHZ:
        file6 << "THOMPSON-COX-HASTINGS (PV) PROFILE, NPROF = " << setw(4) << cntrls->NPROF << endl;
        break;
    }

    if (cntrls->IPREF == 0)
    {
        file6 << "IPREF=0, UDA-RIETVELD PREFERRED ORIENTATION FUNCTION" << endl;
    }
    else
    {
        file6 << "IPREF=1, MARCH-DOLLASE PREFERRED ORIENTATION FUNCTION" << endl;
    }

    //  SURFACE ROUGHNESS
    switch (cntrls->IABSR) {
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
    if (cntrls->IASYM == 0) file6 << "IASYM=0, Usual Rietveld Asymmetry" << endl;
    if (cntrls->IASYM == 1) file6 << "IASYM=1, Asymmetry by Riello et al.:',' Powder Diffraction,10,204-206,1995" << endl;
    //-----FONDO=0 :BKG EVALUATED USING STANDARD METHODS
    //-----FONDO=1 :BKG EVALUATED USING THE ISOTROPIC THERMAL FACTORS
    //-----FONDO=2 :BKG EVALUATED USING THE OVERAL THERMAL FACTORS
    //-----AIR SCATTERING,IAXX= 1 AIR SCATTERING ADDED TO THE BACKGROUND if SCAIR<>0
    //                    IAXX= -1 AIR SCATT. EVALUATED AND SUBTRACTED FROM DATA
    //                    IAXX= 0 AIR SCATT NOT EVALUATED and ibgd = 1
    //-----LINEAR ABSORP.CORR., IAS= 1  DATA CORRECTED FOR LINEAR ABSORPTION
    //
    //-----FI = AIR FRACTION INSIDE THE SAMPLE
    if(codebck->IBCKCODE == -1)
    {
        // line 2.2
        getline(file5,s);
        stringstream(s.substr( 0,4)) >> cntrls->IAS;
        stringstream(s.substr( 4,4)) >> cntrls->FONDO;
    }
    //L5333:
    if (cntrls->IBGD == 1)
    {
        file6 << "NO AMOPHOUS AND COMPTON CORRECTION TO THE BGD" << endl;
        goto L7448;
    }
    if(cntrls->IAS == 0) file6 << "IAS = 0, NO LINEAR ABSORPTION CORRECTION" << endl;
    if(cntrls->IAS == 1) file6 << "IAS = 1, LINEAR ABSORPTION CORRECTION IS APPLIED" << endl;
    if(cntrls->FONDO == 1) file6 << "FONDO = 1,ISOTROPIC B FACTOR USED FOR BKG EVALUATION" << endl;
    if(cntrls->FONDO == 2) file6 << "FONDO = 2,OVERALL Q USED FOR BKG EVALUATION" << endl;


    // line 3
L7448:
    getline (file5,s);
    stringstream(s.substr( 0,1)) >> cntrls->IOT;
    stringstream(s.substr( 1,1)) >> cntrls->IPL;
    stringstream(s.substr( 2,1)) >> cntrls->IPC;
    stringstream(s.substr( 3,1)) >> cntrls->MAT;
    stringstream(s.substr( 4,1)) >> cntrls->NXT;

    stringstream(s.substr( 6,1)) >> cntrls->LST1;
    stringstream(s.substr( 7,1)) >> cntrls->LST2;
    stringstream(s.substr( 8,1)) >> cntrls->LST3;
    stringstream(s.substr( 9,1)) >> cntrls->IPL1;
    stringstream(s.substr(10,1)) >> cntrls->IPL2;

    stringstream(s.substr(12,1)) >> cntrls->IPLST;
    stringstream(s.substr(13,1)) >> cntrls->IPLOSS;
    stringstream(s.substr(14,1)) >> cntrls->IPLCAL;
    stringstream(s.substr(15,1)) >> cntrls->IPLPOL;
    stringstream(s.substr(16,1)) >> cntrls->IPLCOM;

    stringstream(s.substr(18,1)) >> cntrls->IPLDIS;
    stringstream(s.substr(19,1)) >> cntrls->IPLAM;
    stringstream(s.substr(20,1)) >> cntrls->IPBIG;
    if(codebck->IBCKCODE != -1)
    {
        cntrls->IPLOSS=0;
        cntrls->IPLCAL=0;
        cntrls->IPLPOL=0;
        cntrls->IPLCOM=0;
        cntrls->IPLDIS=0;
        cntrls->IPLAM =0;
        cntrls->IPBIG =0;
    }
    simoper->ISIMOP=0;
    if (cntrls->IOT != 0) file6 << "OUTPUT OBSERVED AND CALCULATED INTENSITIES ON LAST CYCLE" << endl;
    if (cntrls->IPL != 0) file6 << "GENERATE LINE PRINTER PLOT" << endl;
    if (cntrls->IPL2 != 0) file6 << "Generate file with obs, calc, dif intensiies, weighted diference, and Bragg peak positions" << endl;
    if (cntrls->IPLST == 1) file6 << "GENERATE PARAMETER LIST" << endl;
    if (cntrls->IPLST == 2) file6 << "GENERATE PARAMETER LIST" << endl << "GENERATE PARAMETERS AND STD.DEV. IN FINAL CYCLE FOR DATA BASE" << endl;
    if (cntrls->IPC == 1) file6 << "OUTPUT INTENSITIES" << endl;
    if (cntrls->IPC == 2) file6 << "OUTPUT ABSOLUTE VALUES OF STRUCTURE FACTORS + PHASE ANGLE" << endl;
    if (cntrls->IPC == 3) file6 << "OUTPUT A AND B (CALC + OBS) STRUCTURE FACTORS" << endl;
    if (cntrls->MAT != 0 && cntrls->JOBTYP < 3) file6 << "OUTPUT CORRELATION MATRIX" << endl;
    if (cntrls->NXT != 0 && cntrls->JOBTYP < 3) file6 << "GENERATE NEW INPUT FILE" << endl;
    if (cntrls->LST1 != 0) file6 << "PRINT REFLECTION LIST" << endl;
    if (cntrls->NPHASE == 1)cntrls->LST3=0;
    if (cntrls->LST3 != 0) file6 << "Print merged reflection list" << endl;
    if (cntrls->LST2 != 0) file6 << "Print corrected data" << endl;

    // line 4
    getline (file5,s);
    stringstream(s.substr( 0*8, 8)) >> g1->LAMDA[1];
    stringstream(s.substr( 1*8, 8)) >> g1->LAMDA[2];
    stringstream(s.substr( 2*8, 8)) >> params->RATIO[2];
    stringstream(s.substr( 3*8, 8)) >> g1->BKPOS;
    stringstream(s.substr( 4*8, 8)) >> g1->WDT;
    stringstream(s.substr( 5*8, 8)) >> g1->CTHM;
    stringstream(s.substr( 6*8, 8)) >> g1->TMV;
    stringstream(s.substr( 7*8, 8)) >> g1->RLIM;
    stringstream(s.substr( 8*8, 8)) >> cntrls->SW;
    params->RATIO[1] = 1.0;
    LAMDAM=(params->RATIO[1]*g1->LAMDA[1]+params->RATIO[2]*g1->LAMDA[2])/(params->RATIO[1]+params->RATIO[2]);
    file6
        << " WAVELENGTHS= "
        << fixed
        << setw(9) << setprecision(6) << g1->LAMDA[1]
        << setw(9) << setprecision(6) << g1->LAMDA[2]
        << " LAMDA MEAN = "
        << setw(9) << setprecision(6) << LAMDAM <<endl;
    //  stop !cp ap 21 97
    file6 << "ALPHA2:ALPHA1 RATIO = " << setw(9) << setprecision(5) << params->RATIO[2] << endl;
    for (IX=1; IX <= 10; ++IX) if(1.03*g1->LAMDA[1] > XRYZ[IX])goto L119;
L119:
    IXRAY=IX;
    file6 << "BASE OF PEAK = 2.0*HW*" << setw(8) << setprecision(2) << g1->WDT << endl;
    file6 << "MONOCHROMATOR CORRECTION =" << setw(8) << setprecision(4) << g1->CTHM << endl;
    file6 << "ABSORPTION CORRECTION COEFFICIENT = " << setw(8) << setprecision(4) << g1->TMV << " CM-1" << endl
        << "SLAB-WIDTH = " << setw(8) << setprecision(4) << cntrls->SW << " CM." << endl;
    if (cntrls->IASYM == 0)
    {
        file6 << "RIETVELD ASYMMETRY CORRECTION FOR ANGLES LESS THAN " << setw(8) << setprecision(3) << g1->RLIM << " DEGREES" << endl;
    }
    else
    {
        file6 << "ASYMMETRY CORRECTION FOR ANGLES LESS THAN "
            << setw(8) << setprecision(3) << 90.0 - g1->RLIM << " DEGREES" << endl
            << "                         AND GREATER THAN "
            << setw(8) << setprecision(3) << 90.0 + g1->RLIM << " DEGREES" << endl;
    }

    // line 5
    getline(file5,s);
    stringstream(s.substr(0,4)) >> cntrls->MCYCLE;
    stringstream(s.substr(4,4)) >> cntrls->EPS;
    stringstream(s.substr(8,4)) >> params->RELAX[1];
    stringstream(s.substr(12,4)) >> params->RELAX[2];
    stringstream(s.substr(16,4)) >> params->RELAX[3];
    stringstream(s.substr(20,4)) >> params->RELAX[4];

    if(cntrls->JOBTYP > 2)
    {
        cntrls->MCYCLE=1;
        stringstream(s.substr(24,8)) >> g1->THMIN;
        stringstream(s.substr(32,8)) >> g1->STEP;
        stringstream(s.substr(40,8)) >> g1->THMAX;
    }
    cntrls->ICYRUN = cntrls->MCYCLE;
    file6 << "NUMBER OF CYCLES = " << setw(4) << cntrls->MCYCLE << endl;
    file6 << "RELAXATION FACTORS" << endl
        << "FOR COORDINATES= " << setw(5) << setprecision(2) << params->RELAX[1] << endl
        << "FOR ANISOTRPIC TEMPERATURE FACTORS= " << setw(5) << setprecision(2) << params->RELAX[2] << endl
        << "FOR FWHM PARAMETERS= " << setw(5) << setprecision(2) << params->RELAX[3] << endl
        << "FOR LATTICE CONSTANTS= " << setw(5) << setprecision(2) << params->RELAX[4] << endl;
    file6 << "EPS-VALUE= " << setw(6) << setprecision(1) << cntrls->EPS << endl;

    if(jnk->NBCKGD >= 2)
    {
        // line 6(*)
        for (I=1; I <= jnk->NBCKGD; ++I)
        {
            getline(file5,s);
            stringstream(s.substr(0,8)) >> jnk->POS[I];
            stringstream(s.substr(8,8)) >> jnk->BCK[I];
        }
        file6 << "BACKGROUND" << endl
            << "POSITION    INTENSITY" << endl;
        for (I=1; I <= jnk->NBCKGD; ++I)
        {
            file6 << setw(9) << setprecision(4) << jnk->POS[I]
            << setw(9) << setprecision(4) << jnk->BCK[I] << endl;
        }
    }

    if(jnk->NEXCRG > 0)
    {
        // line 7(*)
        for (I=1; I <= jnk->NEXCRG; ++I)
        {
            getline(file5,s);
            stringstream(s.substr(0,8)) >> jnk->ALOW[I];
            stringstream(s.substr(8,8)) >> jnk->AHIGH[I];
        }
        file6 << "EXCLUDED REGIONS" << endl
            << "FROM     TO" << endl;
        for (I=1; I <= jnk->NEXCRG; ++I)
        {
            file6 << setw(9) << setprecision(4) << jnk->ALOW[I]
            << setw(9) << setprecision(4) << jnk->AHIGH[I] << endl;
        }
    }

    if(jnk->NSCAT > 0)
    {
        for (I=1; I <= jnk->NSCAT; ++I)
        {
            if(cntrls->JOBTYP == 2 || cntrls->JOBTYP == 4)
            {
                // line 8.1 XRD(*)
                getline(file5,s);
                coefc->NAM[I] = s.substr(0,4);
                stringstream(s.substr(4,8)) >> coeff->DFP[I];
                stringstream(s.substr(12,8)) >> coeff->XMAS[I];
                goto L125;
            }
            // line 8.1 ND(*)
            getline(file5,s);
            coefc->NAM[I] = s.substr(0,4);
            stringstream(s.substr(4,8)) >> coeff->DFP[I];
            stringstream(s.substr(12,8)) >> coeff->DFPP[I];
            stringstream(s.substr(20,8)) >> coeff->XMAS[I];
            K=0;
            // line 8.2 XRD(*)
            L126:
            getline(file5,s);
            for (J=1; J <= 9; ++J) stringstream(s.substr((J-1)*8,8)) >> coeff->AC[J][I];
            if (coeff->AC[1][I] == -100.0) coeff->COEF(&I,&K);
            if(coeff->AC[3][I] != 0.)goto L125;
            K=K+1;
            coeff->POSI[K]=coeff->AC[1][I];
            coeff->SCAT[K]=coeff->AC[2][I];
            if (K <= 29) goto L126;
            file6 << "TOO MANY SCATTERING TABLE ENTRIES" << endl;
            DBWSException("7700");
            L125:;
        }
        if(cntrls->JOBTYP == 1 || cntrls->JOBTYP == 3)goto L129;
        file6 << "SCATTERING LENGTHS" << endl;
        for (I=1; I <= jnk->NSCAT; ++I)
        {
            file6 << "FOR " << setw(4) << coefc->NAM[I] << "     "
                << setw(10) << setprecision(6) << coeff->DFP[I] << endl;
        }
        for (I=1; I <= jnk->NSCAT; ++I)
        {
            for (J=1; J <= 9; ++J) coeff->AC[J][I]=0.0;
            coeff->DFPP[I]=0.0;
        }
        goto L124;
        L129:
        file6 << "FORMFACTORS" << endl;
        for (I=1; I <= jnk->NSCAT; ++I)
        {
            file6 << "FOR " << setw(4) << coefc->NAM[I]
            << " DFP=" << setw(10) << setprecision(6) << coeff->DFP[I]
            << " DFPP=" << setw(10) << setprecision(6) << coeff->DFPP[I] << endl
                << "COEFFICIENTS= ";
            for (J=1; J <= 9; ++J) file6 << setw(10) << setprecision(6) << coeff->AC[J][I];
            file6 << endl;
        }
        L124:;
    }


    // line 9
    getline(file5,s);
    stringstream(s.substr(0,8)) >> cntrls->MAXS;
    if(cntrls->MAXS > MSZ)
    {
        file6 << "* YOU HAVE DECLARED MORE CODEWORDS THAN WILL FIT INTO *" << endl
            << "* THE -MSZ- ARRAY.  EITHER DECREASE THE # OF CODEWORD *" << endl
            << "* OR  INCREASE  THE  -MSZ-  ARRAY  SIZE AND RECOMPILE *" << endl;
        cout  << "* YOU HAVE DECLARED MORE CODEWORDS THAN WILL FIT INTO *" << endl
            << "* THE -MSZ- ARRAY.  EITHER DECREASE THE # OF CODEWORD *" << endl
            << "* OR  INCREASE  THE  -MSZ-  ARRAY  SIZE AND RECOMPILE *" << endl;
        goto L99995;
    }
    cout << "INPUT:    CYCLES =" << setw(4) << cntrls->MCYCLE << "     REFINABLE PARAMETERS =" << setw(4) <<cntrls->MAXS << endl;

    //     CHECK DIMENSIONING FOR SOME EQUIVALENCED ARRAYS
    if (IDSZ < MSZ*cntrls->MAXS)
    {
        file6 << "CHANGE IDSZ OR MSZ SO THAT IDSZ IS EQUAL TO OR GREATER THAN MSZ*MAXS" << endl
            << "IDSZ, MAXIMUM NO. OF DATA POINTS = " << setw(6) << IDSZ << endl
            << "MSZ,  MATRIX SIZE                = " << setw(6) << MSZ << endl
            << "MAXS, NO. OF PARAMETERS VARIED   = " << setw(6) << endl
            << "*** JUST A DIMENSIONING ERROR ***" << endl;
        DBWSException("IDSZ IS LESS THAN MSZ*MAXS");
    }
    if (cntrls->JOBTYP > 2) cntrls->MAXS=0;
    file6 << "NUMBER OF PARAMETERS VARIED= " << setw(5) << cntrls->MAXS << endl;

    // line 10.1
    getline(file5,s);
    params->GLB_[1-1] = s.substr(0,8);
    params->GLB_[10-1] = s.substr(1*8,8);
    params->GLB_[11-1] = s.substr(2*8,8);
    params->GLB_[8-1] = s.substr(3*8,8);
    params->GLB_[9-1] = s.substr(4*8,8);
    params->GLB_[12-1] = s.substr(5*8,8);
    params->GLB_[13-1] = s.substr(6*8,8);

    // line 10.11
    getline(file5,s);
    stringstream(s.substr(0,8))   >> params->GLB_[1-1].codeword;
    stringstream(s.substr(1*8,8)) >> params->GLB_[10-1].codeword;
    stringstream(s.substr(2*8,8)) >> params->GLB_[11-1].codeword;
    stringstream(s.substr(3*8,8)) >> params->GLB_[8-1].codeword;
    stringstream(s.substr(4*8,8)) >> params->GLB_[9-1].codeword;
    stringstream(s.substr(5*8,8)) >> params->GLB_[12-1].codeword;
    stringstream(s.substr(6*8,8)) >> params->GLB_[13-1].codeword;

    file6 << "GLOBAL PARAMETERS AND CODEWORDS" << endl
        << "ZEROPOINT= "
        << setw(8) << setprecision(2) << params->GLB_[1-1]
    << setw(8) << setprecision(2) << params->GLB_[1-1].codeword << endl;

    //-----PMON1,PMON2=PARAMETER OF THE MONOCHROMATOR
    //-----IT READS PARAMETER OF THE MONOCHROMATOR AND AIR SCALE
    //-----if THE MONOCHROMATOR WORKS ON THE INCIDENT BEAM PUT :
    //-----PMON1=1;PMON2=0;
    // next 2 'READ' are for amorphous bkg codes
    // SCAIR=glb(17), FLSCAIR=aglb(17), SCAM =params->GLB_[20-1], FLSCAM=params->GLB_[20-1].A
    // PMON1=params->GLB_[18-1], FLMON1 =params->GLB_[18-1].A, PMON2=params->GLB_[19-1], FLMON2=params->GLB_[19-1].A
    //      READ(5,455,END=99999)SCAIR,FLSCAIR,SCAM,FLSCAM,
    //     *  PMON1,FLMON1,PMON2,FLMON2
    //
    //      READ(5,455,END=99999)glb(17),aglb(17),params->GLB_[20-1], params->GLB_[20-1].A,
    //     *  params->GLB_[18-1],params->GLB_[18-1].A,params->GLB_[19-1], params->GLB_[19-1].A
    //455     FORMAT(BZ,8F8.0)
    //  !cp jun 95 start
    if (cntrls->IBGD != 1)
    {
        // line 10.2 and line 10.21
        getline(file5,s);
        params->GLB_[20-1] = s.substr(0*8,8);
        params->GLB_[18-1] = s.substr(1*8,8);
        params->GLB_[19-1] = s.substr(2*8,8);

        getline(file5,s);
        stringstream(s.substr(0*8,8)) >> params->GLB_[20-1].codeword;
        stringstream(s.substr(1*8,8)) >> params->GLB_[18-1].codeword;
        stringstream(s.substr(2*8,8)) >> params->GLB_[19-1].codeword;

        file6 << "AMORPHOUS SCALE and CODEWORD= "
            << setw(14) << setprecision(4) << params->GLB_[20-1]
        << setw(14) << setprecision(4) << params->GLB_[20-1].codeword << endl;
        file6 << "MONOCROMATOR BANDPASS PARAMETERS AND CODEWORDS" << endl
            << "PARAMETERS MONOC="
            << setw(8) << setprecision(4) << params->GLB_[18-1]
        << setw(8) << setprecision(4) << params->GLB_[18-1].codeword
        << "                  "
            << setw(8) << setprecision(4) << params->GLB_[19-1]
        << setw(8) << setprecision(4) << params->GLB_[19-1].codeword << endl;
    }


    if(jnk->NBCKGD == 0)
    {
        // line 10.4(*)
        // line 10.3 and line 10.31
        getline(file5,s);
        params->GLB_[2-1] = s.substr(0*9,9);
        params->GLB_[3-1] = s.substr(1*9,9);
        params->GLB_[4-1] = s.substr(2*9,9);
        params->GLB_[5-1] = s.substr(3*9,9);
        params->GLB_[6-1] = s.substr(4*9,9);
        params->GLB_[7-1] = s.substr(5*9,9);

        getline(file5,s);
        stringstream(s.substr(0*9,9)) >> params->GLB_[2-1].codeword;
        stringstream(s.substr(1*9,9)) >> params->GLB_[3-1].codeword;
        stringstream(s.substr(2*9,9)) >> params->GLB_[4-1].codeword;
        stringstream(s.substr(3*9,9)) >> params->GLB_[5-1].codeword;
        stringstream(s.substr(4*9,9)) >> params->GLB_[6-1].codeword;
        stringstream(s.substr(5*9,9)) >> params->GLB_[7-1].codeword;


        file6 << "BACKGROUND PARAMETERS AND CODEWORDS" << endl
            << "ORIGIN OF BACKGROUND POLYNOMIAL AT TWO-THETA = "
            << setw(8) << setprecision(3) << g1->BKPOS << "DEGREES" << endl
            << setw(12) << setprecision(4) << params->GLB_[2-1]
            << setw(12) << setprecision(4) << params->GLB_[3-1]
            << setw(12) << setprecision(4) << params->GLB_[4-1]
            << setw(12) << setprecision(4) << params->GLB_[5-1]
            << setw(12) << setprecision(4) << params->GLB_[6-1]
            << setw(12) << setprecision(4) << params->GLB_[7-1] << endl
            << setw(12) << setprecision(3) << params->GLB_[2-1].codeword << "    "
            << setw(12) << setprecision(3) << params->GLB_[3-1].codeword << "    "
            << setw(12) << setprecision(3) << params->GLB_[4-1].codeword << "    "
            << setw(12) << setprecision(3) << params->GLB_[5-1].codeword << "    "
            << setw(12) << setprecision(3) << params->GLB_[6-1].codeword << "    "
            << setw(12) << setprecision(3) << params->GLB_[7-1].codeword << "    " << endl;
    }

    file6 << "DISPLACEMENT PEAKSHIFT PARAMETER AND CODEWORD"
        << setw(8) << setprecision(2) << params->GLB_[10-1]
    << setw(8) << setprecision(2) << params->GLB_[10-1].codeword << endl
        << "TRANSPARENCY PEAKSHIFT PARAMETER AND CODEWORD"
        << setw(8) << setprecision(2) << params->GLB_[11-1]
    << setw(8) << setprecision(2) << params->GLB_[11-1].codeword << endl;
    file6 << "SURFACE ROUGHNESS P PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << params->GLB_[8-1]
    << setw(9) << setprecision(4) << params->GLB_[8-1].codeword << endl
        << "SURFACE ROUGHNESS Q PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << params->GLB_[9-1]
    << setw(9) << setprecision(4) << params->GLB_[9-1].codeword << endl
        << "SURFACE ROUGHNESS R PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << params->GLB_[12-1]
    << setw(9) << setprecision(4) << params->GLB_[12-1].codeword << endl
        << "SURFACE ROUGHNESS T PARAMETER AND CODEWORD"
        << setw(9) << setprecision(4) << params->GLB_[13-1]
    << setw(9) << setprecision(4) << params->GLB_[13-1].codeword << endl;


    if(cntrls->JOBTYP <= 2)
    {

        if(cntrls->JOBTYP == 1 && cntrls->INSTRM == 2)           //-----if PATTERN CALCULATION ONLY FOR SYNCHROTRON X-RAY DATA
        {
            //       READ DATA FROM SYNCHROTRON SOURCE AND CORRECT DATA FOR DEAD TIME,
            //      CALCULATE VARIANCE FOR EACH OF THE DATA POINTS
            //     DATE IS OCT85,FEB86,AUG86,SRS83,SRS91
            //     NRANGE IS THE NO OF BLOCKS IN WHICH THE DATA ARE GIVEN
            //     OFSTI0 - DARK BEAM CURRENT, READING WITH NO ELECTRONS IN THE CHAMBER
            //     OFSTI1 - DETECTOR DARK BEAM CURRENT
            //     CHMBRI - ALL CURRENTS ARE  NORMALISED TO CHMBRI
            getline(file4,s);
            DATE = s.substr(0,5);
            DATAID = s.substr(5,56);
            file6 << "    DATA ID " << DATAID << endl;
            getline(file4,s);
            stringstream(s.substr(0*8,8)) >> NRANGE;
            stringstream(s.substr(1*8,8)) >> CHMBRI;
            stringstream(s.substr(2*8,10)) >> TAUK;

            //     NPTS IS THE COUNTER FOR TOTAL NO OF POINTS IN ALL THE RANGES
            datax->NPTS = 0;
            for (J=1; J <= NRANGE; ++J)
            {
                //     READ INFORMATION ABOUT EACH RANGE
                getline(file4,s);
                stringstream(s.substr(0*8,8)) >> ANGMIN;
                stringstream(s.substr(1*8,8)) >> g1->STEP;
                stringstream(s.substr(2*8,8)) >> ANGMAX;
                stringstream(s.substr(3*8,8)) >> STPTIM;
                stringstream(s.substr(4*8,8)) >> OFSTI0;
                stringstream(s.substr(5*8,8)) >> OFSTI1;
                IPTS = static_cast<int>((ANGMAX-ANGMIN)/g1->STEP +1.5);
                //     FIND MAXIMUM AND MINIMUM TWO THETA IN ALL RANGES
                if (J == 1) g1->THMIN = ANGMIN;
                if (J == NRANGE) g1->THMAX = ANGMAX;
                if (J > 1)
                {
                    getline(file4,s);
                    IPTS = IPTS - 1;
                }
                if (datax->NPTS+IPTS > IDSZ)
                {
                    file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
                        << setw(5) << (datax->NPTS+IPTS) << " POINTS WERE INPUT" << endl
                        << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
                    DBWSException("TOO MANY DATA POINTS");
                }
                //     VAR(I) IS USED FOR TWO QUANTITIES, JUST TO SAVE SOME SPACE.
                if (DATE == "OCT85")
                {
                    for (I=1; I <= IPTS; ++I)
                    {
                        getline(file4,s);
                        stringstream(s.substr(24,7)) >> datax->Y[I+datax->NPTS];
                        stringstream(s.substr(32,7)) >> datax->VAR[I+datax->NPTS];
                    }
                }
                else if (DATE == "FEB86")
                {
                    for (I=1; I <= IPTS; ++I)
                    {
                        getline(file4,s);
                        stringstream(s.substr(24,7)) >> datax->Y[I+datax->NPTS];
                        stringstream(s.substr(32,7)) >> datax->VAR[I+datax->NPTS];
                    }
                }
                else if (DATE == "AUG86")
                {
                    for (I=1; I <= IPTS; ++I)
                    {
                        getline(file4,s);
                        stringstream(s.substr(51,10)) >> datax->VAR[I+datax->NPTS];
                        stringstream(s.substr(62,10)) >> datax->Y[I+datax->NPTS];
                    }
                }
                else if (DATE == "SRS83")
                {
                    for (I=1; I <= IPTS; ++I)
                    {
                        getline(file4,s);
                        stringstream(s.substr(37,9)) >> datax->VAR[I+datax->NPTS];
                        stringstream(s.substr(47,9)) >> datax->Y[I+datax->NPTS];
                    }
                }
                else if (DATE == "SRS91")
                {
                    for (I=1; I <= IPTS; ++I)
                    {
                        getline(file4,s);
                        stringstream(s.substr(29,9)) >> datax->VAR[I+datax->NPTS];
                        stringstream(s.substr(48,9)) >> datax->Y[I+datax->NPTS];
                    }
                }
                else
                {
                    file6 << "    WHEN AND WHERE WERE THESE DATA TAKEN. OCT85,FEB86,AUG86,SRS83,SRS91" << endl;
                }
                for (I=datax->NPTS+1; I <= IPTS+datax->NPTS; ++I)
                {
                    datax->Y[I] = datax->Y[I]/STPTIM;
                    FLEET1 = datax->VAR[I]-OFSTI0;
                    FLEET2= FLEET1/CHMBRI;
                    FLEET2= pow((FLEET2*(1-TAUK*datax->Y[I])) , 2.0);
                    datax->VAR[I] = (STPTIM*datax->Y[I])/FLEET2;
                    datax->Y[I] = STPTIM*(datax->Y[I]/(1-TAUK*datax->Y[I])-OFSTI1)*CHMBRI/FLEET1;
                    if (abs(datax->Y[I]) < 0.01) datax->Y[I]=1.0;
                }
                //L530:
                datax->NPTS=datax->NPTS+IPTS;
            }
            file6 << "DATA RANGE (2THETA):  START = " << setw(8) << setprecision(2) << g1->THMIN
                << ", STOP =" << setw(8) << setprecision(2) << g1->THMAX
                << ", STEP =" << setw(8) << setprecision(2) << g1->STEP << endl;
            cout  << "DATA RANGE (2THETA):  START = " << setw(8) << setprecision(2) << g1->THMIN
                << ", STOP =" << setw(8) << setprecision(2) << g1->THMAX
                << ", STEP =" << setw(8) << setprecision(2) << g1->STEP << endl;
            goto L484;
            //       ENDS SYNCHROTRON DATA MODifICATION
        }

        if(cntrls->JOBTYP == 2 && cntrls->INSTRM == 2)           //-----if PATTERN CALCULATION ONLY FOR MULTIPLE NEUTRON DATA
        {
            //     BEGIN VARIANCE CALCULATION FOR VARYING NO. OF COUNTERS AT EACH STEP
            getline(file4,s);
            stringstream(s.substr(0*8,8)) >> g1->THMIN;
            stringstream(s.substr(1*8,8)) >> g1->STEP;
            stringstream(s.substr(2*8,8)) >> g1->THMAX;
            DATAID = s.substr(3*8,56);
            file6 << "    DATA ID " << DATAID << endl;
            file6 << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(2) << g1->THMIN
                << ", STOP =" << setw(8) << setprecision(2) << g1->THMAX
                << ", STEP =" << setw(8) << setprecision(2) << g1->STEP << endl;
            datax->NPTS = static_cast<int>((g1->THMAX-g1->THMIN)/g1->STEP+1.5);
            if (datax->NPTS > IDSZ)
            {
                file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
                    << datax->NPTS+IPTS << " POINTS WERE INPUT" << endl
                    << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
                DBWSException("TOO MANY DATA POINTS");
            }
            I=0;
            while ( I > datax->NPTS)
            {
                getline(file4,s);
                for (J = 1; J <= 10; ++J)
                {
                    s1 = s.substr((J-1)*8,2);
                    if ( !s1.empty() )
                    {
                        ++I;
                        stringstream(s1) >> datax->VAR[I];
                    }
                    s1 = s.substr((J-1)*8+2,6);
                    if ( !s1.empty() )
                    {
                        stringstream(s1) >> datax->Y[I];
                    }
                }
            }
            for (I=1; I <= datax->NPTS; ++I) datax->VAR[I]=datax->Y[I]/datax->VAR[I];
            goto L484;
        }


        switch (cntrls->IDATA) {
        case 0:
            // read DBWS formated
            getline(file4,s);
            stringstream(s.substr(0*8,8)) >> g1->THMIN;
            stringstream(s.substr(1*8,8)) >> g1->STEP;
            stringstream(s.substr(2*8,8)) >> g1->THMAX;
            DATAID = s.substr(3*8,56);
            cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << g1->THMIN
                << ", STOP =" << setw(8) << setprecision(3) << g1->THMAX
                << ", STEP =" << setw(8) << setprecision(3) << g1->STEP << endl;
            file6 << "    DATA ID " << DATAID << endl;
            cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << g1->THMIN
                << ", STOP =" << setw(8) << setprecision(3) << g1->THMAX
                << ", STEP =" << setw(8) << setprecision(3) << g1->STEP << endl;
            datax->NPTS=static_cast<int>((g1->THMAX-g1->THMIN)/g1->STEP+1.5);
            if (datax->NPTS > IDSZ)
            {
                file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
                    << setw(5) << datax->NPTS << " POINTS WERE INPUT" << endl
                    << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
                cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
                    << setw(5) << datax->NPTS << " POINTS WERE INPUT" << endl
                    << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
                DBWSException("TOO MANY DATA POINTS");
            }
            I=0;
            while ( I < datax->NPTS)
            {
                getline(file4,s);
                for (J = 1; J <= 8; ++J)
                {
                    s1 = s.substr((J-1)*8,7);
                    if ( !s1.empty() )
                    {
                        ++I;
                        stringstream(s1) >> datax->Y[I];
                    }
                }
            }
            maxint->XMAXINT = 0.0;            //     CHECK FOR NON-ZERO COUNTS AT ALL POINTS
            break;
        case 1:
            // read DBWS formated
            getline(file4,s);
            stringstream(s.substr(0*8,8)) >> g1->THMIN;
            stringstream(s.substr(1*8,8)) >> g1->STEP;
            stringstream(s.substr(2*8,8)) >> g1->THMAX;
            DATAID = s.substr(3*8,56);
            cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << g1->THMIN
                << ", STOP =" << setw(8) << setprecision(3) << g1->THMAX
                << ", STEP =" << setw(8) << setprecision(3) << g1->STEP << endl;
            file6 << "    DATA ID " << DATAID << endl;
            cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << g1->THMIN
                << ", STOP =" << setw(8) << setprecision(3) << g1->THMAX
                << ", STEP =" << setw(8) << setprecision(3) << g1->STEP << endl;
            datax->NPTS=static_cast<int>((g1->THMAX-g1->THMIN)/g1->STEP+1.5);
            if (datax->NPTS > IDSZ)
            {
                file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
                    << setw(5) << datax->NPTS << " POINTS WERE INPUT" << endl
                    << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
                cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
                    << setw(5) << datax->NPTS << " POINTS WERE INPUT" << endl
                    << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
                DBWSException("TOO MANY DATA POINTS");
            }
            for (I=1; I <= datax->NPTS; ++I) file4 >> datax->Y[I];
            maxint->XMAXINT = 0.0;            //     CHECK FOR NON-ZERO COUNTS AT ALL POINTS
            break;
        case 2:
            GSASREAD();          // read GSAS data file (read start,stop,step and data)
            break;
        case 3:
            PHILIPSREAD();      // read Philips data file (read start,stop,step and data)
            break;
        case 4:
            scintag();          // read SCINTAG data file (read start,stop,step and data)
            break;
        case 5:
            SIEMENSREAD();     // read  SIEMENS UXD data file (read start,stepsize,stepcount and data)
            break;
        case 6:
            rigakuread();        // read  RIGAKU data file (read start,stepsize,stop,stepcount and data)
            break;
        }

        for (I=1; I <= datax->NPTS; ++I)
        {
            if(datax->Y[I] <= 1.0E-6) datax->Y[I]=1.0;
            if(datax->Y[I] > maxint->XMAXINT) maxint->XMAXINT=datax->Y[I];
        }
        for (I=1; I <= datax->NPTS; ++I)
        {
            //     BUILD UP AMORPHOUS-VECTOR
            //     if REQUIRED MAKE ABSORPTION CORRECTION
            //-----COMPUTE TWO-THETA
            TH = g1->THMIN + static_cast<double>(I-1) * g1->STEP;
            datax->VAR[I]=datax->Y[I];
            //-----COMPUTE SAMPLE ABSORPTION CORRECTION
            if(cntrls->IAS == 1)
            {
                ABSORP(g1->TMV,cntrls->SW,TH,&ABC);
                datax->Y[I] = datax->Y[I] / ABC;
                datax->VAR[I]=datax->VAR[I]/ABC;
            }
        }
        goto L484;


    }



    datax->NPTS = static_cast<int>((g1->THMAX-g1->THMIN)/g1->STEP+1.5);
    if (datax->NPTS > IDSZ)
    {
        file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << datax->NPTS+IPTS << " POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }

L484:
    if (jnk->NBCKGD-1 < 0)
    {
        TH=g1->THMIN-g1->STEP;
        for (I=1; I <= datax->NPTS; ++I)
        {
            TH=TH+g1->STEP;
            THX=TH/g1->BKPOS-1.0;
            datax->BK[I]=params->GLB_[2-1];
            for (J=2; J <= 6; ++J) datax->BK[I]=datax->BK[I] + params->GLB_[J+1-1] * pow(THX,J-1);
        }
        goto L477;
    }
    else if (jnk->NBCKGD-1 == 0)
    {
        file3.open(file3name.data());
        I=0;
        while ( I > datax->NPTS)
        {
            getline(file3,s);
            for (J = 1; J <= 8; ++J)
            {
                s1 = s.substr((J-1)*8,7);
                if ( !s1.empty() )
                {
                    ++I;
                    stringstream(s1) >> datax->BK[I];
                }
            }
        }
        goto L477;
    }
    else
    {
        goto L476;
    }



L476:
    DIFB=jnk->POS[1]-g1->THMIN;
    if (DIFB < 0.0)
    {
        NBX=static_cast<int>((jnk->POS[2]-g1->THMIN)/g1->STEP+1.5);
        datax->BK[1]=jnk->BCK[1]-DIFB/(jnk->POS[2]-jnk->POS[1])*(jnk->BCK[2]-jnk->BCK[1]);
        BSTEP=g1->STEP/(jnk->POS[2]-jnk->POS[1])*(jnk->BCK[2]-jnk->BCK[1]);
        for (I=2; I <= NBX; ++I) datax->BK[I]=datax->BK[I-1]+BSTEP;
        NXX=2;
    }
    else
    {
        NBX=static_cast<int>(DIFB/g1->STEP+2.5);
        for (I=1; I <= NBX; ++I) datax->BK[I]=jnk->BCK[1];
        NXX=1;
    }





    NBC=jnk->NBCKGD-1;
    if(jnk->POS[jnk->NBCKGD] >= g1->THMAX)goto L4767;
    jnk->POS[jnk->NBCKGD+1]=g1->THMAX;
    jnk->BCK[jnk->NBCKGD+1]=jnk->BCK[jnk->NBCKGD];
    NBC=NBC+1;
L4767:
    for (J=NXX; J <= NBC; ++J)
    {
        BSTEP=g1->STEP*(jnk->BCK[J+1]-jnk->BCK[J])/(jnk->POS[J+1]-jnk->POS[J]);
        NINC=static_cast<int>((jnk->POS[J+1]-jnk->POS[J])/g1->STEP+1.5);
        N2X=min(datax->NPTS,NBX+NINC);
        for (I=NBX; I <= N2X; ++I) datax->BK[I]=datax->BK[I-1]+BSTEP;
        NBX=N2X;
        if(NBX == datax->NPTS)goto L477;
    }

L477:
    for (K=1; K <= cntrls->NPHASE; ++K) refls->ICR[K]=0;


    // start loop on phase
    for (K=1; K <= cntrls->NPHASE; ++K)
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
        ISTEST = cntrls->ISPHASE;
        if (K == cntrls->ISPHASE)
        {
            if (phases[K].WTIS == 0.0)
            {
                cntrls->ISPHASE = 0;
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
        if (K == cntrls->ISPHASE)
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
        simoper->ISIMOP=1;
        RTMT(&MULTX_,Y_, &cntrls->IPL1,NCTR_,&K);
        multip->XMLTP[K]=multip->MLTPHASE;
        file6 << "The multiplicity of the general site is " << setw(3) << multip->MLTPHASE << endl;
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
            parac->ATEXT[I+IOF] = s.substr(0,4);
            stringstream(s.substr(5,4)) >> multip->MURT[I+IOF];
            parac->NTYP[I+IOF] = s.substr(10,4);
            stringstream(s.substr(16,8)) >> params->XL[I+IOF][1];
            stringstream(s.substr(24,8)) >> params->XL[I+IOF][2];
            stringstream(s.substr(32,8)) >> params->XL[I+IOF][3];
            stringstream(s.substr(40,8)) >> params->XL[I+IOF][4];
            stringstream(s.substr(48,8)) >> params->XL[I+IOF][5];

            getline(file5,s);
            stringstream(s.substr(16,8)) >> params->A[I+IOF][1];
            stringstream(s.substr(24,8)) >> params->A[I+IOF][2];
            stringstream(s.substr(32,8)) >> params->A[I+IOF][3];
            stringstream(s.substr(40,8)) >> params->A[I+IOF][4];
            stringstream(s.substr(48,8)) >> params->A[I+IOF][5];

            getline(file5,s);
            stringstream(s.substr(0,8)) >> params->XL[I+IOF][6];
            stringstream(s.substr(8,8)) >> params->XL[I+IOF][7];
            stringstream(s.substr(16,8)) >> params->XL[I+IOF][8];
            stringstream(s.substr(24,8)) >> params->XL[I+IOF][9];
            stringstream(s.substr(32,8)) >> params->XL[I+IOF][10];
            stringstream(s.substr(40,8)) >> params->XL[I+IOF][11];

            getline(file5,s);
            stringstream(s.substr(0,8)) >> params->A[I+IOF][6];
            stringstream(s.substr(8,8)) >> params->A[I+IOF][7];
            stringstream(s.substr(16,8)) >> params->A[I+IOF][8];
            stringstream(s.substr(24,8)) >> params->A[I+IOF][9];
            stringstream(s.substr(32,8)) >> params->A[I+IOF][10];
            stringstream(s.substr(40,8)) >> params->A[I+IOF][11];
        }

        // convert lower case to upper case
        for (ITN = 1; ITN <= N; ++ITN)
        {
            for (ITIPO = 1; ITIPO <= 2; ++ITIPO)
            {
                for (IK=1; IK <= 26; ++IK) if (parac->NTYP[ITN+IOF][ITIPO] == LOWER[IK]) parac->NTYP[ITN+IOF][ITIPO]=UPPER[IK];
            }
        }        
        // finish conversion

        for (I=1; I <= N; ++I)
        {
            file6
                << setw(4) << parac->ATEXT[I+IOF] << " " << setw(4) << multip->MURT[I+IOF]
                << " " << parac->NTYP[I+IOF] << "        "
                << setw(10) << setprecision(5) << params->XL[I+IOF][1]
                << setw(10) << setprecision(5) << params->XL[I+IOF][2]
                << setw(10) << setprecision(5) << params->XL[I+IOF][3]
                << setw(10) << setprecision(5) << params->XL[I+IOF][4]
                << setw(10) << setprecision(5) << params->XL[I+IOF][5] << endl
                << "                      "
                << setw(10) << setprecision(5) << params->XL[I+IOF][6]
                << setw(10) << setprecision(5) << params->XL[I+IOF][7]
                << setw(10) << setprecision(5) << params->XL[I+IOF][8]
                << setw(10) << setprecision(5) << params->XL[I+IOF][9]
                << setw(10) << setprecision(5) << params->XL[I+IOF][10]
                << setw(10) << setprecision(5) << params->XL[I+IOF][11] << endl;
        }
        //  !cp jun 96 ... (CONVERT sof MULTIPLICITY)(also changed in OUTPTR)
        for (I=1; I <= N; ++I)
        {
            if(static_cast<int>(params->A[I+IOF][5]/10) != 0 && params->XL[I+IOF][5] == 0)
            {
                params->XL[I+IOF][5]=1e-6;
            }
            params->XL[I+IOF][5] = params->XL[I+IOF][5] * multip->MURT[I+IOF] / multip->XMLTP[K];
        }

        // params->PAR[K][21] introduced below. It is for the term cot**2 in the pv-5 FWHM !cp Aug 95
        // line 11-5, line 11-6, line 11-7, line 11-8 and line 11-9
        getline(file5,s);			// S  O_B (line 11-5)
        phases[K-1].PAR[0] = s.substr(0,8);
        phases[K-1].PAR[1] = s.substr(8,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K-1].PAR[0].codeword;
        stringstream(s.substr(8,8)) >> phases[K-1].PAR[1].codeword;

        getline(file5,s);			// FWHM (line 11-6)
        phases[K-1].PAR[2] = s.substr(0,8);
        phases[K-1].PAR[3] = s.substr(8,8);
        phases[K-1].PAR[4] = s.substr(16,8);
        phases[K-1].PAR[20] = s.substr(24,8);
        phases[K-1].PAR[19] = s.substr(32,8);
        phases[K-1].PAR[14] = s.substr(40,8);
        phases[K-1].PAR[15] = s.substr(48,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K-1].PAR[2].codeword;
        stringstream(s.substr(8,8)) >> phases[K-1].PAR[3].codeword;
        stringstream(s.substr(16,8)) >> phases[K-1].PAR[4].codeword;
        stringstream(s.substr(24,8)) >> phases[K-1].PAR[20].codeword;
        stringstream(s.substr(32,8)) >> phases[K-1].PAR[19].codeword;
        stringstream(s.substr(40,8)) >> phases[K-1].PAR[14].codeword;
        stringstream(s.substr(48,8)) >> phases[K-1].PAR[15].codeword;

        getline(file5,s);
        phases[K-1].PAR[5] = s.substr(0,8);
        phases[K-1].PAR[6] = s.substr(8,8);
        phases[K-1].PAR[7] = s.substr(16,8);
        phases[K-1].PAR[8] = s.substr(24,8);
        phases[K-1].PAR[9] = s.substr(32,8);
        phases[K-1].PAR[10] = s.substr(40,8);

        getline(file5,s);			// Unit cell (line 11-7)
        stringstream(s.substr(0,8)) >> phases[K-1].PAR[5].codeword;
        stringstream(s.substr(8,8)) >> phases[K-1].PAR[6].codeword;
        stringstream(s.substr(16,8)) >> phases[K-1].PAR[7].codeword;
        stringstream(s.substr(24,8)) >> phases[K-1].PAR[8].codeword;
        stringstream(s.substr(32,8)) >> phases[K-1].PAR[9].codeword;
        stringstream(s.substr(40,8)) >> phases[K-1].PAR[10].codeword;

        getline(file5,s);			// G1 G2 P (line 11-8)
        phases[K-1].PAR[11] = s.substr(0,8);
        phases[K-1].PAR[12] = s.substr(8,8);
        phases[K-1].PAR[13] = s.substr(16,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K-1].PAR[11].codeword;
        stringstream(s.substr(8,8)) >> phases[K-1].PAR[12].codeword;
        stringstream(s.substr(16,8)) >> phases[K-1].PAR[13].codeword;

        getline(file5,s);			// NA NB NC (line 11-91)
        phases[K-1].PAR[16] = s.substr(0,8);
        phases[K-1].PAR[17] = s.substr(8,8);
        phases[K-1].PAR[18] = s.substr(16,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K-1].PAR[16].codeword;
        stringstream(s.substr(8,8)) >> phases[K-1].PAR[17].codeword;
        stringstream(s.substr(16,8)) >> phases[K-1].PAR[18].codeword;

        getline(file5,s);			// NA NB NC HS (line 11-93)
        phases[K-1].PAR[23] = s.substr(0,8);
        phases[K-1].PAR[24] = s.substr(8,8);
        phases[K-1].PAR[25] = s.substr(16,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K-1].PAR[23].codeword;
        stringstream(s.substr(8,8)) >> phases[K-1].PAR[24].codeword;
        stringstream(s.substr(16,8)) >> phases[K-1].PAR[25].codeword;

        getline(file5,s);			// s-PVII (line 11-95);
        phases[K-1].PAR[26] = s.substr(0,8);

        getline(file5,s);
        stringstream(s.substr(0,8)) >> phases[K-1].PAR[26].codeword;

        // checking for zeros if TCH-PV is being used
        if(cntrls->NPROF == _TCHZ)
        {
            if(phases[K-1].PAR[14] == 0)phases[K-1].PAR[14]=1e-8;
            if(phases[K-1].PAR[15] == 0)phases[K-1].PAR[15]=1e-8;
            if(phases[K-1].PAR[2] == phases[K-1].PAR[3] && phases[K-1].PAR[3] == phases[K-1].PAR[4] && phases[K-1].PAR[4] == phases[K-1].PAR[19] && phases[K-1].PAR[2] == 0) phases[K-1].PAR[4]=1e-8;
        }
        // checking for zeros in PV #5
        if(cntrls->NPROF == _pseudoVoigt)
        {
            if(phases[K-1].PAR[16] == 0)phases[K-1].PAR[16]=1e-6;
        }
        //      if(int(phases[K-1].PAR[19].APAR/10) != 0 && params->PAR[K][20] == 0)params->PAR[K][20]=1e-9
        // !cp Aug 95 introducing params->PAR[K][21]=ct
        // CHECKING FOR NON-REFINABLE PARAMETERS
        //CCC                        FOR  CT  IN TCHZ AND SPVII FUNCTIONS
        if(cntrls->NPROF == _TCHZ || cntrls->NPROF == _SplitPearsonVII)
        {
            phases[K-1].PAR[20]=0.0;
            if(phases[K-1].PAR[20].codeword != 0.0)
            {
                cout << "NON-REFINABLE PARAMETER TURNED ON (CT) WITH SPLIT PEARSON VII OR TCHZ PROFILE FUNCTION" << endl;
                file6 << "NON-REFINABLE PARAMETER TURNED ON (CT) WITH SPLIT PEARSON VII OR TCHZ PROFILE FUNCTION" << endl;
                DBWSException("");
            }
        }
        //CCC                             FOR  X,Y,Z IN NON-TCHZ FUNCTION
        if(cntrls->NPROF != _TCHZ)
        {
            if(phases[K-1].PAR[19] != 0 || phases[K-1].PAR[15] != 0 || phases[K-1].PAR[14].codeword != 0)
            {
                file6 << "     NON-REFINABLE PARAMETER RESET TO ZERO (Z,X,Y) FOR" << endl
                    << "     NON TCHZ PROFILE FUNCTION" << endl;
                phases[K-1].PAR[19]=0.0;
                phases[K-1].PAR[15]=0.0;
                phases[K-1].PAR[14]=0.0;
            }
            if(phases[K-1].PAR[19].codeword != 0 || phases[K-1].PAR[15].codeword != 0 || phases[K-1].PAR[14].codeword != 0)
            {
                cout << "NON-REFINABLE PARAMETER TURNED ON (Z,X,Y) WITH NOT TCHZ PROFILE FUNCTION" << endl;
                file6 << "NON-REFINABLE PARAMETER TURNED ON (Z,X,Y) WITH NOT TCHZ PROFILE FUNCTION" << endl;
                DBWSException("");
            }
        }
        //CCCC                            FOR  RIET_ASYM,X,Y,Z,CT IN SPVII
        if (cntrls->NPROF == _SplitPearsonVII)
        {
            for (KKS=1; KKS <= 3; ++KKS)
            {
                if(phases[K-1].PAR[13+KKS-1] != 0.0)
                {
                    phases[K-1].PAR[13+KKS-1]=0.0;
                    file6 << " RIET_ASYM X Y RESET TO ZERO FOR SPVII FUNCTION" << endl;
                }
                if(phases[K-1].PAR[13+KKS-1].codeword != 0.0)
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
                if ( phases[K-1].PAR[19+KKS-1] != 0.0 )
                {
                    phases[K-1].PAR[19+KKS-1]=0.0;
                    file6 << "CT,Z  RESET TO ZERO FOR SPVII FUNCTION" << endl;
                }
                if ( phases[K-1].PAR[19+KKS-1].codeword != 0.0 )
                {
                    cout << "NON-REFINEABLE PARAMETER TURNED ON (CT,Z), WITH SPLIT PEARSON VII PROFILE FUNCTION" << endl;
                    file6 << "NON-REFINEABLE PARAMETER TURNED ON (CT,Z), WITH SPLIT PEARSON VII PROFILE FUNCTION" << endl;
                    DBWSException("");
                }
            }
        }
        //CCCCC             FOR HIGH SIDE PARAMETERS IN NON-PVII PROFILE FUNCTION
        if (cntrls->NPROF != _SplitPearsonVII)
        {
            for (KK=1; KK <= 4; ++KK)
            {
                if ( phases[K-1].PAR[23+KK-1] != 0.0 )
                {
                    phases[K-1].PAR[23+KK-1]=0.0;
                    file6 << "NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY RESET TO ZERO FOR A, NON-SPLIT PEARSON VII FUNCTION" << endl;
                }
                if(phases[K-1].PAR[23+KK-1] != 0.0)
                {
                    cout << "NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY ONLY REFINABLE FOR A SPLIT PEARSON VII FUNCTION" << endl;
                    file6 << "NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY ONLY REFINABLE FOR A SPLIT PEARSON VII FUNCTION" << endl;
                }
            }
        }
        //CCCCC END OF CHECKING NON-REFINABLE PARAMETERS
        file6 << "OVERALL SCALE FACTOR=" << scientific << setw(12) << setprecision(6) << phases[K-1].PAR[0] << fixed << endl
            << "OVERALL TEMP. FACTOR=" << setw(12) << setprecision(5) << phases[K-1].PAR[1] << endl;

        file6 << "DIRECT CELL PARAMETERS:" << endl
        << "a = " << setw(9) << setprecision(4) << phases[K-1].PAR[5] << endl
        << "b = " << setw(9) << setprecision(4) << phases[K-1].PAR[6] << endl
        << "c = " << setw(9) << setprecision(4) << phases[K-1].PAR[7] << endl
        << "α = " << setw(9) << setprecision(4) << phases[K-1].PAR[8] << endl
        << "β = " << setw(9) << setprecision(4) << phases[K-1].PAR[9] << endl
        << "γ = " << setw(9) << setprecision(4) << phases[K-1].PAR[10] << endl;

        file6 << "PREFERRED ORIENTATION PARAMETERS="
            << setw(7) << setprecision(3) << phases[K-1].PAR[11]
        << setw(7) << setprecision(3) << phases[K-1].PAR[12] << endl
        << "ASYMMETRY PARAMETER="
            << setw(8) << setprecision(4) << phases[K-1].PAR[13] << endl;

        //-----CHECK FOR SPLIT PEARSON PROFILE
        if (cntrls->NPROF == _SplitPearsonVII)
        {
            file6 << "LOW SIDE EXPONENT COEFFICIENTS="
                << setw(12) << setprecision(4) << phases[K-1].PAR[16]
            << setw(12) << setprecision(4) << phases[K-1].PAR[17]
            << setw(12) << setprecision(4) << phases[K-1].PAR[18]
            << "HIGH SIDE EXPONENT COEFFICIENTS="
                << setw(12) << setprecision(4) << phases[K-1].PAR[23]
            << setw(12) << setprecision(4) << phases[K-1].PAR[24]
            << setw(12) << setprecision(4) << phases[K-1].PAR[25] << endl
            << "SPLIT PEARSON VII ASYMMETRY PARAMETER="
                << setw(8) << setprecision(4) << phases[K-1].PAR[26] << endl;
        }
        else
        {
            //-----if NOT THE SPLIT PEARSON VII PROFILE
            file6 << "MIXING PARAMETERS = "
                << scientific
                << setw(10) << setprecision(3) << phases[K-1].PAR[16] << " "
                << setw(10) << setprecision(3) << phases[K-1].PAR[17] << " "
                << setw(10) << setprecision(3) << phases[K-1].PAR[18] << fixed << endl;
        }
        file6 << "FWHM PARAMETERS (U,V,W,CT,Z,X,Y)="
            << setw(9) << setprecision(4) << phases[K-1].PAR[2]
        << setw(9) << setprecision(4) << phases[K-1].PAR[3]
        << setw(9) << setprecision(4) << phases[K-1].PAR[4]
        << setw(9) << setprecision(4) << phases[K-1].PAR[20]
        << setw(9) << setprecision(4) << phases[K-1].PAR[19]
        << setw(9) << setprecision(4) << phases[K-1].PAR[14]
        << setw(9) << setprecision(4) << phases[K-1].PAR[15] << endl;
        cellx->A=phases[K-1].PAR[5];
        cellx->B=phases[K-1].PAR[6];
        cellx->C=phases[K-1].PAR[7];
        cellx->ALPHA=phases[K-1].PAR[8];
        cellx->BETA=phases[K-1].PAR[9];
        cellx->GAMMA=phases[K-1].PAR[10];
        CELL2(K,LAMDAM);
        // ************************************** !cp ap 97 (from It code)
        if(cntrls->FONDO == 1 || cntrls->FONDO == 2)
        {
            bkgscale->SCABKG[K] = volume->GCOM[K] * phases[K-1].PAR[0];
        }
        // ************************************************************
        file6 << "CELL VOLUME PHASE(" << setw(2) << K << " ) = " << setw(12) << setprecision(4) << volume->VOLI[K] << endl;
        for (I=1; I <= 6; ++I) dc->SAVE[K][I]=phases[K-1].PAR[I+5-1];

        // change to reci
        for (I=1; I <= 3; ++I) phases[K-1].PAR[I+5-1]=cellx->AL[I][I];			// a,b and c cell parameters
        phases[K-1].PAR[8]=cellx->AL[2][3];			// alpha cell
        phases[K-1].PAR[9]=cellx->AL[1][3];			// beta cell
        phases[K-1].PAR[10]=cellx->AL[1][2];			// gamma cell
        file6 << "***Coding of variables***" << endl
            << "ATOM                        X         Y         Z         B         So" << endl
            << "                          B11       B22       B33       B12       B13       B23" << endl;
        //-----PRINT CODEWORDS FOR ATOMIC PARAMETERS
        for (I=1; I <= N; ++I)
        {
            file6 << setw(4) << parac->ATEXT[I+IOF]
            << "                 "
                << setw(10) << setprecision(2) << params->A[I+IOF][1]
            << setw(10) << setprecision(2) << params->A[I+IOF][2]
            << setw(10) << setprecision(2) << params->A[I+IOF][3]
            << setw(10) << setprecision(2) << params->A[I+IOF][4]
            << setw(10) << setprecision(2) << params->A[I+IOF][5] << endl
                << "                     "
                << setw(10) << setprecision(2) << params->A[I+IOF][6]
            << setw(10) << setprecision(2) << params->A[I+IOF][7]
            << setw(10) << setprecision(2) << params->A[I+IOF][8]
            << setw(10) << setprecision(2) << params->A[I+IOF][9]
            << setw(10) << setprecision(2) << params->A[I+IOF][10]
            << setw(10) << setprecision(2) << params->A[I+IOF][11] << endl;
        }
        for (I=1; I <= N; ++I)
        {
            for (J=1; J <= 11; ++J)
            {
                X=params->A[I+IOF][J];
                IYY=static_cast<int>(abs(X)/10.0);
                if(IYY > MSZ) goto L99996;
                params->LP[I+IOF][J]=IYY;
                params->A[I+IOF][J]=(abs(X)-10.*static_cast<double>(IYY))*sign(X);
            }
        }
        //-----PRINT CODEWORDS FOR PROFILE PARAMETERS
        file6 << "OVERALL SCALE FACTOR=" << setw(8) << setprecision(2) << phases[K-1].PAR[0].codeword << endl
            << "OVERALL TEMP. FACTOR=" << setw(8) << setprecision(2) << phases[K-1].PAR[1].codeword << endl
            << "DIRECT CELL PARAMETERS="
            << setw(8) << setprecision(2) << phases[K-1].PAR[5].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[6].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[7].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[8].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[9].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[10].codeword << endl
            << "PREFERRED ORIENTATION PARAMETERS="
            << setw(8) << setprecision(2) << phases[K-1].PAR[11].codeword
        << setw(8) << setprecision(2) << phases[K-1].PAR[12].codeword << endl
            << "ASYMMETRY PARAMETER="
            << setw(8) << setprecision(4) << phases[K-1].PAR[13].codeword << endl;

        // !cp ap 97 (from It code)
        if(cntrls->FONDO == 1 && (phases[K-1].PAR[1] != 0.0 || phases[K-1].PAR[1].codeword != 0.0)) goto L88888;
        if(cntrls->FONDO == 2 && phases[K-1].PAR[1] == 0.0 && phases[K-1].PAR[1].codeword == 0.0) goto L88889;
        //-----CHECK FOR SPLIT PEARSON PROFILE
        if (cntrls->NPROF == _SplitPearsonVII)
        {
            file6 << "LOW SIDE EXPONENT COEFFICIENTS="
                << setw(8) << setprecision(2) << phases[K-1].PAR[16].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[17].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[18].codeword << endl
                << "HIGH SIDE EXPONENT COEFFICIENTS="
                << setw(8) << setprecision(2) << phases[K-1].PAR[23].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[24].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[25].codeword << endl
                << "SPLIT PEARSON VII ASSYMETRY PARAMETER="
                << setw(8) << setprecision(2) << phases[K-1].PAR[26].codeword << endl;
        }
        else
        {
            //-----if NOT THE SPLIT PEARSON VII PROFILE
            file6
                << "MIXING PARAMETERS = "
                << setw(10) << setprecision(3) << phases[K-1].PAR[16].codeword
                << setw(10) << setprecision(3) << phases[K-1].PAR[17].codeword
                << setw(10) << setprecision(3) << phases[K-1].PAR[18].codeword << endl;
        }
        file6
            << "FWHM PARAMETERS (U,V,W,CT,Z,X,Y)="
            << setw(8) << setprecision(2) << phases[K-1].PAR[2].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[3].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[4].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[20].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[19].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[14].codeword
            << setw(8) << setprecision(2) << phases[K-1].PAR[15].codeword << endl;
        for (I=1; I <= 27; ++I)
        {
            X=phases[K-1].PAR[I-1].codeword;
            IYY=static_cast<int>(abs(X)/10.0);
            phases[K-1].PAR[I-1].L=IYY;
            phases[K-1].PAR[I-1].codeword=(abs(X)-10.*static_cast<double>(IYY))*sign(X);
        }
        //L151:
        LOOKUP(K,phases[K].AtomCount,jnk->NSCAT,IXRAY,cntrls->JOBTYP);
        if(cntrls->FONDO == 1 || cntrls->FONDO == 2) FINDC(K,jnk->NSCAT);
        g1->U=phases[K-1].PAR[2];
        g1->V=phases[K-1].PAR[3];
        g1->W=phases[K-1].PAR[4];
        g1->ZZZ = phases[K-1].PAR[19];
        g1->UC = phases[K-1].PAR[20];
        if (cntrls->NPROF == _TCHZ)
        {
            g1->ULOR=phases[K-1].PAR[14];
            g1->VLOR=phases[K-1].PAR[15];
        }
        REFGEN(K,params->GLB_[1-1],params->GLB_[10-1],params->GLB_[11-1],phases[K-1].PAR[11],NCTR_);
        simoper->ISIMOP=0;
        RTMT(&MULTX_,Y_,&cntrls->IPL1,NCTR_,&K);
        ICY=1;
        ICZ=0;
        for (IIPHAS=1; IIPHAS <= K; ++IIPHAS) ICZ = ICZ + refls->ICR[IIPHAS];
        if(K >= 2) ICY=1+ICZ-refls->ICR[K];
        if(cntrls->LST1 != 1)goto L479;
        IXDEL=0;
        for (IXX=ICY; IXX <= ICZ; ++IXX)
        {
            IX=IXX-IXDEL;
            IRL=(refls->IREFS[IXX] % 256)-128;
            IRK=((refls->IREFS[IXX]/256) % 256)-128;
            IRH=((refls->IREFS[IXX]/(256*256)) % 256)-128;
            IRC=(refls->IREFS[IXX]/(256*256*256)) % 8;
            MLTT = static_cast<int>(refls->FMGNTD[IXX]);
            if (cntrls->NPROF == _TCHZ)
            {
                if((IX-1 % 60) == 0) file6 << "NO.  CODE    H   K   L  MULT   HW     POSN      FACTOR       HWL       HWG     ETA" << endl;
                file6 << setw(4) << IX
                    << setw(4) << IRC << "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(6) << MLTT
                    << setw(8) << setprecision(3) << refls->REFS[IXX][1]
                << setw(8) << setprecision(3) << refls->REFS[IXX][2]
                << setw(10) << setprecision(6) << refls->REFS[IXX][3]
                << setw(8) << setprecision(3) << refls->HALFL[IXX]
                << setw(8) << setprecision(3) << refls->HALFG[IXX]
                << setw(8) << setprecision(3) << refls->GAM[IXX] << endl;
            }
            else if (cntrls->NPROF == _SplitPearsonVII)
            {
                refls->FWHM[IXX][1]=2.0*(refls->REFS[IXX][1])*phases[K-1].PAR[26]/(1.0+phases[K-1].PAR[26]);
                refls->FWHM[IXX][2]=2.0*(refls->REFS[IXX][1])/(1.0+phases[K-1].PAR[26]);
                if((IX-1 % 60) == 0) file6 << "NO.  CODE    H   K   L  MULT      HWL    HWH     FWHM   POSN    FACTOR" << endl;
                file6 << setw(4) << IX
                    << setw(4) << IRC << "   "
                    << setw(4) << IRH
                    << setw(4) << IRK
                    << setw(4) << IRL
                    << setw(6) << MLTT
                    << setw(8) << setprecision(3) << refls->FWHM[IXX][1]
                << setw(8) << setprecision(3) << refls->FWHM[IXX][2]
                << setw(8) << setprecision(3) << refls->REFS[IXX][1]
                << setw(8) << setprecision(3) << refls->REFS[IXX][2]
                << setw(10) << setprecision(6) << refls->REFS[IXX][3] << endl;
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
                    << setw(8) << setprecision(3) << refls->REFS[IXX][1]
                << setw(8) << setprecision(3) << refls->REFS[IXX][2]
                << setw(10) << setprecision(6) << refls->REFS[IXX][3] << endl;
            }
        }
L479:
        for (IX=ICY; IX <= ICZ; ++IX) refls->REFS[IX][3]=refls->REFS[IX][3]*refls->FMGNTD[IX]; // double(MLTT(IX))
        //L81:;
    }
    // end of great loop on phases

    // FIND CODEWORDS FOR GLOBAL PARAMETERS: LGLB(I) & AGLB(I)
    for (I=1; I <= 20; ++I)
    {
        X=params->GLB_[I-1].codeword;
        IYY=static_cast<int>(abs(X)/10.0);
        params->GLB_[I-1].L=IYY;
        params->GLB_[I-1].codeword=(abs(X)-10.*static_cast<double>(IYY))*sign(X);
    }
    qpainit();
    if (cntrls->JOBTYP > 2) return;
    NATOMS = 0;
    for (K=1; K <= cntrls->NPHASE; ++K) NATOMS = NATOMS + phases[K].AtomCount;
    //     COUNT NO. OF USES OF SAME LOCATION IN NORMAL MATRIX
    for (I=1; I <= cntrls->MAXS; ++I)
    {
        LCOUNT = 0;
        for (J=1; J <= 20; ++J) if ( I == params->GLB_[J-1].L ) LCOUNT = LCOUNT + 1;
        for (K=1; K <= cntrls->NPHASE; ++K)
        {
            for (J=1; J <= 27; ++J) if ( I == phases[K-1].PAR[J-1].L )  LCOUNT = LCOUNT + 1;
        }
        for (K=1; K <= NATOMS; ++K)
        {
            for (J=1; J <= 11; ++J) if(I == params->LP[K][J]) LCOUNT = LCOUNT + 1;
        }
        if (LCOUNT >= 2) file6 << "******* CODEWORD " << setw(3) << I << " is used " << setw(3) << LCOUNT << " times **" << endl;
    }
    //     END COUNT NO. OF USES OF SAME LOCATION IN NORMAL MATRIX
    //     CHECK FOR HOLES IN THE NORMAL MATRIX
    ISTOP=0;
    for (I=1; I <= cntrls->MAXS; ++I)
    {
        for (J=1; J <= 20; ++J) if (I == params->GLB_[J-1].L) goto L7200;
        //  HOLES IN PHASE PARAMETER
        for (K=1; K <= cntrls->NPHASE; ++K)
        {
            for (J=1; J <= 27; ++J) if(I == phases[K-1].PAR[J-1].L) goto L7200;
        }
        // HOLES IN ATOMS PARAMETERS
        for (K=1; K <= NATOMS; ++K)
        {
            for (J=1; J <= 11; ++J) if(I == params->LP[K][J]) goto L7200;
        }
        file6 << "     ***** HOLE IN THE MATRIX ******. ELEMENT " << setw(3) << I << " IN THE NORMAL MATRIX IS MISSING" << endl;
        cout << "     ***** HOLE IN THE MATRIX ******. ELEMENT " << setw(3) << I << " IN THE NORMAL MATRIX IS MISSING" << endl;
        ISTOP = 1;
L7200:;
    }
    if (ISTOP == 1) DBWSException("");
    return;
    //L99999:
    //	DBWSException("END OF FILE TAPE5");
    //L99998:
    //	DBWSException("END OF FILE TAPE4");
    //L99997:
    //	DBWSException("END OF FILE TAPE3");
L99996:
    DBWSException("MATRIX SIZE IS TOO SMALL");
L99995:
    DBWSException("MAXS > MSZ");
L88888:
    DBWSException("WITH FONDO=1 YOU MUST USE ISOTROPIC THERMAL FACTORS");
L88889:
    DBWSException("if FONDO=2 YOU MUST USE THE OVERALL THERMAL FACTOR");

//L88088:
//    if (file4.is_open()) file4.close();
//    if (file5.is_open()) file5.close();
//    if (file6.is_open()) file6.close();
//    if (file8o.is_open()) file8o.close(); //	close(8,status='DELETE');
//    if (file8i.is_open()) file8i.close();
//    DBWSException("");
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
    if (THMIN1 != g1->THMIN)
    {
        file6 << "AMORPHOUS THMIN="
            << setw(8) << setprecision(2) << THMIN1
            << "IS DIFFERENT FROM THE DATA THMIN="
            << setw(8) << setprecision(2) << g1->THMIN << endl;
    }
    if (STEP1 != g1->STEP)
    {
        file6 << "AMORPHOUS STEP="
            << setw(8) << setprecision(2) << STEP1
            << "IS DIFFERENT FROM THE DATA STEP="
            << setw(8) << setprecision(2) << g1->STEP << endl;
    }
    if (THMAX1 != g1->THMAX)
    {
        file6 << "AMORPHOUS THMAX="
            << setw(8) << setprecision(2) << THMAX1
            << "IS DIFFERENT FROM THE DATA THMAX="
            << setw(8) << setprecision(2) << g1->THMAX << endl;
    }
    //      read the rest of the file in free format
    for (I=1; I <= datax->NPTS; ++I) file11 >> datax->AMORPHOUS[I];
    for (I = 1; I <= datax->NPTS; ++I)
    {
        TH=g1->THMIN+(I-1)*g1->STEP;
        //------AMORPHOUS CORRECTION FOR AIR SCATTERING
        //       CALL ARIA(TMV1,SW1,FI1,TH,SCA)
        //       AMORPHOUS(I)=0.0
        if (cntrls->IAS == 1)
        {
            ABSORP(TMV1,SW1,TH,&ABC);
            datax->AMORPHOUS[I] = datax->AMORPHOUS[I] / ABC;
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
        app.file4name = string(argv[1]);
        app.file4.open(app.file4name.data());

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



                if ( app.file4.is_open() ) app.file4.close();
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
