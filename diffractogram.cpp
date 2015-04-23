#include "diffractogram.h"

Diffractogram::Diffractogram()
{
}

void Diffractogram::CalcLamdaM(void)
{
    LAMDAM=(RATIO[1]*LAMDA[1]+RATIO[2]*LAMDA[2])/(RATIO[1]+RATIO[2]);
}

int Diffractogram::SearchWavelength(void)
{
    int i,tmp = 0;

    for (i=1; i <= 10; ++i)
    {
        if(1.03*LAMDA[1] > XRYZ[i])
        {
            tmp=i;
            break;
        }
    }
    return tmp;
}

void Diffractogram::ReadDBWS(void)
{
    string s,s1;
    int i,j;

    // read DBWS formated
    getline(file4,s);
    stringstream(s.substr(0*8,8)) >> THMIN;
    stringstream(s.substr(1*8,8)) >> STEP;
    stringstream(s.substr(2*8,8)) >> THMAX;
    DATAID = s.substr(3*8,56);
    cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << THMIN
        << ", STOP =" << setw(8) << setprecision(3) << THMAX
        << ", STEP =" << setw(8) << setprecision(3) << STEP << endl;
    *file6 << "    DATA ID " << DATAID << endl;
    cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << THMIN
        << ", STOP =" << setw(8) << setprecision(3) << THMAX
        << ", STEP =" << setw(8) << setprecision(3) << STEP << endl;
    NPTS=static_cast<int>((THMAX-THMIN)/STEP+1.5);
    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
            << setw(5) << NPTS << " POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
            << setw(5) << NPTS << " POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        //DBWSException("TOO MANY DATA POINTS");
    }
    i=0;
    while ( i < NPTS)
    {
        getline(file4,s);
        for (j = 1; j <= 8; ++j)
        {
            s1 = s.substr((j-1)*8,7);
            if ( !s1.empty() )
            {
                ++i;
                stringstream(s1) >> Y[i];
            }
        }
    }
    XMAXINT = 0.0;            //     CHECK FOR NON-ZERO COUNTS AT ALL POINTS
}

void Diffractogram::ReadDBWSF(void)
{
    string s;
    int i;

    getline(file4,s);
    stringstream(s.substr(0*8,8)) >> THMIN;
    stringstream(s.substr(1*8,8)) >> STEP;
    stringstream(s.substr(2*8,8)) >> THMAX;
    DATAID = s.substr(3*8,56);
    cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << THMIN
        << ", STOP =" << setw(8) << setprecision(3) << THMAX
        << ", STEP =" << setw(8) << setprecision(3) << STEP << endl;
    *file6 << "    DATA ID " << DATAID << endl;
    cout << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(3) << THMIN
        << ", STOP =" << setw(8) << setprecision(3) << THMAX
        << ", STEP =" << setw(8) << setprecision(3) << STEP << endl;
    NPTS=static_cast<int>((THMAX-THMIN)/STEP+1.5);
    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
            << setw(5) << NPTS << " POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << " POINTS" << endl
            << setw(5) << NPTS << " POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (i=1; i <= NPTS; ++i)
    {
        file4 >> Y[i];
    }
    XMAXINT = 0.0;            //     CHECK FOR NON-ZERO COUNTS AT ALL POINTS
}

void Diffractogram::GSASREAD(void)
{
    int I, J, K, L4, L5, L10, L11, L12, L13;
    int LAB[20+1];
    int NCTR[IDSZ+1];
    string s,TEST;

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
    stringstream(TEST.substr(L4,L5-L4+1)) >> NPTS;			// npts
    stringstream(TEST.substr(L10,L11-L10+1)) >> THMIN;			// start*100
    stringstream(TEST.substr(L12,L13-L12+1)) >> STEP;			// step*100
    THMIN=THMIN/100;
    STEP=STEP/100;
    THMAX =THMIN+NPTS*STEP;
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    *file6 << "    DATA ID " << setw(56) << DATAID << endl;
    *file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;

    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= NPTS; ++I)
    {
        file4 >> NCTR[I] >> Y[I];			//format(10(i2,f6.0))
    }
    for (I=1; I <= NPTS; ++I)
    {
        if (NCTR[I] == 0) NCTR[I]=1;
        if (Y[I] == 0.0) Y[I]=0.0;		// TODO: ???????
    }
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
}

void Diffractogram::PHILIPSREAD(void)
{
    int I,LB1,LB2,LB3;
    string s,TEST,ATHMIN,ATHMAX;

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
    stringstream(s.substr(13,8)) >> STEP;
L30:
    getline(file4,s);
    TEST = " " + s.substr(0,56);
    for (I=1; I <= 56; ++I)
    {
        if( TEST[I] == ',') LB1=I;
    }

    if( TEST.substr(1,LB1-1-1+1) != "RawScan") goto L30;
    //                            RawScan
    stringstream(ATHMIN) >> THMIN;
    stringstream(ATHMAX) >> THMAX;

    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    *file6 << "    DATA ID " << setw(56) << DATAID << endl;
    *file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    NPTS=static_cast<int>((THMAX-THMIN)/STEP+1.5);
    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= NPTS; ++I) file4 >> Y[I];
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
}

void Diffractogram::scintag(void)
{
    int I;
    double X1,X2,X3;
    string s,ASTEP,  ATHMIN, ATHMAX;

    getline(file4,s);
    DATAID = " "+s.substr(0,56);
L10:
    getline(file4,s);
    TEST_ = " "+s.substr(0,56);
    if( TEST_.substr(1,12) != "Start Angle:") goto L10;
    //                         Start Angle:
    readasc();
    ATHMIN=TEST_.substr(LB2_,LB3_-LB2_+1);
L20:
    getline(file4,s);
    TEST_ = " "+s.substr(0,56);
    if(TEST_.substr(1,11) != "Stop Angle:")goto L20;
    //                         Stop Angle:
    readasc();
    ATHMAX=TEST_.substr(LB2_,LB3_-LB2_+1);
L30:
    getline(file4,s);
    TEST_ = " "+s.substr(0,56);
    if(TEST_.substr(1,10) != "Step Size:")goto L30;

    //                         Step Size:
    readasc();
    ASTEP=TEST_.substr(LB2_,LB3_-LB2_+1);

L70:
    getline(file4,s);
    if( TEST_.substr(1,5) != "Range")goto L70;
    stringstream(ATHMIN) >> THMIN;
    stringstream(ATHMAX) >> THMAX;
    stringstream(ASTEP) >> STEP;
    NPTS=static_cast<int>((THMAX-THMIN)/STEP+1.5);
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    *file6 << "    DATA ID " << setw(56) << DATAID << endl;
    *file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;

    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= NPTS; ++I) file4 >> X1 >> Y[I] >> X2 >> X3;
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
}

void Diffractogram::SIEMENSREAD(void)
{
    int I;
    string s,TEST, ASTEP, ANPTS, ATHMIN;


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
    stringstream(ATHMIN) >> THMIN;
    stringstream(ASTEP) >> STEP;
    stringstream(ANPTS) >> NPTS;
    THMAX=THMIN + STEP*NPTS;
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    *file6 << "    DATA ID " << setw(56) << DATAID << endl;
    *file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    for (I=1; I <= NPTS; ++I) file4 >> Y[I];
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
    //L99999:
    //	DBWSException("IS THE FILE NOT UXD SIEMENS FORMAT?");
}

void  Diffractogram::rigakuread(void)
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
    stringstream(ATHMIN) >> THMIN;
    stringstream(ASTEP) >> STEP;
    stringstream(ATHMAX) >> THMAX;
    stringstream(ANPTS) >> NPTS;
    cout << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    //file6 << "    DATA ID " << setw(56) << DATAID << endl;
    *file6 << "DATA RANGE (2THETA):  START ="
        << setw(8) << setprecision(3) << THMIN
        << ", STOP ="
        << setw(8) << setprecision(3) << THMAX
        << ", STEP ="
        << setw(8) << setprecision(3) << STEP << endl;
    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        cout  << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << setw(5) << NPTS << "POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    // read the data
    for (I=1; I <= NPTS; ++I) file4 >> Y[I];
    return;
    //L99998:
    //	DBWSException("END OF FILE unit=4");
    //L99999:
    //	DBWSException("IS THE FILE NOT RIGAKU FORMATED?");
}

void Diffractogram::ReadSynchrotronData(void)
{
    string s,DATE;
    int I,J,NRANGE, IPTS;
    double CHMBRI, ANGMAX, ANGMIN, TAUK, STPTIM, OFSTI0, OFSTI1,FLEET1, FLEET2;

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
    *file6 << "    DATA ID " << DATAID << endl;
    getline(file4,s);
    stringstream(s.substr(0*8,8)) >> NRANGE;
    stringstream(s.substr(1*8,8)) >> CHMBRI;
    stringstream(s.substr(2*8,10)) >> TAUK;

    //     NPTS IS THE COUNTER FOR TOTAL NO OF POINTS IN ALL THE RANGES
    NPTS = 0;
    for (J=1; J <= NRANGE; ++J)
    {
        //     READ INFORMATION ABOUT EACH RANGE
        getline(file4,s);
        stringstream(s.substr(0*8,8)) >> ANGMIN;
        stringstream(s.substr(1*8,8)) >> STEP;
        stringstream(s.substr(2*8,8)) >> ANGMAX;
        stringstream(s.substr(3*8,8)) >> STPTIM;
        stringstream(s.substr(4*8,8)) >> OFSTI0;
        stringstream(s.substr(5*8,8)) >> OFSTI1;
        IPTS = static_cast<int>((ANGMAX-ANGMIN)/STEP +1.5);
        //     FIND MAXIMUM AND MINIMUM TWO THETA IN ALL RANGES
        if (J == 1) THMIN = ANGMIN;
        if (J == NRANGE) THMAX = ANGMAX;
        if (J > 1)
        {
            getline(file4,s);
            IPTS = IPTS - 1;
        }
        if (NPTS+IPTS > IDSZ)
        {
            *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
                << setw(5) << (NPTS+IPTS) << " POINTS WERE INPUT" << endl
                << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
            DBWSException("TOO MANY DATA POINTS");
        }
        //     VAR(I) IS USED FOR TWO QUANTITIES, JUST TO SAVE SOME SPACE.
        if (DATE == "OCT85")
        {
            for (I=1; I <= IPTS; ++I)
            {
                getline(file4,s);
                stringstream(s.substr(24,7)) >> Y[I+NPTS];
                stringstream(s.substr(32,7)) >> VAR[I+NPTS];
            }
        }
        else if (DATE == "FEB86")
        {
            for (I=1; I <= IPTS; ++I)
            {
                getline(file4,s);
                stringstream(s.substr(24,7)) >> Y[I+NPTS];
                stringstream(s.substr(32,7)) >> VAR[I+NPTS];
            }
        }
        else if (DATE == "AUG86")
        {
            for (I=1; I <= IPTS; ++I)
            {
                getline(file4,s);
                stringstream(s.substr(51,10)) >> VAR[I+NPTS];
                stringstream(s.substr(62,10)) >> Y[I+NPTS];
            }
        }
        else if (DATE == "SRS83")
        {
            for (I=1; I <= IPTS; ++I)
            {
                getline(file4,s);
                stringstream(s.substr(37,9)) >> VAR[I+NPTS];
                stringstream(s.substr(47,9)) >> Y[I+NPTS];
            }
        }
        else if (DATE == "SRS91")
        {
            for (I=1; I <= IPTS; ++I)
            {
                getline(file4,s);
                stringstream(s.substr(29,9)) >> VAR[I+NPTS];
                stringstream(s.substr(48,9)) >> Y[I+NPTS];
            }
        }
        else
        {
            *file6 << "    WHEN AND WHERE WERE THESE DATA TAKEN. OCT85,FEB86,AUG86,SRS83,SRS91" << endl;
        }
        for (I=NPTS+1; I <= IPTS+NPTS; ++I)
        {
            Y[I] = Y[I]/STPTIM;
            FLEET1 = VAR[I]-OFSTI0;
            FLEET2= FLEET1/CHMBRI;
            FLEET2= pow((FLEET2*(1-TAUK*Y[I])) , 2.0);
            VAR[I] = (STPTIM*Y[I])/FLEET2;
            Y[I] = STPTIM*(Y[I]/(1-TAUK*Y[I])-OFSTI1)*CHMBRI/FLEET1;
            if (abs(Y[I]) < 0.01) Y[I]=1.0;
        }
        //L530:
        NPTS=NPTS+IPTS;
    }
    *file6 << "DATA RANGE (2THETA):  START = " << setw(8) << setprecision(2) << THMIN
        << ", STOP =" << setw(8) << setprecision(2) << THMAX
        << ", STEP =" << setw(8) << setprecision(2) << STEP << endl;
    cout  << "DATA RANGE (2THETA):  START = " << setw(8) << setprecision(2) << THMIN
        << ", STOP =" << setw(8) << setprecision(2) << THMAX
        << ", STEP =" << setw(8) << setprecision(2) << STEP << endl;
}

void Diffractogram::ReadNeutronData(void)
{
    int I,J, IPTS = 0;
    string s,s1;

    //     BEGIN VARIANCE CALCULATION FOR VARYING NO. OF COUNTERS AT EACH STEP
    getline(file4,s);
    stringstream(s.substr(0*8,8)) >> THMIN;
    stringstream(s.substr(1*8,8)) >> STEP;
    stringstream(s.substr(2*8,8)) >> THMAX;
    DATAID = s.substr(3*8,56);
    *file6 << "    DATA ID " << DATAID << endl;
    *file6 << "DATA RANGE (2THETA):  START =" << setw(8) << setprecision(2) << THMIN
        << ", STOP =" << setw(8) << setprecision(2) << THMAX
        << ", STEP =" << setw(8) << setprecision(2) << STEP << endl;
    NPTS = static_cast<int>((THMAX-THMIN)/STEP+1.5);
    if (NPTS > IDSZ)
    {
        *file6 << "PROGRAM CAN HANDLE " << setw(5) << IDSZ << "POINTS" << endl
            << NPTS+IPTS << " POINTS WERE INPUT" << endl
            << "INCREASE IDSZ IN PARAMETER STATEMENT" << endl;
        DBWSException("TOO MANY DATA POINTS");
    }
    I=0;
    while ( I > NPTS)
    {
        getline(file4,s);
        for (J = 1; J <= 10; ++J)
        {
            s1 = s.substr((J-1)*8,2);
            if ( !s1.empty() )
            {
                ++I;
                stringstream(s1) >> VAR[I];
            }
            s1 = s.substr((J-1)*8+2,6);
            if ( !s1.empty() )
            {
                stringstream(s1) >> Y[I];
            }
        }
    }
    for (I=1; I <= NPTS; ++I) VAR[I]=Y[I]/VAR[I];
}

void Diffractogram::ReadFileBackground()
{
    int I,J;
    string s,s1;

    file3.open(file3name.data());
    I=0;
    while ( I > NPTS)
    {
        getline(file3,s);
        for (J = 1; J <= 8; ++J)
        {
            s1 = s.substr((J-1)*8,7);
            if ( !s1.empty() )
            {
                ++I;
                stringstream(s1) >> BK[I];
            }
        }
    }
}

void Diffractogram::InitPolynomialBackground()
{
    int I,J, NBX, NBC, NXX, NINC, N2X;
    double DIFB,BSTEP;

    DIFB=POS[1]-THMIN;
    if (DIFB < 0.0)
    {
        NBX=static_cast<int>((POS[2]-THMIN)/STEP+1.5);
        BK[1]=BCK[1]-DIFB/(POS[2]-POS[1])*(BCK[2]-BCK[1]);
        BSTEP=STEP/(POS[2]-POS[1])*(BCK[2]-BCK[1]);
        for (I=2; I <= NBX; ++I) BK[I]=BK[I-1]+BSTEP;
        NXX=2;
    }
    else
    {
        NBX=static_cast<int>(DIFB/STEP+2.5);
        for (I=1; I <= NBX; ++I) BK[I]=BCK[1];
        NXX=1;
    }
    NBC=NBCKGD-1;
    if (POS[NBCKGD] < THMAX)
    {
        POS[NBCKGD+1]=THMAX;
        BCK[NBCKGD+1]=BCK[NBCKGD];
        NBC=NBC+1;
    }
    for (J=NXX; J <= NBC; ++J)
    {
        BSTEP=STEP*(BCK[J+1]-BCK[J])/(POS[J+1]-POS[J]);
        NINC=static_cast<int>((POS[J+1]-POS[J])/STEP+1.5);
        N2X=min(NPTS,NBX+NINC);
        for (I=NBX; I <= N2X; ++I) BK[I]=BK[I-1]+BSTEP;
        NBX=N2X;
        if(NBX == NPTS) break;
    }
}

void Diffractogram::readasc(void)
{
    int I;

    for (I=1; I <= 56; ++I)
    {
        if( TEST_[I] == ':')
        {
            LB1_=I;
            goto L31;
        }
    }
L31:
    for (I=LB1_+1; I <= 56; ++I)
    {
        if( TEST_[I] != ' ')
        {
            LB2_=I;
            goto L41;
        }
    }
L41:
    for (I=LB2_; I <= 56; ++I)
    {
        if( TEST_[I] == ' ')
        {
            LB3_=I-1;
            goto L51;
        }
    }
L51:;
}
