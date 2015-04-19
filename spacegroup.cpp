#include "spacegroup.h"

SpaceGroup::SpaceGroup()
{
    value = "";
    NSPGRP = 0;
}

void SpaceGroup::set(string p)
{
    value = p;
}

string SpaceGroup::get(void)
{
    return value;
}

SpaceGroup& SpaceGroup::operator=(string p)
{
    set(p);

    return *this;
}

SpaceGroup::operator string()
{
    return value;
}

ostream& operator<<(ostream& os, SpaceGroup& p)
{
    os << setw(20) << p.get();
    return os;
}

void SpaceGroup::SPGP(string SPG)
{
    //
    //                             THIS SR INTERPRETS THE SPACE GROUP SYMBOL
    //                           AND GENERATES OPERATORS WHICH SYNSPGRPOPR USES TO
    //                           GENERATE THE FULL SET OF EQUIVALENT
    //                           POSITIONS AND OP1 AND SMTRY2 USE TO GENERATE
    //                           THE FULL SET OF EQUIVALENT INDICES.

    const int NDAT[7+1] = { 0,3,1,2,0,2,1,3};
    const string CHR = "  CBAPFIRMND123456-/H.";
    int I,J,K,M,M1,M2,M3,N,LD,IJ,KJ,MJ,IX,IN1,IN2,NI4,NI5,NI6,NK4,NK5,NK6,NM4,
        NM5,NM6,NCI,NCK,NCM,LOD,NIJ,NXC,NYC,NZC,NJX,ICO4,ICO5,ICO6,KCO4,KCO5,KCO6,
        MCO4,MCO5,IKO4,IKO5,IKO6,MCO6,LOD1;
    int L[4+1][4+1];
    int ICO[6+1][15+1];


    for (I=0; I <= 6; ++I)
    {
        for (J=0; J <= 15; ++J)
        {
            ICO[I][J]=0;
        }
    }
    for (I=0; I <= 4; ++I)
    {
        for (J=0; J <= 4; ++J)
        {
            L[I][J]=0;
        }
    }
    M2 = 0;
    K = 1;
    M = 1;
    NSPGRP = 0;
    NAXIS = 0;
    N = 0;


    for (J=1; J <= 20; ++J)
    {
        for (I=1; I <= 21; ++I)
        {
            if ( SPG[J] == CHR[I] ) goto L101;
        }
        goto L110;
        L101:
        if ( K+M+I == 3 ) goto L110;
        if ( I == 1 ) goto L108;
        L[M][K] = I;
        N = 0;
        M = M+1;
        if ( M-4 <= 0 ) goto L110; else goto L108;
        L108:
        if ( N == 1 ) goto L110;
        N = 1;
        M = 1;
        K = K+1;
        if ( K > 4 ) goto L200;
        L110:;
    }

    L200:
    if ( K <= 2 ) goto L500;
    //DB 3
    //     PRINT 3,N,M,K,J,L
    //   3 FORMAT (*0N =*I3*  M =*I3*  K =*I3*  J =*I3/
    //    1  * L =*4I4/(4X4I4))
    if ( L[1][1] > 8 ) goto L500;
    J = 1;
    I = L[1][1]-1;
    switch (I) {
    case 1:
        goto L203;
        break;
    case 2:
        goto L202;
        break;
    case 3:
        goto L201;
        break;
    case 4:
        goto L119;
        break;
    case 5:
        goto L204;
        break;
    case 6:
        goto L205;
        break;
    case 7:
        goto L206;
        break;
    }
    GOTOER();
    L201:
    NCONT[J] = 576;
    goto L209;
    L202:
    NCONT[J] = 516;
    goto L209;
    L203:
    NCONT[J] = 68;
    goto L209;
    L204:
    NCONT[J] = 576;
    J = J+1;
    goto L202;
    L205:
    NCONT[J] = 580;
    goto L209;
    L207:
    K = K-1;
    NSPGRP = K+3;
    goto L210;
    L206:
    if ( L[1][3] == 8  || L[1][4] == 8 ) goto L207;
    NCONT[J] = 8192;
    if ( L[1][K-1] == 20 ) K=K-1;
    NSPGRP = K+5;
    J = J+1;
    goto L210;
    L209:
    J = J+1;
    L119:
    if ( K-4 < 0 )
    {
        goto L120;
    }
    else if ( K-4 == 0)
    {
        goto L130;
    }
    else
    {
        goto L140;
    }
    L120:
    if ( L[1][2] == 18 ) goto L121;
    if ( L[1][2] == 17 ) goto L122;
    if ( L[1][2] == 14 ) goto L123;
    if ( L[1][2] == 15 ) goto L124;
    if ( L[1][2] == 12 ) goto L125;
    NSPGRP = 2;
    NAXIS = 2;
    goto L210;
    L121:
    if ( L[2][2] == 17 ) goto L122;
    if ( L[2][2] == 14 ) goto L123;
    if ( L[2][2] == 15 ) goto L124;
    if ( L[2][2] == 12 ) goto L125;
    goto L500;
    L122:
    NSPGRP = 11;
    goto L210;
    L123:
    NSPGRP = 8;
    goto L210;
    L124:
    NSPGRP = 4;
    goto L210;
    L125:
    NSPGRP = 1;
    goto L210;
    L130:
    if ( L[1][3] != 14 ) goto L500;
    NSPGRP = 13;
    goto L210;
    L140:
    if ( L[1][3] == 14 ) goto L151;
    if ( L[1][2] == 18 ) goto L152;
    if ( L[1][2] == 17 ) goto L153;
    if ( L[1][2] == 14 ) goto L154;
    if ( L[1][2] == 15 ) goto L155;
    if ( L[1][2] == 12 ) goto L141;
    if ( L[1][3] == 12 ) goto L142;
    if ( L[1][4] == 12 ) goto L500;
    NSPGRP = 3;
    goto L210;
    L141:
    if ( L[1][3] == 12 ) goto L143;
    if ( L[1][4] != 12 ) goto L500;
    NSPGRP = 2;
    NAXIS = 2;
    goto L210;
    L142:
    if ( L[1][4] != 12 ) goto L500;
    NSPGRP = 2;
    NAXIS = 1;
    goto L210;
    L143:
    if ( L[1][4] == 12 ) goto L500;
    NSPGRP = 2;
    NAXIS = 3;
    goto L210;
    L151:
    NSPGRP = 14;
    goto L210;
    L152:
    if ( L[2][2] == 17 ) goto L153;
    if ( L[2][2] == 14 ) goto L154;
    if ( L[2][2] == 15 ) goto L155;
    goto L500;
    L153:
    NSPGRP = 12;
    goto L210;
    L154:
    if ( L[1][3] == 12 ) goto L156;
    if ( L[1][4] != 12 ) goto L500;
    NSPGRP = 9;
    goto L210;
    L155:
    NSPGRP = 5;
    goto L210;
    L156:
    if ( L[1][4] == 12 ) goto L123;
    NSPGRP = 10;
    L210:
    K = K-1;

    //DB 2
    //     PRINT 4,K,NSPGRP,NAXIS,J,NCONT
    //   4 FORMAT ('0K =',I3,'  NSPGRP =',2I4,'  J=',I3/' NCONT =',10I4)
    N = 1;
    for (M=2; M <= K; ++M)
    {
        if ( L[1][M] == 0 ) goto L500;
        L218:
        I = abs(L[1][M]-5);
        L219:
        if ( I <= 0 || I > 15 ) goto L500;
        switch (I) {
        case 1:
            goto L220;  // A
            break;
        case 2:
            goto L220;  // B
            break;
        case 3:
            goto L220;  // C
            break;
        case 4:
            goto L220;  // M
            break;
        case 5:
            goto L220;  // N
            break;
        case 6:
            goto L245;  // D
            break;
        case 7:
            goto L250;  // 1
            break;
        case 8:
            goto L255;  // 2
            break;
        case 9:
            goto L260;  // 3
            break;
        case 10:
            goto L265;  // 4
            break;
        case 11:
            goto L500;  // 5
            break;
        case 12:
            goto L270;  // 6
            break;
        case 13:
            goto L275;  // -
            break;
        case 14:
            goto L300;  // /
            break;
        case 15:
            goto L300;
            break;
        }
        GOTOER();
        //        H
        L220:
        switch (M) {
        case 1:
            goto L500;
            break;
        case 2:
            goto L221;
            break;
        case 3:
            goto L222;
            break;
        case 4:
            goto L223;
            break;
        }
        GOTOER();
        L221:
        if ( L[1][3] == 14 ) goto L2230;
        if ( L[1][2] == 15 ) goto L2230;
        if ( L[1][2] == 17 ) goto L2230;
        if ( K == 2 ) goto L2220;
        L2210:
        ICO[1][N] = 2;
        ICO[2][N] = 0;
        ICO[3][N] = 0;
        ICO[4][N] = 4;
        if ( I == 2 || I == 5 ) ICO[5][N] = 1;
        if ( I == 3 || I == 5 ) ICO[6][N] = 1;
        L2211:
        N = N+1;
        if ( L[1][2] != 15 ) goto L300;
        ICO[1][N] = ICO[2][N-1];
        ICO[2][N] = ICO[1][N-1];
        ICO[3][N] = ICO[3][N-1];
        ICO[4][N] = ICO[5][N-1];
        ICO[5][N] = ICO[4][N-1];
        ICO[6][N] = ICO[6][N-1];
        N= N+1;
        goto L300;
        L222:
        if ( L[1][2] == 14 || L[1][2] == 17 ) goto L225;
        if ( L[1][2] == 15 ) goto L2210;
        L2220:
        ICO[1][N] = 0;
        ICO[2][N] = 2;
        ICO[3][N] = 0;
        ICO[5][N] = 4;
        if ( I == 1 || I == 5 ) ICO[4][N] = 1;
        if ( I == 3 || I == 5 ) ICO[6][N] = 1;
        N = N+1;
        goto L300;
        L223:
        if ( L[1][3] == 14 || L[1][2] == 15 ) goto L224;
        if ( L[1][2] == 14 || L[1][2] == 17 ) goto L224;
        L2230:
        ICO[1][N] = 0;
        ICO[2][N] = 0;
        ICO[3][N] = 2;
        ICO[6][N] = 4;
        if ( I == 1 || I == 5 ) ICO[4][N] = 1;
        if ( I == 2 || I == 5 ) ICO[5][N] = 1;
        N = N+1;
        if ( M == 2 && L[1][2] == 15 && I == 5 ) ICO[4][N-3]=1;
        if ( M !=  2 || L[1][2] != 17) goto L300;
        if ( L[1][3] == 2 ) ICO[6][N-1]=1;
        if ( L[1][4] == 2 ) ICO[6][N-1]= (ICO[6][N-1]+1) % 2;
        goto L300;
        L224:
        ICO[1][N] = 1;
        ICO[2][N] = 1;
        ICO[3][N] = 0;
        ICO[4][N] = 4;
        ICO[5][N] = 4;
        if ( I == 3 || I == 5 ) ICO[6][N] = 1;
        if ( NSPGRP == 7 && I == 3 ) goto L2240;
        if ( I != 5 ) goto L226;
        L2240:
        ICO[4][N] = 5;
        ICO[5][N] = 5;
        N = N+1;
        goto L300;
        L225:
        if ( NSPGRP == 7 ) goto L224;
        ICO[1][N] = 3;
        ICO[2][N] = 3;
        ICO[3][N] = 0;
        ICO[4][N] = 4;
        ICO[5][N] = 4;
        if ( I == 3 ) ICO[6][N] = 1;
        L226:
        N = N+1;
        goto L300;
        L245:
        //    D TYPE MIRROR
        switch (M) {
        case 1:
            goto L500;
            break;
        case 2:
            goto L246;
            break;
        case 3:
            goto L247;
            break;
        case 4:
            goto L248;
            break;
        }
        GOTOER();
        L246:
        if ( K == 2 ) goto L247;
        ICO[1][N] = 2;
        ICO[2][N] = 0;
        ICO[3][N] = 0;
        ICO[4][N] = 6;
        ICO[5][N] = 6;
        ICO[6][N] = 6;
        N = N+1;
        goto L300;
        L247:
        ICO[1][N] = 0;
        ICO[2][N] = 2;
        ICO[3][N] = 0;
        ICO[4][N] = 6;
        ICO[5][N] = 6;
        ICO[6][N] = 6;
        N = N+1;
        goto L300;
        L248:
        if ( L[1][2] == 15 || L[1][3] == 14 ) goto L249;
        ICO[1][N] = 0;
        ICO[2][N] = 0;
        ICO[3][N] = 2;
        ICO[4][N] = 6;
        ICO[5][N] = 6;
        ICO[6][N] = 6;
        N = N+1;
        goto L300;
        L249:
        ICO[1][N] = 1;
        ICO[2][N] = 1;
        ICO[3][N] = 0;
        ICO[4][N] = 6;
        ICO[5][N] = 6;
        ICO[6][N] = 6;
        if (L[1][3] != 13) goto L226;
        ICO[4][N] = 0;
        ICO[5][N] = 1;
        N = N+1;
        goto L300;
        L250:
        //    1 FOLD ROTATION
        if ( L[2][M] != 3 ) goto L300;
        ICO[1][N] = 2;
        ICO[2][N] = 2;
        ICO[3][N] = 2;
        ICO[4][N] = 4;
        ICO[5][N] = 4;
        ICO[6][N] = 4;
        N = N+1;
        goto L300;
        L255:
        switch (M) {
        case 1:
            goto L500;
            break;
        case 2:
            goto L256;
            break;
        case 3:
            goto L257;
            break;
        case 4:
            goto L258;
            break;
        }
        GOTOER();
        //    2 FOLD ROTATION
        L256:
        if ( K == 2 ) goto L2571;
        L2561:
        ICO[1][N] = 0;
        ICO[2][N] = 2;
        ICO[3][N] = 2;
        ICO[4][N] = 0;
        ICO[5][N] = 4;
        ICO[6][N] = 4;
        if ( abs(L[2][M]-13) == 1 ) ICO[4][N] = 1;
        L2560:
        N = N+1;
        for (I=2; I <= 4; ++I)
        {
            if ( L[I][M] == 19 ) goto L2563;
        }
        goto L300;
        L2563:
        if ( L[I+1][M] <= 1 ) goto L500;
        I = abs(L[I+1][M]-5);
        goto L219;
        L257:
        if ( L[1][2] == 14 || L[1][2] == 17 ) goto L259;
        L2571:
        ICO[1][N] = 2;
        ICO[2][N] = 0;
        ICO[3][N] = 2;
        ICO[4][N] = 4;
        ICO[5][N] = 0;
        ICO[6][N] = 4;
        if ( L[2][M] == 12 ) ICO[5][N]=1;
        if ( L[1][2] == 15 ) goto L2211;
        goto L2560;
        L258:
        if ( L[1][2] >= 14   ) goto L2595;
        if ( L[1][3] == 14 ) goto L259;
        L2581:
        ICO[1][N] = 2;
        ICO[2][N] = 2;
        ICO[3][N] = 0;
        ICO[4][N] = 4;
        ICO[5][N] = 4;
        ICO[6][N] = 0;
        if ( abs(L[2][M]-13) == 1 ) ICO[6][N] = 1;
        if ( L[2][M] == 16 ) ICO[6][N] = 1;
        goto L2560;
        L259:
        if ( L[1][3] == 8 || L[1][4] == 8 ) goto L2595;
        ICO[1][N] = 1;
        ICO[2][N] = 1;
        ICO[3][N] = 2;
        ICO[4][N] = 0;
        ICO[5][N] = 0;
        ICO[6][N] = 4;
        goto L2560;
        L2595:
        if ( L[1][2] == 15)goto L259;
        ICO[1][N] = 3;
        ICO[2][N] = 3;
        ICO[3][N] = 2;
        ICO[4][N] = 4;
        ICO[5][N] = 4;
        ICO[6][N] = 4;
        goto L2560;
        //    3 FOLD ROTATION
        L260:
        switch (M) {
        case 1:
            goto L500;
            break;
        case 2:
            goto L261;
            break;
        case 3:
            goto L262;
            break;
        case 4:
            goto L500;
            break;
        }
        GOTOER();
        L261:
        if ( L[1][1] == 8 && J == 1 ) goto L262;
        IX = 0;
        if ( L[2][M] == 12 ) goto L2611;
        if ( L[2][M] == 13 ) goto L2612;
        L2610:
        NCONT[J] = 8196;
        goto L2613;
        L2611:
        NCONT[J] = 8200;
        if ( L[1][4] != 13 ) goto L2613;
        J = J+1;
        NCONT[J] = 2321;
        if ( L[1][2] == 14 ) NCONT[J]=4403;
        L[1][3] = 12;
        L[1][4] = 12;
        goto L2613;
        L2612:
        NCONT[J] = 8204;
        if ( L[1][4] != 13 ) goto L2613;
        J = J+1;
        NCONT[J] = 4369;
        if ( L[1][2] == 14 ) NCONT[J]=2355;
        L[1][3] = 12;
        L[1][4] = 12;
        L2613:
        J = J+1;
        if ( IX-1 < 0 )
        {
            goto L250;
        } else if ( IX-1 == 0)
        {
            goto L2581;
        }
        else
        {
            goto L2230;
        }
        //    CUBIC OR RHOMBOHEDERAL
        L262:
        NCONT[J] = 16384;
        if (L[2][M] != 3 ) goto L2620;
        J = J+1;
        NCONT[J] = 290;
        L2620:
        J = J+1;
        if ( N == 1 ) goto L300;
        I = N-1;
        IN1 = N;
        IN2 = N+1;
        ICO[1][IN1] = ICO[3][I];
        ICO[2][IN1] = ICO[1][I];
        ICO[3][IN1] = ICO[2][I];
        ICO[4][IN1] = ICO[6][I];
        ICO[5][IN1] = ICO[4][I];
        ICO[6][IN1] = ICO[5][I];
        ICO[1][IN2] = ICO[2][I];
        ICO[2][IN2] = ICO[3][I];
        ICO[3][IN2] = ICO[1][I];
        ICO[4][IN2] = ICO[5][I];
        ICO[5][IN2] = ICO[6][I];
        //L263:
        ICO[6][IN2] = ICO[4][I];
        N = IN2+1;
        goto L300;
        L265:
        if ( M != 2 ) goto L500;
        ICO[1][N] = 3;
        ICO[2][N] = 1;
        ICO[3][N] = 0;
        ICO[4][N] = 4;
        ICO[5][N] = 4;
        ICO[6][N] = 0;
        if ( L[2][2] == 12 ) ICO[6][N] = 2;
        if ( L[2][2] == 13 ) ICO[6][N] = 1;
        if ( L[2][2] == 14 ) ICO[6][N] = 3;
        if ( L[2][2] ==  3 ) ICO[3][N] = 2;
        N = N+1;
        if ( L[2][2] == 3 && L[1][3] == 14 && L[1][4] == 11 ) goto L266;
        if ( L[1][3] == 14 ) goto L2561;
        if ( K > 2 || L[1][1] != 7 ) goto L2581;
        if ( L[2][2]+L[3][2] != 12 ) goto L2581;
        L[2][2] = 0;
        ICO[5][N-1] = 1;
        goto L2581;
        L266:
        if ( L[1][1] != 7 || L[1][2] != 15 ) goto L500;
        NCONT[J] = 16384;
        NCONT[J+1] = 356;
        NCONT[J+2] = 834;
        NCONT[J+3] = 1177;
        J = J+3;
        goto L410;
        L270:
        if ( M != 2 ) goto L500;
        IX = 1;
        if ( L[2][2] == 12 ) goto L2611;
        if ( L[2][2] == 13 ) goto L2612;
        if ( L[2][2] == 14 ) goto L2610;
        if ( L[2][2] == 15 ) goto L2611;
        if ( L[2][2] == 16 ) goto L2612;
        if ( L[2][2] ==  3 ) goto L271;
        goto L2610;
        L271:
        IX = 2;
        goto L2610;
        L275:
        if ( M != 2 ) goto L500;
        L[1][2] = L[2][2];
        L[2][2] = 3;
        goto L218;
        L300:;
    }
    N = N-1;
    IJ = J;
    if ( N > 0 ) goto L301;
    J = J-1;
    goto L410;
    L301:
    for (I=1; I <= N; ++I)
    {
        NCONT[J] = ICO[1][I]+16*ICO[2][I]+128*ICO[3][I]+4*(ICO[4][I] % 4) +64*(ICO[5][I] % 4)+512*(ICO[6][I] % 4);
        J = J+1;
    }
    L306:
    J = J-1;
    //DB 4
    //     PRINT 5,J,IJ,N,NSPGRP,NAXIS,NCONT
    //    1 ,ICO
    //   5 FORMAT ('0J =',I3,'  IJ =',I3,'  N =',I3,'  NSP =',2I5/
    //    1 ' NCONT =',10I4/' ICO =',6I3/(6X6I3))
    //L310:
    for (I=IJ; I <= J; ++I)
    {
        NCI = I-IJ+1;
        ICO4 = ICO[4][NCI] % 4;
        ICO5 = ICO[5][NCI] % 4;
        ICO6 = ICO[6][NCI] % 4;
        NI4 = NDAT[ICO4+4];
        NI5 = NDAT[ICO5+4];
        NI6 = NDAT[ICO6+4];
        KJ = I+1;
        if ( (NCONT[I] & 307)-290 == 0) goto L315; else goto L320;
        L315:
        if ( (NCONT[I] & 7884) == 0 ) goto L318; else goto L316;
        //   A CENTER IS PRESENT, NOT AT 0,0,0
        L316:;
        //*outputfile << endl << " A 1-bar site is present, But NOT at 0,0,0" << endl;
        goto L320;
        //   A CENTER AT 0,0,0 IS PRESENT
        L318:;
        L320:;
        if ( KJ > J ) goto L400;
        for (K=KJ; K <= J; ++K)
        {
            NCK = K-IJ+1;
            KCO4 = ICO[4][NCK] % 4;
            KCO5 = ICO[5][NCK] % 4;
            KCO6 = ICO[6][NCK] % 4;
            NK4 = NDAT[KCO4+4];
            NK5 = NDAT[KCO5+4];
            NK6 = NDAT[KCO6+4];
            NXC = NI4+NK4 % 4;
            NYC = NI5+NK5 % 4;
            NZC = NI6+NK6 % 4;
            MJ = K+1;
            LOD = (NCONT[K] ^ NCONT[I]) & 307;
            if ( LOD != 290 ) goto L330;
            //   A CENTER IS GENERATED
            if ( NXC+NYC+NZC == 0 ) goto L330;
            if ( L[1][2] == 15 && I == IJ+1 ) goto L324;
            L321:
            if ( NXC == 0 ) goto L322;
            if ( ICO[4][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NXC]*4;
            if ( ICO[4][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NXC]*4;
            L322:
            if ( NYC == 0 ) goto L323;
            if ( ICO[5][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NYC]*64;
            if ( ICO[5][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NYC]*64;
            L323:
            if ( NZC == 0 ) goto L325;
            if ( ICO[6][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NZC]*512;
            if ( ICO[6][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NZC]*512;
            goto L325;
            L324:
            L[1][2] = 30;
            NIJ = ICO[6][1];
            NIJ = NDAT[NIJ+4];
            NJX = NIJ+NXC % 4;
            NCONT[IJ] = (NCONT[IJ] & 2035)+NDAT[NJX+4]*4;
            NCONT[IJ] = (NCONT[IJ] & 1855)+NDAT[NIJ+4]*64;
            goto L321;
            L325:
            ICO[4][NCI] = 0;
            ICO[4][NCK] = 0;
            ICO[5][NCK] = 0;
            ICO[5][NCI] = 0;
            ICO[6][NCI] = 0;
            ICO[6][NCK] = 0;
            L330:
            ICO4 = (NCONT[I]/4) % 4;
            IKO4 = (NCONT[K]/4) % 4;
            ICO5 = (NCONT[I]/64) % 4;
            IKO5 = (NCONT[K]/64) % 4;
            ICO6 = (NCONT[I]/512) % 4;
            IKO6 = (NCONT[K]/512) % 4;
            NI4 = NDAT[ICO4+4];
            NK4 = NDAT[IKO4+4];
            NI5 = NDAT[ICO5+4];
            NK5 = NDAT[IKO5+4];
            NI6 = NDAT[ICO6+4];
            NK6 = NDAT[IKO6+4];
            if ( MJ > J ) goto L390;
            for (M=MJ; M <= J; ++M)
            {
                NCM = M-IJ+1;
                MCO4 = ICO[4][NCM] % 4;
                MCO5 = ICO[5][NCM] % 4;
                MCO6 = ICO[6][NCM] % 4;
                NM4 = NDAT[MCO4+4];
                NM5 = NDAT[MCO5+4];
                NM6 = NDAT[MCO6+4];
                LOD1 = (LOD ^ NCONT[M]) & 307;
                NXC = NI4+NK4+NM4 % 4;
                NYC = NI5+NK5+NM5 % 4;
                NZC = NI6+NK6+NM6 % 4;
                if ( LOD1 != 290 ) goto L350;
                if  ( L[1][2] == 11 ) goto L334;
                if ( NXC+NYC+NZC == 0 ) goto L380;
                if ( IJ == 1 ) goto L331;
                if(NXC +NYC  == 0 )goto L333;
                if ( L[1][3] == 14 ) goto L331;
                if ( L[1][2] == 15 ) goto L331;
                if ( L[1][2] == 14 ) goto L331;
                if ( L[1][2] == 17 ) goto L331;
                if ( L[1][2]+L[1][3]+L[1][4] < 18 && L[1][2] != 9 ) goto L331;
                NXC = (NCONT[1]/2+NXC) % 4;
                NYC = (NCONT[1]/32+NYC) % 4;
                NZC = (NCONT[1]/256+NZC) % 4;
                L331:
                if ( NXC == 0 ) goto L332;
                if ( ICO[4][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NXC]*4;
                if ( ICO[4][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NXC]*4;
                if ( ICO[4][NCM] == 4 ) NCONT[M] = NCONT[M]+NDAT[NXC]*4;
                L332:
                if ( NYC == 0 ) goto L333;
                if ( ICO[5][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NYC]*64;
                if ( ICO[5][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NYC]*64;
                if ( ICO[5][NCM] == 4 ) NCONT[M] = NCONT[M]+NDAT[NYC]*64;
                L333:
                if ( NZC == 0 ) goto L335;
                if ( ICO[6][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NZC]*512;
                if ( ICO[6][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NZC]*512;
                if ( ICO[6][NCM] == 4 ) NCONT[M] = NCONT[M]+NDAT[NZC]*512;
                goto L335;
                L334:
                NCONT[I] = (NCONT[I] & 2035);
                NCONT[K] = (NCONT[K] & 1855);
                NCONT[M] = (NCONT[M] & 511);
                L335:
                ICO[4][NCI] = 0;
                ICO[5][NCI] = 0;
                ICO[6][NCI] = 0;
                ICO[4][NCK] = 0;
                ICO[5][NCK] = 0;
                ICO[6][NCK] = 0;
                ICO[4][NCM] = 0;
                ICO[5][NCM] = 0;
                ICO[6][NCM] = 0;
                goto L380;
                L350:
                if ( LOD1 != 0 ) goto L380;
                //DB 3
                //     PRINT 9,LOD1,NXC,NYC,NZC,NI4,NK4,NI5,NK5,NI6,NK6
                //   9 FORMAT (* LOD1 = *O10* NXC =*I4* NYC =*I4* NZC =*I4/
                //    1 * NI4 =*I3* NK4 =*I3* NI5 =*I3* NK5 =*I3* NI6 =*I3* NK6 =*I3)
                if ( NXC+NYC+NZC == 0 ) goto L359;
                if ( L[1][2] == 30 ) goto L3500;
                if ( L[1][2] == 15 ) goto L355;
                if ( L[1][2] == 17 ) goto L359;
                if ( L[1][2] == 13 ) goto L351;
                if ( L[1][3] == 13 ) goto L352;
                if ( L[1][4] == 13 ) goto L353;
                goto L500;
                L3500:
                if ( NXC == 0 ) goto L359;
                if ( L[4][2] == 4 ) goto L359;
                if ( NXC == NK4 ) goto L359;
                NCONT[K] = NCONT[K] ^ (NDAT[NXC]*4);
                if ( L[1][4] == 13 && L[2][4] == 0 ) goto L3580;
                goto L359;
                L351:
                NYC = (NK5+NM5) % 4;
                NZC = (NK6+NM6) % 4;
                if ( L[1][3] == 13 ) goto L354;
                if ( L[1][3] == 14 ) goto L354;
                if ( L[1][4] == 13 ) goto L354;
                M2 = I+1;
                //L3511:
                if(L[2][2] == 12 && (L[1][3] == 9  || L[1][4] == 9))  goto L360;
                if ( ICO[5][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NYC]*64;
                if ( ICO[5][NCM] == 4 ) NCONT[M] = NCONT[M]+NDAT[NYC]*64;
                //L3510:
                if ( NZC <= 0 ) goto L360;
                if ( ICO[6][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NZC]*512;
                if ( ICO[6][NCM] == 4 ) NCONT[M] = NCONT[M]+NDAT[NZC]*512;
                goto L360;
                L352:
                NXC = (NI4+NM4) % 4;
                NZC = (NI6+NM6) % 4;
                M2 = K+1;
                if(L[2][3] == 12  && (L[1][2] == 9  ||  L[1][4] == 9)) goto L360;
                if ( NXC <= 0 ) goto L3520;
                if ( ICO[4][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NXC]*4;
                if ( ICO[4][NCM] == 4 ) NCONT[M] = NCONT[M]+NDAT[NXC]*4;
                L3520:
                if ( NZC <= 0 ) goto L360;
                if ( ICO[6][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NZC]*512;
                if ( ICO[6][NCM] == 4 ) NCONT[M] = NCONT[M]+NDAT[NZC]*512;
                goto L360;
                L353:
                NXC = (NI4+NK4) % 4;
                NYC = (NI5+NK5) % 4;
                if(L[2][4] == 12  &&  ( L[1][2] == 9 || L[1][3] == 9)) goto L359;
                if ( NXC <= 0 ) goto L3530;
                if ( ICO[4][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NXC]*4;
                if ( ICO[4][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NXC]*4;
                L3530:
                if ( NYC <= 0 ) goto L359;
                if ( ICO[5][NCI] == 4 ) NCONT[I] = NCONT[I]+NDAT[NYC]*64;
                if ( ICO[5][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NYC]*64;
                goto L359;
                L354:
                if ( NYC <= 0 ) goto L3540;
                if ( abs(L[2][2]-13) == 1 ) NCONT[I]=NCONT[I]+NDAT[NYC]*64;
                L3540:
                if ( NZC <= 0 ) goto L3541;
                if ( ICO[6][NCK] == 4 ) NCONT[K] = NCONT[K]+NDAT[NZC]*512;
                L3541:
                if ( NXC <= 0 ) goto L359;
                if ( L[1][4] == 0 || L[1][3] == 14 ) goto L359;
                if ( L[2][4] != 12 ) NCONT[K]=NCONT[K]+NDAT[NXC]*4;
                goto L359;
                L355:
                L[1][2] = 30;
                if ( L[1][3] == 14 ) goto L370;
                if ( L[1][3] == 13 && L[2][3] == 12 ) goto L3561;
                if ( L[1][1] ==  7 && L[2][2] == 12 ) goto L357;
                if ( L[1][1] ==  7 && L[2][2] == 14 ) goto L357;
                if ( L[1][1] ==  6 && L[2][2] == 12 ) goto L358;
                if ( L[1][1] ==  6 && L[2][2] == 14 ) goto L358;
                if ( L[1][1] ==  7 && L[1][3] == 13 ) goto L3581;
                if ( (NCONT[I] % 1024) == 19 ) goto L350;
                if ( NXC+NYC > 0 ) goto L356;
                if ( NZC != NM6 ) goto L359;
                if ( NZC <= 0 ) goto L359;
                NCONT[K] = NCONT[K]+NDAT[NZC]*512;
                goto L359;
                L356:
                if ( L[3][2] == 10 || L[4][2] == 10 ) goto L359;
                if ( L[1][3] != 10 ) goto L359;
                L3561:
                if ( L[2][2] == 3 ) goto L3560;
                NCONT[I] = NCONT[I] ^ 68;
                goto L3580;
                L3560:
                if ( L[1][4] == 2 ) NCONT[K]= NCONT[K] ^ 512;
                goto L359;
                L357:
                NCONT[I] = NCONT[I]+64;
                if ( L[3][2] != 19 ) NCONT[I+1] = NCONT[I+1]-512;
                if ( L[1][4] == 11 && L[3][2] == 19 )NCONT[I] = NCONT[I]+588;
                goto L3580;
                L358:
                NCONT[I] = NCONT[I]+136;
                L3580:
                if ( L[1][3] != 13 ) goto L359;
                L3581:
                M2 = K+1;
                goto L360;
                L359:
                M2 = M+1;
                L360:
                if ( M2 > J ) goto L306;
                //DB 2
                //     PRINT 6,M2,I,J,K,M
                //   6 FORMAT (*0M2 =*I3*  I =*I3*  J =*I3* K =*I3*  M =*I3)
                for (M1=M2; M1 <= J; ++M1)
                {
                    NCONT[M1-1] = NCONT[M1];
                    for (M3=1; M3 <= 6; ++M3) ICO[M3][M1-1] = ICO[M3][M1];
                }
                goto L306;
                L370:
                if ( L[1][4] > 12 ) goto L371;
                M2 = I+1;
                goto L360;
                L371:
                L[1][2] = 13;
                if ( L[2][2] == 12 && L[1][1] == 5 ) NCONT[I] = NCONT[I]+200;
                if ( L[2][2] == 12 && L[1][1] == 6 ) NCONT[I] = NCONT[I]+136;
                if ( L[2][2] == 12 && L[1][1] == 7 ) NCONT[I] = NCONT[I]+652;
                if ( L[2][2] == 13 )                 NCONT[I] = NCONT[I]+68;
                if ( L[2][2] == 14 )                 NCONT[I] = NCONT[I]+140;
                goto L359;
                L380:;
            }
            L390:;
        }
        L400:;
    }
    LD = 0;
    for (I=1; I <= J; ++I) LD = (LD ^ NCONT[I]) & 307;
    if ( LD != 0 ) goto L410;
    NCONT[J-1] = NCONT[J];
    J = J-1;
    L410:
    NCONT[J+1] = 0;
    //DB 1
    //     PRINT 4,K,NSPGRP,NAXIS,J,NCONT
    L450:
    return;
    L500:
    //*outputfile << " Did YOU enter the space group correctly?" << endl
    //    << "  ERRORS of some sort were detected." << endl
    //    << "  The card was read as:" << SPG << endl;
    NAXIS = 4;
    goto L450;
}

void SpaceGroup::GOTOER(void)
{
    //*outputfile << "subroutine gotoer was called" << endl;
    throw;
    //throw DBWSException("subroutine gotoer was called");
}
