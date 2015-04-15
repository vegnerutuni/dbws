      PROGRAM DBWS9807         
C     RELEASE 27.02.99
C     PROGRAM PFSR(TAPE3=O"1001",TAPE4=O"1001",TAPE5=O"1001",TAPE6,OUTPU
C    * T=O"301",TAPE7=OUTPUT,TAPE8=O"10001")
C   PROGRAM FOR RIETVELD ANALYSIS OF X-RAY AND NEUTRON (NUCLEAR SCATTERING
C       ONLY) POWDER DIFFRACTION PATTERNS
C   SCHOOL OF PHYSICS (R. A. YOUNG)
C   GEORGIA INSTITUTE OF TECHNOLOGY, ATLANTA, GA. 30332  USA
C   PRECURSOR PROGRAM, DBW3.2, WRITTEN BY D. B. WILES
C   MUCH OF SUBROUTINE SPGP AND RELATED SYMMETRY AND MULTIPLICITY 
C     CALCULATIONS TAKEN FROM `SPGP' BY ALLEN C. LARSON while at LOS ALAMOS 
C     NATIONAL LABORATORY
C   DBW3.2S CREATED BY A. SAKTHIVEL BY MODIFICATION OF DBW3.2
C   DURING JAN 86 - AUG 88. 
C   DBWS-9006PC CREATED BY A. SAKTHIVEL APR-JUN 90. MARCH-DOLLASE PREFERRED 
C       ORIENTATION ADDED, ANSI FILE OPEN CHANGES MADE, CALPLT DROPPED 
C   INPTR MOD 12.6.90 TO USE SRS DATA 
C   NEW TYPO RE IOF IN CALCUL AFTER LINE 4340 CORRECTED  17.8.90
C   BREAK IN LINE 3387 REMOVED OCT. 90
C   FUNCTION MULT CORRECTED FOR MULTIPLICITY PER RJ HILL 28.11.90
C   SIGN ERRORS IN DELTA F' TABLE CORRECTED 8.4.91
C   OUTPUT OF BRAGG MARKER LOCATIONS TO PLOTINFO MODIFIED TO INCLUDE ZERO
C      ERROR AND SAMPLE DISPLACEMENT CORRECTION 4.7.91
C   OUTPUT INFORMATION ON NO. OF PHASES AND NO.OF REFLECTION IN EACH PHASE
C      TO PLOTINFO
C   OUTPUT PHASE NAME TO PLOTINFO ON 22.7.91
C   THE FOLLOWING CHANGES MADE BY T.S. MOSS, JULY 1994 - APRIL 1995:
C      CHANGES MADE TO ALLOW FOR SURFACE ROUGHNESS MODELING USING EQUATIONS
C         FROM PITSCHKE, ET AL., AND SPARKS, ET AL.
C      CORRECTION MADE TO THE MULT FUNCTION:  REDUCED CALC'D MULTIPLICITY BY
C         ONE-HALF IN THE -3, -3M1, 6/M, AND 6/MMM LAUE GROUPS
C      INCORPORATED QUANTITATIVE ANALYSIS ROUTINE PROVIDED BY CARLOS O. 
C        PAIVA-SANTOS (2 STEPS) WHICH NOW (30.3.95 RELEASE) FULLY UTILIZES
C        REFINED VALUES OF LATTICE AND OCCUPANCY PARAMETERS PLUS ADDED 
C        TABLE OF ATOMIC WEIGHTS
C      ADDITIONAL SURFACE ROUGHNESS MODEL BY THE MODIFIED EQUATION
C         OF SUORRTI AND BY THE EQUATION OF YOUNG
C      SPLIT PEARSON VII PROFILE ADDED BASED ON THE CODING OF TORAYA
C      ADDITIONAL ANNOTATION TO THE INPUT FILE ADDED
C      CORRECTION OF R-EXP ERROR DURING THE FINAL CYCLE
C      FURTHER CORRECTIONS TO THE COMBINATION SURFACE ROUGHNESS MODEL 
C      CORRECTION OF DURBIN-WATSON STATISTIC TO USE UNWEIGHTED DATA
C      CORRECTION OF ESD'S OF MOLAR PERCENTAGE IN QUANTITATIVE ROUTINE
C   DR. H. MARCINIAK, AUTHOR OF DMPLOT, DURING 1993-94 CONTRIBUTED TO THIS 
C       SOURCE CODE:   
C     (1) THE DYNAMIC OUTPUT TO SCREEN, (2) THE CHANGE TO USE THE `PARAM.INC' 
C     FILE  TO FACILITATE RESETTING SOME ARRAY SIZES, AND (3) THE CODE CHANGES

C     NEEDED TO MODIFY THE `PLOTINFO' AND `PLOTINFO.BIN' FILES TO TAKE 
C     ADVANTAGE OF SOME ADVANCED CAPABILITIES OF DMPLOT. 
C   Line marker 475 in subroutine INPTR moved 20.1.95 to the command
C            line  'open(3,file=' ',status='unknown')'                  
C   THE FOLLOWING CHANGES MADE BY C.O.PAIVA-SANTOS, 
C    MARCH 1995:
C    INCLUSION OF ATOMIC WEIGHTS TABLE, CHANGES IN THE icf FORMAT 
C    (line 8 and line 11)and changes to allow the refined values of cell
C    and occupancy parameters the computation of QPA. The unit cell density
C    of each phase is also given.
C    SEPTEMBER AND OCTOBER 1996: CHANGES in the occupancy code, to make use of
C     the site multiplicity (M) instead of the occupancy (N: site multiplicity
C     divided by general multiplicity). A refined site occupation (so)
C     parameter for each atom is now substituting N.
C   COPIES OF THE DISTRIBUTION PACKAGE, WHICH CONTAINS THIS SOURCE CODE, A 
C    USER'S GUIDE, TEST CASES AND COPIES OF THE SHAREWARE PLOT PROGRAMS ARE 
C    AVAILABLE FROM PROF. R. A. YOUNG AT THE ABOVE ADDRESS.
      include 'param.inc'
C      PARAMETER (idsz=2048,irs=512,NATS=64,msz=32 ,nov=128)
C     CHANGE THESE VARIABLES TO SUIT THE SIZE OF THE COMPUTER AND THE PROBLEM
C     IDSZ = MAXIMUM NUMBER OF DATA POINTS IN TAPE 4
C     IRS  = MAXIMUM NUBMER OF REFLECTIONS IN THE PATTERN
C     NATS = MAXIMUM NUMBER OF ATOMS IN THE PROBLEM
C     MSZ  = MAXIMUM NUMBER OF REFINABLE PARAMETERS (MATRIX SIZE)
C     NOV  = MAXIMUM NUMBER OF BRAGG REFLECTIONS CONTRIBUTING TO ONE STEP
      REAL LAMDA
      INTEGER PTR
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
C     COMMON changed  !cp ap 97
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
C     *IPLCOM,IPLDIS,IPLAM,IPLAIR,IPBIG
C     * ,IPREF,IABSR,FONDO,IAX,IAS,CTIME,FI,iasym,sw,iaxx,ibgd
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
C      COMMON/ALLP/FINAL(8*MSZ,2),ILOC
      COMMON/ALLP/FINAL(nfinal,2),ILOC
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),xmas(16)
C      COMMON/AIRPAR/AI1,AI2,AI3,AI4,AI5,AI6,AI7,AI8,AI9,AI10
      CHARACTER*4 NAM(16)
      COMMON/BKGSCALE/SCABKG(99)
      COMMON/COEFC/NAM
C      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ),BK(IDSZ),NPTS
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON /RTMTX/ IVEC(99,192),MLTPHS(99),ICNTPHS(99)
      common/sizestrain/sizeG(15),strainG(15),sizeL(15),strainL(15)
     *,siz(15),strain(15),NsizeStrain
      common/convert/icnvt
c      DATA UPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c      DATA LOWER/'abcdefghijklmnopqrstuvwxyz'/
c
c
c
C unit 4=data file
      open (4,file=' ',status='unknown')
C unit 5 Input control File      
      open (5,file=' ',status='OLD')
C unit 6 = output file      
      open (6,file=' ',status='unknown')
C unit 7 = r indexes      
      open (7,file='unit7',status='unknown')
C unit 8 = scratch      
      open (8,file='unit8',form='unformatted',status='unknown')  
      CALL INPTR  
C      if(icnvt.eq.1) goto 800
C Canton et all code starts here !cp may 03 97 
C-----OPEN,IF NECESSARY, FILE CONTAINING AMORPHOUS SCATTERING
      IF(GLB(20).NE.0.0.OR.AGLB(20).NE.0.0)call inpam
C Canton et all code stops here
      CALL ITER
      IF (MAXS.GT.0) THEN
      MCYCLX = MCYCLE
      MAXSX  = MAXS
      MCYCLE = 1
      MAXS   = 0
C-------LAST CALL TO ITER IF MCYCLE = 1 & MAXS = 0 
      CALL ITER
      MCYCLE = MCYCLX
      MAXS   = MAXSX
      END IF
      CALL EXPUT
      close (4)
      close (5)
      close (6)
      close (7)
      close (8,status='DELETE') 
      close (11)     
      close (20)
800   STOP   
      END 
c
c
c
c
      SUBROUTINE COEF(J,K)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),XMAS(16)
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/CNV/N,SINLAM(30),F(30)
      DIMENSION AB(10)
      EXTERNAL CV1,CV2
c
      DO 10 I=1,K
      SINLAM(I)=POSI(I)*POSI(I)
10      F(I)=SCAT(I)
      N=K 
      NN=MIN0(N,9)
      DO 15 I=1,NN
15      AB(I)=I
      NN1=NN+1
      DO 16 I=NN1,10
16      AB(I)=0
      CALL STEEP(AB,3,10,CV2,CV1)
      IF(AB(3).EQ.0.)AB(3)=-1E-6
      DO 20 I=1,9
20      AC(I,J)=AB(I)
      RETURN
      END 
c
c
c
c
      SUBROUTINE CV1(AB,V)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON/CNV/N,SINLAM(30),F(30)
      DIMENSION AB(10),V(10)
      MAXX=MIN0(9,N)
      DO 30 K=1,MAXX,2
      V(K)=0
      DO 20 I=1,N 
        V(K)=V(K)+F(I)
        DO 20 J=1,MAXX,2
20          V(K)=V(K)-AB(J)*EXP(-AB(J+1)*SINLAM(I))
      V(K)=-2.*V(K)*EXP(-AB(K+1)*SINLAM(K))
30      V(K+1)=-V(K)*SINLAM(K)*AB(K)
      NN=MAXX+1
      DO 10 I=NN,10 
10      V(I)=0
      RETURN
      END 
c
c
c
c
      SUBROUTINE CV2(AB,G)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON/CNV/N,SINLAM(30),F(30)
      DIMENSION AB(10)
c
	print *,'teste'
      MAXX=MIN0(9,N)
      G=0 
      DO 30 I=1,N
      T=F(I)
      DO 20 J=1,MAXX,2
20        T=T-AB(J)*EXP(-AB(J+1)*SINLAM(I))
30      G=G+T*T
      RETURN
      END 
c
c
c
c 
      SUBROUTINE SPGP(SPG)
C
C                             THIS SR INTERPRETS THE SPACE GROUP SYMBOL
C                           AND GENERATES OPERATORS WHICH SYMOPR USES TO
C                           GENERATE THE FULL SET OF EQUIVALENT
C                           POSITIONS AND OP1 AND SMTRY2 USE TO GENERATE
C                           THE FULL SET OF EQUIVALENT INDICES.
      CHARACTER*1 SPG(20),CHR(25)
      DIMENSION L(4,4),ICO(6,15)
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      DIMENSION NDAT(7)
      common/convert/icnvt
C     INTEGER AND,OR
C     AND(I,J) = IAND(I,J)
C     OR(I,J) = IOR(I,J)
      LODGDF(I,J) = iand( ieor(I,J),307)
1     FORMAT (/'0THE SPACE GROUP IS ',80A1)
2     FORMAT (' Did YOU enter the space group correctly?',/,
     *   '   ERRORS of some sort were detected.',/, 
     *   /'  The card was read as:',80A1)
24    FORMAT (/,' A 1-bar site is present, But NOT at 0,0,0')
      DATA NDAT(1),NDAT(2),NDAT(3),NDAT(4),NDAT(5),NDAT(6),NDAT(7)
     *       /  3,      1,      2,      0,      2,      1,      3/
      DATA CHR(1),CHR(2),CHR(3),CHR(4),CHR(5),CHR(6),CHR(7),CHR(8)
     *       /' ',   'C',   'B',   'A',   'P',   'F',   'I',   'R'/
      DATA CHR(9),CHR(10),CHR(11),CHR(12),CHR(13),CHR(14),CHR(15)
     *       /'M',    'N',    'D',    '1',    '2',    '3',    '4'/
      DATA CHR(16),CHR(17),CHR(18),CHR(19),CHR(20),CHR(21)
     *       / '5',    '6',    '-',    '/',    'H',    '.'/
CDB 2
C     PRINT 8,CHR
C   8 FORMAT (1X25A1)
      if(icnvt.ne.1) WRITE(6,1)SPG
      DO 97 I=1,6
      DO 97 J=1,15
97        ICO(I,J)=0
      DO 98 I=1,4
      DO 98 J=1,4 
98        L(I,J) = 0
      M2 = 0
      K = 1
      M = 1
      NSPGRP = 0
      NAXIS = 0
      N = 0
      DO 110 J=1,20 
99      DO 100 I=1,21
        IF ( SPG(J).EQ.CHR(I) ) GO TO 101
100       CONTINUE
CDB 2
C     PRINT 7,SPG(J),J
C   7 FORMAT (1XA1* IS AN UNRECOGNIZABLE CHARACTER IN COL. *I2)
      GO TO 110
101     IF ( K+M+I.EQ.3 ) GO TO 110
      IF ( I.EQ.1 ) GO TO 108
      L(M,K) = I
      N = 0
      M = M+1
      IF ( M-4 ) 110,110,108
108     IF ( N.EQ.1 ) GO TO 110
      N = 1
      M = 1
      K = K+1
      IF ( K.GT.4 ) GO TO 200
110     CONTINUE
200   IF ( K.LE.2 ) GO TO 500 
CDB 3
C     PRINT 3,N,M,K,J,L
C   3 FORMAT (*0N =*I3*  M =*I3*  K =*I3*  J =*I3/
C    1  * L =*4I4/(4X4I4))
      IF ( L(1,1).GT.8 ) GO TO 500
      J = 1
      I = L(1,1)-1
      GO TO (203,202,201,119,204,205,206),I
      CALL GOTOER
201   NCONT(J) = 576
      GO TO 209
202   NCONT(J) = 516
      GO TO 209
203   NCONT(J) = 68 
      GO TO 209
204   NCONT(J) = 576
      J = J+1
      GO TO 202
205   NCONT(J) = 580
      GO TO 209
207   K = K-1
      NSPGRP = K+3
      GO TO 210
206   IF ( L(1,3).EQ.8.OR.L(1,4).EQ.8 ) GO TO 207 
      NCONT(J) = 8192
      IF ( L(1,K-1).EQ.20 ) K=K-1
      NSPGRP = K+5
      J = J+1
      GO TO 210
209   J = J+1
119   IF ( K-4 ) 120,130,140
120   IF ( L(1,2).EQ.18 ) GO TO 121
      IF ( L(1,2).EQ.17 ) GO TO 122
      IF ( L(1,2).EQ.14 ) GO TO 123
      IF ( L(1,2).EQ.15 ) GO TO 124
      IF ( L(1,2).EQ.12 ) GO TO 125
      NSPGRP = 2
      NAXIS = 2
      GO TO 210
121   IF ( L(2,2).EQ.17 ) GO TO 122
      IF ( L(2,2).EQ.14 ) GO TO 123
      IF ( L(2,2).EQ.15 ) GO TO 124
      IF ( L(2,2).EQ.12 ) GO TO 125
      GO TO 500
122   NSPGRP = 11
      GO TO 210
123   NSPGRP = 8
      GO TO 210
124   NSPGRP = 4
      GO TO 210
125   NSPGRP = 1
      GO TO 210
130   IF ( L(1,3).NE.14 ) GO TO 500
      NSPGRP = 13
      GO TO 210
140   IF ( L(1,3).EQ.14 ) GO TO 151
      IF ( L(1,2).EQ.18 ) GO TO 152
      IF ( L(1,2).EQ.17 ) GO TO 153
      IF ( L(1,2).EQ.14 ) GO TO 154
      IF ( L(1,2).EQ.15 ) GO TO 155
      IF ( L(1,2).EQ.12 ) GO TO 141
      IF ( L(1,3).EQ.12 ) GO TO 142
      IF ( L(1,4).EQ.12 ) GO TO 500
      NSPGRP = 3
      GO TO 210
141   IF ( L(1,3).EQ.12 ) GO TO 143
      IF ( L(1,4).NE.12 ) GO TO 500
      NSPGRP = 2
      NAXIS = 2
      GO TO 210
142   IF ( L(1,4).NE.12 ) GO TO 500
      NSPGRP = 2
      NAXIS = 1
      GO TO 210
143   IF ( L(1,4).EQ.12 ) GO TO 500
      NSPGRP = 2
      NAXIS = 3
      GO TO 210
151   NSPGRP = 14
      GO TO 210
152   IF ( L(2,2).EQ.17 ) GO TO 153
      IF ( L(2,2).EQ.14 ) GO TO 154
      IF ( L(2,2).EQ.15 ) GO TO 155
      GO TO 500
153   NSPGRP = 12
      GO TO 210
154   IF ( L(1,3).EQ.12 ) GO TO 156
      IF ( L(1,4).NE.12 ) GO TO 500
      NSPGRP = 9
      GO TO 210
155   NSPGRP = 5
      GO TO 210
156   IF ( L(1,4).EQ.12 ) GO TO 123
      NSPGRP = 10
210   K = K-1
CDB 2
C     PRINT 4,K,NSPGRP,NAXIS,J,NCONT
C   4 FORMAT ('0K =',I3,'  NSPGRP =',2I4,'  J=',I3/' NCONT =',10I4)
      N = 1
      DO 300 M=2,K
      IF ( L(1,M).EQ.0 ) GO TO 500
218     I = IABS(L(1,M)-5)
219     IF ( I.LE.0.OR.I.GT.15 ) GO TO 500
C              A   B   C   M   N   D   1   2   3   4   5   6   -   /
      GO TO (220,220,220,220,220,245,250,255,260,265,500,270,275,300,
     *      300),I
      CALL GOTOER 
C        H
220     GO TO (500,221,222,223),M
      CALL GOTOER 
221   IF ( L(1,3).EQ.14 ) GO TO 2230
      IF ( L(1,2).EQ.15 ) GO TO 2230
      IF ( L(1,2).EQ.17 ) GO TO 2230
      IF ( K.EQ.2 ) GO TO 2220
2210  ICO(1,N) = 2
      ICO(2,N) = 0
      ICO(3,N) = 0
      ICO(4,N) = 4
      IF ( I.EQ.2.OR.I.EQ.5 ) ICO(5,N) = 1
      IF ( I.EQ.3.OR.I.EQ.5 ) ICO(6,N) = 1
2211    N = N+1
      IF ( L(1,2).NE.15 ) GO TO 300
      ICO(1,N) = ICO(2,N-1) 
      ICO(2,N) = ICO(1,N-1) 
      ICO(3,N) = ICO(3,N-1) 
      ICO(4,N) = ICO(5,N-1) 
      ICO(5,N) = ICO(4,N-1) 
      ICO(6,N) = ICO(6,N-1) 
      N= N+1
      GO TO 300
222     IF ( L(1,2).EQ.14.OR.L(1,2).EQ.17 ) GO TO 225
      IF ( L(1,2).EQ.15 ) GO TO 2210
2220  ICO(1,N) = 0
      ICO(2,N) = 2
      ICO(3,N) = 0
      ICO(5,N) = 4
      IF ( I.EQ.1.OR.I.EQ.5 ) ICO(4,N) = 1
      IF ( I.EQ.3.OR.I.EQ.5 ) ICO(6,N) = 1
      N = N+1
      GO TO 300
223     IF ( L(1,3).EQ.14.OR.L(1,2).EQ.15 ) GO TO 224
        IF ( L(1,2).EQ.14.OR.L(1,2).EQ.17 ) GO TO 224
2230    ICO(1,N) = 0
      ICO(2,N) = 0
      ICO(3,N) = 2
      ICO(6,N) = 4
      IF ( I.EQ.1.OR.I.EQ.5 ) ICO(4,N) = 1
      IF ( I.EQ.2.OR.I.EQ.5 ) ICO(5,N) = 1
      N = N+1
      IF ( M.EQ.2.AND.L(1,2).EQ.15.AND.I.EQ.5 ) ICO(4,N-3)=1
      IF ( M.NE. 2.OR.L(1,2).NE.17) GO TO 300
      IF ( L(1,3).EQ.2 ) ICO(6,N-1)=1 
      IF ( L(1,4).EQ.2 ) ICO(6,N-1)=MOD(ICO(6,N-1)+1,2)
      GO TO 300
224   ICO(1,N) = 1
      ICO(2,N) = 1
      ICO(3,N) = 0
      ICO(4,N) = 4
      ICO(5,N) = 4
      IF ( I.EQ.3.OR.I.EQ.5 ) ICO(6,N) = 1
      IF ( NSPGRP.EQ.7.AND.I.EQ.3 ) GO TO 2240
      IF ( I.NE.5 ) GO TO 226
2240    ICO(4,N) = 5
      ICO(5,N) = 5
      N = N+1
      GO TO 300
225     IF ( NSPGRP.EQ.7 ) GO TO 224
      ICO(1,N) = 3
      ICO(2,N) = 3
      ICO(3,N) = 0
      ICO(4,N) = 4
      ICO(5,N) = 4
      IF ( I.EQ.3 ) ICO(6,N) = 1
226     N = N+1
      GO TO 300
C    D TYPE MIRROR
245     GO TO (500,246,247,248),M
      CALL GOTOER 
246     IF ( K.EQ.2 ) GO TO 247
      ICO(1,N) = 2
      ICO(2,N) = 0
      ICO(3,N) = 0
      ICO(4,N) = 6
      ICO(5,N) = 6
      ICO(6,N) = 6
      N = N+1
      GO TO 300
247   ICO(1,N) = 0
      ICO(2,N) = 2
      ICO(3,N) = 0
      ICO(4,N) = 6
      ICO(5,N) = 6
      ICO(6,N) = 6
      N = N+1
      GO TO 300
248   IF ( L(1,2).EQ.15.OR.L(1,3).EQ.14 ) GO TO 249
      ICO(1,N) = 0
      ICO(2,N) = 0
      ICO(3,N) = 2
      ICO(4,N) = 6
      ICO(5,N) = 6
      ICO(6,N) = 6
      N = N+1
      GO TO 300
249     ICO(1,N) = 1
      ICO(2,N) = 1
      ICO(3,N) = 0
      ICO(4,N) = 6
      ICO(5,N) = 6
      ICO(6,N) = 6
      IF (L(1,3).NE.13) GO TO 226
      ICO(4,N) = 0
      ICO(5,N) = 1
      N = N+1
      GO TO 300
C    1 FOLD ROTATION
250     IF ( L(2,M).NE.3 ) GO TO 300
      ICO(1,N) = 2
      ICO(2,N) = 2
      ICO(3,N) = 2
      ICO(4,N) = 4
      ICO(5,N) = 4
      ICO(6,N) = 4
      N = N+1
      GO TO 300
255     GO TO (500,256,257,258),M
      CALL GOTOER 
C    2 FOLD ROTATION
256     IF ( K.EQ.2 ) GO TO 2571
2561    ICO(1,N) = 0
      ICO(2,N) = 2
      ICO(3,N) = 2
      ICO(4,N) = 0
      ICO(5,N) = 4
      ICO(6,N) = 4
      IF ( IABS(L(2,M)-13).EQ.1 ) ICO(4,N) = 1
2560    N = N+1
      DO 2562 I=2,4
        IF ( L(I,M).EQ.19 ) GO TO 2563
2562      CONTINUE
      GO TO 300
2563    IF ( L(I+1,M).LE.1 ) GO TO 500
      I = IABS(L(I+1,M)-5)
      GO TO 219
257     IF ( L(1,2).EQ.14.OR.L(1,2).EQ.17 ) GO TO 259
2571    ICO(1,N) = 2
      ICO(2,N) = 0
      ICO(3,N) = 2
      ICO(4,N) = 4
      ICO(5,N) = 0
      ICO(6,N) = 4
      IF ( L(2,M).EQ.12 ) ICO(5,N)=1
      IF ( L(1,2).EQ.15 ) GO TO 2211
      GO TO 2560
258   IF ( L(1,2).GE.14   ) GO TO 2595
      IF ( L(1,3).EQ.14 ) GO TO 259
2581  ICO(1,N) = 2
      ICO(2,N) = 2
      ICO(3,N) = 0
      ICO(4,N) = 4
      ICO(5,N) = 4
      ICO(6,N) = 0
      IF ( IABS(L(2,M)-13).EQ.1 ) ICO(6,N) = 1
      IF (      L(2,M).EQ.16    ) ICO(6,N) = 1
      GO TO 2560
259     IF ( L(1,3).EQ.8.OR.L(1,4).EQ.8 ) GO TO 2595
      ICO(1,N) = 1
      ICO(2,N) = 1
      ICO(3,N) = 2
      ICO(4,N) = 0
      ICO(5,N) = 0
      ICO(6,N) = 4
      GO TO 2560
2595    IF ( L(1,2).EQ.15)GO TO 259
      ICO(1,N) = 3
      ICO(2,N) = 3
      ICO(3,N) = 2
      ICO(4,N) = 4
      ICO(5,N) = 4
      ICO(6,N) = 4
      GO TO 2560
C    3 FOLD ROTATION
260     GO TO (500,261,262,500),M
      CALL GOTOER 
261     IF ( L(1,1).EQ.8.AND.J.EQ.1 ) GO TO 262
      IX = 0
      IF ( L(2,M).EQ.12 ) GO TO 2611
      IF ( L(2,M).EQ.13 ) GO TO 2612
2610    NCONT(J) = 8196
      GO TO 2613
2611    NCONT(J) = 8200
      IF ( L(1,4).NE.13 ) GO TO 2613
      J = J+1
      NCONT(J) = 2321
      IF ( L(1,2).EQ.14 ) NCONT(J)=4403
      L(1,3) = 12 
      L(1,4) = 12 
      GO TO 2613
2612    NCONT(J) = 8204
      IF ( L(1,4).NE.13 ) GO TO 2613
      J = J+1
      NCONT(J) = 4369
      IF ( L(1,2).EQ.14 ) NCONT(J)=2355
      L(1,3) = 12 
      L(1,4) = 12 
2613    J = J+1
      IF ( IX-1 ) 250,2581,2230
C    CUBIC OR RHOMBOHEDERAL
262     NCONT(J) = 16384
      IF (L(2,M).NE.3 ) GO TO 2620
      J = J+1
      NCONT(J) = 290
2620    J = J+1
      IF ( N.EQ.1 ) GO TO 300
      I = N-1
      IN1 = N
      IN2 = N+1
      ICO(1,IN1) = ICO(3,I) 
      ICO(2,IN1) = ICO(1,I) 
      ICO(3,IN1) = ICO(2,I) 
      ICO(4,IN1) = ICO(6,I) 
      ICO(5,IN1) = ICO(4,I) 
      ICO(6,IN1) = ICO(5,I) 
      ICO(1,IN2) = ICO(2,I) 
      ICO(2,IN2) = ICO(3,I) 
      ICO(3,IN2) = ICO(1,I) 
      ICO(4,IN2) = ICO(5,I) 
      ICO(5,IN2) = ICO(6,I) 
263     ICO(6,IN2) = ICO(4,I) 
      N = IN2+1
      GO TO 300
265     IF ( M.NE.2 ) GO TO 500
      ICO(1,N) = 3
      ICO(2,N) = 1
      ICO(3,N) = 0
      ICO(4,N) = 4
      ICO(5,N) = 4
      ICO(6,N) = 0
      IF ( L(2,2).EQ.12 ) ICO(6,N) = 2
      IF ( L(2,2).EQ.13 ) ICO(6,N) = 1
      IF ( L(2,2).EQ.14 ) ICO(6,N) = 3
      IF ( L(2,2).EQ. 3 ) ICO(3,N) = 2
      N = N+1
      IF ( L(2,2).EQ.3.AND.L(1,3).EQ.14.AND.L(1,4).EQ.11 ) GO TO 266
      IF ( L(1,3).EQ.14 ) GO TO 2561
      IF ( K.GT.2.OR.L(1,1).NE.7 ) GO TO 2581
      IF ( L(2,2)+L(3,2).NE.12 ) GO TO 2581
      L(2,2) = 0
      ICO(5,N-1) = 1
      GO TO 2581
266     IF ( L(1,1).NE.7.OR.L(1,2).NE.15 ) GO TO 500
      NCONT(J  ) = 16384
      NCONT(J+1) = 356
      NCONT(J+2) = 834
      NCONT(J+3) = 1177
      J = J+3
      GO TO 410
270     IF ( M.NE.2 ) GO TO 500
      IX = 1
      IF ( L(2,2).EQ.12 ) GO TO 2611
      IF ( L(2,2).EQ.13 ) GO TO 2612
      IF ( L(2,2).EQ.14 ) GO TO 2610
      IF ( L(2,2).EQ.15 ) GO TO 2611
      IF ( L(2,2).EQ.16 ) GO TO 2612
      IF ( L(2,2).EQ. 3 ) GO TO 271
      GO TO 2610
271     IX = 2
      GO TO 2610
275     IF ( M.NE.2 ) GO TO 500
      L(1,2) = L(2,2)
      L(2,2) = 3
      GO TO 218
300     CONTINUE
      N = N-1
      IJ = J
      IF ( N.GT.0 ) GO TO 301 
      J = J-1
      GO TO 410
301   DO 305 I=1,N
      NCONT(J) = ICO(1,I)+16*ICO(2,I)+128*ICO(3,I)+4*MOD(ICO(4,I),4)
     *      +64*MOD(ICO(5,I),4)+ 512*MOD(ICO(6,I),4)
305     J = J+1
306   J = J-1
CDB 4
C     PRINT 5,J,IJ,N,NSPGRP,NAXIS,NCONT
C    1 ,ICO
C   5 FORMAT ('0J =',I3,'  IJ =',I3,'  N =',I3,'  NSP =',2I5/ 
C    1 ' NCONT =',10I4/' ICO =',6I3/(6X6I3))
310   DO 400 I=IJ,J 
      NCI = I-IJ+1
      ICO4 = MOD(ICO(4,NCI),4)
      ICO5 = MOD(ICO(5,NCI),4)
      ICO6 = MOD(ICO(6,NCI),4)
      NI4 = NDAT(ICO4+4)
      NI5 = NDAT(ICO5+4)
      NI6 = NDAT(ICO6+4)
      KJ = I+1
      IF ( iand(NCONT(I),307)-290 ) 320,315,320
315     IF ( iand(NCONT(I),7884) ) 316,318,316
C   A CENTER IS PRESENT, NOT AT 0,0,0
316     CONTINUE
      WRITE(6,24) 
      GO TO 320
C   A CENTER AT 0,0,0 IS PRESENT
318     CONTINUE
320     IF ( KJ.GT.J ) GO TO 400
      DO 390 K=KJ,J
        NCK = K-IJ+1
        KCO4 = MOD(ICO(4,NCK),4)
        KCO5 = MOD(ICO(5,NCK),4)
        KCO6 = MOD(ICO(6,NCK),4)
        NK4 = NDAT(KCO4+4)
        NK5 = NDAT(KCO5+4)
        NK6 = NDAT(KCO6+4)
        NXC = MOD(NI4+NK4,4)
        NYC = MOD(NI5+NK5,4)
        NZC = MOD(NI6+NK6,4)
        MJ = K+1
        LOD = LODGDF(NCONT(K),NCONT(I))
        IF ( LOD.NE.290 ) GO TO 330
C   A CENTER IS GENERATED
        IF ( NXC+NYC+NZC.EQ.0 ) GO TO 330
        IF ( L(1,2).EQ.15.AND.I.EQ.IJ+1 ) GO TO 324
321       IF ( NXC.EQ.0 ) GO TO 322
        IF ( ICO(4,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NXC)*4
        IF ( ICO(4,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NXC)*4
322       IF ( NYC.EQ.0 ) GO TO 323
        IF ( ICO(5,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NYC)*64
        IF ( ICO(5,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NYC)*64
323       IF ( NZC.EQ.0 ) GO TO 325
        IF ( ICO(6,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NZC)*512
        IF ( ICO(6,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NZC)*512
        GO TO 325 
324       L(1,2) = 30
        NIJ = ICO(6,1)
        NIJ = NDAT(NIJ+4)
        NJX = MOD(NIJ+NXC,4)
        NCONT(IJ) = iand(NCONT(IJ),2035)+NDAT(NJX+4)*4
        NCONT(IJ) = iand(NCONT(IJ),1855)+NDAT(NIJ+4)*64
        GO TO 321 
325       CONTINUE
        ICO(4,NCI) = 0
        ICO(4,NCK) = 0
        ICO(5,NCK) = 0
        ICO(5,NCI) = 0
        ICO(6,NCI) = 0
        ICO(6,NCK) = 0
330       CONTINUE
        ICO4 = MOD(NCONT(I)/4  ,4)
        IKO4 = MOD(NCONT(K)/4  ,4)
        ICO5 = MOD(NCONT(I)/64 ,4)
        IKO5 = MOD(NCONT(K)/64 ,4)
        ICO6 = MOD(NCONT(I)/512,4)
        IKO6 = MOD(NCONT(K)/512,4)
        NI4 = NDAT(ICO4+4)
        NK4 = NDAT(IKO4+4)
        NI5 = NDAT(ICO5+4)
        NK5 = NDAT(IKO5+4)
        NI6 = NDAT(ICO6+4)
        NK6 = NDAT(IKO6+4)
        IF ( MJ.GT.J ) GO TO 390
        DO 380 M=MJ,J
          NCM = M-IJ+1
          MCO4 = MOD(ICO(4,NCM),4)
          MCO5 = MOD(ICO(5,NCM),4)
          MCO6 = MOD(ICO(6,NCM),4)
          NM4 = NDAT(MCO4+4)
          NM5 = NDAT(MCO5+4)
          NM6 = NDAT(MCO6+4)
          LOD1 = LODGDF(LOD,NCONT(M)) 
          NXC = MOD(NI4+NK4+NM4,4)
          NYC = MOD(NI5+NK5+NM5,4)
          NZC = MOD(NI6+NK6+NM6,4)
          IF ( LOD1.NE.290 ) GO TO 350
          IF  ( L(1,2).EQ.11 ) GO TO 334
          IF ( NXC+NYC+NZC.EQ.0 ) GO TO 380
          IF ( IJ.EQ.1 ) GO TO 331
          IF(NXC +NYC .EQ.0 )GO TO 333
          IF ( L(1,3).EQ.14 ) GO TO 331
          IF ( L(1,2).EQ.15 ) GO TO 331
          IF ( L(1,2).EQ.14 ) GO TO 331
          IF ( L(1,2).EQ.17 ) GO TO 331
          IF ( L(1,2)+L(1,3)+L(1,4).LT.18.AND.L(1,2).NE.9 ) GO TO 331
          NXC = MOD(NCONT(1)/2  +NXC,4)
          NYC = MOD(NCONT(1)/32 +NYC,4)
          NZC = MOD(NCONT(1)/256+NZC,4)
331         CONTINUE
          IF ( NXC.EQ.0 ) GO TO 332
          IF ( ICO(4,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NXC)*4
          IF ( ICO(4,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NXC)*4
          IF ( ICO(4,NCM).EQ.4 ) NCONT(M) = NCONT(M)+NDAT(NXC)*4
332         IF ( NYC.EQ.0 ) GO TO 333
          IF ( ICO(5,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NYC)*64
          IF ( ICO(5,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NYC)*64
          IF ( ICO(5,NCM).EQ.4 ) NCONT(M) = NCONT(M)+NDAT(NYC)*64
333         IF ( NZC.EQ.0 ) GO TO 335
          IF ( ICO(6,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NZC)*512
          IF ( ICO(6,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NZC)*512
          IF ( ICO(6,NCM).EQ.4 ) NCONT(M) = NCONT(M)+NDAT(NZC)*512
          GO TO 335
334         NCONT(I) = iand(NCONT(I),2035)
          NCONT(K) = iand(NCONT(K),1855)
          NCONT(M) = iand(NCONT(M),511)
335       ICO(4,NCI) = 0
          ICO(5,NCI) = 0
          ICO(6,NCI) = 0
          ICO(4,NCK) = 0
          ICO(5,NCK) = 0
          ICO(6,NCK) = 0
          ICO(4,NCM) = 0
          ICO(5,NCM) = 0
          ICO(6,NCM) = 0
          GO TO 380
350       IF ( LOD1.NE.0 ) GO TO 380
CDB 3
C     PRINT 9,LOD1,NXC,NYC,NZC,NI4,NK4,NI5,NK5,NI6,NK6
C   9 FORMAT (* LOD1 = *O10* NXC =*I4* NYC =*I4* NZC =*I4/
C    1 * NI4 =*I3* NK4 =*I3* NI5 =*I3* NK5 =*I3* NI6 =*I3* NK6 =*I3)
          IF ( NXC+NYC+NZC.EQ.0 ) GO TO 359
          IF ( L(1,2).EQ.30 ) GO TO 3500
          IF ( L(1,2).EQ.15 ) GO TO 355
          IF ( L(1,2).EQ.17 ) GO TO 359
          IF ( L(1,2).EQ.13 ) GO TO 351
          IF ( L(1,3).EQ.13 ) GO TO 352
          IF ( L(1,4).EQ.13 ) GO TO 353
          GO TO 500
3500        IF ( NXC.EQ.0 ) GO TO 359
          IF ( L(4,2).EQ.4 ) GO TO 359
          IF ( NXC.EQ.NK4 ) GO TO 359 
          NCONT(K) = ior(NCONT(K),NDAT(NXC)*4)
          IF ( L(1,4).EQ.13.AND.L(2,4).EQ.0 ) GO TO 3580
          GO TO 359
351         NYC = MOD(NK5+NM5,4)
          NZC = MOD(NK6+NM6,4)
          IF ( L(1,3).EQ.13 ) GO TO 354
          IF ( L(1,3).EQ.14 ) GO TO 354
          IF ( L(1,4).EQ.13 ) GO TO 354
          M2 = I+1
3511        CONTINUE
          IF(L(2,2).EQ.12.AND.(L(1,3).EQ.9 .OR.L(1,4).EQ.9))  GO TO 36
     *       0
          IF ( ICO(5,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NYC)*64
          IF ( ICO(5,NCM).EQ.4 ) NCONT(M) = NCONT(M)+NDAT(NYC)*64
3510        IF ( NZC.LE.0 ) GO TO 360
          IF ( ICO(6,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NZC)*512
          IF ( ICO(6,NCM).EQ.4 ) NCONT(M) = NCONT(M)+NDAT(NZC)*512
          GO TO 360
352         NXC = MOD(NI4+NM4,4)
          NZC = MOD(NI6+NM6,4)
          M2 = K+1
          IF(L(2,3).EQ.12 .AND.(L(1,2).EQ.9 .OR. L(1,4).EQ.9)) GO TO 3
     *       60
          IF ( NXC.LE.0 ) GO TO 3520
          IF ( ICO(4,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NXC)*4
          IF ( ICO(4,NCM).EQ.4 ) NCONT(M) = NCONT(M)+NDAT(NXC)*4
3520        IF ( NZC.LE.0 ) GO TO 360
          IF ( ICO(6,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NZC)*512
          IF ( ICO(6,NCM).EQ.4 ) NCONT(M) = NCONT(M)+NDAT(NZC)*512
          GO TO 360
353         NXC = MOD(NI4+NK4,4)
          NYC = MOD(NI5+NK5,4)
          IF(L(2,4).EQ.12 .AND. ( L(1,2).EQ.9.OR.L(1,3).EQ.9)) GO TO 3
     *       59
          IF ( NXC.LE.0 ) GO TO 3530
          IF ( ICO(4,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NXC)*4
          IF ( ICO(4,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NXC)*4
3530      IF ( NYC.LE.0 ) GO TO 359
          IF ( ICO(5,NCI).EQ.4 ) NCONT(I) = NCONT(I)+NDAT(NYC)*64
          IF ( ICO(5,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NYC)*64
          GO TO 359
354       IF ( NYC.LE.0 ) GO TO 3540
          IF ( IABS(L(2,2)-13).EQ.1 ) NCONT(I)=NCONT(I)+NDAT(NYC)*64
3540      IF ( NZC.LE.0 ) GO TO 3541
          IF ( ICO(6,NCK).EQ.4 ) NCONT(K) = NCONT(K)+NDAT(NZC)*512
3541      IF ( NXC.LE.0 ) GO TO 359
          IF ( L(1,4).EQ.0.OR.L(1,3).EQ.14 ) GO TO 359
          IF ( L(2,4).NE.12 ) NCONT(K)=NCONT(K)+NDAT(NXC)*4
          GO TO 359
355       L(1,2) = 30
          IF ( L(1,3).EQ.14 ) GO TO 370
          IF ( L(1,3).EQ.13.AND.L(2,3).EQ.12 ) GO TO 3561 
          IF ( L(1,1).EQ. 7.AND.L(2,2).EQ.12 ) GO TO 357
          IF ( L(1,1).EQ. 7.AND.L(2,2).EQ.14 ) GO TO 357
          IF ( L(1,1).EQ. 6.AND.L(2,2).EQ.12 ) GO TO 358
          IF ( L(1,1).EQ. 6.AND.L(2,2).EQ.14 ) GO TO 358
          IF ( L(1,1).EQ.7.AND.L(1,3).EQ.13 ) GO TO 3581
          IF ( MOD(NCONT(I),1024).EQ.19 ) GO TO 350
          IF ( NXC+NYC.GT.0 ) GO TO 356
          IF ( NZC.NE.NM6 ) GO TO 359 
          IF ( NZC.LE.0 ) GO TO 359
          NCONT(K) = NCONT(K)+NDAT(NZC)*512
          GO TO 359
356       IF ( L(3,2).EQ.10.OR.L(4,2).EQ.10 ) GO TO 359
          IF ( L(1,3).NE.10 ) GO TO 359
3561      IF ( L(2,2).EQ.3 ) GO TO 3560
          NCONT(I) = ior(NCONT(I),68)
          GO TO 3580
3560      IF ( L(1,4).EQ.2 ) NCONT(K)= ior(NCONT(K),512)
          GO TO 359
357       NCONT(I) = NCONT(I)+64
          IF ( L(3,2).NE.19 ) NCONT(I+1) = NCONT(I+1)-512 
          IF ( L(1,4).EQ.11.AND.L(3,2).EQ.19 )NCONT(I) = NCONT(I)+588
          GO TO 3580
358       NCONT(I) = NCONT(I)+136
3580      IF ( L(1,3).NE.13 ) GO TO 359
3581      M2 = K+1
          GO TO 360
359       M2 = M+1
360       CONTINUE
          IF ( M2.GT.J ) GO TO 306
CDB 2
C     PRINT 6,M2,I,J,K,M
C   6 FORMAT (*0M2 =*I3*  I =*I3*  J =*I3* K =*I3*  M =*I3) 
          DO 365 M1=M2,J
            NCONT(M1-1) = NCONT(M1)
            DO 365 M3=1,6
365             ICO(M3,M1-1) = ICO(M3,M1)
          GO TO 306
370         IF ( L(1,4).GT.12 ) GO TO 371
          M2 = I+1
          GO TO 360
371         CONTINUE
          L(1,2) = 13
          IF ( L(2,2).EQ.12.AND.L(1,1).EQ.5 ) NCONT(I) = NCONT(I)+200
          IF ( L(2,2).EQ.12.AND.L(1,1).EQ.6 ) NCONT(I) = NCONT(I)+136
          IF ( L(2,2).EQ.12.AND.L(1,1).EQ.7 ) NCONT(I) = NCONT(I)+652
          IF ( L(2,2).EQ.13 )                 NCONT(I) = NCONT(I)+68
          IF ( L(2,2).EQ.14 )                 NCONT(I) = NCONT(I)+140
          GO TO 359
380         CONTINUE
390       CONTINUE
400     CONTINUE
      LD = 0
      DO 405 I=1,J
405     LD = LODGDF(LD,NCONT(I))
      IF ( LD.NE.0 ) GO TO 410
      NCONT(J-1) = NCONT(J)
      J = J-1
410   CONTINUE
      NCONT(J+1) = 0
CDB 1
C     PRINT 4,K,NSPGRP,NAXIS,J,NCONT
450   RETURN
500   WRITE(6,2)SPG 
      NAXIS = 4
      GO TO 450
      END 
      SUBROUTINE ASSIGN
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      REAL LAMDA
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      DIMENSION KRR(2,IDSZ)
      ICX=0
      DO 4100 J =1,NPHASE
4100    ICX = ICX + ICR(J)
      IF(NPHASE.GT.1)CALL SORT(1000)
      DO 2 J=1,2
      DO 2 K=1,IDSZ
2         KRR(J,K)=0
      IF(LST3.EQ.0)GOTO 7
      DO 8 I=1,ICX
      IRL=MOD(IREFS(I),256)-128
      IRK=MOD(IREFS(I)/256,256)-128
      IRH=MOD(IREFS(I)/(256*256),256)-128
      IRC=MOD(IREFS(I)/(256*256*256),8) 
      IPHAS=IREFS(I)/(256*256*256*8)
      IF(MOD(I-1,60).EQ.0)WRITE(6,9)
8       WRITE(6,10)I,IRC,IRH,IRK,IRL,IPHAS,(REFS(I,JX),JX=1,2)
9     FORMAT(38H1NO.  CODE    H   K   L  PHASE  HW       ,
     * 16HPOSN              /)
10    FORMAT(1X,2I4,3X,3I4,I6,2F8.3)
7     DO 479 I=1,ICX
      if (nprof.eq.5) then
        RMIN=REFS(I,2)-WDT*fwhm(I,1)
        RMAX=REFS(I,2)+WDT*fwhm(I,2)
      else
        RMIN=REFS(I,2)-WDT*REFS(I,1)
        RMAX=REFS(I,2)+WDT*REFS(I,1)
      end if
      MIN=(RMIN-THMIN)/STEP+1.5
      MAX=(RMAX-THMIN)/STEP+1.5
      MIN=MAX0(MIN,1)
      MIN=MIN0(MIN,NPTS)
      MAX=MIN0(MAX,NPTS)
      DO 477 K=MIN,MAX
        KRR(2,K)=I
477   IF(KRR(1,K).EQ.0)KRR(1,K)=I
479   CONTINUE
      DO 482 J=1,NEXCRG
      IF(AHIGH(J).LE.THMIN)GOTO 482
      IN1=(ALOW(J)-THMIN)/STEP+1.5
      IN1=MAX0(IN1,1)
      IN2=(AHIGH(J)-THMIN)/STEP+1.5
      IN2=MIN0(IN2,NPTS)
      DO 483 I=IN1,IN2
        KRR(2,I)=1
483     KRR(1,I)=0
      IF(IN2.EQ.NPTS)GOTO 484
482     CONTINUE
484   IF(LST2.EQ.0)GOTO 530
      WRITE(6,589)THMIN,THMAX,STEP
589   FORMAT(13H1PATTERN FROM,F8.4,1X,2HTO,F8.4,1X,11HIN STEPS OF,
     * F8.4,1X,7HDEGREES)
      WRITE(6,544)
544   FORMAT(6H0 POSN,6X,3HI+B,5X,1HB,5X,1HI,7X,5H100*W,3X,3HK11,2X,
     * 3HK21)
      DO 520 J=1,NPTS
      PX=THMIN+FLOAT(J-1)*STEP
      YUN=Y(J)-BK(J)
      if(jobtyp.lt.3) WT=1./VAR(J)
      WRITE(6,559)PX,Y(J),BK(J),YUN,WT,(KRR(J2,J),J2=1,2) 
559     FORMAT(1X,F8.4,3F7.0,2PF9.4,2I5)
520     CONTINUE
530   CONTINUE
      DO 600 I=1,NPTS
      IF(KRR(2,I)-KRR(1,I).GT.NOV)THEN
        WRITE(7,601) I,KRR(2,I)-KRR(1,I),NOV
        WRITE(6,601) I,KRR(2,I)-KRR(1,I),NOV
601       FORMAT(' EXCESSIVE PEAK OVERLAP',/,
     *          5X,'AT THE',I5,'TH STEP THERE ARE ',I5,' REFLECTIONS',/,
     *               5X,'INCREASE THE VALUE OF *NOV* WHICH IS NOW',I5/)
        STOP 'EXCESSIVE PEAK OVERLAP' 
      ENDIF
600     CONTINUE
      DO 531 I=1,NPTS
      KR(I)=KRR(1,I)+IRS *KRR(2,I)
531     CONTINUE
      DO 602 I=1,NPTS
      IF(KR(I).NE.0.AND.KR(I).NE.IRS )GOTO 603
602     CONTINUE
      STOP 'NO REFLECTIONS FOUND '
C test for detecting the asymmetry model required    !cp ap 16 97 
C new code included in line 2 of ICF
C iasym = 0 (usual Rietveld asymmetry)
C iasym = 1 (new asymmetry. Riello, Canton & Fagherazzi.PD 10,3,204-206,1997)
603   if (iasym.eq.0)then
C THE FOLLOWING TEST IS REALLY ONLY VALID FOR THE SINGLE PHASE CASE
        IF(REFS(KRR(1,I),2).GE.RLIM.AND.LPAR(1,14).NE.0)THEN
          WRITE(6,604)
          WRITE(7,604)
         STOP 'ASYMMETRY PARAMETER USAGE INVALID'
        ENDIF
      else     
        IF(REFS(KRR(1,I),2).GE.(90.0-RLIM).AND.LPAR(1,14).NE.0)THEN
          WRITE(6,604)
          WRITE(7,604)
          STOP 'ASYMMETRY PARAMETER USAGE INVALID'
         endif
      endif
604     FORMAT(1X,'ASYMMETRY PARAMETER USAGE INVALID')
      DO 606 I=1,NPHASE
      IF(PAR(I,12).EQ.0..AND.PAR(I,13).EQ.0..AND.LPAR(I,13).NE.0)THEN
        WRITE(6,605)
        WRITE(7,605)
605     FORMAT(1X,'PREFERRED ORIENTATION USAGE INVALID')
        STOP 'PREFERRED ORIENTATION USAGE INVALID'
      ENDIF
606     CONTINUE
      RETURN
      END 
c
      SUBROUTINE CHISQ
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/G2/S1,S2,ss2,S3,ss4,D1,D2,D4,R1,R2,R3,r2nobk,r3nobk
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
C     COMMON/JNK/FILL(155),NBCKGD,NEXCRG,NSCAT
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
C      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
C     * LST3,IPL1,IPL2,IPLST,MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN
C     * ,IPREF,IABSR
C changing to italian codes
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      S1=0. 
      s1nobk=0.
      S2=0.    
      s2nobk = 0.
      ss2=0.
      S3=0.
      s3nobk = 0.
C     S4=0.
      ss4=0.
      D1=0.
      D2=0.
      D4=0.
C     DELW =  1.0E+25
C     DELWP=  1.0E+25
      sdelw = 1.0E+25
      sdelwp= 1.0E+25 
      d1nobk = 0.
      d2nobk = 0.
      DO 10 I=1,NPTS
      IF((NBCKGD.NE.0.AND.KR(I).EQ.0).OR.KR(I).EQ.IRS )GOTO 10
      IF(NEXCRG.GT.0) THEN
        DO 1400 IEXC=1,NEXCRG
          NPTLOW = (ALOW(IEXC)-THMIN)/STEP + 1.5
          NPTHI  = (AHIGH(IEXC)-THMIN)/STEP + 1.5
          IF (I.GE.NPTLOW.AND.I.LE.NPTHI) GO TO 10
1400        CONTINUE
      END IF
C     IF (DELW.LT.1.0E+24) DELWP = DELW
      if (sdelw.lt.1.0e+24) sdelwp = sdelw
      DEL=Y(I)-BK(I)-YC(I)
      DELW = DEL/SQRT(VAR(I))
      sdelw = del
      S1=S1+ABS(DEL) 
      s1nobk = s1nobk + del*(y(i)-bk(i))/y(i)                              !numerador para calcular Rp(-bcg)
      S2=S2+DEL*DEL/VAR(I)    
c      varbk = abs(var(i)-bk(i))
c      if(varbk.eq.0.) varbk = 1e-06                                       !variancia removendo o background
c      s2nobk = s2nobk+del*del/varbk                                       !nao deu certo  03nov00
      s2nobk=s2nobk+del*del*(y(i)-bk(i))*(y(i)-bk(i))/(y(i)*y(i))/var(i)   !numerador para calcular Rwp(-bcg)
      ss2 = ss2 + del*del
C     IF (DELWP.LT.1.0E+24) S4=S4+(DELW-DELWP)*(DELW-DELWP)
      if (sdelwp.lt.1.0e+24) ss4=ss4+(sdelw-sdelwp)*(sdelw-sdelwp)
10      CONTINUE 
      DO 15 I=1,NPTS
      IF((NBCKGD.NE.0.AND.KR(I).EQ.0).OR.KR(I).EQ.IRS )GOTO 15
      IF(NEXCRG.GT.0) THEN
        DO 1401 IEXC=1,NEXCRG
          NPTLOW = (ALOW(IEXC)-THMIN)/STEP + 1.5
          NPTHI  = (AHIGH(IEXC)-THMIN)/STEP + 1.5
          IF (I.GE.NPTLOW.AND.I.LE.NPTHI) GO TO 15
1401        CONTINUE
      END IF
      D1=D1+Y(I) 
      d1nobk = d1nobk+(y(i)-bk(i))                                !denominador para calcular Rp(-bcg)
      D2=D2+Y(I)*Y(I)/VAR(I)
c      if(varbk.eq.0.) varbk = 1e-06                             
      d2nobk = d2nobk + (y(i)-bk(i))*(y(i)-bk(i))/var(i)          !denominador para calcular Rwp(-bcg)
c                                                                 !Paiva-Santos 03nov00  (calcula Rwp sem o bck)
      S3=S3+BK(I)+YC(I)
c      s3nobk = s3nobk + yc(i)
15      CONTINUE
      R1=0
      R2=100.*S1/D1  
      r2nobk = 100.*s1nobk/d1nobk                                 !Rp without background 03nov00 (no writable)
      R3=100.*SQRT(S2/D2)    
      r3nobk = 100.*sqrt(s2nobk/d2nobk)                           !Rwp without background 03nov00 (no writable)
      RETURN
      END 
      SUBROUTINE ITER
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)   
      real lamda
      COMMON/F1/RJAC(MSZ,MSZ),VX(MSZ)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/G2/S1,S2,ss2,S3,ss4,D1,D2,D4,R1,R2,R3,r2nobk,r3nobk
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/BKGSCALE/SCABKG(99)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
C      COMMON/ALLP/FINAL(8*MSZ,2),ILOC
      COMMON/ALLP/FINAL(nfinal,2),ILOC
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),XMAS(16)
C      COMMON/AIRPAR/AI1,AI2,AI3,AI4,AI5,AI6,AI7,AI8,AI9,AI10
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/G3/COND,IORD1,IORD2,TH,NUM 
C   common !cp ap 16 97
      COMMON/FONDI/BKCOM(IDSZ),BKDIS(IDSZ),BKAM(IDSZ),BKPOL(IDSZ)
C      DIMENSION RJ(MSZ,MSZ),VLAST(MSZ),KS(50)   ! RJ(MSZ,MSZ) is not used !!!
      dimension vlast(msz),ks(50)               
C     EQUIVALENCE (RJ(1,1),KR(1))
      DIMENSION DYCDD(99),CSK(99),DISK(99),DERISO(NATS)
      DIMENSION ISODER(NATS)
      INTEGER FONDO
      REAL ISODER     

      DO 70 I=1,MSZ 
70      VLAST(I)=0. 
C     **** start cycles ****
      DO 10 IX=1,MCYCLE
        CALL ASSIGN 
        DO 1 I=1,MSZ
          VX(I)=0.
          DO 1 J=1,MSZ
1           RJAC(I,J)=0. 
        DO 2 I=1,NATS
          if (ibgd.eq.1) then
            isoder(i) = 1
          else
            ISODER(I)=0
          end if
2         CONTINUE
C
        NUM=0
        IPREV=1
        IF(NBCKGD.NE.0)GOTO 30
        TH=THMIN-STEP
        DO 31 I=1,NPTS
          TH=TH+STEP
          THX=TH/BKPOS-1.
          BK(I)=GLB(2)
          DO 29 J=2,6
29          BK(I)=BK(I)+GLB(J+1)*THX**(J-1)
31        IF(MCYCLE.EQ.1.AND.MAXS.EQ.0.AND.IPLPOL.EQ.1)BKPOL(I)=BK(I)
30      TH=THMIN-STEP
C           ****** START GREAT LOOP ON NPTS STEP-POINTS ******
        DO 9 IPM=1,NPTS
          if ( mod(ipm,npts/80+1) .eq. 0) print '(a\)', '.'
          TH=TH+STEP
C !cp jun 97 test for bck option
          if (ibgd.eq.1) then
            csk(1)    = 1
            disk(1)   = 1
            dycdd(1)  = 1
C            isoder = 1
            totcs  = 1
            goto 4334 
          end if
C     begin of italian code   !cp ap 20 97
C-----COMPUTE FPOL = POLARIZATION FACTOR FOR THIS IPM POINT
C     TH  = 2.0 * THETA
C     STH = SIN(THETA)
          STH = SIN(TH * 0.008726646)
          FPOL = (1.0 + ( (1.0 - 2.0 * STH * STH)**2 ) * CTHM ) / 2.0
          IF(KR(IPM).EQ.IRS) THEN
C ----IDERIV = 1 MEANS ONLY BACKGROUND (FOR THE STEP-POINTS THAT
C                                       ARE IN EXCLUDED REGIONS)
            IDERIV = 1
          ELSE
C-----IDERIV = 2 MEANS BACKGROUND + DERIVATIVES (FOR STEP-POINTS THAT ARE
C                                                IN INCLUDED REGIONS)
            IDERIV = 2
          ENDIF
C-----IF NECESSARY COMPUTE TOTCS :
C     TOTCS = TOTAL COMPTON INTENSITY SCATTERED BY ALL CRYSTALLINE
C             PHASES AT THE IPM-TH POINT OF THE X-RAY PATTERN.
          TOTCS = 0.0
          IF ( FONDO.EQ.1 .OR. FONDO.EQ.2 ) THEN
            DO 120 K = 1,NPHASE
              CALL COMPTON(K,STH,CISK)
              CSK(K) = CISK * FPOL
              TOTCS = TOTCS + SCABKG(K) * CSK(K)
  120         CONTINUE 
C-----COMPTON UPDATE
            IF ( MCYCLE.EQ.1 .AND. MAXS.EQ.0 .AND. IPLCOM.EQ.1 )
     1          BKCOM(IPM) = TOTCS
            BK(IPM) = BK(IPM) + TOTCS
          END IF
C-----IF NECESSARY CALL DIS = SUBROUTINE TO COMPUTE THERMAL AND LATTICE
C                             DISORDER SCATTERING SDK  AND  DERIVATIVES
C                             DYC IN THE K-TH PHASE
C     TOTDS = TOTAL DISORDER SCATTERING
          TOTDS = 0.0
          IF (FONDO.EQ.1.OR.FONDO.EQ.2) THEN
            DO 150 K = 1,NPHASE
              IOF = 0
              IF(K.GT.1) THEN
                DO 111 I = 2,K
111               IOF = IOF + NATOM(I-1)
              ENDIF
              DO 112 I=1,NATS
                DERISO(I)=0
112             CONTINUE
              CALL DISORDER(K,STH,IDERIV,SDK,DYC,FONDO,DERISO)
              DISK(K)  =   SDK * FPOL
             TOTDS    = TOTDS + SCABKG(K) * DISK(K)
C     UPDATING DERIVATE OF ISOTROPIC THERMAL FACTORS
              IF(FONDO.EQ.1) THEN
                DO 121 I=1,NATOM(K)
                  ISODER(IOF+I)= ISODER(IOF+I)+
     1                    SCABKG(K)*FPOL*DERISO(IOF+I)
121               CONTINUE
              END IF
C     UPDATING DERIVATE OF OVERALL THERMAL FACTOR
              IF(FONDO.EQ.2) DYCDD(K)=DYC*SCABKG(K)*FPOL
150           CONTINUE
C-----DISORDER UPDATE
            IF(MCYCLE.EQ.1 .AND. MAXS.EQ.0 .AND. IPLDIS.EQ.1 )
     1          BKDIS(IPM) = TOTDS
            BK(IPM) = BK(IPM) + TOTDS
          END IF
C------AMORPHOUS EVALUATIONS
          TOTAS =0.0
          TOTAS = GLB(20) * AMORPHOUS(IPM)
C-----AMORPHOUS UPDATE
          IF ( MCYCLE.EQ.1 .AND. MAXS.EQ.0 .AND. IPLAM.EQ.1 ) 
     1        BKAM(IPM) = TOTAS
C      IF(MCYCLE.EQ.1.AND.MAXS.EQ.0.AND.IPLAIR.EQ.1) BKAIR(IPM)=TOTAIR
C      BK(IPM) = BK(IPM) + TOTAS + TOTAIR
          BK(IPM) = BK(IPM) + TOTAS
C  end of italian code     !cp ap 20 97
4334      IF(KR(IPM).EQ.IRS )GOTO 9
          IORD1=MOD(KR(IPM),IRS )
          IORD2=MOD(KR(IPM)/IRS ,IRS )
          IF ( (IORD2.EQ.0 .OR. IORD1.EQ.0) .AND. NBCKGD.NE.0 ) GOTO 9
C     IF (IX.EQ.1) GOTO 15
C     IF ((IORD2-IORD1).LT.NOV) GOTO 15 
C     WRITE (6,1234) IORD2-IORD1,NOV,IPM
C1234 FORMAT(1X,'NO.OF REFLECTIONS IS ',I4,2X,'AT STEP NO.',I6)
C     STOP
15        NUM=NUM+1 
          IORDLIM2=IORD2
          IF(IPREV.GT.IORDLIM2) GO TO 40
          DO 21 J=IPREV,IORDLIM2 
21          CALL CALCUL(J)
c40        IPREV=MAX0(IPREV,IORD2+1)    !cp may 01 97
40        IPREV=MAX0(IPREV,IORDlim2+1)  
108       CALL SUMMAT(IPM,CSK,DISK,DYCDD,ISODER,TOTCS)
9         CONTINUE
        print '(a)', ':'
        IF ( JOBTYP.GT.2 ) RETURN 
        DO 50 I=1,MAXS
          DO 50 J=I,MAXS
50          RJAC(J,I) = RJAC(I,J)
        COND = DPINV(RJAC(1,1),VX(1),MAXS)
        CALL CHISQ
        OOM=S2/FLOAT(NUM-MAXS)
        DO 55 I=1,MAXS
          DO 55 J=1,MAXS
55          RJAC(I,J) = RJAC(I,J)*OOM
C              CODE TO ATTEMPT TO STABILIZE OSCILLATIONS
        DO 60 I=1,MAXS
          IF(SIGN(1.,VX(I)).EQ.SIGN(1.,VLAST(I)))GOTO 60
          IF(ABS(VX(I)).GT.1.2*ABS(VLAST(I)))GOTO 60
          IF(ABS(VX(I)).LT..8*ABS(VLAST(I)))GOTO 60
          VX(I)=VX(I)/2.
60        VLAST(I)=VX(I)
        CALL OUTPTR(IX)
        DO 5 I=1,MAXS
          X1=SQRT(ABS(RJAC(I,I)))*EPS
          IF(ABS(VX(I)).GT.X1)GOTO 10
5         CONTINUE
        IF (MAXS.GT.0) WRITE (6,1111)
1111    FORMAT(/,10X,'***** EPSED OUT *****'/)
        IF (MAXS.GT.0) ICYRUN = IX
        GOTO 20
10      CONTINUE 
C                 ***** END CYCLES *****
20    continue
C        DO 6 I=1,MAXS 
C        DO 6 J=1,MAXS
c6         RJ(I,J)=RJAC(I,J)
C     CODE FOR PRINTING CORRELATION MATRIX MOVED FROM EXPUT (65-85)
      IF(MAT.EQ.0.OR.MAXS.EQ.0) RETURN
      WRITE(6,7)
      IA=1
      LIM=19
38    LIM=MIN0(MAXS,LIM)
      WRITE(6,8)(I,I=IA,LIM)
      DO 39 I=1,MAXS
      L=0
      X=RJAC(I,I) 
      DO 400 J=IA,LIM
        L=L+1
        YY=RJAC(J,J)*X
400       KS(L)=100.*RJAC(I,J)/(SQRT(ABS(YY))*SIGN(1.,YY))+0.5
39      WRITE(6,12)I,(KS(J),J=1,L)
      IF(LIM.GE.MAXS)GOTO 37
      IA=LIM+1
      LIM=LIM+19
      GOTO 38
7     FORMAT(//,20H0CORRELATION MATRIX=)
8     FORMAT(5H0    ,19I6)
12    FORMAT(I5,19I6)
37    CONTINUE
      RETURN
      END 
      FUNCTION DPINV(A1,B1,N) 
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      DOUBLE PRECISION A,B,C,T,U
      DIMENSION A(MSZ,MSZ),A1(MSZ,MSZ),B(MSZ),B1(MSZ),U(MSZ),M(MSZ)
      IF ( N.NE.1 ) GOTO 50
      A1(1,1) = 1./A1(1,1)
      B1(1) = B1(1)*A1(1,1)
      DPINV = 1
      GOTO 99
50    DPINV = 0
      DO 1 I=1,N
        B(I)=B1(I)
        DO 1 J=1,N
1         A(I,J)=A1(I,J)  
      DO 2 I=1,N
        U(I)=DABS(A(I,1))
        C=U(I)
        DO 3 J=2,N
          U(I)=DMAX1(U(I),DABS(A(I,J))) 
3         C=C+DABS(A(I,J))
        DPINV=AMAX1(DPINV,real(C))
        M(I)=I
        B(I)=B(I)/U(I)
        DO 2 J=1,N
2         A(I,J)=A(I,J)/U(I)
      DO 10 I=1,N
      IF(I.EQ.N)GOTO 15
      J=I
      I1=I+1
      MI=M(I)
      S=DABS(A(MI,I))
      DO 4 K=I1,N 
        MI=M(K)
        IF(DABS(A(MI,I)).LE.S)GOTO 4
        S=DABS(A(MI,I))
        J=K
4         CONTINUE
      MK1=M(J)
      M(J)=M(I)
      M(I)=MK1
15      MI1=M(I)
      A(MI1,I)=1./A(MI1,I)
      B(MI1)=B(MI1)*A(MI1,I)
      DO 20 J=1,N 
        IF(I.EQ.J)GOTO 20
        A(MI1,J)=A(MI1,J)*A(MI1,I)
20        CONTINUE
      DO 30 J=1,N 
        IF(J.EQ.I)GOTO 30
        MJ=M(J)
        T=A(MJ,I)
        MI1=M(I)
        B(MJ)=B(MJ)-T*B(MI1)
        A(MJ,I)=-T*A(MI1,I)
        DO 40 K=1,N
          IF(K.EQ.I)GOTO 40 
          A(MJ,K)=A(MJ,K)-T*A(MI1,K)
40          CONTINUE
30        CONTINUE
10      CONTINUE
      DO 93 I=1,N
      MI1=M(I)
      B1(I)=B(MI1)
      DO 94 J=1,N 
        MJ=M(J)
94        A1(I,MJ)=A(MI1,J)/U(MJ)
93      CONTINUE
      C=0 
      DO 91 I=1,N
      D=ABS(A1(I,1))
      DO 92 J=2,N 
92        D=D+ABS(A1(I,J))
91      C=AMAX1(real(C),real(D))
      DPINV=DPINV*C 
99    RETURN
      END 
      SUBROUTINE DIRECT(SM,V,IPH)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON/F1/SMM(MSZ,MSZ),VX(MSZ)
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      DIMENSION V(6),SM(6,6)
      DO 1 J=1,6
        V(J)=PAR(IPH,J+5)
        K=LPAR(IPH,J+5)
        X=APAR(IPH,J+5)
        IF(K.EQ.0)X=0.
        DO 1 L=J,6
          M=LPAR(IPH,L+5)
          SM(L,J)=0.
          IF((M.EQ.0).OR.(K.EQ.0))GOTO 1
          SM(L,J)=SMM(M,K)*X*APAR(IPH,L+5)
1         SM(J,L)=SM(L,J)
      CALL ESD(SM,V,1.)
      RETURN
      END 
      SUBROUTINE ESD(SM,V,SUM)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      DIMENSION SM(6,6),V(6),S(6),R(6),E(6)
      I=0 
      DO 1 IA=1,3
        I=I+1
        IB=IA+1
        IF(IB.GT.3)IB=IB-3
        IC=IB+1
        IF(IC.GT.3)IC=IC-3
        ID=IA+3
        IE=IB+3
        IP=IC+3
        FNUM=4.*V(IB)*V(IC)-V(ID)*V(ID) 
        DEN=4.*V(IA)*V(IB)*V(IC)-V(IA)*V(ID)*V(ID)-V(IB)*V(IE)*V(IE)-V(I
     *   C)*V(IP)*V(IP)+V(ID)*V(IE)*V(IP)
        S(I)=SQRT(FNUM/DEN)
        R(IA)=-S(I)**4
        DENL=1./(DEN*DEN)
        R(IB)=DENL*(4.*V(IC)*DEN-(4.*V(IA)*V(IC)-V(IE)*V(IE))*FNUM)
        R(IC)=DENL*(4.*V(IB)*DEN-(4.*V(IA)*V(IB)-V(IP)*V(IP))*FNUM)
        R(ID)=-DENL*(2.*V(ID)*DEN-(2.*V(IA)*V(ID)-V(IE)*V(IP))*FNUM)
        R(IE)=DENL*FNUM*(2.*V(IB)*V(IE)-V(ID)*V(IP))
        R(IP)=DENL*FNUM*(2.*V(IC)*V(IP)-V(ID)*V(IE))
        E(I)=ERROR(SM(1,1),R(1),SUM)
        E(I)=E(I)/(2.*S(I))
        FNUM=V(IE)*V(IP)-2.*V(IA)*V(ID) 
        DEN=16.*V(IA)*V(IA)*V(IB)*V(IC)+V(IE)*V(IE)*V(IP)*V(IP)-4.*(V(IA
     *   )*V(IB)*V(IE)*V(IE)+V(IC)*V(IA)*V(IP)*V(IP))
        DENL=1./(DEN*DEN)
        FNUML=FNUM*FNUM
        S(I+3)=FNUML/DEN
        R(IA)=-4.*FNUM*DENL*(V(ID)*DEN+FNUM*(8.*V(IA)*V(IB)*V(IC)-V(IB)*
     *   V(IE)*V(IE)-V(IC)*V(IP)*V(IP)))
        R(IB)=-4.*FNUML*DENL*(4.*V(IA)*V(IA)*V(IC)-V(IA)*V(IE)*V(IE)) 
        R(IC)=-4.*FNUML*DENL*(4.*V(IA)*V(IA)*V(IB)-V(IA)*V(IP)*V(IP)) 
        R(ID)=-4.*FNUM/DEN*V(IA)
        R(IE)=2.*FNUM*DENL*(V(IP)*DEN-FNUM*(V(IE)*V(IP)*V(IP)-4.*V(IA)*V
     *   (IB)*V(IE)))
        R(IP)=2.*FNUM*DENL*(V(IE)*DEN-FNUM*(V(IE)*V(IE)*V(IP)-4.*V(IC)*V
     *   (IA)*V(IP)))
        E(I+3)=ERROR(SM(1,1),R(1),SUM)
        IF(S(I+3).EQ.0.0) GOTO 3
        E(I+3)=E(I+3)/(2.*SQRT(S(I+3)*(1.-S(I+3))))*180./3.14159265359
        S(I+3)=ATAN(SQRT((1.-S(I+3))/S(I+3)))*180./3.14159265359
        GOTO 4
3       S(I+3)=90.
        E(I+3)=0.0
c4       CONTINUE
C       IF(FNUM.LT.0.)S(I+3)=180.-S(I+3)   
C                     !cp ap 23 97
4       IF(FNUM.LT.0.)S(I+3)=180.-S(I+3)
        IF(SM(ID,ID).EQ.0.)E(ID)=0.
1       CONTINUE
      IF(ABS(S(1)-S(2)).LT..00008)E(1)=E(2)
      DO 2 I=1,6
      SM(I,I)=E(I)
2       V(I)=S(I)
      RETURN
      END 
      FUNCTION ERROR(A,B,OMEGA)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      DIMENSION A(6,6),B(6)
      SUM=0.
      DO 1 I=1,6
      DO 1 J=I,6
        X=2.
        IF(I.EQ.J)X=1.
1         SUM=SUM+A(I,J)*B(I)*B(J)*X
      IF(SUM.LT.0.)SUM=0.
      ERROR=SQRT(SUM*OMEGA)
      RETURN
      END 
      SUBROUTINE STEEP(X,N,M,CHISQ,GRAD)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      DIMENSION Z(16),X(16),T(16),ZNT(16),SIGT(16),SIG(16)
      DIMENSION Y(16,16),XI(16),ZN(16)
      IF(N.EQ.0)RETURN
      CALL GRAD(X,Z)
      CALL CHISQ(X,G)
      DO 10 I=1,M
      DO 10 J=1,M 
        Y(I,J)=0.0
        IF(I.EQ.J)Y(I,J)=1.0
10        CONTINUE
      DO 200 K=1,N
      U=0
      DO 20 I=1,M 
        T1=0
        DO 30 J=1,M
30          T1=T1+Z(J)*Y(I,J) 
        SIG(I)=T1 
20        U=U-SIG(I)*Z(I)
      ETA=AMIN1(1.0,ABS(2.*G/U))
      DO 40 I=1,M 
40        T(I)=X(I)-ETA*SIG(I)
      CALL GRAD(T,ZN)
      V=0
      DO 50 I=1,M 
50        V=V-ZN(I)*SIG(I)
      CALL CHISQ(T,G1)
      D=3.*ABS(G-G1)/ETA+U+V
      W=SQRT(D*D-U*V)
      ALPHA=ETA*(1.-(V+W-D)/(V-U+2.*W))
      DO 60 I=1,M 
        SIG(I)=-1.0*ALPHA*SIG(I)
        X(I)=X(I)+SIG(I)
60        ZN(I)=Z(I)
      CALL GRAD(X,Z)
      CALL CHISQ(X,G)
      T1=0
      DO 70 I=1,M 
        XI(I)=Z(I)-ZN(I)
70        T1=T1+XI(I)*SIG(I)
      DO 80 I=1,M 
80        SIGT(I)=SIG(I)/T1
      T1=0
      DO 90 I=1,M 
        T2=0
        DO 100 J=1,M
100         T2=T2+XI(J)*Y(I,J)
        ZN(I)=T2
90        T1=T1+ZN(I)*XI(I)
      DO 110 I=1,M
        T2=0
        DO 120 J=1,M
120         T2=T2+XI(J)*Y(J,I)
110       ZNT(I)=T2/T1
      DO 131 I=1,M
        DO 130 J=1,M
130         Y(I,J)=Y(I,J)+SIG(I)*SIGT(J)-ZN(I)*ZNT(J)
131       CONTINUE
200     CONTINUE
      RETURN
      END 
      SUBROUTINE RTMT(IPRT,IPHASE)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
C    THE MATRICES ARE IN THE FIRST THREE ROWS OF VCTR AND THE VECTORS 
C    ARE IN THE FOURTH ROW OF VCTR.
C    THE MATRICES ARE PACKED INTO IVEC. 
      COMMON/BLNK1/N,NMAX,MULTX,Y(192,3),XLT,IXB,FILL(192)
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      COMMON /RTMTX/ IVEC(99,192),MLTPHS(99),ICNTPHS(99)
      COMMON/MULTIP/TMASSA(99),MLTPHASE,xmltp(99),MURT(NATS),SAQF(99)
     *,wtis(99)
      common/simoper/isimop
      DIMENSION VCTR(3,4)
6     FORMAT(3H0 (,3F3.0,18H )     ( X )     (,F6.3,13H )     ( X2 ) /
     *  3H ( ,3F3.0,18H  ) * (  Y  ) + ( ,F6.3,14H  ) = (  Y2  )/
     *   3H  (,3F3.0,18H )     ( Z )     (,F6.3,13H )     ( Z2 )
     *    ,10X,5HIVEC(,I3,3H) =,I15)
10    FORMAT ('1The operations of the space group are')
      Y(1,1) = 2.0/108.0
      Y(1,2) = 3.0/108.0
      Y(1,3) = 4.0/108.0
      CALL SYMOPR
      MULTP = MULTX-1
      mltphase=MULTP
      DO 310 I=1,MULTP
        IVEC(IPHASE,I) = 0
        DO 300 K=1,3
          IX = Y(I,K)*108.0+0.1
          VCTR(K,4) = FLOAT((IX+4)/9)/12.0
          VCTR(K,4) = AMOD(VCTR(K,4),1.0)
          IY = MOD(IX+4,9)-4
          IZ = IABS(IY)
          GO TO ( 220,230,240,250),IZ
          CALL GOTOER
220       VCTR(K,1) = -1.0*FLOAT(ISIGN(1,IY))
          VCTR(K,2) =  1.0*FLOAT(ISIGN(1,IY))
          VCTR(K,3) = 0.0
          GO TO 300 
230       VCTR(K,1) =  1.0*FLOAT(ISIGN(1,IY))
          VCTR(K,2) = 0.0
          VCTR(K,3) = 0.0
          GO TO 300 
240       VCTR(K,1) = 0.0
          VCTR(K,2) =  1.0*FLOAT(ISIGN(1,IY))
          VCTR(K,3) = 0.00
          GO TO 300 
250       VCTR(K,1) = 0.0
          VCTR(K,2) = 0.0
          VCTR(K,3) =  1.0*FLOAT(ISIGN(1,IY))
300       IVEC(IPHASE,I) = IVEC(IPHASE,I)+(9*(INT(VCTR(K,1)+1.5))
     *                    +3*(INT(VCTR(K,2)+1.5)) 
     *                    +INT(VCTR(K,3)+1.5))*32768*32**(3-K)
     *                    +(INT(VCTR(K,4)*12.0+.5)+16)*32**(3-K)
          IF ( IPRT.LE.0 ) GOTO 310
      if(isimop.eq.1)goto 310
          IF ( MOD(I-1,15).EQ.0)WRITE(6,10)
          WRITE(6,6)((VCTR(K,J),J=1,4),K=1,3),I,IVEC(IPHASE,I)
310   CONTINUE
      MLTPHS(IPHASE)=MULTP
      ICNTPHS(IPHASE)=1
      IF ( NC.EQ.0 ) ICNTPHS(IPHASE)=2
      RETURN
      END 
      SUBROUTINE OPERTR(I,L2) 
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
C
C                             THIS SR GENERATES THE SYMMETRY EQUIVALENT
C                           POSITIONS. IT IS CALLED BY SYMOPR.
C
      COMMON/BLNK1/ N(2),L,X(192),Y(192),Z(192),XLT,IXB,FILL2(192)
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP

      I1 = MOD(NCONT(I),4)+1
      I2 = MOD(NCONT(I)/4,4)+1
      I3 = MOD(NCONT(I)/16,4)+1
      I4 = MOD(NCONT(I)/64,4)+1
      I5 = MOD(NCONT(I)/256,2)
      I6 = MOD(NCONT(I)/512,4)+1
      I7 = MOD(NCONT(I)/2048,4)+1 
      I8 = NCONT(I)/8192+1
      GO TO (100,200,300),I8
      CALL GOTOER
100   GO TO (101,102,103,104),I1
      CALL GOTOER
101   X(L) = X(L2)
      GO TO 105
102   X(L) = Y(L2)
      GO TO 105
103   X(L) = -X(L2) 
      GOTO 105
104   X(L) = -Y(L2) 
105   GO TO (110,107,108,106),I2
      CALL GOTOER
106   X(L) = X(L)+0.25
107   X(L) = X(L)+0.25
108   X(L) = X(L)+0.25
110   GO TO (111,112,113,114),I3
      CALL GOTOER
111   Y(L) = Y(L2)
      GO TO 115
112   Y(L) = X(L2)
      GO TO 115
113   Y(L) = -Y(L2) 
      GOTO 115
114   Y(L) = -X(L2) 
115   GO TO (120,117,118,116),I4
      CALL GOTOER
116   Y(L) = Y(L)+0.25
117   Y(L) = Y(L)+0.25
118   Y(L) = Y(L)+0.25
120   Z(L) = Z(L2)*FLOAT(1-2*I5)
      GO TO (125,122,123,121),I6
      CALL GOTOER
121   Z(L) = Z(L)+0.25
122   Z(L) = Z(L)+0.25
123   Z(L) = Z(L)+0.25
125   GO TO (130,127,126),I7
      CALL GOTOER
126   Z(L) = Z(L)+0.333333333 
127   Z(L) = Z(L)+0.333333333 
130   L = L+1
      CALL CEL000
      RETURN
200   X(L) = -Y(L2) 
      Y(L) = X(L2)-Y(L2)
      Z(L) = Z(L2)+0.333333333*FLOAT(I2-2)
      L = L+1
      CALL CEL000
      X(L) = -Y(L-1)
      Y(L) = -X(L2) 
      Z(L) = Z(L2)-0.333333333*FLOAT(I2-2)
      GO TO 130
300   I1 = L2
301   X(L) = Z(I1)
      Y(L) = X(I1)
      Z(L) = Y(I1)
      IF ( I1-L2 ) 302,302,130
302   I1 = L
      L = L+1
      CALL CEL000
      GO TO 301
      END 
      SUBROUTINE SYMOPR
C
C                             THIS SR TAKES A POINT XYZ, DETERMINES THE
C                           SITE MULTIPLICITY, THE SITE SYMMETRY, AND 
C                           THE SYMMETRY CONSTRAINTS ON THE POSITIONAL
C                           AND THERMAL PARAMETERS.
C
      COMMON/BLNK1/ N(2),L,X(192),Y(192),Z(192),XLT,KXB,NCTR(192)
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP

C     INTEGER AND,OR
C     AND(I,J) = IAND(I,J)
C     OR(I,J) = IOR(I,J)

      ISCK = 0
103   L = 2
      ISRH = 0
      KXB = 0
      XLT = 1.0
      J = 1
      CALL CEL000
      L1=J
      DO 100 I=1,8
        IOP = 0
        IF ( NCONT(I).EQ.0 ) GO TO 104
1       DO 99 L2=J,L1
          L3=L-1
          IF ( NCONT(I)-8192 ) 4,5,4
5         X(L)=X(L2)+1.0/3.0
          Y(L)=Y(L2)+2.0/3.0
          Z(L)=Z(L2)+2.0/3.0
          X(L+1)=X(L2)+2.0/3.0
          Y(L+1)=Y(L2)+1.0/3.0
          Z(L+1)=Z(L2)+1.0/3.0
          L=L+1
          CALL CEL000
          L=L+1
          CALL CEL000
          GO TO 6
4         CALL OPERTR(I,L2)
6         L5=L
          L6=NCONT(I)
7         L = L3+1
8         L4 = 1
          IF ( ABS(X(L)-X(L4))-0.0001 ) 2,2,90
2         IF ( ABS(Y(L)-Y(L4))-0.0001 ) 3,3,90
3         IF ( ABS(Z(L)-Z(L4))-0.0001 ) 91,91,90
91        CONTINUE
          IF ( NCTR(L)-8192 ) 92,94,93
92        KX = MOD(NCTR(L),4)+1
          KY = MOD(NCTR(L)/16,4)+1
          KZ = MOD(NCTR(L)/256,2)
          IF ( IOP.EQ.0 ) XLT=XLT/2.0
          IOP = 1
          GO TO (920,930,940,950),KX
          CALL GOTOER
920       GO TO (921,1000,922,1000),KY
          CALL GOTOER
921       IF ( KZ.EQ.0 ) GO TO 1000
        KXB = ior(KXB,2621450)
        GO TO 90
922       CONTINUE
        IF ( KZ.EQ.0 ) KXB= ior(KXB,2228292)
        IF ( KZ.EQ.1 ) KXB= ior(KXB, 655432)
        GO TO 90
930       KXB = ior(KXB,32768)
        GO TO (1000,931,1000,932),KY
        CALL GOTOER
931       CONTINUE
        IF ( KZ.EQ.0 ) KXB= ior(KXB, 4194464)
        IF ( KZ.EQ.1 ) KXB= ior(KXB,12583048)
        GO TO 90
932       KXB = ior(KXB,2753088)
        IF ( KZ.EQ.1 ) KXB = ior(KXB,8)
        GO TO 90
940       KXB = ior(KXB,512)
        GO TO (941,1000,942,1000),KY
        CALL GOTOER
941       KXB = ior(KXB,131072)
        IF ( KZ.EQ.1 ) GO TO 943
        KXB = ior(KXB,525312)
        GO TO 90
942       KXB = ior(KXB,64)
        IF ( KZ.EQ.1 ) GO TO 944
        KXB = ior(KXB,2623488)
        GO TO 90
943       KXB = ior(KXB,2097152)
944       KXB = ior(KXB,8)
        GO TO 90
950       KXB = ior(KXB,32768)
        GO TO (1000,951,1000,952),KY
        CALL GOTOER
951       KXB = ior(KXB,2753088)
        IF ( KZ.EQ.1 ) KXB = ior(KXB,8)
        GO TO 90
952       IF ( KZ.EQ.0 ) KXB = ior(KXB,8388864)
        IF ( KZ.EQ.1 ) KXB= ior(KXB,4194568)
        GO TO 90
93        CONTINUE
        IF ( MOD(NCTR(L),16384).EQ.0 ) GO TO 95 
        ISRH = 1
        IF ( iand(307,NCTR(L)-NCONT(I)).EQ.0 ) GO TO 92
        GO TO 90
94        KXB = ior(KXB,2916928)
        XLT = XLT/3.0
        GO TO 99
95        KXB = ior(KXB,17924240)
        XLT = XLT/3.0
        GO TO 99
90        CONTINUE
        IF ( L6-8192 ) 99, 98,98
98        L6=1
        L = L+1
        GO TO 8
99        L=L5
100     L1=L-1
104   CONTINUE
      IF ( ISRH.EQ.0 ) GO TO 200
      IF ( ISCK.GT.1 ) GO TO 2000
      ISCK = ISCK+1 
      DO 1010 J=1,I 
      IF ( NCONT(J).EQ.16384 ) GO TO 1011 
1010    CONTINUE
      DO 101 L2 = 1,L1
      IF ( ABS(X(L2)-Y(L2)).LE.0.0001 ) GO TO 102
      IF ( ABS(X(L2)+Y(L2)-1.0).LE.0.0002 ) GO TO 102
101     CONTINUE
      ISCK = 2
      WRITE(6,110) X(1),Y(1),Z(1),XLT
      GO TO 200
1011  CONTINUE
      IF ( ABS(X(1)-Y(1)).LT.0.0001 ) GO TO 120
      WRITE(6,111)X(1),Y(1),Z(1),Y(1),Z(1),X(1)
      X3 = X(1)
      X(1) = Y(1)
      Y(1) = Z(1)
      Z(1) = X3
      GO TO 103
110   FORMAT (37H0NO X,X,Z OR X,-X,Z SET WAS FOUND FOR,4F8.5)
102   IF ( L2.EQ.1 ) GO TO 120
      WRITE(6,111)X(1),Y(1),Z(1),X(L2),Y(L2),Z(L2)
111   FORMAT (12H0THE ATOM AT ,3F8.5,11H WAS PUT AT,3F8.5)
      X(1) = X(L2)
      Y(1) = Y(L2)
      Z(1) =Z(L2)
      GO TO 103
1000  WRITE(6,112) NCTR(L)
112   FORMAT (22H0THE SYMMETRY OPERATOR ,I10,9H IS WRONG)
      STOP 7701
2000  CONTINUE
      WRITE(6,113)X(1),Y(1),Z(1)
113   FORMAT (  '0THERE ARE UNRESOLVABLE PROBLEMS WITH THE POSITION SET'
     *  ,3F8.5,/,' ATTEMPTS TO DEFINE THE TRUE SITE SYMMETRY WERE ',
     * 'ABANDONED.')
120   IF ( ABS(X(1)-Y(1))+ABS(X(1)-Z(1)).LE.0.0002 ) GO TO 200
      DO 130 L2=2,L1
      IF ( ABS(X(L2)-Y(L2))+ABS(X(L2)-Z(L2)).LE.0.0002 ) GO TO 102
130     CONTINUE
200   CONTINUE
      KXB = MAX0( iand(KXB,18874367), iand(KXB,16777215))
      IF (  iand(KXB,384).EQ.384 ) KXB= ior(KXB,512)
      RETURN
      END 
      SUBROUTINE CEL000

C                             THIS SR IS USED TO FORCE A POINT XYZ TO 
C                           LIE IN THE UNIT CELL AT 000.
C
      COMMON/BLNK1/ N(2),L,X(192),Y(192),Z(192),XLT,IXB,FILL(192)
      X(L-1) = AMOD(X(L-1)+8.0,1.0)
      IF ( X(L-1)-.99999 ) 2,1,1
1     X(L-1) = 0.0
2     Y(L-1) = AMOD(Y(L-1)+8.0,1.0)
      IF ( Y(L-1)-.99999 ) 4,3,3
3     Y(L-1) = 0.0
4     Z(L-1) = AMOD(Z(L-1)+8.0,1.0)
      IF ( Z(L-1)-.99999 ) 6,5,5
5     Z(L-1) = 0.0
6     RETURN
      END 
      SUBROUTINE OP1(IPHASE)
        include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
C
C                             THIS SR GENERATES THE OPERATORS FOR
C                           SMTRY2, DETERMINES THE POLARITY OF THE SPACE
C                           GROUP, NOTES THE PRESENCE OF A CENTER(1BAR),
C                           AND GENERATES A SET OF PSEUDO-OPERATORS FOR
C                           USE IN SYMOPR TO DEFINE SITE SYMMETRY.
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      COMMON/BLNK1/N(581),NCTR(192)
      COMMON /HKLCTL/ IHKL(3,48),AZ(48),NC1(10,7,99),ICHKL(99),
     1    N1HKL(99),IER
C     INTEGER AND
C     AND(I,J) = IAND(I,J)
      K = 2
      NPOL = 7
      DO 40 I=1,7
40      NC1(1,I,IPHASE) = 0
      NCTR(1) = 0
      DO 100 I=1,8
        IF ( NCONT(I).LE.0 ) GO TO 101
        IF ( NCONT(I).LT.8192 ) GO TO 80
        IF ( NCONT(I).LE.08192 ) GO TO 60
        IF ( NCONT(I).EQ.16384 ) GO TO 70
        NPOL = 4
        NC1(I,1,IPHASE) = 2
50      DO 51 J=2,7 
51        NC1(I,J,IPHASE) = 0
        NC1(I,6,IPHASE) = MOD(NCONT(I)/4,32)
        GO TO 90
60      NC1 (I,1,IPHASE) = 3
        NPOL = 1
        GO TO 50
70      NC1(I,1,IPHASE) = 1
        GO TO 50
80      NC1(I,1,IPHASE) = 0
        NC1(I,2,IPHASE) = MOD(NCONT(I),4)
        NC1(I,3,IPHASE) = MOD(NCONT(I)/4,4)
        NC1(I,4,IPHASE) = MOD(NCONT(I)/16,4)
        NC1(I,5,IPHASE) = MOD(NCONT(I)/64,4)
        NC1(I,6,IPHASE) = MOD(NCONT(I)/256,2)
        NC1(I,7,IPHASE) = MOD(NCONT(I)/512,16) 
90      CONTINUE
        L = K-1
        DO 100 J=1,L
          IF ( NCONT(I).EQ.8192 ) GO TO 95
          NCTR(K) = iand( ieor(NCONT(I),NCTR (J )),57651)
          IF( iand(NCONT(I),17).EQ.0) GO TO 92
          IF( iand(NCTR(J),34).EQ.0) GO TO 92
          IF( iand(NCTR(J),34).EQ.34) GO TO 92
          NCTR(K)= ieor(NCTR(K),34)
92        CONTINUE
          K = K+1
          IF ( NCONT(I).LT.8192 ) GO TO 100
          NCTR(K) = NCTR(K-1)+32768
          K = K+1
          GO TO 100 
95        CONTINUE
          NCTR(K  ) = 0
          NCTR(K+1) = 0
          K = K+2
100       CONTINUE
101   MULTP = K-1
      N1HKL(IPHASE) = I-1
      IF ( N1HKL(IPHASE).LT.2 ) GO TO 105
      DO 102 N2=2,N1HKL(IPHASE)
        IF ( NC1(N2-1,1,IPHASE).EQ.1 ) GO TO 103
102     CONTINUE
      GO TO 105
103   DO 104 I=N2,N1HKL(IPHASE)
        DO 104 J=1,7
104       NC1(I-1,J,IPHASE) = NC1(I,J,IPHASE)
      NC1(N1HKL(IPHASE),1,IPHASE) = 1 
105   CONTINUE
      NC = 1
      DO 200 I=1,MULTP
        IF (  ieor(NCTR(I),290).EQ.0 ) NC = 0
        IF ( iand(NCTR(I),3).GT.0 ) NPOL = iand(NPOL,6)
        IF ( iand(NCTR(I),48).GT.0 ) NPOL = iand(NPOL,5)
200     IF ( iand(NCTR(I),256).GT.0 ) NPOL = iand(NPOL,3)
      RETURN
      END 
      SUBROUTINE SMTRY2(IPHASE)
        include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
C                             THIS SR GENERATES THE FULL SET OF
C                           EQUIVALENT HKLS FROM ANY MEMBER OF THE SET.
C                           IT ALSO DETERMINES THE PHASE SHIFTS RELATIVE
C                           TO THE INPUT HKL PHASE.
      COMMON /HKLCTL/ IHKL(3,48),AZ(48),NC1(10,7,99),ICHKL(99),
     1    N1HKL(99),IER
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      DIMENSION KD(6),KE(3),KF(3)
      DIMENSION NCH(48),NCK(48),NCL(48) 
      DIMENSION LD(3,48)
      DIMENSION IPH(48)
      DATA KD(1),KD(2),KD(3),KD(4),KD(5),KD(6)/0,6,3,9,4,8/ 
      DATA KE(1),KE(2),KE(3),KF(1),KF(2),KF(3)/0,4,8,0,8,4/ 
      
      IF ( NC.GT.0 ) GO TO 14 
      IF ( IHKL(1,1) ) 11,9,14
9     IF ( IHKL(2,1) ) 12,10,14
10    IF ( IHKL(3,1) ) 13,14,14
11    IHKL(1,1) = -IHKL(1,1)
12    IHKL(2,1) = -IHKL(2,1)
13    IHKL(3,1) = -IHKL(3,1)
14    CONTINUE
      IPH(1) = IHKL(3,1)+512*(IHKL(2,1)+512*IHKL(1,1))
      CI = 1.0
      ICHKL(IPHASE) = 1
      AZ(1) = 0.0
      LD(1,1) = 0
      LD(2,1) = 0
      LD(3,1) = 0
      NCH(1) = 0
      NCK(1) = 0
      NCL(1) = 1
      IF ( N1HKL(IPHASE).EQ.0 ) GO TO 1000
      L = 2
      DO 900 I=1,N1HKL(IPHASE) 
        CI = CI*2.0 
        IF ( NC1(I,1,IPHASE).GT.0 ) CI = CI*1.5
        DO 890 J=1,ICHKL(IPHASE)
          NCO = NC1(I,1,IPHASE)+1
          GO TO (20,700,800,15),NCO
          CALL GOTOER
15        IF(MOD(IHKL(1,1)-IHKL(2,1)-IHKL(3,1),3).NE.0 ) GO TO 1002
          GO TO 900 
20        NCH(L) =  ieor(NCH(J),NC1(I,2,IPHASE))
          IF ( MOD(NC1(I,2,IPHASE),2).NE.0 ) THEN
            NCH(L) = ieor(NCK(J),NC1(I,2,IPHASE))
          END IF
          M = 1+MOD(NCH(L),2) 
          MS = 1-2*MOD(NCH(L)/2,2)
          N = IABS(NCL(J))
          IHKL(M,L) = IHKL(1,N)*MS
          M = 1+MOD(NC1(I,2,IPHASE),2)
          MS = 1-2*MOD(NC1(I,2,IPHASE)/2,2)
          NCX = NC1(I,3,IPHASE)+1
          LD(1,L) = KD(NCX)+LD(M,J)*MS
          NCK(L) =  ieor(NCK(J),NC1(I,4,IPHASE))
          IF ( MOD(NC1(I,4,IPHASE),2).NE.0 ) THEN
            NCK(L)= ieor(NCH(J),NC1(I,4,IPHASE))
          END IF
          M = 2-MOD(NCK(L),2) 
          MS = 1-2*MOD(NCK(L)/2,2)
          IHKL(M,L) = IHKL(2,N)*MS
          M = 2-MOD(NC1(I,4,IPHASE),2)
          MS = 1-2*MOD(NC1(I,4,IPHASE)/2,2)
          NCX = NC1(I,5,IPHASE)+1
          LD(2,L) = KD(NCX)+LD(M,J)*MS
          MS = (1-2*NC1(I,6,IPHASE))*ISIGN(1,NCL(J))
          NCL(L) = MS*NCL(N)
          IHKL(3,L) = IHKL(3,N)*MS
70        NCO = NC1(I,7,IPHASE)+1
          IF ( NCO.GT.5 ) NCO=6
          MS = 1-2*NC1(I,6,IPHASE) 
          IF ( N.GT.1.AND.NCK(J).GE.4.AND.M.EQ.1 ) MS=-MS
          LD(3,L) = KD(NCO)+LD(3,J)*MS
80        IXIT = 0
81        CONTINUE
          IPH(L) = IHKL(3,L)+512*(IHKL(2,L)+512*IHKL(1,L)) 
87        LA = LD(1,L)*IHKL(1,1)+LD(2,L)*IHKL(2,1)+LD(3,L)*IHKL(3,1) 
          AZ(L) = FLOAT(LA)/12.0
CDB 2
C     IF ( J1(20).GT.0 ) PRINT 1,L,(IHKL(M,L),M=1,3),(LD(M,L),M=1,3),
C    Z AZ(L),NCH(L),NCK(L),NCL(L),(NC1(I,M,IPHASE),M=1,7)
          IF ( NC.NE.0 ) GO TO 89
          IF ( NSPGRP.LT.12 ) GO TO 89
          IF ( NC1(I,1,IPHASE).GT.0 ) GO TO 89
          IF ( MOD(NC1(I,2,IPHASE),2).EQ.1 ) GO TO 89
          IF ( IPH(L).LE.0 ) GO TO 680
89        CONTINUE
          DO 90 M=2,L
            IF ( IPH(M-1).EQ.IPH(L) ) GO TO 690
            IF ( NC.NE.0 ) GO TO 90
            IF ( IPH(M-1).EQ.-IPH(L) ) GO TO 690
90          CONTINUE
          L = L+1
680       IF ( IXIT-1 ) 890,750,850
690       IF ( AMOD(ABS(AZ(L)-AZ(M-1)),1.).NE.0. ) GO TO 1002 
          GO TO 680 
700       IXIT = 2
          JD = J
750       IXIT = IXIT-1
          LD(1,L) = LD(1,JD)
          LD(2,L) = LD(2,JD)
          LD(3,L) = LD(3,JD)
          NCH(L) = 0
          NCK(L) = 0
          NCL(L) = L
          K1JD = IHKL(1,JD)
          IHKL(1,L) = IHKL(2,JD)
          IHKL(2,L) = IHKL(3,JD)
          IHKL(3,L) = K1JD
          JD = L
          GO TO 81
800       IXIT = 4
          JD = J
          NCO = NC1(I,6,IPHASE)
850       IXIT = IXIT-2
          NCH(L) = 0
          NCK(L) = 4
          NCL(L) = L
          LD(1,L) = 0
          LD(2,L) = 0
          IHKL(1,L) = -IHKL(1,JD)-IHKL(2,JD)
          IHKL(2,L) = IHKL(1,JD)
          IHKL(3,L) = IHKL(3,JD)
          LD(3,L) = KF(NCO)+LD(3,JD)
          IF ( MOD(NCH(J),2).NE.0 ) LD(3,L)=KE(NCO)+LD(3,JD)
          JD = L
          GO TO 81
890       CONTINUE
900     ICHKL(IPHASE) = L-1
1000  IER = 0
      IF ( NC.NE.0 ) GO TO 1001
      DO 1003 M=2,ICHKL(IPHASE)
        IF ( IPH(M).GT.0 ) GO TO 1003
        IHKL(1,M) = -IHKL(1,M)
        IHKL(2,M) = -IHKL(2,M)
        IHKL(3,M) = -IHKL(3,M)
1003    CONTINUE
1001  RETURN
1002  IER = 1
      GO TO 1001
      END 
       SUBROUTINE LOOKUP(K,N,NSCAT,IXRAY,JOB)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      real lamda
      INTEGER PTR,TBXPTR
      CHARACTER*4 TABNC(128),TBXC(256)
      COMMON/TABCHR/TABNC,TBXC
      COMMON/TABLES/TBX(256,10),TBD(128,20),TABN(128),tbm(128)
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),XMAS(16)
      CHARACTER*4 NAM(16)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/COEFC/NAM
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCATX
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      IOF=0
      NS=0
      IF(K.GT.1) THEN
      DO 4320 IIPHAS=2,K
4320      IOF = IOF + NATOM(IIPHAS-1)
      END IF
      IF(K.GT.1)NS=NSAVE
      IF(JOB.EQ.2.OR.JOB.EQ.4)GOTO 70
      DO 10 I=1,N
      NS=MAX0(NSCAT,NS)
      DO 20 J=1,NS
        IF(NTYP(I+IOF).EQ.NAM(J))GOTO 30
20        CONTINUE
      DO 40 J=1,212
        IF(NTYP(I+IOF).EQ.TBXC(J))GOTO 50
40        CONTINUE
      WRITE(6,41)NTYP(I+IOF)
41      FORMAT(' SCATTERING COEFFICIENTS NOT FOUND FOR ',A4)
      STOP ' SCATTERING DATA MISSING' 
30      PTR(I+IOF)=J
      GOTO 10
50      NS=NS+1
      PTR(I+IOF)=NS
      DO 60 L=1,9 
60        AC(L,NS)=TBX(J,L)
      NAM(NS)=TBXC(J)
      TBXPTR=TBX(J,10)+0.5
      DFP(NS)=TBD(TBXPTR,IXRAY)
C      DFPP(NS)=TBD(TBXPTR,IXRAY+5) !To 10 wavelenghts (june 98)
      DFPP(NS)=TBD(TBXPTR,IXRAY+10)
C      XMAS(NS)=TBM(TBXPTR+1)       !To 10 wavelenghts (june 98)
      XMAS(NS)=TBM(TBXPTR)
C TO CHECK THE  DFP, DFPP AND MASS VALUES (drx), UNCOMMENT THE LINES BELOW
C      WRITE(*,9752)NAm(NS),DFP(NS),DFPP(NS), XMAS(NS),lamda(1)
C      WRITE (6,9752)NAm(NS),DFP(NS),DFPP(NS), XMAS(NS),lamda(1)
c9752  FORMAT(' ATOM =',a4,1X,'dfp,dfpp=',2(F10.5,1X),' MASS=',
C     *F8.3, '  Wlenght:',f8.6)
C      WRITE(*,9751) (AC(IJI,NS),IJI=1,9) 
C      WRITE(6,9751) (AC(IJI,NS),IJI=1,9) 
c9751  FORMAT(' Scat Coef.=', 9(F9.5,1X)) 
C END OF CHECKING THE VALUES (drx)
10      CONTINUE
      NSAVE=NS
      RETURN
70    DO 80 I=1,N
      NS=MAX0(NSCAT,NS)
      DO 90 J=1,NS
        IF(NTYP(I+IOF).EQ.NAM(J))GOTO 100
90        CONTINUE
      DO 110 J=1,85
        IF(NTYP(I+IOF).EQ.TABNC(J))GOTO 120
110       CONTINUE
      WRITE(6,111)NTYP(I+IOF)
111     FORMAT(34H SCATTERING LENGTHS NOT FOUND FOR ,A4)
      STOP ' SCATTERING DATA MISSING' 
100     PTR(I+IOF)=J
      GOTO 80
120     NS=NS+1
      PTR(I+IOF)=NS
      DFP(NS)=TABN(J)
      NAM(NS)=TABNC(J)
         if (j .ge. 61 .and. j .le. 81) goto 81
         if (j .eq. 82) goto 82
         if (j .ge. 83 .and. j .le. 85) goto 83
           goto 84
81          nsl = j+2
           goto 180
82           nsl = 90
           goto 180
83           nsl = j+9
           goto 180
84          nsl = j
180      XMAS(NS)=TBM(nsl)
80      CONTINUE
      NSAVE=NS
      RETURN
      END 
      SUBROUTINE CELL2(NPHASE,LAMDAM)
      include 'param.inc'
      REAL LAMDAM
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
     COMMON/CELLX/A,B,C,ALPHA,BETA,GAMMA,AL(3,3) 
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      EQUIVALENCE (AL(1,1),ASTAR),(AL(2,2),BSTAR),(AL(3,3),CSTAR),
     * (AL(1,2),BACOSC),(AL(1,3),CACOSB),(AL(2,3),BCCOSA)
35    FORMAT (22H0THE LAUE SYMMETRY IS ,A8)
      CHARACTER*8 LAU(14)
      DATA LAU(1),LAU(2),LAU(3),LAU(4),LAU(5),LAU(6),LAU(7),LAU(8),
     * LAU(9),LAU(10),LAU(11),LAU(12),LAU(13),LAU(14)
     *  /'1BAR','2/M','MMM','4/M','4/MMM','3BAR   R','3BAR M R','3BAR',
     *   '3BAR M 1','3BAR 1 M','6/M','6/MMM','M3','M3M'/
C      SIN1(A) = SIN(A*6.28318531)              ! not used
c      COS1(A) = COS(A*6.28318531)
      WRITE(6,35)LAU(NSPGRP)
      IF ( NAXIS.GT.3 ) STOP 5001 
      GO TO (1002,1071,1072,1073,1073,1074,1074,1075,1075,1075,1075,
     *  1075,1076,1076),NSPGRP
      CALL GOTOER
1071  IF ( NAXIS.NE.1 ) ALPHA=90.0
      IF ( NAXIS.NE.2 ) BETA=90.0 
      IF ( NAXIS.NE.3) GAMMA=90.0 
      GO TO 1002
1073  B = A
1072  ALPHA = 90.0
      BETA = 90.0
      GAMMA = 90.0
      GO TO 1002
1074  B = A
      C = A
      BETA = ALPHA
      GAMMA = ALPHA 
      GO TO 1002
1075  B = A
      ALPHA = 90.
      BETA = 90.0
      GAMMA = 120.0 
      GO TO 1002
1076  B = A
      C = A
      GO TO 1072
1002  COSA = COS(6.28318531*ALPHA/360.0)
105   COSB = COS(6.28318531*BETA/360.0) 
1052  COSC = COS(6.28318531*GAMMA/360.0)
1053  SINA=SQRT(1.0-COSA**2)
      SINB=SQRT(1.0-COSB**2)
      SINC=SQRT(1.0-COSC**2)
C
C-----THE VOLUME OF THE CELL IN THE DIRECT SPACE OF THE REAL CELL OF
C-----THE K-TH PHASE=NPHASE IS CALCULATED
C-----AND THE RAY GMAX OF THE SPHERE WITH THE VOLUME EQUIVALENT TO THE
C-----BRILLOUIN CELL 4/3*PI*GMAX**3=1/VOL
C
      VOLi(NPHASE) = A * B * C * SQRT(1.0 - COSA**2 - COSB**2 - COSC**2
     & + 2.0 * COSA * COSB * COSC )
C
      GCOM(NPHASE) = .877298169 * VOLi(NPHASE) / LAMDAM**3
C
      COSASR=(COSB*COSC-COSA)/(SINB*SINC)
      COSBSR=(COSA*COSC-COSB)/(SINA*SINC)
      COSCSR=(COSA*COSB-COSC)/(SINA*SINB)
      SINASR=SQRT(1.0-COSASR**2)
      SINBSR=SQRT(1.0-COSBSR**2)
      ASTAR=1.0/(A*SINC*SINBSR)
      BSTAR=1.0/(B *SINC*SINASR)
      CSTAR=1.0/(C*SINA*SINBSR)
      BCCOSA=BSTAR*CSTAR*COSASR*2.0
      CACOSB=CSTAR*ASTAR*COSBSR*2.0
      BACOSC=BSTAR*ASTAR*COSCSR*2.0
      ASTAR=ASTAR*ASTAR
      BSTAR=BSTAR*BSTAR
      CSTAR=CSTAR*CSTAR
      RETURN
      END 
      FUNCTION MULT(IH,IK,IL,KXIS)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      P = 2.0
      GO TO (220,100,110,120,130,150,180,140,170,160,170,190,200,210),
     * NSPGRP
      CALL GOTOER
100   CONTINUE
      IF ( IH ) 220, 220, 101
101   CONTINUE
      IF ( KXIS.EQ.3 ) THEN
      IF ( IL.GT.0 .AND. IH+IABS(IK).GT.0 ) P=2.0*P
      ELSE
      IF ( IK.GT.0 .AND. IH+IABS(IL).GT.0 ) P=2.0*P
      END IF
      GO TO 220
110   P = 1.0
      IF ( IH ) 112,112,111
111   P = 2.0*P
112   IF ( IK ) 114,114,113
113   P = 2.0*P
114   IF ( IL ) 220,220,115
115   P = 2.0*P
      GO TO 220
120   IF ( IK ) 220,220,121
121   P = 2.0*P
      IF ( IL ) 220,220,122
122   P = 2.0*P
      GO TO 220
130   IF ( IK ) 220,220,131
131   P = 2.0*P
      IF ( IH ) 134,134,132
132   IF ( IH-IK ) 133,134,134
133   P = 2.0*P
134   IF ( IL ) 220,220,135
135   P = 2.0*P
      GO TO 220
140   IF ( IH+IABS(IK) ) 142,142,141
141   P = 3.0*P
142   IF ( IL ) 220,220,143
143   P = 2.0*P
      GO TO 220
150   IF ( IH+IK+IL ) 151,152,151
151   P = 2.0*P
152   IF ( IH-IK ) 153,154,153
153   P = 3.0*P
      GO TO 220
154   IF ( IH-IL ) 153,220,153
160   IF ( IK ) 220,220,161
161   P = 3.0*P
      IF ( IK-IH ) 220,220,162
162   IF ( IH ) 220,220,135
170   IF ( IL ) 172,172,171
171   P = 2.0*P
172   IF ( IH ) 174,174,173
173   P = 3.0*P
      IF ( IK ) 220,220,135
174   IF ( IK ) 220,220,153
180   IF ( IH-IK ) 182,181,182
181   IF ( IH-IL ) 183,220,183
182   IF ( IK-IL ) 184,183,184
183   P = P/2.0
184   P = 3.0*P
      IF ( IH+IK+IL ) 135,220,135
190   IF ( IK ) 194,194,191
191   P = 6.0*P
      IF ( IH ) 192,194,192
192   IF ( IK-IH ) 193,194,193
193   P = 2.0*P
194   IF ( IL ) 135,220,135
200   IF ( IH-IK ) 203,201,203
201   IF ( IK-IL ) 203,202,203
202   P = 4.0*P
      GO TO 220
203   P = 3.0*P
      IF ( IK ) 205,205,204
204   P = 2.0*P
205   IF ( IH ) 220,220,135
210   IF ( IH ) 211,214,211
211   P = 4.0*P
      IF ( IH-IL ) 212,220,212
212   P = 3.0*P
      IF ( IH-IK ) 213,220,213
213   IF ( IK-IL ) 135,220,135
214   P = 3.0*P
      IF ( IK ) 215,220,215
215   P = 2.0*P
      IF ( IK-IL ) 135,220,135
C-----THIS LINE ADDED TO CORRECT THE PROBLEMS WITH THE -3, -3M1, 6/M,
C-----AND 6/MMM SPACE GROUPS THAT WERE OFF BY A FACTOR OF TWO
220   if (nspgrp.eq.8 .or. nspgrp.eq.9 .or. nspgrp.eq.11 .or.
     1 nspgrp.eq.12) p=0.5*p
      MULT = P
      RETURN
      END 
      BLOCK DATA FTABLE
      CHARACTER*4 TABNC(128),TBXC(256)
      COMMON/TABCHR/TABNC,TBXC
      COMMON/TABLES/TBX(256,10),TBD(128,20),TABN(128),tbm(128)
      DATA TABNC(1),TABN(1)/'H',-0.372/ 
      DATA TABNC(2),TABN(2)/'HE',0.3/
      DATA TABNC(3),TABN(3)/'LI',-0.194/
      DATA TABNC(4),TABN(4)/'BE',0.774/ 
      DATA TABNC(5),TABN(5)/'B',0.54/
      DATA TABNC(6),TABN(6)/'C',0.663/
      DATA TABNC(7),TABN(7)/'N',0.94/
      DATA TABNC(8),TABN(8)/'O',0.575/
      DATA TABNC(9),TABN(9)/'F',0.57/
      DATA TABNC(10),TABN(10)/'NE',0.46/
      DATA TABNC(11),TABN(11)/'NA',0.351/
      DATA TABNC(12),TABN(12)/'MG',0.52/
      DATA TABNC(13),TABN(13)/'AL',0.345/
      DATA TABNC(14),TABN(14)/'SI',0.415/
      DATA TABNC(15),TABN(15)/'P',0.51/ 
      DATA TABNC(16),TABN(16)/'S',0.285/
      DATA TABNC(17),TABN(17)/'CL',0.958/
      DATA TABNC(18),TABN(18)/'AR',0.20/
      DATA TABNC(19),TABN(19)/'K',0.37/ 
      DATA TABNC(20),TABN(20)/'CA',0.47/
      DATA TABNC(21),TABN(21)/'SC',1.18/
      DATA TABNC(22),TABN(22)/'TI',-0.33/
      DATA TABNC(23),TABN(23)/'V',-0.05/
      DATA TABNC(24),TABN(24)/'CR',0.352/
      DATA TABNC(25),TABN(25)/'MN',-0.36/
      DATA TABNC(26),TABN(26)/'FE',0.95/
      DATA TABNC(27),TABN(27)/'CO',0.25/
      DATA TABNC(28),TABN(28)/'NI',1.03/
      DATA TABNC(29),TABN(29)/'CU',0.76/
      DATA TABNC(30),TABN(30)/'ZN',0.57/
      DATA TABNC(31),TABN(31)/'GA',0.72/
      DATA TABNC(32),TABN(32)/'GE',0.84/
      DATA TABNC(33),TABN(33)/'AS',0.64/
      DATA TABNC(34),TABN(34)/'SE',0.80/
      DATA TABNC(35),TABN(35)/'BR',0.68/
      DATA TABNC(36),TABN(36)/'KR',0.74/
      DATA TABNC(37),TABN(37)/'RB',0.70/
      DATA TABNC(38),TABN(38)/'SR',0.69/
      DATA TABNC(39),TABN(39)/'Y',0.79/ 
      DATA TABNC(40),TABN(40)/'ZR',0.69/
      DATA TABNC(41),TABN(41)/'NB',0.71/
      DATA TABNC(42),TABN(42)/'MO',0.66/
      DATA TABNC(43),TABN(43)/'TC',0.68/
      DATA TABNC(44),TABN(44)/'RU',0.73/
      DATA TABNC(45),TABN(45)/'RH',0.59/
      DATA TABNC(46),TABN(46)/'PD',0.60/
      DATA TABNC(47),TABN(47)/'AG',0.61/
      DATA TABNC(48),TABN(48)/'CD',0.37/
      DATA TABNC(49),TABN(49)/'IN',0.39/
      DATA TABNC(50),TABN(50)/'SN',0.61/
      DATA TABNC(51),TABN(51)/'SB',0.56/
      DATA TABNC(52),TABN(52)/'TE',0.54/
      DATA TABNC(53),TABN(53)/'I',0.52/ 
      DATA TABNC(54),TABN(54)/'XE',0.47/
      DATA TABNC(55),TABN(55)/'CS',0.55/
      DATA TABNC(56),TABN(56)/'BA',0.52/
      DATA TABNC(57),TABN(57)/'LA',0.83/
      DATA TABNC(58),TABN(58)/'CE',0.55/
      DATA TABNC(59),TABN(59)/'PR',0.44/
      DATA TABNC(60),TABN(60)/'ND',0.72/
      DATA TABNC(61),TABN(61)/'EU',0.55/
      DATA TABNC(62),TABN(62)/'GD',1.5/ 
      DATA TABNC(63),TABN(63)/'TB',0.76/
      DATA TABNC(64),TABN(64)/'DY',1.69/
      DATA TABNC(65),TABN(65)/'HO',0.85/
      DATA TABNC(66),TABN(66)/'ER',0.79/
      DATA TABNC(67),TABN(67)/'TM',0.69/
      DATA TABNC(68),TABN(68)/'YB',1.26/
      DATA TABNC(69),TABN(69)/'LU',0.73/
      DATA TABNC(70),TABN(70)/'HF',0.78/
      DATA TABNC(71),TABN(71)/'TA',0.70/
      DATA TABNC(72),TABN(72)/'W',0.48/ 
      DATA TABNC(73),TABN(73)/'RE',0.92/
      DATA TABNC(74),TABN(74)/'OS',1.07/
      DATA TABNC(75),TABN(75)/'IR',1.06/
      DATA TABNC(76),TABN(76)/'PT',0.95/
      DATA TABNC(77),TABN(77)/'AU',0.76/
      DATA TABNC(78),TABN(78)/'HG',1.27/
      DATA TABNC(79),TABN(79)/'TL',0.89/
      DATA TABNC(80),TABN(80)/'PB',0.94/
      DATA TABNC(81),TABN(81)/'BI',0.85/
      DATA TABNC(82),TABN(82)/'TH',0.99/
      DATA TABNC(83),TABN(83)/'U',0.84/ 
      DATA TABNC(84),TABN(84)/'NP',1.055/
      DATA TABNC(85),TABN(85)/'PU',0.75/
C  Dispersion Coefficients from International Tables
C  Volume C page 219-222 for wavelengths
C    2.748510   2.289620   1.935970   1.788965   1.540520
C      Ti         Cr         Fe         Co         Cu
C    0.709260   0.559360   0.215947   0.209010   0.180195
C      Mo         Ag         Ta         W          Au
C     first  2 lines = 10 fp  values
C     second 2 lines = 10 fpp values
C  From the codes of Ian Madsen , added in DBWS on 11 Feb 98
C  Dispersion Coefficients for H
      DATA (TBD( 1,i),i=1,10)
     */    0.0000,    0.0000,    0.0000,    0.0000,    0.0000,
     *     0.0000,    0.0000,    0.0000,    0.0000,    0.0000/
      DATA (TBD( 1,i),i=11,20)
     */    0.0000,    0.0000,    0.0000,    0.0000,    0.0000,
     *     0.0000,    0.0000,    0.0000,    0.0000,    0.0000/
C  Dispersion Coefficients for He
      DATA (TBD( 2,i),i=1,10)
     */    0.0000,    0.0000,    0.0000,    0.0000,    0.0000,
     *     0.0000,    0.0000,    0.0000,    0.0000,    0.0000/
      DATA (TBD( 2,i),i=11,20)
     */    0.0000,    0.0000,    0.0000,    0.0000,    0.0000,
     *     0.0000,    0.0000,    0.0000,    0.0000,    0.0000/
C  Dispersion Coefficients for Li
      DATA (TBD( 3,i),i=1,10)
     */    0.0035,    0.0023,    0.0015,    0.0013,    0.0008,
     *    -0.0003,   -0.0004,   -0.0006,   -0.0006,   -0.0006/
      DATA (TBD( 3,i),i=11,20)
     */    0.0013,    0.0008,    0.0006,    0.0005,    0.0003,
     *     0.0001,    0.0000,    0.0000,    0.0000,    0.0000/
C  Dispersion Coefficients for Be
      DATA (TBD( 4,i),i=1,10)
     */    0.0117,    0.0083,    0.0060,    0.0052,    0.0038,
     *     0.0005,    0.0001,   -0.0005,   -0.0005,   -0.0005/
      DATA (TBD( 4,i),i=11,20)
     */    0.0050,    0.0033,    0.0023,    0.0019,    0.0014,
     *     0.0002,    0.0001,    0.0000,    0.0000,    0.0000/
C  Dispersion Coefficients for B
      DATA (TBD( 5,i),i=1,10)
     */    0.0263,    0.0190,    0.0140,    0.0121,    0.0090,
     *     0.0013,    0.0004,   -0.0009,   -0.0009,   -0.0010/
      DATA (TBD( 5,i),i=11,20)
     */    0.0139,    0.0094,    0.0065,    0.0055,    0.0039,
     *     0.0007,    0.0004,    0.0000,    0.0000,    0.0000/
C  Dispersion Coefficients for C
      DATA (TBD( 6,i),i=1,10)
     */    0.0490,    0.0364,    0.0273,    0.0237,    0.0181,
     *     0.0033,    0.0015,   -0.0012,   -0.0013,   -0.0014/
      DATA (TBD( 6,i),i=11,20)
     */    0.0313,    0.0213,    0.0148,    0.0125,    0.0091,
     *     0.0016,    0.0009,    0.0001,    0.0001,    0.0001/
C  Dispersion Coefficients for N
      DATA (TBD( 7,i),i=1,10)
     */    0.0807,    0.0606,    0.0461,    0.0403,    0.0311,
     *     0.0061,    0.0030,   -0.0020,   -0.0020,   -0.0023/
      DATA (TBD( 7,i),i=11,20)
     */    0.0606,    0.0416,    0.0293,    0.0248,    0.0180,
     *     0.0033,    0.0019,    0.0002,    0.0002,    0.0001/
C  Dispersion Coefficients for O
      DATA (TBD( 8,i),i=1,10)
     */    0.1213,    0.0928,    0.0716,    0.0630,    0.0492,
     *     0.0106,    0.0056,   -0.0025,   -0.0026,   -0.0030/
      DATA (TBD( 8,i),i=11,20)
     */    0.1057,    0.0731,    0.0518,    0.0440,    0.0322,
     *     0.0060,    0.0036,    0.0004,    0.0004,    0.0003/
C  Dispersion Coefficients for F
      DATA (TBD( 9,i),i=1,10)
     */    0.1700,    0.1324,    0.1037,    0.0920,    0.0727,
     *     0.0171,    0.0096,   -0.0027,   -0.0028,   -0.0034/
      DATA (TBD( 9,i),i=11,20)
     */    0.1710,    0.1192,    0.0851,    0.0725,    0.0534,
     *     0.0103,    0.0061,    0.0007,    0.0007,    0.0005/
C  Dispersion Coefficients for Ne
      DATA (TBD( 10,i),i=1,10)
     */    0.2257,    0.1793,    0.1426,    0.1273,    0.1019,
     *     0.0259,    0.0152,   -0.0025,   -0.0028,   -0.0037/
      DATA (TBD( 10,i),i=11,20)
     */    0.2621,    0.1837,    0.1318,    0.1126,    0.0833,
     *     0.0164,    0.0098,    0.0012,    0.0011,    0.0008/
C  Dispersion Coefficients for Na
      DATA (TBD(11,i),i=1,10)
     */    0.2801,    0.2295,    0.1857,    0.1670,    0.1353,
     *     0.0362,    0.0218,   -0.0028,   -0.0031,   -0.0044/
      DATA (TBD(11,i),i=11,20)
     */    0.3829,    0.2699,    0.1957,    0.1667,    0.1239,
     *     0.0249,    0.0150,    0.0019,    0.0017,    0.0012/
C  Dispersion Coefficients for Mg
      DATA (TBD(12,i),i=1,10)
     */    0.3299,    0.2778,    0.2309,    0.2094,    0.1719,
     *     0.0486,    0.0298,   -0.0030,   -0.0034,   -0.0052/
      DATA (TBD(12,i),i=11,20)
     */    0.5365,    0.3812,    0.2765,    0.2373,    0.1771,
     *     0.0363,    0.0220,    0.0028,    0.0026,    0.0018/
C  Dispersion Coefficients for Al
      DATA (TBD(13,i),i=1,10)
     */    0.3760,    0.3260,    0.2774,    0.2551,    0.2130,
     *     0.0645,    0.0406,   -0.0020,   -0.0026,   -0.0050/
      DATA (TBD(13,i),i=11,20)
     */    0.7287,    0.5212,    0.3807,    0.3276,    0.2455,
     *     0.0514,    0.0313,    0.0040,    0.0037,    0.0027/
C  Dispersion Coefficients for Si
      DATA (TBD(14,i),i=1,10)
     */    0.3921,    0.3647,    0.3209,    0.2979,    0.2541,
     *     0.0817,    0.0522,   -0.0017,   -0.0025,   -0.0055/
      DATA (TBD(14,i),i=11,20)
     */    0.9619,    0.6921,    0.5081,    0.4384,    0.3302,
     *     0.0704,    0.0431,    0.0056,    0.0052,    0.0038/
C  Dispersion Coefficients for P
      DATA (TBD(15,i),i=1,10)
     */    0.3821,    0.3898,    0.3592,    0.3388,    0.2955,
     *     0.1023,    0.0667,   -0.0002,   -0.0012,   -0.0050/
      DATA (TBD(15,i),i=11,20)
     */    1.2423,    0.8984,    0.6628,    0.5731,    0.4335,
     *     0.0942,    0.0580,    0.0077,    0.0071,    0.0052/
C  Dispersion Coefficients for S
      DATA (TBD(16,i),i=1,10)
     */    0.3167,    0.3899,    0.3848,    0.3706,    0.3331,
     *     0.1246,    0.0826,    0.0015,    0.0003,   -0.0045/
      DATA (TBD(16,i),i=11,20)
     */    1.5665,    1.1410,    0.8457,    0.7329,    0.5567,
     *     0.1234,    0.0763,    0.0103,    0.0096,    0.0069/
C  Dispersion Coefficients for Cl
      DATA (TBD(17,i),i=1,10)
     */    0.1832,    0.3508,    0.3920,    0.3892,    0.3639,
     *     0.1484,    0.0998,    0.0032,    0.0017,   -0.0042/
      DATA (TBD(17,i),i=11,20)
     */    1.9384,    1.4222,    1.0596,    0.9202,    0.7018,
     *     0.1585,    0.0984,    0.0134,    0.0125,    0.0091/
C  Dispersion Coefficients for Ar
      DATA (TBD(18,i),i=1,10)
     */   -0.0656,    0.2609,    0.3696,    0.3880,    0.3843,
     *     0.1743,    0.1191,    0.0059,    0.0041,   -0.0030/
      DATA (TBD(18,i),i=11,20)
     */    2.3670,    1.7458,    1.3087,    1.1388,    0.8717,
     *     0.2003,    0.1249,    0.0174,    0.0162,    0.0118/
C  Dispersion Coefficients for K
      DATA (TBD(19,i),i=1,10)
     */   -0.5083,    0.0914,    0.3068,    0.3532,    0.3868,
     *     0.2009,    0.1399,    0.0089,    0.0067,   -0.0017/
      DATA (TBD(19,i),i=11,20)
     */    2.8437,    2.1089,    1.5888,    1.3865,    1.0657,
     *     0.2494,    0.1562,    0.0219,    0.0204,    0.0149/
C  Dispersion Coefficients for Ca
      DATA (TBD(20,i),i=1,10)
     */   -1.3666,   -0.1987,    0.1867,    0.2782,    0.3641,
     *     0.2262,    0.1611,    0.0122,    0.0097,   -0.0002/
      DATA (TBD(20,i),i=11,20)
     */    3.3694,    2.5138,    1.9032,    0.6648,    1.2855,
     *     0.3064,    0.1926,    0.0273,    0.0255,    0.0187/
C  Dispersion Coefficients for Sc
      DATA (TBD(21,i),i=1,10)
     */   -5.4265,   -0.6935,   -0.0120,    0.1474,    0.3119,
     *     0.2519,    0.1829,    0.0159,    0.0130,    0.0015/
      DATA (TBD(21,i),i=11,20)
     */    4.0017,    2.9646,    2.2557,    1.9774,    1.5331,
     *     0.3716,    0.2348,    0.0338,    0.0315,    0.0231/
C  Dispersion Coefficients for Ti
      DATA (TBD(22,i),i=1,10)
     */   -2.2250,   -1.6394,   -0.3318,   -0.0617,    0.2191,
     *     0.2776,    0.2060,    0.0212,    0.0179,    0.0047/
      DATA (TBD(22,i),i=11,20)
     */    0.5264,    3.4538,    2.6425,    2.3213,    1.8069,
     *     0.4457,    0.2830,    0.0414,    0.0387,    0.0284/
C  Dispersion Coefficients for V
      DATA (TBD(23,i),i=1,10)
     */   -1.6269,   -4.4818,   -0.8645,   -0.3871,    0.0687,
     *     0.3005,    0.2276,    0.0259,    0.0221,    0.0070/
      DATA (TBD(23,i),i=11,20)
     */    0.6340,    0.4575,    3.0644,    2.6994,    2.1097,
     *     0.5294,    0.3376,    0.0500,    0.0468,    0.0344/
C  Dispersion Coefficients for Cr
      DATA (TBD(24,i),i=1,10)
     */   -1.2999,   -2.1308,   -1.9210,   -0.9524,   -0.1635,
     *     0.3209,    0.2496,    0.0314,    0.0272,    0.0101/
      DATA (TBD(24,i),i=11,20)
     */    0.7569,    0.5468,    3.5251,    3.1130,    2.4439,
     *     0.6236,    0.3992,    0.0599,    0.0561,    0.0413/
C  Dispersion Coefficients for Mn
      DATA (TBD(25,i),i=1,10)
     */   -1.0732,   -1.5980,   -3.5716,   -2.0793,   -0.5299,
     *     0.3368,    0.2704,    0.0377,    0.0330,    0.0139/
      DATA (TBD(25,i),i=11,20)
     */    0.8956,    0.6479,    0.4798,    3.5546,    2.8052,
     *     0.7283,    0.4681,    0.0712,    0.0666,    0.0492/
C  Dispersion Coefficients for Fe
      DATA (TBD(26,i),i=1,10)
     */   -0.8901,   -1.2935,   -2.0554,   -3.3307,   -1.1336,
     *     0.3463,    0.2886,    0.0438,    0.0386,    0.0173/
      DATA (TBD(26,i),i=11,20)
     */    1.0521,    0.7620,    0.5649,    0.4901,    3.1974,
     *     0.8444,    0.5448,    0.0840,    0.0787,    0.0582/
C  Dispersion Coefficients for Co
      DATA (TBD(27,i),i=1,10)
     */   -0.7307,   -1.0738,   -1.5743,   -2.0230,   -2.3653,
     *     0.3494,    0.3050,    0.0512,    0.0454,    0.0219/
      DATA (TBD(27,i),i=11,20)
     */    1.2272,    0.8897,    0.6602,    0.5731,    3.6143,
     *     0.9721,    0.6296,    0.0984,    0.0921,    0.0682/
C  Dispersion Coefficients for Ni
      DATA (TBD(28,i),i=1,10)
     */   -0.5921,   -0.9005,   -1.2894,   -1.5664,   -3.0029,
     *     0.3393,    0.3147,    0.0563,    0.0500,    0.0244/
      DATA (TBD(28,i),i=11,20)
     */    1.4240,    1.0331,    0.7671,    0.6662,    0.5091,
     *     1.1124,    0.7232,    0.1146,    0.1074,    0.0796/
C  Dispersion Coefficients for Cu
      DATA (TBD(29,i),i=1,10)
     */   -0.4430,   -0.7338,   -1.0699,   -1.2789,   -1.9646,
     *     0.3201,    0.3240,    0.0647,    0.0579,    0.0298/
      DATA (TBD(29,i),i=11,20)
     */    1.6427,    1.1930,    0.8864,    0.7700,    0.5888,
     *     1.2651,    0.8257,    0.1326,    0.1242,    0.0922/
C  Dispersion Coefficients for Zn
      DATA (TBD(30,i),i=1,10)
     */   -0.3524,   -0.6166,   -0.9134,   -1.0843,   -1.5491,
     *     0.2839,    0.3242,    0.0722,    0.0648,    0.0344/
      DATA (TBD(30,i),i=11,20)
     */    1.8861,    1.3712,    1.0193,    0.8857,    0.6778,
     *     1.4301,    0.9375,    0.1526,    0.1430,    0.1063/
C  Dispersion Coefficients for Ga
      DATA (TBD(31,i),i=1,10)
     */   -0.2524,   -0.4989,   -0.7701,   -0.9200,   -1.2846,
     *     0.2307,    0.3179,    0.0800,    0.0721,    0.0393/
      DATA (TBD(31,i),i=11,20)
     */    2.1518,    1.5674,    1.1663,    1.0138,    0.7763,
     *     1.6083,    1.0589,    0.1745,    0.1636,    0.1218/
C  Dispersion Coefficients for Ge
      DATA (TBD(32,i),i=1,10)
     */   -0.1549,   -0.3858,   -0.6412,   -0.7781,   -1.0885,
     *     0.1547,    0.3016,    0.0880,    0.0796,    0.0445/
      DATA (TBD(32,i),i=11,20)
     */    2.4445,    1.7841,    1.3291,    1.1557,    0.8855,
     *     1.8001,    1.1903,    0.1987,    0.1863,    0.1389/
C  Dispersion Coefficients for As
      DATA (TBD(33,i),i=1,10)
     */   -0.0687,   -0.2871,   -0.5260,   -0.6523,   -0.9300,
     *     0.0499,    0.2758,    0.0962,    0.0873,    0.0501/
      DATA (TBD(33,i),i=11,20)
     */    2.7627,    2.0194,    1.5069,    1.3109,    1.0051,
     *     2.0058,    1.3314,    0.2252,    0.2112,    0.1576/
C  Dispersion Coefficients for Se
      DATA (TBD(34,i),i=1,10)
     */    0.0052,   -0.1919,   -0.4179,   -0.5390,   -0.7943,
     *    -0.0929,    0.2367,    0.1047,    0.0954,    0.0560/
      DATA (TBD(34,i),i=11,20)
     */    3.1131,    2.2784,    1.7027,    1.4821,    1.1372,
     *     2.2259,    1.4831,    0.2543,    0.2386,    0.1782/
C  Dispersion Coefficients for Br
      DATA (TBD(35,i),i=1,10)
     */    0.0592,   -0.1095,   -0.3244,   -0.4363,   -0.6763,
     *    -0.2901,    0.1811,    0.1106,    0.1026,    0.0613/
      DATA (TBD(35,i),i=11,20)
     */    3.4901,    2.5578,    1.9140,    1.6673,    1.2805,
     *     2.4595,    1.6452,    0.2858,    0.2682,    0.2006/
C  Dispersion Coefficients for Kr
      DATA (TBD(36,i),i=1,10)
     */    0.1009,   -0.0316,   -0.2303,   -0.3390,   -0.5657,
     *    -0.5574,    0.1067,    0.1180,    0.1082,    0.0668/
      DATA (TBD(36,i),i=11,20)
     */    3.9083,    2.8669,    2.1472,    1.8713,    1.4385,
     *     2.7079,    1.8192,    0.3197,    0.3003,    0.2251/
C  Dispersion Coefficients for Rb
      DATA (TBD(37,i),i=1,10)
     */    0.1056,    0.0247,   -0.1516,   -0.2535,   -0.4688,
     *    -0.9393,    0.0068,    0.1247,    0.1146,    0.0717/
      DATA (TBD(37,i),i=11,20)
     */    4.3505,    3.1954,    2.3960,    2.0893,    1.6079,
     *     2.9676,    2.0025,    0.3561,    0.3346,    0.2514/
C  Dispersion Coefficients for Sr
      DATA (TBD(38,i),i=1,10)
     */    0.1220,    0.1037,   -0.0489,   -0.1448,   -0.3528,
     *    -1.5307,   -0.1172,    0.1321,    0.1219,    0.0769/
      DATA (TBD(38,i),i=11,20)
     */    4.8946,    3.6029,    2.7060,    2.3614,    1.8200,
     *     3.2498,    2.2025,    0.3964,    0.3726,    0.2805/
C  Dispersion Coefficients for Y
      DATA (TBD(39,i),i=1,10)
     */    0.0654,    0.1263,    0.0138,   -0.0720,   -0.2670,
     *    -2.7962,   -0.2879,    0.1380,    0.1278,    0.0819/
      DATA (TBD(39,i),i=11,20)
     */    5.4198,    3.9964,    3.0054,    2.6241,    2.0244,
     *     3.5667,    2.4099,    0.4390,    0.4128,    0.3112/
C  Dispersion Coefficients for Zr
      DATA (TBD(40,i),i=1,10)
     */   -0.0304,    0.1338,    0.0659,   -0.0066,   -0.1862,
     *    -2.9673,   -0.5364,    0.1431,    0.1329,    0.0863/
      DATA (TBD(40,i),i=11,20)
     */    5.9818,    4.4226,    3.3301,    2.9086,    2.2449,
     *     0.5597,    2.6141,    0.4852,    0.4562,    0.3443/
C  Dispersion Coefficients for Nb
      DATA (TBD(41,i),i=1,10)
     */   -0.1659,    0.1211,    0.1072,    0.0496,   -0.1121,
     *    -2.0727,   -0.8282,    0.1471,    0.1371,    0.0905/
      DATA (TBD(41,i),i=11,20)
     */    6.5803,    4.8761,    3.6768,    3.2133,    2.4826,
     *     0.6215,    2.8404,    0.5342,    0.5025,    0.3797/
C  Dispersion Coefficients for Mo
      DATA (TBD(42,i),i=1,10)
     */   -0.3487,    0.0801,    0.1301,    0.0904,   -0.0483,
     *    -1.6832,   -1.2703,    0.1487,    0.1391,    0.0934/
      DATA (TBD(42,i),i=11,20)
     */    7.2047,    5.3484,    4.0388,    3.5326,    2.7339,
     *     0.6857,    3.0978,    0.5862,    0.5517,    0.4177/
C  Dispersion Coefficients for Tc
      DATA (TBD(43,i),i=1,10)
     */   -0.6073,   -0.0025,    0.1314,    0.1164,    0.0057,
     *    -1.4390,   -2.0087,    0.1496,    0.1406,    0.0960/
      DATA (TBD(43,i),i=11,20)
     */    7.8739,    5.8597,    4.4331,    3.8799,    3.0049,
     *     0.7593,    3.3490,    0.6424,    0.6047,    0.4582/
C  Dispersion Coefficients for Ru
      DATA (TBD(44,i),i=1,10)
     */   -0.9294,   -0.1091,    0.1220,    0.1331,    0.0552,
     *    -1.2594,   -5.3630,    0.1491,    0.1409,    0.0981/
      DATA (TBD(44,i),i=11,20)
     */    8.5988,    6.4069,    4.8540,    4.2509,    3.2960,
     *     0.8363,    3.6506,    0.7016,    0.6607,    0.5014/
C  Dispersion Coefficients for Rh
      DATA (TBD(45,i),i=1,10)
     */   -1.3551,   -0.2630,    0.0861,    0.1305,    0.0927,
     *    -1.1178,   -2.5280,    0.1445,    0.1373,    0.0970/
      DATA (TBD(45,i),i=11,20)
     */    9.3504,    6.9820,    5.2985,    4.6432,    3.6045,
     *     0.9187,    0.5964,    0.7639,    0.7195,    0.5469/
C  Dispersion Coefficients for Pd
      DATA (TBD(46,i),i=1,10)
     */   -1.9086,   -0.4640,    0.0279,    0.1128,    0.1215,
     *    -0.9988,   -1.9556,    0.1387,    0.1327,    0.0959/
      DATA (TBD(46,i),i=11,20)
     */   10.1441,    7.5938,    5.7719,    5.0613,    3.9337,
     *     1.0072,    0.6546,    0.8302,    0.7822,    0.5955/
C  Dispersion Coefficients for Ag
      DATA (TBD(47,i),i=1,10)
     */   -2.5003,   -0.7387,   -0.0700,    0.0634,    0.1306,
     *    -0.8971,   -1.6473,    0.1295,    0.1251,    0.0928/
      DATA (TBD(47,i),i=11,20)
     */   10.9916,    8.2358,    6.2709,    5.5027,    4.2820,
     *     1.1015,    0.7167,    0.9001,    0.8484,    0.6469/
C  Dispersion Coefficients for Cd
      DATA (TBD(48,i),i=1,10)
     */   -3.5070,   -1.1086,   -0.2163,   -0.0214,    0.1185,
     *    -0.8075,   -1.4396,    0.1171,    0.1147,    0.0881/
      DATA (TBD(48,i),i=11,20)
     */   11.9019,    8.9174,    6.8017,    5.9728,    4.6533,
     *     1.2024,    0.7832,    0.9741,    0.9185,    0.7013/
C  Dispersion Coefficients for In
      DATA (TBD(49,i),i=1,10)
     */   -5.1325,   -1.5975,   -0.4165,   -0.1473,    0.0822,
     *    -0.7276,   -1.2843,    0.1013,    0.1012,    0.0816/
      DATA (TBD(49,i),i=11,20)
     */   12.6310,    9.6290,    7.3594,    6.4674,    5.0449,
     *     1.3100,    0.8542,    1.0519,    0.9922,    0.7587/
C  Dispersion Coefficients for Sn
      DATA (TBD(50,i),i=1,10)
     */   -7.5862,   -2.2019,   -0.6686,   -0.3097,    0.0259,
     *    -0.6537,   -1.1587,    0.0809,    0.0839,    0.0728/
      DATA (TBD(50,i),i=11,20)
     */   13.5168,   10.3742,    7.9473,    6.9896,    5.4591,
     *     1.4246,    0.9299,    1.1337,    1.0697,    0.8192/
C  Dispersion Coefficients for Sb
      DATA (TBD(51,i),i=1,10)
     */   -9.2145,   -3.0637,   -0.9868,   -0.5189,   -0.0562,
     *    -0.5866,   -1.0547,    0.0559,    0.0619,    0.0613/
      DATA (TBD(51,i),i=11,20)
     */   12.7661,   11.1026,    8.5620,    7.5367,    5.8946,
     *     1.5461,    1.0104,    1.2196,    1.1512,    0.8830/
C  Dispersion Coefficients for Te
      DATA (TBD(52,i),i=1,10)
     */  -11.6068,   -4.2407,   -1.4022,   -0.7914,   -0.1759,
     *    -0.5308,   -0.9710,    0.0216,    0.0316,    0.0435/
      DATA (TBD(52,i),i=11,20)
     */   10.1013,   11.8079,    9.2067,    8.1113,    6.3531,
     *     1.6751,    1.0960,    1.3095,    1.2366,    0.9499/
C  Dispersion Coefficients for I
      DATA (TBD(53,i),i=1,10)
     */  -13.9940,   -5.6353,   -1.9032,   -1.1275,   -0.3257,
     *    -0.4742,   -0.8919,   -0.0146,   -0.0001,    0.0259/
      DATA (TBD(53,i),i=11,20)
     */    3.4071,   12.6156,    9.8852,    8.7159,    6.8362,
     *     1.8119,    1.1868,    1.4037,    1.3259,    1.0201/
C  Dispersion Coefficients for Xe
      DATA (TBD(54,i),i=1,10)
     */   -9.6593,   -8.1899,   -2.6313,   -1.5532,   -0.5179,
     *    -0.4205,   -0.8200,   -0.0565,   -0.0367,    0.0057/
      DATA (TBD(54,i),i=11,20)
     */    3.7063,   11.7407,   10.5776,    9.3585,    7.3500,
     *     1.9578,    1.2838,    1.5023,    1.4195,    1.0938/
C  Dispersion Coefficients for Cs
      DATA (TBD(55,i),i=1,10)
     */   -8.1342,  -10.3310,   -3.5831,   -2.1433,   -0.7457,
     *    -0.3680,   -0.7527,   -0.1070,   -0.0809,   -0.0194/
      DATA (TBD(55,i),i=11,20)
     */    4.0732,   12.8551,   11.2902,   10.0454,    7.9052,
     *     2.1192,    1.3916,    1.6058,    1.5179,    1.1714/
C  Dispersion Coefficients for Ba
      DATA (TBD(56,i),i=1,10)
     */   -7.2079,  -11.0454,   -4.6472,   -2.7946,   -1.0456,
     *    -0.3244,   -0.6940,   -0.1670,   -0.1335,   -0.0494/
      DATA (TBD(56,i),i=11,20)
     */    4.4110,   10.0919,   12.0003,   10.7091,    8.4617,
     *     2.2819,    1.5004,    1.7127,    1.6194,    1.2517/
C  Dispersion Coefficients for La
      DATA (TBD(57,i),i=1,10)
     */   -6.5722,  -12.8190,   -6.3557,   -3.6566,   -1.4094,
     *    -0.2871,   -0.6411,   -0.2363,   -0.1940,   -0.0835/
      DATA (TBD(57,i),i=11,20)
     */    4.7587,    3.5648,   12.8927,   11.4336,    9.0376,
     *     2.4523,    1.6148,    1.8238,    1.7250,    1.3353/
C  Dispersion Coefficients for Ce
      DATA (TBD(58,i),i=1,10)
     */   -6.0641,   -9.3304,   -8.0962,   -4.8792,   -1.8482,
     *    -0.2486,   -0.5890,   -0.3159,   -0.2633,   -0.1222/
      DATA (TBD(58,i),i=11,20)
     */    5.1301,    3.8433,   11.8734,   12.1350,    9.6596,
     *     2.6331,    1.7358,    1.9398,    1.8353,    1.4227/
C  Dispersion Coefficients for Pr
      DATA (TBD(59,i),i=1,10)
     */   -5.6727,   -7.9841,  -10.9279,   -6.7923,   -2.4164,
     *    -0.2180,   -0.5424,   -0.4096,   -0.3443,   -0.1666/
      DATA (TBD(59,i),i=11,20)
     */    5.5091,    4.1304,    9.2394,   12.8653,   10.2820,
     *     2.8214,    1.8624,    2.0599,    1.9496,    1.5136/
C  Dispersion Coefficients for Nd
      DATA (TBD(60,i),i=1,10)
     */   -5.3510,   -7.1451,  -10.5249,   -8.1618,   -3.1807,
     *    -0.1943,   -0.5012,   -0.5194,   -0.4389,   -0.2183/
      DATA (TBD(60,i),i=11,20)
     */    5.9005,    4.4278,    9.9814,   11.9121,   10.9079,
     *     3.0179,    1.9950,    2.1843,    2.0679,    1.6077/
C  Dispersion Coefficients for Pm
      DATA (TBD(61,i),i=1,10)
     */   -5.0783,   -6.5334,  -13.2062,  -10.0720,   -4.0598,
     *    -0.1753,   -0.4626,   -0.6447,   -0.5499,   -0.2776/
      DATA (TBD(61,i),i=11,20)
     */    6.3144,    4.7422,    3.6278,    9.2324,   11.5523,
     *     3.2249,    2.1347,    2.3143,    2.1906,    1.7056/
C  Dispersion Coefficients for Sm
      DATA (TBD(62,i),i=1,10)
     */   -4.8443,   -6.0570,   -9.3497,  -10.2609,   -5.3236,
     *    -0.1638,   -0.4287,   -0.7989,   -0.6734,   -0.3455/
      DATA (TBD(62,i),i=11,20)
     */    6.7524,    5.0744,    3.8839,    9.9412,   12.2178,
     *     3.4418,    2.2815,    2.4510,    2.3197,    1.8069/
C  Dispersion Coefficients for Eu
      DATA (TBD(63,i),i=1,10)
     */   -4.6288,   -5.6630,   -7.9854,  -13.5405,   -8.9294,
     *    -0.1578,   -0.3977,   -0.9903,   -0.8137,   -0.4235/
      DATA (TBD(63,i),i=11,20)
     */    7.2035,    5.4178,    4.1498,    3.6550,   11.1857,
     *     3.6682,    2.4351,    2.5896,    2.4526,    1.9120/
C  Dispersion Coefficients for Gd
      DATA (TBD(64,i),i=1,10)
     */   -4.5094,   -5.3778,   -7.1681,   -9.3863,   -8.8380,
     *    -0.1653,   -0.3741,   -1.2279,   -1.0234,   -0.5140/
      DATA (TBD(64,i),i=11,20)
     */    7.6708,    5.7756,    4.4280,    3.9016,   11.9157,
     *     3.9035,    2.5954,    2.7304,    2.5878,    2.0202/
C  Dispersion Coefficients for Tb
      DATA (TBD(65,i),i=1,10)
     */   -4.3489,   -5.0951,   -6.5583,   -8.0413,   -9.1472,
     *    -0.1723,   -0.3496,   -1.5334,   -1.2583,   -0.6165/
      DATA (TBD(65,i),i=11,20)
     */    8.1882,    6.1667,    4.7292,    4.1674,    9.1891,
     *     4.1537,    2.7654,    2.8797,    2.7310,    2.1330/
C  Dispersion Coefficients for Dy
      DATA (TBD(66,i),i=1,10)
     */   -4.1616,   -4.8149,   -6.0597,   -7.1503,   -9.8046,
     *    -0.1892,   -0.3302,   -1.9594,   -1.5632,   -0.7322/
      DATA (TBD(66,i),i=11,20)
     */    8.6945,    6.5527,    5.0280,    4.4320,    9.8477,
     *     4.4098,    2.9404,    3.0274,    2.8733,    2.2494/
C  Dispersion Coefficients for Ho
      DATA (TBD(67,i),i=1,10)
     */   -4.0280,   -4.5887,   -5.6628,   -6.5338,  -14.9734,
     *    -0.2175,   -0.3168,   -2.6705,   -1.9886,   -0.8709/
      DATA (TBD(67,i),i=11,20)
     */    9.2302,    6.9619,    5.3451,    4.7129,    3.7046,
     *     4.6783,    3.1241,    3.1799,    3.0218,    2.3711/
C  Dispersion Coefficients for Er
      DATA (TBD(68,i),i=1,10)
     */   -3.9471,   -4.4106,   -5.3448,   -6.0673,   -9.4367,
     *    -0.2586,   -0.3091,   -5.5645,   -2.6932,   -1.0386/
      DATA (TBD(68,i),i=11,20)
     */    9.7921,    7.3910,    5.6776,    5.0074,    3.9380,
     *     4.9576,    3.3158,    0.6167,    3.1695,    2.4949/
C  Dispersion Coefficients for Tm
      DATA (TBD(69,i),i=1,10)
     */   -3.9079,   -4.2698,   -5.0823,   -5.6969,   -8.0393,
     *    -0.3139,   -0.3084,   -2.8957,   -5.6057,   -1.2397/
      DATA (TBD(69,i),i=11,20)
     */   10.3763,    7.8385,    6.0249,    5.3151,    4.1821,
     *     5.2483,    3.5155,    0.6569,    0.6192,    2.6240/
C  Dispersion Coefficients for Yb
      DATA (TBD(70,i),i=1,10)
     */   -3.8890,   -4.1523,   -4.8591,   -5.3940,   -7.2108,
     *    -0.3850,   -0.3157,   -2.4144,   -2.9190,   -1.4909/
      DATA (TBD(70,i),i=11,20)
     */   10.9742,    8.2969,    6.3813,    5.6309,    4.4329,
     *     5.5486,    3.7229,    0.6994,    0.6592,    2.7538/
C  Dispersion Coefficients for Lu
      DATA (TBD(71,i),i=1,10)
     */   -3.9056,   -4.0630,   -4.6707,   -5.1360,   -6.6179,
     *    -0.4720,   -0.3299,   -2.1535,   -2.4402,   -1.8184/
      DATA (TBD(71,i),i=11,20)
     */   11.5787,    8.7649,    6.7484,    5.9574,    4.6937,
     *     5.8584,    3.9377,    0.7436,    0.7010,    2.8890/
C  Dispersion Coefficients for Hf
      DATA (TBD(72,i),i=1,10)
     */   -4.0452,   -4.0564,   -4.4593,   -4.9466,   -6.1794,
     *    -0.5830,   -0.3548,   -1.9785,   -2.1778,   -2.2909/
      DATA (TBD(72,i),i=11,20)
     */   12.2546,    9.2832,    7.1518,    6.3150,    4.9776,
     *     6.1852,    4.1643,    0.7905,    0.7454,    3.0246/
C  Dispersion Coefficients for Ta
      DATA (TBD(73,i),i=1,10)
     */   -4.0905,   -3.9860,   -4.3912,   -4.7389,   -5.7959,
     *    -0.7052,   -0.3831,   -1.8534,   -2.0068,   -3.1639/
      DATA (TBD(73,i),i=11,20)
     */   12.9479,    9.8171,    7.5686,    6.6850,    5.2718,
     *     6.5227,    4.3992,    0.8392,    0.7915,    3.1610/
C  Dispersion Coefficients for W
      DATA (TBD(74,i),i=1,10)
     */   -4.1530,   -3.9270,   -4.2486,   -4.5529,   -5.4734,
     *    -0.8490,   -0.4201,   -1.7565,   -1.8819,   -3.8673/
      DATA (TBD(74,i),i=11,20)
     */   13.6643,   10.3696,    8.0005,    7.0688,    5.5774,
     *     6.8722,    4.6430,    0.8905,    0.8388,    0.6433/
C  Dispersion Coefficients for Re
      DATA (TBD(75,i),i=1,10)
     */   -4.2681,   -3.9052,   -4.1390,   -4.4020,   -5.2083,
     *    -1.0185,   -0.4693,   -1.6799,   -1.7868,   -2.8429/
      DATA (TBD(75,i),i=11,20)
     */   14.3931,   10.9346,    8.4435,    7.4631,    5.8923,
     *     7.2310,    4.8944,    0.9441,    0.8907,    0.6827/
C  Dispersion Coefficients for Os
      DATA (TBD(76,i),i=1,10)
     */   -4.4183,   -3.9016,   -4.0478,   -4.2711,   -4.9801,
     *    -1.2165,   -0.5280,   -1.6170,   -1.7107,   -2.4688/
      DATA (TBD(76,i),i=11,20)
     */   15.1553,   11.5251,    8.9067,    7.8753,    6.2216,
     *     7.6030,    5.1558,    1.0001,    0.9437,    0.7238/
C  Dispersion Coefficients for Ir
      DATA (TBD(77,i),i=1,10)
     */   -4.5860,   -3.9049,   -3.9606,   -4.1463,   -4.7710,
     *    -1.4442,   -0.5977,   -1.5648,   -1.6486,   -2.2499/
      DATA (TBD(77,i),i=11,20)
     */   15.9558,   12.1453,    9.3923,    8.3074,    6.5667,
     *     7.9887,    5.4269,    1.0589,    0.9993,    0.7669/
C  Dispersion Coefficients for Pt
      DATA (TBD(78,i),i=1,10)
     */   -4.8057,   -3.9435,   -3.8977,   -4.0461,   -4.5932,
     *    -1.7033,   -0.6812,   -1.5228,   -1.5998,   -2.1036/
      DATA (TBD(78,i),i=11,20)
     */   16.7870,   12.7910,    9.8985,    8.7578,    6.9264,
     *     8.3905,    5.7081,    1.1193,    1.0565,    0.8116/
C  Dispersion Coefficients for Au
      DATA (TBD(79,i),i=1,10)
     */   -5.0625,   -3.9908,   -3.8356,   -3.9461,   -4.4197,
     *    -2.0133,   -0.7638,   -1.4693,   -1.5404,   -1.9775/
      DATA (TBD(79,i),i=11,20)
     */   17.6400,   13.4551,   10.4202,    9.2222,    7.2980,
     *     8.8022,    5.9978,    1.1833,    1.1171,    0.8589/
C  Dispersion Coefficients for Hg
      DATA (TBD(80,i),i=1,10)
     */   -5.4327,   -4.1029,   -3.8228,   -3.8921,   -4.2923,
     *    -2.3894,   -0.8801,   -1.4389,   -1.5055,   -1.8958/
      DATA (TBD(80,i),i=11,20)
     */   18.5241,   14.1473,   10.9650,    9.7076,    7.6849,
     *     9.2266,    6.2989,    1.2483,    1.1796,    0.9080/
C  Dispersion Coefficients for Tl
      DATA (TBD(81,i),i=1,10)
     */   -5.8163,   -4.2233,   -3.8103,   -3.8340,   -4.1627,
     *    -2.8358,   -1.0117,   -1.4111,   -1.4740,   -1.8288/
      DATA (TBD(81,i),i=11,20)
     */   19.4378,   14.8643,   11.5300,   10.2108,    8.0900,
     *     9.6688,    6.6090,    1.3189,    1.2456,    0.9594/
C  Dispersion Coefficients for Pb
      DATA (TBD(82,i),i=1,10)
     */   -6.4779,   -4.4167,   -3.8519,   -3.8236,   -4.0753,
     *    -3.3944,   -1.1676,   -1.3897,   -1.4497,   -1.7773/
      DATA (TBD(82,i),i=11,20)
     */   20.3336,   15.5987,   12.1106,   10.7292,    8.5060,
     *    10.1111,    6.9287,    1.3909,    1.3137,    1.0127/
C  Dispersion Coefficients for Bi
      DATA (TBD(83,i),i=1,10)
     */   -7.0419,   -4.6533,   -3.9228,   -3.8408,   -4.0111,
     *    -4.1077,   -1.3494,   -1.3721,   -1.4290,   -1.7346/
      DATA (TBD(83,i),i=11,20)
     */   21.2196,   16.3448,   12.7017,   11.2575,    8.9310,
     *    10.2566,    7.2566,    1.4661,    1.3851,    1.0685/
C  Dispersion Coefficients for Po
      DATA (TBD(84,i),i=1,10)
     */   -7.7195,   -4.9604,   -4.0267,   -3.8855,   -3.9670,
     *    -5.1210,   -1.5613,   -1.3584,   -1.4133,   -1.7005/
      DATA (TBD(84,i),i=11,20)
     */   22.1974,   17.1410,   13.3329,   11.8209,    9.3834,
     *    11.0496,    7.5986,    1.5443,    1.4592,    1.1266/
C  Dispersion Coefficients for At
      DATA (TBD(85,i),i=1,10)
     */   -8.5994,   -5.3399,   -4.1781,   -3.9706,   -3.9588,
     *    -7.9122,   -1.8039,   -1.3540,   -1.4066,   -1.6784/
      DATA (TBD(85,i),i=11,20)
     */   23.2213,   17.9390,   13.9709,   12.3915,    9.8433,
     *     9.9777,    7.9509,    1.6260,    1.5367,    1.1876/
C  Dispersion Coefficients for Rn
      DATA (TBD(86,i),i=1,10)
     */  -10.2749,   -5.7275,   -4.3331,   -4.0549,   -3.9487,
     *    -8.0659,   -2.0847,   -1.3475,   -1.3982,   -1.6571/
      DATA (TBD(86,i),i=11,20)
     */   24.2613,   18.7720,   14.6313,   12.9815,   10.3181,
     *    10.4580,    8.3112,    1.7103,    1.6167,    1.2504/
C  Dispersion Coefficients for Fr
      DATA (TBD(87,i),i=1,10)
     */  -10.8938,   -6.2180,   -4.5387,   -4.1818,   -3.9689,
     *    -7.2224,   -2.4129,   -1.3404,   -1.3892,   -1.6367/
      DATA (TBD(87,i),i=11,20)
     */   24.3041,   19.6009,   15.3016,   13.5825,   10.8038,
     *     7.7847,    8.6839,    1.7986,    1.7004,    1.3162/
C  Dispersion Coefficients for Ra
      DATA (TBD(88,i),i=1,10)
     */  -12.3462,   -6.7502,   -4.7764,   -4.3309,   -4.0088,
     *    -6.7704,   -2.8081,   -1.3462,   -1.3931,   -1.6299/
      DATA (TBD(88,i),i=11,20)
     */   25.5374,   20.4389,   15.9778,   14.1902,   11.2969,
     *     8.1435,    9.0614,    1.8891,    1.7863,    1.3840/
C  Dispersion Coefficients for Ac
      DATA (TBD(89,i),i=1,10)
     */  -12.3496,   -7.4161,   -5.0617,   -4.5270,   -4.0794,
     *    -6.8494,   -3.2784,   -1.3473,   -1.3922,   -1.6190/
      DATA (TBD(89,i),i=11,20)
     */   25.1363,   21.3053,   16.6687,   14.8096,   11.7994,
     *     8.5178,    9.4502,    1.9845,    1.8770,    1.4553/
C  Dispersion Coefficients for Th
      DATA (TBD(90,i),i=1,10)
     */  -13.6049,   -8.2118,   -5.3692,   -4.7310,   -4.1491,
     *    -7.2400,   -3.8533,   -1.3524,   -1.3955,   -1.6136/
      DATA (TBD(90,i),i=11,20)
     */   26.2511,   22.2248,   17.4018,   15.4642,   12.3296,
     *     8.8979,    9.8403,    2.0819,    1.9695,    1.5284/
C  Dispersion Coefficients for Pa
      DATA (TBD(91,i),i=1,10)
     */  -14.4639,   -9.4459,   -5.7337,   -4.9639,   -4.2473,
     *    -8.0334,   -4.6067,   -1.3672,   -1.4083,   -1.6170/
      DATA (TBD(91,i),i=11,20)
     */   27.4475,   23.1548,   18.1406,   16.1295,   12.8681,
     *     9.2807,   10.2413,    2.1835,    2.0661,    1.6047/
C  Dispersion Coefficients for U
      DATA (TBD(92,i),i=1,10)
     */  -12.3528,   -9.9362,   -6.1485,   -5.2392,   -4.3638,
     *    -9.6767,   -5.7225,   -1.3792,   -1.4184,   -1.6188/
      DATA (TBD(92,i),i=11,20)
     */   30.1725,   23.1239,   18.8728,   16.7952,   13.4090,
     *     9.6646,   10.6428,    2.2876,    1.1650,    1.6831/
C  Dispersion Coefficients for Np
      DATA (TBD(93,i),i=1,10)
     */  -17.4143,  -11.1080,   -6.6136,   -5.5633,   -4.5053,
     *   -11.4937,   -6.9995,   -1.3941,   -1.4312,   -1.6231/
      DATA (TBD(93,i),i=11,20)
     */   31.7405,   24.1168,   19.6379,   17.4837,   13.9666,
     *     4.1493,    9.5876,    2.3958,    2.2679,    1.7648/
C  Dispersion Coefficients for Pu
      DATA (TBD(94,i),i=1,10)
     */  -18.0862,  -11.4073,   -6.9721,   -5.8130,   -4.6563,
     *    -9.4100,  -13.5905,   -1.4180,   -1.4527,   -1.6351/
      DATA (TBD(94,i),i=11,20)
     */   33.8963,   23.2960,   20.1548,   17.9579,   14.3729,
     *     4.3056,    6.9468,    2.4979,    2.3652,    1.8430/
C  Dispersion Coefficients for Am
      DATA (TBD(95,i),i=1,10)
     */  -19.7042,  -11.7097,   -7.7881,   -6.2920,   -4.8483,
     *    -7.8986,   -6.7022,   -1.4359,   -1.4684,   -1.6424/
      DATA (TBD(95,i),i=11,20)
     */   37.3716,   24.5715,   21.1738,   18.8618,   15.0877,
     *     4.5125,    7.3108,    2.6218,    2.4829,    1.9358/
C  Dispersion Coefficients for Cm
      DATA (TBD(96,i),i=1,10)
     */  -24.9307,  -10.4100,   -8.6102,   -6.7506,   -5.0611,
     *    -7.3248,   -6.2891,   -1.4655,   -1.4952,   -1.6592/
      DATA (TBD(96,i),i=11,20)
     */   41.4852,   25.8115,   21.8880,   19.5119,   15.6355,
     *     4.6980,    7.6044,    2.7421,    2.5974,    2.0271/
C  Dispersion Coefficients for Bk
      DATA (TBD(97,i),i=1,10)
     */  -32.8492,   -9.2185,   -9.3381,   -7.4293,   -5.3481,
     *    -6.8498,   -6.3438,   -1.4932,   -1.5203,   -1.6746/
      DATA (TBD(97,i),i=11,20)
     */   32.5421,   29.3028,   21.9514,   20.3581,   16.3190,
     *     4.9086,    7.9477,    2.8653,    2.7147,    2.1208/
C  Dispersion Coefficients for Cf
      DATA (TBD(98,i),i=1,10)
     */  -23.6520,  -23.5202,   -9.7799,   -7.8616,   -5.5545,
     *    -6.6561,   -6.4144,   -1.5323,   -1.5562,   -1.6984/
      DATA (TBD(98,i),i=11,20)
     */   21.9334,   31.2999,   22.4858,   20.8536,   16.7428,
     *     5.0785,    8.1930,    2.9807,    2.8250,    2.2102/
C Coefficients for analytical approximation to the scattering factors
C International Tables for Crystallography, Vol C, 501-503 [blue book]
      DATA TBXC( 1),TBX(  1,1),TBX(  1,2),TBX(  1,3),TBX(  1,4)
     *    /'H   '   ,.493002   ,10.5109   ,.322912   ,26.1257   /
      DATA TBX(  1,5),TBX(  1,6),TBX(  1,7),TBX(  1,8),TBX(  1,9)
     *    / .140191   ,3.14236   ,.040810   ,57.7997   ,.003038   /
      DATA TBX(  1,10)/1./
      DATA TBXC(  2),TBX(  2,1),TBX(  2,2),TBX(  2,3),TBX(  2,4)
     *    /'H   '   ,.489918   ,20.6593   ,.262003   ,7.74039   /
      DATA TBX(  2,5),TBX(  2,6),TBX(  2,7),TBX(  2,8),TBX(  2,9)
     *    / .196767   ,49.5519   ,.049879   ,2.20159   ,.001305   /
      DATA TBX(  2,10)/1./
      DATA TBXC(  3),TBX(  3,1),TBX(  3,2),TBX(  3,3),TBX(  3,4)
     *    /'H-1 '   ,.897661   ,53.1368   ,.565616   ,15.1870   /
      DATA TBX(  3,5),TBX(  3,6),TBX(  3,7),TBX(  3,8),TBX(  3,9)
     *    / .415815   ,186.576   ,.116973   ,3.56709   ,.002389   /
      DATA TBX(  3,10)/1./
      DATA TBXC(  4),TBX(  4,1),TBX(  4,2),TBX(  4,3),TBX(  4,4)
     *    /'HE  '   ,  .8734   , 9.1037   ,  .6309   , 3.3568   /
      DATA TBX(  4,5),TBX(  4,6),TBX(  4,7),TBX(  4,8),TBX(  4,9)
     *    /   .3112   ,22.9276   ,  .1780   ,  .9821   ,  .0064   /
      DATA TBX(  4,10)/2./
      DATA TBXC(  5),TBX(  5,1),TBX(  5,2),TBX(  5,3),TBX(  5,4)
     *    /'LI  '   , 1.1282   , 3.9546   ,  .7508   , 1.0524   /
      DATA TBX(  5,5),TBX(  5,6),TBX(  5,7),TBX(  5,8),TBX(  5,9)
     *    /   .6175   ,85.3905   ,  .4653   ,168.261   ,  .0377   /
      DATA TBX(  5,10)/3./
      DATA TBXC(  6),TBX(  6,1),TBX(  6,2),TBX(  6,3),TBX(  6,4)
     *    /'LI+1'   ,  .6968   , 4.6237   ,  .7888   , 1.9557   /
      DATA TBX(  6,5),TBX(  6,6),TBX(  6,7),TBX(  6,8),TBX(  6,9)
     *    /   .3414   ,  .6316   ,  .1563   ,10.0953   ,  .0167   /
      DATA TBX(  6,10)/3./
      DATA TBXC(  7),TBX(  7,1),TBX(  7,2),TBX(  7,3),TBX(  7,4)
     *    /'BE  '   , 1.5919   ,43.6427   , 1.1278   , 1.8623   /
      DATA TBX(  7,5),TBX(  7,6),TBX(  7,7),TBX(  7,8),TBX(  7,9)
     *    /   .5391   ,103.483   ,  .7029   ,  .5420   ,  .0385   /
      DATA TBX(  7,10)/4./
      DATA TBXC(  8),TBX(  8,1),TBX(  8,2),TBX(  8,3),TBX(  8,4)
     *    /'BE+2'   , 6.2603   ,  .0027   ,  .8849   ,  .8313   /
      DATA TBX(  8,5),TBX(  8,6),TBX(  8,7),TBX(  8,8),TBX(  8,9)
     *    /   .7993   , 2.2758   ,  .1647   , 5.1146   ,-6.1092   /
      DATA TBX(  8,10)/4./
      DATA TBXC(  9),TBX(  9,1),TBX(  9,2),TBX(  9,3),TBX(  9,4)
     *    /'B   '   , 2.0545   ,23.2185   , 1.3326   , 1.0210   /
      DATA TBX(  9,5),TBX(  9,6),TBX(  9,7),TBX(  9,8),TBX(  9,9)
     *    /  1.0979   ,60.3498   ,  .7068   ,  .1403   , -.1932   /
      DATA TBX(  9,10)/5./
      DATA TBXC( 10),TBX( 10,1),TBX( 10,2),TBX( 10,3),TBX( 10,4)
     *    /'C   '   , 2.3100   ,20.8439   , 1.0200   ,10.2075   /
      DATA TBX( 10,5),TBX( 10,6),TBX( 10,7),TBX( 10,8),TBX( 10,9)
     *    /  1.5886   ,  .5687   ,  .8650   ,51.6512   ,  .2156   /
      DATA TBX( 10,10)/6./
      DATA TBXC( 11),TBX( 11,1),TBX( 11,2),TBX( 11,3),TBX( 11,4)
     *    /'CVAL'   ,2.26069   ,22.6907   ,1.56165   ,.656665   /
      DATA TBX( 11,5),TBX( 11,6),TBX( 11,7),TBX( 11,8),TBX( 11,9)
     *    / 1.05075   ,9.75618   ,.839259   ,55.5949   ,.286977   /
      DATA TBX( 11,10)/6./
      DATA TBXC( 12),TBX( 12,1),TBX( 12,2),TBX( 12,3),TBX( 12,4)
     *    /'N   '   ,12.2126   ,  .0057   , 3.1322   , 9.8933   /
      DATA TBX( 12,5),TBX( 12,6),TBX( 12,7),TBX( 12,8),TBX( 12,9)
     *    /  2.0125   ,28.9975   , 1.1663   ,  .5826   ,-11.529   /
      DATA TBX( 12,10)/7./
      DATA TBXC( 13),TBX( 13,1),TBX( 13,2),TBX( 13,3),TBX( 13,4)
     *    /'O   '   , 3.0485   ,13.2771   , 2.2868   , 5.7011   /
      DATA TBX( 13,5),TBX( 13,6),TBX( 13,7),TBX( 13,8),TBX( 13,9)
     *    /  1.5463   ,  .3239   ,  .8670   ,32.9089   ,  .2508   /
      DATA TBX( 13,10)/8./
      DATA TBXC( 14),TBX( 14,1),TBX( 14,2),TBX( 14,3),TBX( 14,4)
     *    /'O-1 '   ,4.19160   ,12.8573   ,1.63969   ,4.17236   /
      DATA TBX( 14,5),TBX( 14,6),TBX( 14,7),TBX( 14,8),TBX( 14,9)
     *    / 1.52673   ,47.0179   ,-20.307   ,-.01404   ,21.9412   /
      DATA TBX( 14,10)/8./
      DATA TBXC( 15),TBX( 15,1),TBX( 15,2),TBX( 15,3),TBX( 15,4)
     *    /'F   '   , 3.5392   ,10.2825   , 2.6412   , 4.2944   /
      DATA TBX( 15,5),TBX( 15,6),TBX( 15,7),TBX( 15,8),TBX( 15,9)
     *    /  1.5170   ,  .2615   , 1.0243   ,26.1476   ,  .2776   /
      DATA TBX( 15,10)/9./
      DATA TBXC( 16),TBX( 16,1),TBX( 16,2),TBX( 16,3),TBX( 16,4)
     *    /'F-1 '   ,3.63220   ,5.27756   ,3.51057   ,14.7353   /
      DATA TBX( 16,5),TBX( 16,6),TBX( 16,7),TBX( 16,8),TBX( 16,9)
     *    / 1.26064   ,.442258   ,.940706   ,47.3437   ,.653396   /
      DATA TBX( 16,10)/9./
      DATA TBXC( 17),TBX( 17,1),TBX( 17,2),TBX( 17,3),TBX( 17,4)
     *    /'NE  '   , 3.9553   , 8.4042   , 3.1125   , 3.4262   /
      DATA TBX( 17,5),TBX( 17,6),TBX( 17,7),TBX( 17,8),TBX( 17,9)
     *    /  1.4546   ,  .2306   , 1.1251   ,21.7184   ,  .3515   /
      DATA TBX( 17,10)/10./
      DATA TBXC( 18),TBX( 18,1),TBX( 18,2),TBX( 18,3),TBX( 18,4)
     *    /'NA  '   , 4.7626   , 3.2850   , 3.1736   , 8.8422   /
      DATA TBX( 18,5),TBX( 18,6),TBX( 18,7),TBX( 18,8),TBX( 18,9)
     *    /  1.2674   ,  .3136   , 1.1128   ,129.424   ,  .6760   /
      DATA TBX( 18,10)/11./
      DATA TBXC( 19),TBX( 19,1),TBX( 19,2),TBX( 19,3),TBX( 19,4)
     *    /'NA+1'   , 3.2565   , 2.6671   , 3.9362   , 6.1153   /
      DATA TBX( 19,5),TBX( 19,6),TBX( 19,7),TBX( 19,8),TBX( 19,9)
     *    /  1.3998   ,  .2001   , 1.0032   ,14.0390   ,  .4040   /
      DATA TBX( 19,10)/11./
      DATA TBXC( 20),TBX( 20,1),TBX( 20,2),TBX( 20,3),TBX( 20,4)
     *    /'MG  '   , 5.4204   , 2.8275   , 2.1735   ,79.2611   /
      DATA TBX( 20,5),TBX( 20,6),TBX( 20,7),TBX( 20,8),TBX( 20,9)
     *    /  1.2269   ,  .3808   , 2.3073   , 7.1937   ,  .8584   /
      DATA TBX( 20,10)/12./
      DATA TBXC( 21),TBX( 21,1),TBX( 21,2),TBX( 21,3),TBX( 21,4)
     *    /'MG+2'   , 3.4988   , 2.1676   , 3.8378   , 4.7542   /
      DATA TBX( 21,5),TBX( 21,6),TBX( 21,7),TBX( 21,8),TBX( 21,9)
     *    /  1.3284   ,  .1850   ,  .8497   ,10.1411   ,  .4853   /
      DATA TBX( 21,10)/12./
      DATA TBXC( 22),TBX( 22,1),TBX( 22,2),TBX( 22,3),TBX( 22,4)
     *    /'AL  '   , 6.4202   , 3.0387   , 1.9002   ,  .7426   /
      DATA TBX( 22,5),TBX( 22,6),TBX( 22,7),TBX( 22,8),TBX( 22,9)
     *    /  1.5936   ,31.5472   , 1.9646   ,85.0886   , 1.1151   /
      DATA TBX( 22,10)/13./
      DATA TBXC( 23),TBX( 23,1),TBX( 23,2),TBX( 23,3),TBX( 23,4)
     *    /'AL+3'   ,4.17448   ,1.93816   ,3.38760   ,4.14553   /
      DATA TBX( 23,5),TBX( 23,6),TBX( 23,7),TBX( 23,8),TBX( 23,9)
     *    / 1.20296   ,.228753   ,.528137   ,8.28524   ,.706786   /
      DATA TBX( 23,10)/13./
      DATA TBXC( 24),TBX( 24,1),TBX( 24,2),TBX( 24,3),TBX( 24,4)
     *    /'SI  '   , 6.2915   , 2.4386   , 3.0353   ,32.3337   /
      DATA TBX( 24,5),TBX( 24,6),TBX( 24,7),TBX( 24,8),TBX( 24,9)
     *    /  1.9891   ,  .6785   , 1.5410   ,81.6937   , 1.1407   /
      DATA TBX( 24,10)/14./
      DATA TBXC( 25),TBX( 25,1),TBX( 25,2),TBX( 25,3),TBX( 25,4)
     *    /'SIV '   ,5.66269   ,2.66520   ,3.07164   ,38.6634   /
      DATA TBX( 25,5),TBX( 25,6),TBX( 25,7),TBX( 25,8),TBX( 25,9)
     *    / 2.62446   ,.916946   ,1.39320   ,93.5458   ,1.24707   /
      DATA TBX( 25,10)/14./
      DATA TBXC( 26),TBX( 26,1),TBX( 26,2),TBX( 26,3),TBX( 26,4)
     *    /'SI+4'   ,4.43918   ,1.64167   ,3.20345   ,3.43757   /
      DATA TBX( 26,5),TBX( 26,6),TBX( 26,7),TBX( 26,8),TBX( 26,9)
     *    / 1.19453   ,.214900   ,.416530   ,6.65365   ,.746297   /
      DATA TBX( 26,10)/14./
      DATA TBXC( 27),TBX( 27,1),TBX( 27,2),TBX( 27,3),TBX( 27,4)
     *    /'P   '   , 6.4345   , 1.9067   , 4.1791   ,27.1570   /
      DATA TBX( 27,5),TBX( 27,6),TBX( 27,7),TBX( 27,8),TBX( 27,9)
     *    /  1.7800   ,  .5260   , 1.4908   ,68.1645   , 1.1149   /
      DATA TBX( 27,10)/15./
      DATA TBXC( 28),TBX( 28,1),TBX( 28,2),TBX( 28,3),TBX( 28,4)
     *    /'S   '   , 6.9053   , 1.4679   , 5.2034   ,22.2151   /
      DATA TBX( 28,5),TBX( 28,6),TBX( 28,7),TBX( 28,8),TBX( 28,9)
     *    /  1.4379   ,  .2536   , 1.5863   ,56.1720   ,  .8669   /
      DATA TBX( 28,10)/16./
      DATA TBXC( 29),TBX( 29,1),TBX( 29,2),TBX( 29,3),TBX( 29,4)
     *    /'CL  '   ,11.4604   ,  .0104   , 7.1964   , 1.1662   /
      DATA TBX( 29,5),TBX( 29,6),TBX( 29,7),TBX( 29,8),TBX( 29,9)
     *    /  6.2556   ,18.5194   , 1.6455   ,47.7784   ,-9.5574   /
      DATA TBX( 29,10)/17./
      DATA TBXC( 30),TBX( 30,1),TBX( 30,2),TBX( 30,3),TBX( 30,4)
     *    /'CL-1'   ,18.2915   ,  .0066   , 7.2084   , 1.1717   /
      DATA TBX( 30,5),TBX( 30,6),TBX( 30,7),TBX( 30,8),TBX( 30,9)
     *    /  6.5337   ,19.5424   , 2.3386   ,60.4486   ,-16.378   /
      DATA TBX( 30,10)/17./
      DATA TBXC( 31),TBX( 31,1),TBX( 31,2),TBX( 31,3),TBX( 31,4)
     *    /'AR  '   , 7.4845   ,  .9072   , 6.7723   ,14.8407   /
      DATA TBX( 31,5),TBX( 31,6),TBX( 31,7),TBX( 31,8),TBX( 31,9)
     *    /   .6539   ,43.8983   , 1.6442   ,33.3929   , 1.4445   /
      DATA TBX( 31,10)/18./
      DATA TBXC( 32),TBX( 32,1),TBX( 32,2),TBX( 32,3),TBX( 32,4)
     *    /'K   '   , 8.2186   ,12.7949   , 7.4398   ,  .7748   /
      DATA TBX( 32,5),TBX( 32,6),TBX( 32,7),TBX( 32,8),TBX( 32,9)
     *    /  1.0519   ,213.187   ,  .8659   ,41.6841   , 1.4228   /
      DATA TBX( 32,10)/19./
      DATA TBXC( 33),TBX( 33,1),TBX( 33,2),TBX( 33,3),TBX( 33,4)
     *    /'K+1 '   , 7.9578   ,12.6331   , 7.4917   ,  .7674   /
      DATA TBX( 33,5),TBX( 33,6),TBX( 33,7),TBX( 33,8),TBX( 33,9)
     *    /  6.3590   , -.0020   , 1.1915   ,31.9128   ,-4.9978   /
      DATA TBX( 33,10)/19./
      DATA TBXC( 34),TBX( 34,1),TBX( 34,2),TBX( 34,3),TBX( 34,4)
     *    /'CA  '   , 8.6266   ,10.4421   , 7.3873   ,  .6599   /
      DATA TBX( 34,5),TBX( 34,6),TBX( 34,7),TBX( 34,8),TBX( 34,9)
     *    /  1.5899   ,85.7484   , 1.0211   ,178.437   , 1.3751   /
      DATA TBX( 34,10)/20./
      DATA TBXC( 35),TBX( 35,1),TBX( 35,2),TBX( 35,3),TBX( 35,4)
     *    /'CA+2'   ,15.6348   , -.0074   , 7.9518   ,  .6089   /
      DATA TBX( 35,5),TBX( 35,6),TBX( 35,7),TBX( 35,8),TBX( 35,9)
     *    /  8.4372   ,10.3116   ,  .8537   ,25.9905   ,-14.875   /
      DATA TBX( 35,10)/20./
      DATA TBXC( 36),TBX( 36,1),TBX( 36,2),TBX( 36,3),TBX( 36,4)
     *    /'SC  '   , 9.1890   , 9.0213   , 7.3679   ,  .5729   /
      DATA TBX( 36,5),TBX( 36,6),TBX( 36,7),TBX( 36,8),TBX( 36,9)
     *    /  1.6409   ,136.108   , 1.4680   ,51.3531   , 1.3329   /
      DATA TBX( 36,10)/21./
      DATA TBXC( 37),TBX( 37,1),TBX( 37,2),TBX( 37,3),TBX( 37,4)
     *    /'SC+3'   ,13.4008   ,.298540   ,8.02730   ,7.96290   /
      DATA TBX( 37,5),TBX( 37,6),TBX( 37,7),TBX( 37,8),TBX( 37,9)
     *    / 1.65943   ,-.28604   ,1.57936   ,16.0662   ,-6.6667   /
      DATA TBX( 37,10)/21./
      DATA TBXC( 38),TBX( 38,1),TBX( 38,2),TBX( 38,3),TBX( 38,4)
     *    /'TI  '   , 9.7595   , 7.8508   , 7.3558   ,  .5000   /
      DATA TBX( 38,5),TBX( 38,6),TBX( 38,7),TBX( 38,8),TBX( 38,9)
     *    /  1.6991   ,35.6338   , 1.9021   ,116.105   , 1.2807   /
      DATA TBX( 38,10)/22./
      DATA TBXC( 39),TBX( 39,1),TBX( 39,2),TBX( 39,3),TBX( 39,4)
     *    /'TI+2'   ,9.11423   ,7.52430   ,7.62174   ,.457585   /
      DATA TBX( 39,5),TBX( 39,6),TBX( 39,7),TBX( 39,8),TBX( 39,9)
     *    / 2.27930   ,19.5361   ,.087899   ,61.6558   ,.897155   /
      DATA TBX( 39,10)/22./
      DATA TBXC( 40),TBX( 40,1),TBX( 40,2),TBX( 40,3),TBX( 40,4)
     *    /'TI+3'   ,17.7344   ,.220610   ,8.73816   ,7.04716   /
      DATA TBX( 40,5),TBX( 40,6),TBX( 40,7),TBX( 40,8),TBX( 40,9)
     *    / 5.25691   ,-.15762   ,1.92134   ,15.9768   ,-14.652   /
      DATA TBX( 40,10)/22./
      DATA TBXC( 41),TBX( 41,1),TBX( 41,2),TBX( 41,3),TBX( 41,4)
     *    /'TI+4'   ,19.5114   ,.178847   ,8.23473   ,6.67018   /
      DATA TBX( 41,5),TBX( 41,6),TBX( 41,7),TBX( 41,8),TBX( 41,9)
     *    / 2.01341   ,-.29263   ,1.52080   ,12.9464   ,-13.280   /
      DATA TBX( 41,10)/22./
      DATA TBXC( 42),TBX( 42,1),TBX( 42,2),TBX( 42,3),TBX( 42,4)
     *    /'V   '   ,10.2971   , 6.8657   , 7.3511   ,  .4385   /
      DATA TBX( 42,5),TBX( 42,6),TBX( 42,7),TBX( 42,8),TBX( 42,9)
     *    /  2.0703   ,26.8938   , 2.0571   ,102.478   , 1.2199   /
      DATA TBX( 42,10)/23./
      DATA TBXC( 43),TBX( 43,1),TBX( 43,2),TBX( 43,3),TBX( 43,4)
     *    /'V+2 '   ,10.1060   , 6.8818   , 7.3541   ,  .4409   /
      DATA TBX( 43,5),TBX( 43,6),TBX( 43,7),TBX( 43,8),TBX( 43,9)
     *    /  2.2884   ,20.3004   ,  .0223   ,115.122   , 1.2298   /
      DATA TBX( 43,10)/23./
      DATA TBXC( 44),TBX( 44,1),TBX( 44,2),TBX( 44,3),TBX( 44,4)
     *    /'V+3 '   ,9.43141   ,6.39535   ,7.74190   ,.383349   /
      DATA TBX( 44,5),TBX( 44,6),TBX( 44,7),TBX( 44,8),TBX( 44,9)
     *    / 2.15343   ,15.1908   ,.016865   ,63.9690   ,.656565   /
      DATA TBX( 44,10)/23./
      DATA TBXC( 45),TBX( 45,1),TBX( 45,2),TBX( 45,3),TBX( 45,4)
     *    /'V+5 '   ,15.6887   ,.679003   ,8.14208   ,5.40135   /
      DATA TBX( 45,5),TBX( 45,6),TBX( 45,7),TBX( 45,8),TBX( 45,9)
     *    / 2.03081   ,9.97278   ,-9.5760   ,.940464   ,1.71430   /
      DATA TBX( 45,10)/23./
      DATA TBXC( 46),TBX( 46,1),TBX( 46,2),TBX( 46,3),TBX( 46,4)
     *    /'CR  '   ,10.6406   , 6.1038   , 7.3537   ,  .3920   /
      DATA TBX( 46,5),TBX( 46,6),TBX( 46,7),TBX( 46,8),TBX( 46,9)
     *    /  3.3240   ,20.2626   , 1.4922   ,98.7399   , 1.1832   /
      DATA TBX( 46,10)/24./
      DATA TBXC( 47),TBX( 47,1),TBX( 47,2),TBX( 47,3),TBX( 47,4)
     *    /'CR+2'   ,9.54034   ,5.66078   ,7.75090   ,.344261   /
      DATA TBX( 47,5),TBX( 47,6),TBX( 47,7),TBX( 47,8),TBX( 47,9)
     *    / 3.58274   ,13.3075   ,.509107   ,32.4224   ,.616898   /
      DATA TBX( 47,10)/24./
      DATA TBXC( 48),TBX( 48,1),TBX( 48,2),TBX( 48,3),TBX( 48,4)
     *    /'CR+3'   ,9.68090   ,5.59463   ,7.81136   ,.334393   /
      DATA TBX( 48,5),TBX( 48,6),TBX( 48,7),TBX( 48,8),TBX( 48,9)
     *    / 2.87603   ,12.8288   ,.113575   ,32.8761   ,.518275   /
      DATA TBX( 48,10)/24./
      DATA TBXC( 49),TBX( 49,1),TBX( 49,2),TBX( 49,3),TBX( 49,4)
     *    /'MN  '   ,11.2819   , 5.3409   , 7.3573   ,  .3432   /
      DATA TBX( 49,5),TBX( 49,6),TBX( 49,7),TBX( 49,8),TBX( 49,9)
     *    /  3.0193   ,17.8674   , 2.2441   ,83.7543   , 1.0896   /
      DATA TBX( 49,10)/25./
      DATA TBXC( 50),TBX( 50,1),TBX( 50,2),TBX( 50,3),TBX( 50,4)
     *    /'MN+2'   ,10.8061   , 5.2796   , 7.3620   ,  .3435   /
      DATA TBX( 50,5),TBX( 50,6),TBX( 50,7),TBX( 50,8),TBX( 50,9)
     *    /  3.5268   ,14.3430   ,  .2184   ,41.3235   , 1.0874   /
      DATA TBX( 50,10)/25./
      DATA TBXC( 51),TBX( 51,1),TBX( 51,2),TBX( 51,3),TBX( 51,4)
     *    /'MN+3'   ,9.84521   ,4.91797   ,7.87194   ,.294393   /
      DATA TBX( 51,5),TBX( 51,6),TBX( 51,7),TBX( 51,8),TBX( 51,9)
     *    / 3.56531   ,10.8171   ,.323613   ,24.1281   ,.393974   /
      DATA TBX( 51,10)/25./
      DATA TBXC( 52),TBX( 52,1),TBX( 52,2),TBX( 52,3),TBX( 52,4)
     *    /'MN+4'   ,9.96253   ,4.84850   ,7.97057   ,.283303   /
      DATA TBX( 52,5),TBX( 52,6),TBX( 52,7),TBX( 52,8),TBX( 52,9)
     *    / 2.76067   ,10.4852   ,.054447   ,27.5730   ,.251877   /
      DATA TBX( 52,10)/25./
      DATA TBXC( 53),TBX( 53,1),TBX( 53,2),TBX( 53,3),TBX( 53,4)
     *    /'FE  '   ,11.7695   , 4.7611   , 7.3573   ,  .3072   /
      DATA TBX( 53,5),TBX( 53,6),TBX( 53,7),TBX( 53,8),TBX( 53,9)
     *    /  3.5222   ,15.3535   , 2.3045   ,76.8805   , 1.0369   /
      DATA TBX( 53,10)/26./
      DATA TBXC( 54),TBX( 54,1),TBX( 54,2),TBX( 54,3),TBX( 54,4)
     *    /'FE+2'   ,11.0424   , 4.6538   , 7.3740   ,  .3053   /
      DATA TBX( 54,5),TBX( 54,6),TBX( 54,7),TBX( 54,8),TBX( 54,9)
     *    /  4.1346   ,12.0546   ,  .4399   ,31.2809   , 1.0097   /
      DATA TBX( 54,10)/26./
      DATA TBXC( 55),TBX( 55,1),TBX( 55,2),TBX( 55,3),TBX( 55,4)
     *    /'FE+3'   ,11.1764   , 4.6147   , 7.3863   ,  .3005   /
      DATA TBX( 55,5),TBX( 55,6),TBX( 55,7),TBX( 55,8),TBX( 55,9)
     *    /  3.3948   ,11.6729   ,  .0724   ,38.5566   ,  .9707   /
      DATA TBX( 55,10)/26./
      DATA TBXC( 56),TBX( 56,1),TBX( 56,2),TBX( 56,3),TBX( 56,4)
     *    /'CO  '   ,12.2841   , 4.2791   , 7.3409   ,  .2784   /
      DATA TBX( 56,5),TBX( 56,6),TBX( 56,7),TBX( 56,8),TBX( 56,9)
     *    /  4.0034   ,13.5359   , 2.3488   ,71.1692   , 1.0118   /
      DATA TBX( 56,10)/27./
      DATA TBXC( 57),TBX( 57,1),TBX( 57,2),TBX( 57,3),TBX( 57,4)
     *    /'CO+2'   ,11.2296   , 4.1231   , 7.3883   ,  .2726   /
      DATA TBX( 57,5),TBX( 57,6),TBX( 57,7),TBX( 57,8),TBX( 57,9)
     *    /  4.7393   ,10.2443   ,  .7108   ,25.6466   ,  .9324   /
      DATA TBX( 57,10)/27./
      DATA TBXC( 58),TBX( 58,1),TBX( 58,2),TBX( 58,3),TBX( 58,4)
     *    /'CO+3'   ,10.3380   ,3.90969   ,7.88173   ,.238668   /
      DATA TBX( 58,5),TBX( 58,6),TBX( 58,7),TBX( 58,8),TBX( 58,9)
     *    / 4.76795   ,8.35583   ,.725591   ,18.3491   ,.286667   /
      DATA TBX( 58,10)/27./
      DATA TBXC( 59),TBX( 59,1),TBX( 59,2),TBX( 59,3),TBX( 59,4)
     *    /'NI  '   ,12.8376   , 3.8785   , 7.2920   ,  .2565   /
      DATA TBX( 59,5),TBX( 59,6),TBX( 59,7),TBX( 59,8),TBX( 59,9)
     *    /  4.4438   ,12.1763   , 2.3800   ,66.3421   , 1.0341   /
      DATA TBX( 59,10)/28./
      DATA TBXC( 60),TBX( 60,1),TBX( 60,2),TBX( 60,3),TBX( 60,4)
     *    /'NI+2'   ,11.4166   , 3.6766   , 7.4005   ,  .2449   /
      DATA TBX( 60,5),TBX( 60,6),TBX( 60,7),TBX( 60,8),TBX( 60,9)
     *    /  5.3442   , 8.8730   ,  .9773   ,22.1626   ,  .8614   /
      DATA TBX( 60,10)/28./
      DATA TBXC( 61),TBX( 61,1),TBX( 61,2),TBX( 61,3),TBX( 61,4)
     *    /'NI+3'   ,10.7806   ,3.54770   ,7.75868   ,.223140   /
      DATA TBX( 61,5),TBX( 61,6),TBX( 61,7),TBX( 61,8),TBX( 61,9)
     *    / 5.22746   ,7.64468   ,.847114   ,16.9673   ,.386044   /
      DATA TBX( 61,10)/28./
      DATA TBXC( 62),TBX( 62,1),TBX( 62,2),TBX( 62,3),TBX( 62,4)
     *    /'CU  '   ,13.3380   , 3.5828   , 7.1676   ,  .2470   /
      DATA TBX( 62,5),TBX( 62,6),TBX( 62,7),TBX( 62,8),TBX( 62,9)
     *    /  5.6158   ,11.3966   , 1.6735   ,64.8126   , 1.1910   /
      DATA TBX( 62,10)/29./
      DATA TBXC( 63),TBX( 63,1),TBX( 63,2),TBX( 63,3),TBX( 63,4)
     *    /'CU+1'   ,11.9475   , 3.3669   , 7.3573   ,  .2274   /
      DATA TBX( 63,5),TBX( 63,6),TBX( 63,7),TBX( 63,8),TBX( 63,9)
     *    /  6.2455   , 8.6625   , 1.5578   ,25.8487   ,  .8900   /
      DATA TBX( 63,10)/29./
      DATA TBXC( 64),TBX( 64,1),TBX( 64,2),TBX( 64,3),TBX( 64,4)
     *    /'CU+2'   ,11.8168   ,3.37484   ,7.11181   ,.244078   /
      DATA TBX( 64,5),TBX( 64,6),TBX( 64,7),TBX( 64,8),TBX( 64,9)
     *    / 5.78135   ,7.98760   ,1.14523   ,19.8970   ,1.14431   /
      DATA TBX( 64,10)/29./
      DATA TBXC( 65),TBX( 65,1),TBX( 65,2),TBX( 65,3),TBX( 65,4)
     *    /'ZN  '   ,14.0743   , 3.2655   , 7.0318   ,  .2333   /
      DATA TBX( 65,5),TBX( 65,6),TBX( 65,7),TBX( 65,8),TBX( 65,9)
     *    /  5.1652   ,10.3163   , 2.4100   ,58.7097   , 1.3041   /
      DATA TBX( 65,10)/30./
      DATA TBXC( 66),TBX( 66,1),TBX( 66,2),TBX( 66,3),TBX( 66,4)
     *    /'ZN+2'   ,11.9719   , 2.9946   , 7.3862   ,  .2031   /
      DATA TBX( 66,5),TBX( 66,6),TBX( 66,7),TBX( 66,8),TBX( 66,9)
     *    /  6.4668   , 7.0826   , 1.3940   ,18.0995   ,  .7807   /
      DATA TBX( 66,10)/30./
      DATA TBXC( 67),TBX( 67,1),TBX( 67,2),TBX( 67,3),TBX( 67,4)
     *    /'GA  '   ,15.2354   , 3.0669   , 6.7006   ,  .2412   /
      DATA TBX( 67,5),TBX( 67,6),TBX( 67,7),TBX( 67,8),TBX( 67,9)
     *    /  4.3591   ,10.7805   , 2.9623   ,61.4135   , 1.7189   /
      DATA TBX( 67,10)/31./
      DATA TBXC( 68),TBX( 68,1),TBX( 68,2),TBX( 68,3),TBX( 68,4)
     *    /'GA+3'   ,12.6920   ,2.81262   ,6.69883   ,.227890   /
      DATA TBX( 68,5),TBX( 68,6),TBX( 68,7),TBX( 68,8),TBX( 68,9)
     *    / 6.06692   ,6.36441   ,1.00660   ,14.4122   ,1.53545   /
      DATA TBX( 68,10)/31./
      DATA TBXC( 69),TBX( 69,1),TBX( 69,2),TBX( 69,3),TBX( 69,4)
     *    /'GE  '   ,16.0816   , 2.8509   , 6.3747   ,  .2516   /
      DATA TBX( 69,5),TBX( 69,6),TBX( 69,7),TBX( 69,8),TBX( 69,9)
     *    /  3.7068   ,11.4468   , 3.6830   ,54.7625   , 2.1313   /
      DATA TBX( 69,10)/32./
      DATA TBXC( 70),TBX( 70,1),TBX( 70,2),TBX( 70,3),TBX( 70,4)
     *    /'GE+4'   ,12.9172   ,2.53718   ,6.70003   ,.205855   /
      DATA TBX( 70,5),TBX( 70,6),TBX( 70,7),TBX( 70,8),TBX( 70,9)
     *    / 6.06791   ,5.47913   ,.859041   ,11.6030   ,1.45572   /
      DATA TBX( 70,10)/32./
      DATA TBXC( 71),TBX( 71,1),TBX( 71,2),TBX( 71,3),TBX( 71,4)
     *    /'AS  '   ,16.6723   , 2.6345   , 6.0701   ,  .2647   /
      DATA TBX( 71,5),TBX( 71,6),TBX( 71,7),TBX( 71,8),TBX( 71,9)
     *    /  3.4313   ,12.9479   , 4.2779   ,47.7972   , 2.5310   /
      DATA TBX( 71,10)/33./
      DATA TBXC( 72),TBX( 72,1),TBX( 72,2),TBX( 72,3),TBX( 72,4)
     *    /'SE  '   ,17.0006   , 2.4098   , 5.8196   ,  .2726   /
      DATA TBX( 72,5),TBX( 72,6),TBX( 72,7),TBX( 72,8),TBX( 72,9)
     *    /  3.9731   ,15.2372   , 4.3543   ,43.8163   , 2.8409   /
      DATA TBX( 72,10)/34./
      DATA TBXC( 73),TBX( 73,1),TBX( 73,2),TBX( 73,3),TBX( 73,4)
     *    /'BR  '   ,17.1789   , 2.1723   , 5.2358   ,16.5796   /
      DATA TBX( 73,5),TBX( 73,6),TBX( 73,7),TBX( 73,8),TBX( 73,9)
     *    /  5.6377   ,  .2609   , 3.9851   ,41.4328   , 2.9557   /
      DATA TBX( 73,10)/35./
      DATA TBXC( 74),TBX( 74,1),TBX( 74,2),TBX( 74,3),TBX( 74,4)
     *    /'BR-1'   ,17.1718   , 2.2059   , 6.3338   ,19.3345   /
      DATA TBX( 74,5),TBX( 74,6),TBX( 74,7),TBX( 74,8),TBX( 74,9)
     *    /  5.5754   ,  .2871   , 3.7272   ,58.1535   , 3.1776   /
      DATA TBX( 74,10)/35./
      DATA TBXC( 75),TBX( 75,1),TBX( 75,2),TBX( 75,3),TBX( 75,4)
     *    /'KR  '   ,17.3555   , 1.9384   , 6.7286   ,16.5623   /
      DATA TBX( 75,5),TBX( 75,6),TBX( 75,7),TBX( 75,8),TBX( 75,9)
     *    /  5.5493   ,  .2261   , 3.5375   ,39.3972   , 2.8250   /
      DATA TBX( 75,10)/36./
      DATA TBXC( 76),TBX( 76,1),TBX( 76,2),TBX( 76,3),TBX( 76,4)
     *    /'RB  '   ,17.1784   , 1.7888   , 9.6435   ,17.3151   /
      DATA TBX( 76,5),TBX( 76,6),TBX( 76,7),TBX( 76,8),TBX( 76,9)
     *    /  5.1399   ,  .2748   , 1.5292   ,164.934   , 3.4873   /
      DATA TBX( 76,10)/37./
      DATA TBXC( 77),TBX( 77,1),TBX( 77,2),TBX( 77,3),TBX( 77,4)
     *    /'RB+1'   ,17.5816   , 1.7139   , 7.6598   ,14.7957   /
      DATA TBX( 77,5),TBX( 77,6),TBX( 77,7),TBX( 77,8),TBX( 77,9)
     *    /  5.8981   ,  .1603   , 2.7817   ,31.2087   , 2.0782   /
      DATA TBX( 77,10)/37./
      DATA TBXC( 78),TBX( 78,1),TBX( 78,2),TBX( 78,3),TBX( 78,4)
     *    /'SR  '   ,17.5663   , 1.5564   , 9.8184   ,14.0988   /
      DATA TBX( 78,5),TBX( 78,6),TBX( 78,7),TBX( 78,8),TBX( 78,9)
     *    /  5.4220   ,  .1664   , 2.6694   ,132.376   , 2.5064   /
      DATA TBX( 78,10)/38./
      DATA TBXC( 79),TBX( 79,1),TBX( 79,2),TBX( 79,3),TBX( 79,4)
     *    /'SR+2'   ,18.0874   , 1.4907   , 8.1373   ,12.6963   /
      DATA TBX( 79,5),TBX( 79,6),TBX( 79,7),TBX( 79,8),TBX( 79,9)
     *    /  2.5654   ,24.5651   ,-34.193   , -.0138   ,41.4025   /
      DATA TBX( 79,10)/38./
      DATA TBXC( 80),TBX( 80,1),TBX( 80,2),TBX( 80,3),TBX( 80,4)
     *    /'Y   '   ,17.7760   ,1.40290   ,10.2946   ,12.8006   /
      DATA TBX( 80,5),TBX( 80,6),TBX( 80,7),TBX( 80,8),TBX( 80,9)
     *    / 5.72629   ,.125599   ,3.26588   ,104.354   ,1.91213   /
      DATA TBX( 80,10)/39./
      DATA TBXC( 81),TBX( 81,1),TBX( 81,2),TBX( 81,3),TBX( 81,4)
     *    /'Y+3 '   ,17.9268   ,1.35417   ,9.15310   ,11.2145   /
      DATA TBX( 81,5),TBX( 81,6),TBX( 81,7),TBX( 81,8),TBX( 81,9)
     *    / 1.76795   ,22.6599   ,-33.108   ,-.01319   ,40.2602   /
      DATA TBX( 81,10)/39./
      DATA TBXC( 82),TBX( 82,1),TBX( 82,2),TBX( 82,3),TBX( 82,4)
     *    /'ZR  '   ,17.8765   ,1.27618   ,10.9480   ,11.9160   /
      DATA TBX( 82,5),TBX( 82,6),TBX( 82,7),TBX( 82,8),TBX( 82,9)
     *    / 5.41732   ,.117622   ,3.65721   ,87.6627   ,2.06929   /
      DATA TBX( 82,10)/40./
      DATA TBXC( 83),TBX( 83,1),TBX( 83,2),TBX( 83,3),TBX( 83,4)
     *    /'ZR+4'   ,18.1668   ,1.21480   ,10.0562   ,10.1483   /
      DATA TBX( 83,5),TBX( 83,6),TBX( 83,7),TBX( 83,8),TBX( 83,9)
     *    / 1.01118   ,21.6054   ,-2.6479   ,-.10276   ,9.41454   /
      DATA TBX( 83,10)/40./
      DATA TBXC( 84),TBX( 84,1),TBX( 84,2),TBX( 84,3),TBX( 84,4)
     *    /'NB  '   ,17.6142   ,1.18865   ,12.0144   ,11.7660   /
      DATA TBX( 84,5),TBX( 84,6),TBX( 84,7),TBX( 84,8),TBX( 84,9)
     *    / 4.04183   ,.204785   ,3.53346   ,69.7957   ,3.75591   /
      DATA TBX( 84,10)/41./
      DATA TBXC( 85),TBX( 85,1),TBX( 85,2),TBX( 85,3),TBX( 85,4)
     *    /'NB+3'   ,19.8812   ,.019175   ,18.0653   ,1.13305   /
      DATA TBX( 85,5),TBX( 85,6),TBX( 85,7),TBX( 85,8),TBX( 85,9)
     *    / 11.0177   ,10.1621   ,1.94715   ,28.3389   ,-12.912   /
      DATA TBX( 85,10)/41./
      DATA TBXC( 86),TBX( 86,1),TBX( 86,2),TBX( 86,3),TBX( 86,4)
     *    /'NB+5'   ,17.9163   ,1.12446   ,13.3417   ,.028781   /
      DATA TBX( 86,5),TBX( 86,6),TBX( 86,7),TBX( 86,8),TBX( 86,9)
     *    / 10.7990   ,9.28206   ,.337905   ,25.7228   ,-6.3934   /
      DATA TBX( 86,10)/41./
      DATA TBXC( 87),TBX( 87,1),TBX( 87,2),TBX( 87,3),TBX( 87,4)
     *    /'MO  '   , 3.7025   ,  .2772   ,17.2356   , 1.0958   /
      DATA TBX( 87,5),TBX( 87,6),TBX( 87,7),TBX( 87,8),TBX( 87,9)
     *    / 12.8876   ,11.0040   , 3.7429   ,61.6584   , 4.3875   /
      DATA TBX( 87,10)/42./
      DATA TBXC( 88),TBX( 88,1),TBX( 88,2),TBX( 88,3),TBX( 88,4)
     *    /'MO+3'   ,21.1664   ,.014734   ,18.2017   ,1.03031   /
      DATA TBX( 88,5),TBX( 88,6),TBX( 88,7),TBX( 88,8),TBX( 88,9)
     *    / 11.7423   ,9.53659   ,2.30951   ,26.6307   ,-14.421   /
      DATA TBX( 88,10)/42./
      DATA TBXC( 89),TBX( 89,1),TBX( 89,2),TBX( 89,3),TBX( 89,4)
     *    /'MO+5'   ,21.0149   ,.014345   ,18.0992   ,1.02238   /
      DATA TBX( 89,5),TBX( 89,6),TBX( 89,7),TBX( 89,8),TBX( 89,9)
     *    / 11.4632   ,8.78809   ,.740625   ,23.3452   ,-14.316   /
      DATA TBX( 89,10)/42./
      DATA TBXC( 90),TBX( 90,1),TBX( 90,2),TBX( 90,3),TBX( 90,4)
     *    /'MO+6'   ,17.8871   ,1.03649   ,11.1750   ,8.48061   /
      DATA TBX( 90,5),TBX( 90,6),TBX( 90,7),TBX( 90,8),TBX( 90,9)
     *    / 6.57891   ,.058881   ,.0000   ,.0000    ,.344941   /
      DATA TBX( 90,10)/42./
      DATA TBXC( 91),TBX( 91,1),TBX( 91,2),TBX( 91,3),TBX( 91,4)
     *    /'TC  '   ,19.1301   ,.864132   ,11.0948   ,8.14487   /
      DATA TBX( 91,5),TBX( 91,6),TBX( 91,7),TBX( 91,8),TBX( 91,9)
     *    / 4.64901   ,21.5707   ,2.71263   ,86.8472   ,5.40428   /
      DATA TBX( 91,10)/43./
      DATA TBXC( 92),TBX( 92,1),TBX( 92,2),TBX( 92,3),TBX( 92,4)
     *    /'RU  '   ,19.2674   ,.808520   ,12.9182   ,8.43467   /
      DATA TBX( 92,5),TBX( 92,6),TBX( 92,7),TBX( 92,8),TBX( 92,9)
     *    / 4.86337   ,24.7997   ,1.56756   ,94.2928   ,5.37874   /
      DATA TBX( 92,10)/44./
      DATA TBXC( 93),TBX( 93,1),TBX( 93,2),TBX( 93,3),TBX( 93,4)
     *    /'RU+3'   ,18.5638   ,.847329   ,13.2885   ,8.37164   /
      DATA TBX( 93,5),TBX( 93,6),TBX( 93,7),TBX( 93,8),TBX( 93,9)
     *    / 9.32602   ,.017662   ,3.00964   ,22.8870   ,-3.1892   /
      DATA TBX( 93,10)/44./
      DATA TBXC( 94),TBX( 94,1),TBX( 94,2),TBX( 94,3),TBX( 94,4)
     *    /'RU+4'   ,18.5003   ,.844582   ,13.1787   ,8.12534   /
      DATA TBX( 94,5),TBX( 94,6),TBX( 94,7),TBX( 94,8),TBX( 94,9)
     *    / 4.71304   ,.036495   ,2.18535   ,20.8504   ,1.42357   /
      DATA TBX( 94,10)/44./
      DATA TBXC( 95),TBX( 95,1),TBX( 95,2),TBX( 95,3),TBX( 95,4)
     *    /'RH  '   ,19.2957   ,.751536   ,14.3501   ,8.21758   /
      DATA TBX( 95,5),TBX( 95,6),TBX( 95,7),TBX( 95,8),TBX( 95,9)
     *    / 4.73425   ,25.8749   ,1.28918   ,98.6062   ,5.32800   /
      DATA TBX( 95,10)/45./
      DATA TBXC( 96),TBX( 96,1),TBX( 96,2),TBX( 96,3),TBX( 96,4)
     *    /'RH+3'   ,18.8785   ,.764252   ,14.1259   ,7.84438   /
      DATA TBX( 96,5),TBX( 96,6),TBX( 96,7),TBX( 96,8),TBX( 96,9)
     *    / 3.32515   ,21.2487   ,-6.1989   ,-.01036   ,11.8678   /
      DATA TBX( 96,10)/45./
      DATA TBXC( 97),TBX( 97,1),TBX( 97,2),TBX( 97,3),TBX( 97,4)
     *    /'RH+4'   ,18.8545   ,.760825   ,13.9806   ,7.62436   /
      DATA TBX( 97,5),TBX( 97,6),TBX( 97,7),TBX( 97,8),TBX( 97,9)
     *    / 2.53464   ,19.3317   ,-5.6526   ,-.01020   ,11.2835   /
      DATA TBX( 97,10)/45./
      DATA TBXC( 98),TBX( 98,1),TBX( 98,2),TBX( 98,3),TBX( 98,4)
     *    /'PD  '   ,19.3319   ,.698655   ,15.5017   ,7.98929   /
      DATA TBX( 98,5),TBX( 98,6),TBX( 98,7),TBX( 98,8),TBX( 98,9)
     *    / 5.29537   ,25.2052   ,.605844   ,76.8986   ,5.26593   /
      DATA TBX( 98,10)/46./
      DATA TBXC( 99),TBX( 99,1),TBX( 99,2),TBX( 99,3),TBX( 99,4)
     *    /'PD+2'   ,19.1701   ,.696219   ,15.2096   ,7.55573   /
      DATA TBX( 99,5),TBX( 99,6),TBX( 99,7),TBX( 99,8),TBX( 99,9)
     *    / 4.32234   ,22.5057   ,.0000   ,.0000   ,5.29160   /
      DATA TBX( 99,10)/46./
      DATA TBXC(100),TBX(100,1),TBX(100,2),TBX(100,3),TBX(100,4)
     *    /'PD+4'   ,19.2493   ,.683839   ,14.7900   ,7.14833   /
      DATA TBX(100,5),TBX(100,6),TBX(100,7),TBX(100,8),TBX(100,9)
     *    / 2.89289   ,17.9144   ,-7.9492   ,.005127   ,13.0174   /
      DATA TBX(100,10)/46./
      DATA TBXC(101),TBX(101,1),TBX(101,2),TBX(101,3),TBX(101,4)
     *    /'AG  '   ,19.2808   ,  .6446   ,16.6885   , 7.4726   /
      DATA TBX(101,5),TBX(101,6),TBX(101,7),TBX(101,8),TBX(101,9)
     *    /  4.8045   ,24.6605   , 1.0463   ,99.8156   , 5.1790   /
      DATA TBX(101,10)/47./
      DATA TBXC(102),TBX(102,1),TBX(102,2),TBX(102,3),TBX(102,4)
     *    /'AG+1'   ,19.1812   ,.646179   ,15.9719   ,7.19123   /
      DATA TBX(102,5),TBX(102,6),TBX(102,7),TBX(102,8),TBX(102,9)
     *    / 5.27475   ,21.7326   ,.357534   ,66.1147   ,5.21572   /
      DATA TBX(102,10)/47./
      DATA TBXC(103),TBX(103,1),TBX(103,2),TBX(103,3),TBX(103,4)
     *    /'AG+2'   ,19.1643   ,.645643   ,16.2456   ,7.18544   /
      DATA TBX(103,5),TBX(103,6),TBX(103,7),TBX(103,8),TBX(103,9)
     *    / 4.37090   ,21.4072   ,.0000   ,.0000   ,5.21404   /
      DATA TBX(103,10)/47./
      DATA TBXC(104),TBX(104,1),TBX(104,2),TBX(104,3),TBX(104,4)
     *    /'CD  '   ,19.2214   ,  .5946   ,17.6444   , 6.9089   /
      DATA TBX(104,5),TBX(104,6),TBX(104,7),TBX(104,8),TBX(104,9)
     *    /  4.4610   ,24.7008   , 1.6029   ,87.4825   , 5.0694   /
      DATA TBX(104,10)/48./
      DATA TBXC(105),TBX(105,1),TBX(105,2),TBX(105,3),TBX(105,4)
     *    /'CD+2'   ,19.1514   ,.597922   ,17.2535   ,6.80639   /
      DATA TBX(105,5),TBX(105,6),TBX(105,7),TBX(105,8),TBX(105,9)
     *    / 4.47128   ,20.2521   ,.0000   ,.0000   ,5.11937   /
      DATA TBX(105,10)/48./
      DATA TBXC(106),TBX(106,1),TBX(106,2),TBX(106,3),TBX(106,4)
     *    /'IN  '   ,19.1624   ,  .5476   ,18.5596   , 6.3776   /
      DATA TBX(106,5),TBX(106,6),TBX(106,7),TBX(106,8),TBX(106,9)
     *    /  4.2948   ,25.8499   , 2.0396   ,92.8029   , 4.9391   /
      DATA TBX(106,10)/49./
      DATA TBXC(107),TBX(107,1),TBX(107,2),TBX(107,3),TBX(107,4)
     *    /'IN+3'   ,19.1045   ,.551522   ,18.1108   ,6.32470   /
      DATA TBX(107,5),TBX(107,6),TBX(107,7),TBX(107,8),TBX(107,9)
     *    / 3.78897   ,17.3595   ,.0000   ,.0000   ,4.99635   /
      DATA TBX(107,10)/49./
      DATA TBXC(108),TBX(108,1),TBX(108,2),TBX(108,3),TBX(108,4)
     *    /'SN  '   ,19.1889   , 5.8303   ,19.1005   ,  .5031   /
      DATA TBX(108,5),TBX(108,6),TBX(108,7),TBX(108,8),TBX(108,9)
     *    /  4.4585   ,26.8909   , 2.4663   ,83.9571   , 4.7821   /
      DATA TBX(108,10)/50./
      DATA TBXC(109),TBX(109,1),TBX(109,2),TBX(109,3),TBX(109,4)
     *    /'SN+2'   ,19.1094   ,  .5036   ,19.0548   , 5.8378   /
      DATA TBX(109,5),TBX(109,6),TBX(109,7),TBX(109,8),TBX(109,9)
     *    /  4.5648   ,23.3752   ,  .4870   ,62.2061   , 4.7861   /
      DATA TBX(109,10)/50./
      DATA TBXC(110),TBX(110,1),TBX(110,2),TBX(110,3),TBX(110,4)
     *    /'SN+4'   ,18.9333   , 5.7640   ,19.7131   ,  .4655   /
      DATA TBX(110,5),TBX(110,6),TBX(110,7),TBX(110,8),TBX(110,9)
     *    /  3.4182   ,14.0049   ,  .0193   , -.7583   , 3.9182   /
      DATA TBX(110,10)/50./
      DATA TBXC(111),TBX(111,1),TBX(111,2),TBX(111,3),TBX(111,4)
     *    /'SB  '   ,19.6418   , 5.3034   ,19.0455   ,  .4607   /
      DATA TBX(111,5),TBX(111,6),TBX(111,7),TBX(111,8),TBX(111,9)
     *    /  5.0371   ,27.9074   , 2.6827   ,75.2825   , 4.5909   /
      DATA TBX(111,10)/51./
      DATA TBXC(112),TBX(112,1),TBX(112,2),TBX(112,3),TBX(112,4)
     *    /'SB+3'   ,18.9755   ,.467196   ,18.9330   ,5.22126   /
      DATA TBX(112,5),TBX(112,6),TBX(112,7),TBX(112,8),TBX(112,9)
     *    / 5.10789   ,19.5902   ,.288753   ,55.5113   ,4.69626   /
      DATA TBX(112,10)/51./
      DATA TBXC(113),TBX(113,1),TBX(113,2),TBX(113,3),TBX(113,4)
     *    /'SB+5'   ,19.8685   ,5.44853   ,19.0302   ,.467973   /
      DATA TBX(113,5),TBX(113,6),TBX(113,7),TBX(113,8),TBX(113,9)
     *    / 2.41253   ,14.1259   ,.0000   ,.0000   ,4.69263   /
      DATA TBX(113,10)/51./
      DATA TBXC(114),TBX(114,1),TBX(114,2),TBX(114,3),TBX(114,4)
     *    /'TE  '   ,19.9644   ,4.81742   ,19.0138   ,.420885   /
      DATA TBX(114,5),TBX(114,6),TBX(114,7),TBX(114,8),TBX(114,9)
     *    / 6.14487   ,28.5284   ,2.52390   ,70.8403   ,4.35200   /
      DATA TBX(114,10)/52./
      DATA TBXC(115),TBX(115,1),TBX(115,2),TBX(115,3),TBX(115,4)
     *    /'I   '   ,20.1472   , 4.3470   ,18.9949   ,  .3814   /
      DATA TBX(115,5),TBX(115,6),TBX(115,7),TBX(115,8),TBX(115,9)
     *    /  7.5138   ,27.7660   , 2.2735   ,66.8776   , 4.0712   /
      DATA TBX(115,10)/53./
      DATA TBXC(116),TBX(116,1),TBX(116,2),TBX(116,3),TBX(116,4)
     *    /'I-1 '   ,20.2332   , 4.3579   ,18.9970   ,  .3815   /
      DATA TBX(116,5),TBX(116,6),TBX(116,7),TBX(116,8),TBX(116,9)
     *    /  7.8069   ,29.5259   , 2.8868   ,84.9304   , 4.0714   /
      DATA TBX(116,10)/53./
      DATA TBXC(117),TBX(117,1),TBX(117,2),TBX(117,3),TBX(117,4)
     *    /'XE  '   ,20.2933   , 3.9282   ,19.0298   ,  .3440   /
      DATA TBX(117,5),TBX(117,6),TBX(117,7),TBX(117,8),TBX(117,9)
     *    /  8.9767   ,26.4659   , 1.9900   ,64.2658   , 3.7118   /
      DATA TBX(117,10)/54./
      DATA TBXC(118),TBX(118,1),TBX(118,2),TBX(118,3),TBX(118,4)
     *    /'CS  '   ,20.3892   , 3.5690   ,19.1062   ,  .3107   /
      DATA TBX(118,5),TBX(118,6),TBX(118,7),TBX(118,8),TBX(118,9)
     *    / 10.6620   ,24.3879   , 1.4953   ,213.904   , 3.3352   /
      DATA TBX(118,10)/55./
      DATA TBXC(119),TBX(119,1),TBX(119,2),TBX(119,3),TBX(119,4)
     *    /'CS+1'   ,20.3524   , 3.5520   ,19.1278   ,  .3086   /
      DATA TBX(119,5),TBX(119,6),TBX(119,7),TBX(119,8),TBX(119,9)
     *    / 10.2821   ,23.7128   ,  .9615   ,59.4565   , 3.2791   /
      DATA TBX(119,10)/55./
      DATA TBXC(120),TBX(120,1),TBX(120,2),TBX(120,3),TBX(120,4)
     *    /'BA  '   ,20.3361   , 3.2160   ,19.2970   ,  .2756   /
      DATA TBX(120,5),TBX(120,6),TBX(120,7),TBX(120,8),TBX(120,9)
     *    / 10.8880   ,20.2073   , 2.6959   ,167.202   , 2.7731   /
      DATA TBX(120,10)/56./
      DATA TBXC(121),TBX(121,1),TBX(121,2),TBX(121,3),TBX(121,4)
     *    /'BA+2'   ,20.1807   ,3.21367   ,19.1136   ,.283310   /
      DATA TBX(121,5),TBX(121,6),TBX(121,7),TBX(121,8),TBX(121,9)
     *    / 10.9054   ,20.0558   ,.773634   ,51.7460   ,3.02902   /
      DATA TBX(121,10)/56./
      DATA TBXC(122),TBX(122,1),TBX(122,2),TBX(122,3),TBX(122,4)
     *    /'LA  '   ,20.5780   ,2.94817   ,19.5990   ,.244475   /
      DATA TBX(122,5),TBX(122,6),TBX(122,7),TBX(122,8),TBX(122,9)
     *    / 11.3727   ,18.7726   ,3.28719   ,133.124   ,2.14678   /
      DATA TBX(122,10)/57./
      DATA TBXC(123),TBX(123,1),TBX(123,2),TBX(123,3),TBX(123,4)
     *    /'LA+3'   ,20.2489   ,2.92070   ,19.3763   ,.250698   /
      DATA TBX(123,5),TBX(123,6),TBX(123,7),TBX(123,8),TBX(123,9)
     *    / 11.6323   ,17.8211   ,.336048   ,54.9453   ,2.40860   /
      DATA TBX(123,10)/57./
      DATA TBXC(124),TBX(124,1),TBX(124,2),TBX(124,3),TBX(124,4)
     *    /'CE  '   ,21.1671   ,2.81219   ,19.7695   ,.226836   /
      DATA TBX(124,5),TBX(124,6),TBX(124,7),TBX(124,8),TBX(124,9)
     *    / 11.8513   ,17.6083   ,3.33049   ,127.113   ,1.86264   /
      DATA TBX(124,10)/58./
      DATA TBXC(125),TBX(125,1),TBX(125,2),TBX(125,3),TBX(125,4)
     *    /'CE+3'   ,20.8036   ,2.77691   ,19.5590   ,.231540   /
      DATA TBX(125,5),TBX(125,6),TBX(125,7),TBX(125,8),TBX(125,9)
     *    / 11.9369   ,16.5408   ,.612376   ,43.1692   ,2.09013   /
      DATA TBX(125,10)/58./
      DATA TBXC(126),TBX(126,1),TBX(126,2),TBX(126,3),TBX(126,4)
     *    /'CE+4'   ,20.3235   ,2.65941   ,19.8186   ,.218850   /
      DATA TBX(126,5),TBX(126,6),TBX(126,7),TBX(126,8),TBX(126,9)
     *    / 12.1233   ,15.7992   ,.144583   ,62.2355   ,1.59180   /
      DATA TBX(126,10)/58./
      DATA TBXC(127),TBX(127,1),TBX(127,2),TBX(127,3),TBX(127,4)
     *    /'PR  '   ,22.0440   ,2.77393   ,19.6697   ,.222087   /
      DATA TBX(127,5),TBX(127,6),TBX(127,7),TBX(127,8),TBX(127,9)
     *    / 12.3856   ,16.7669   ,2.82428   ,143.644   ,2.05830   /
      DATA TBX(127,10)/59./
      DATA TBXC(128),TBX(128,1),TBX(128,2),TBX(128,3),TBX(128,4)
     *    /'PR+3'   ,21.3727   ,2.64520   ,19.7491   ,.214299   /
      DATA TBX(128,5),TBX(128,6),TBX(128,7),TBX(128,8),TBX(128,9)
     *    / 12.1329   ,15.3230   ,.975180   ,36.4065   ,1.77132   /
      DATA TBX(128,10)/59./
      DATA TBXC(129),TBX(129,1),TBX(129,2),TBX(129,3),TBX(129,4)
     *    /'PR+4'   ,20.9413   ,2.54467   ,20.0539   ,.202481   /
      DATA TBX(129,5),TBX(129,6),TBX(129,7),TBX(129,8),TBX(129,9)
     *    / 12.4668   ,14.8137   ,.296689   ,45.4643   ,1.24285   /
      DATA TBX(129,10)/59./
      DATA TBXC(130),TBX(130,1),TBX(130,2),TBX(130,3),TBX(130,4)
     *    /'ND  '   ,22.6845   ,2.66248   ,19.6847   ,.210628   /
      DATA TBX(130,5),TBX(130,6),TBX(130,7),TBX(130,8),TBX(130,9)
     *    / 12.7740   ,15.8850   ,2.85137   ,137.903   ,1.98486   /
      DATA TBX(130,10)/60./
      DATA TBXC(131),TBX(131,1),TBX(131,2),TBX(131,3),TBX(131,4)
     *    /'ND+3'   ,21.9610   ,2.52722   ,19.9339   ,.199237   /
      DATA TBX(131,5),TBX(131,6),TBX(131,7),TBX(131,8),TBX(131,9)
     *    / 12.1200   ,14.1783   ,1.51031   ,30.8717   ,1.47588   /
      DATA TBX(131,10)/60./
      DATA TBXC(132),TBX(132,1),TBX(132,2),TBX(132,3),TBX(132,4)
     *    /'PM  '   ,23.3405   ,2.56270   ,19.6095   ,.202088   /
      DATA TBX(132,5),TBX(132,6),TBX(132,7),TBX(132,8),TBX(132,9)
     *    / 13.1235   ,15.1009   ,2.87516   ,132.721   ,2.02876   /
      DATA TBX(132,10)/61./
      DATA TBXC(133),TBX(133,1),TBX(133,2),TBX(133,3),TBX(133,4)
     *    /'PM+3'   ,22.5527   ,2.41740   ,20.1108   ,.185769   /
      DATA TBX(133,5),TBX(133,6),TBX(133,7),TBX(133,8),TBX(133,9)
     *    / 12.0671   ,13.1275   ,2.07492   ,27.4491   ,1.19499   /
      DATA TBX(133,10)/61./
      DATA TBXC(134),TBX(134,1),TBX(134,2),TBX(134,3),TBX(134,4)
     *    /'SM  '   ,24.0042   ,2.47274   ,19.4258   ,.196451   /
      DATA TBX(134,5),TBX(134,6),TBX(134,7),TBX(134,8),TBX(134,9)
     *    / 13.4396   ,14.3996   ,2.89604   ,128.007   ,2.20963   /
      DATA TBX(134,10)/62./
      DATA TBXC(135),TBX(135,1),TBX(135,2),TBX(135,3),TBX(135,4)
     *    /'SM+3'   ,23.1504   ,2.31641   ,20.2599   ,.174081   /
      DATA TBX(135,5),TBX(135,6),TBX(135,7),TBX(135,8),TBX(135,9)
     *    / 11.9202   ,12.1571   ,2.71488   ,24.8242   ,.954586   /
      DATA TBX(135,10)/62./
      DATA TBXC(136),TBX(136,1),TBX(136,2),TBX(136,3),TBX(136,4)
     *    /'EU  '   ,24.6274   , 2.3879   ,19.0886   ,  .1942   /
      DATA TBX(136,5),TBX(136,6),TBX(136,7),TBX(136,8),TBX(136,9)
     *    / 13.7603   ,13.7546   , 2.9227   ,123.174   , 2.5745   /
      DATA TBX(136,10)/63./
      DATA TBXC(137),TBX(137,1),TBX(137,2),TBX(137,3),TBX(137,4)
     *    /'EU+2'   ,24.0063   ,2.27783   ,19.9504   ,.173530   /
      DATA TBX(137,5),TBX(137,6),TBX(137,7),TBX(137,8),TBX(137,9)
     *    / 11.8034   ,11.6096   ,3.87243   ,26.5156   ,1.36389   /
      DATA TBX(137,10)/63./
      DATA TBXC(138),TBX(138,1),TBX(138,2),TBX(138,3),TBX(138,4)
     *    /'EU+3'   ,23.7497   ,2.22258   ,20.3745   ,.163940   /
      DATA TBX(138,5),TBX(138,6),TBX(138,7),TBX(138,8),TBX(138,9)
     *    / 11.8509   ,11.3110   ,3.26503   ,22.9966   ,.759344   /
      DATA TBX(138,10)/63./
      DATA TBXC(139),TBX(139,1),TBX(139,2),TBX(139,3),TBX(139,4)
     *    /'GD  '   ,25.0709   ,2.25341   ,19.0798   ,.181951   /
      DATA TBX(139,5),TBX(139,6),TBX(139,7),TBX(139,8),TBX(139,9)
     *    / 13.8518   ,12.9331   ,3.54545   ,101.398   ,2.41960   /
      DATA TBX(139,10)/64./
      DATA TBXC(140),TBX(140,1),TBX(140,2),TBX(140,3),TBX(140,4)
     *    /'GD+3'   ,24.3466   ,2.13553   ,20.4208   ,.155525   /
      DATA TBX(140,5),TBX(140,6),TBX(140,7),TBX(140,8),TBX(140,9)
     *    / 11.8708   ,10.5782   ,3.71490   ,21.7029   ,.645089   /
      DATA TBX(140,10)/64./
      DATA TBXC(141),TBX(141,1),TBX(141,2),TBX(141,3),TBX(141,4)
     *    /'TB  '   ,25.8976   ,2.24256   ,18.2185   ,.196143   /
      DATA TBX(141,5),TBX(141,6),TBX(141,7),TBX(141,8),TBX(141,9)
     *    / 14.3167   ,12.6648   ,2.95354   ,115.362   ,3.58324   /
      DATA TBX(141,10)/65./
      DATA TBXC(142),TBX(142,1),TBX(142,2),TBX(142,3),TBX(142,4)
     *    /'TB+3'   ,24.9559   ,2.05601   ,20.3271   ,.149525   /
      DATA TBX(142,5),TBX(142,6),TBX(142,7),TBX(142,8),TBX(142,9)
     *    / 12.2471   ,10.0499   ,3.77300   ,21.2773   ,.691967   /
      DATA TBX(142,10)/65./
      DATA TBXC(143),TBX(143,1),TBX(143,2),TBX(143,3),TBX(143,4)
     *    /'DY  '   ,26.5070   ,2.18020   ,17.6383   ,.202172   /
      DATA TBX(143,5),TBX(143,6),TBX(143,7),TBX(143,8),TBX(143,9)
     *    / 14.5596   ,12.1899   ,2.96577   ,111.874   ,4.29728   /
      DATA TBX(143,10)/66./
      DATA TBXC(144),TBX(144,1),TBX(144,2),TBX(144,3),TBX(144,4)
     *    /'DY+3'   ,25.5395   ,1.98040   ,20.2861   ,.143384   /
      DATA TBX(144,5),TBX(144,6),TBX(144,7),TBX(144,8),TBX(144,9)
     *    / 11.9812   ,9.34972   ,4.50073   ,19.5810   ,.689690   /
      DATA TBX(144,10)/66./
      DATA TBXC(145),TBX(145,1),TBX(145,2),TBX(145,3),TBX(145,4)
     *    /'HO  '   ,26.9049   ,2.07051   ,17.2940   ,.197940   /
      DATA TBX(145,5),TBX(145,6),TBX(145,7),TBX(145,8),TBX(145,9)
     *    / 14.5583   ,11.4407   ,3.63837   ,92.6566   ,4.56796   /
      DATA TBX(145,10)/67./
      DATA TBXC(146),TBX(146,1),TBX(146,2),TBX(146,3),TBX(146,4)
     *    /'HO+3'   ,26.1296   ,1.91072   ,20.0994   ,.139358   /
      DATA TBX(146,5),TBX(146,6),TBX(146,7),TBX(146,8),TBX(146,9)
     *    / 11.9788   ,8.80018   ,4.93676   ,18.5908   ,.852795   /
      DATA TBX(146,10)/67./
      DATA TBXC(147),TBX(147,1),TBX(147,2),TBX(147,3),TBX(147,4)
     *    /'ER  '   ,27.6563   ,2.07356   ,16.4285   ,.223545   /
      DATA TBX(147,5),TBX(147,6),TBX(147,7),TBX(147,8),TBX(147,9)
     *    / 14.9779   ,11.3604   ,2.98233   ,105.703   ,5.92046   /
      DATA TBX(147,10)/68./
      DATA TBXC(148),TBX(148,1),TBX(148,2),TBX(148,3),TBX(148,4)
     *    /'ER+3'   ,26.7220   ,1.84659   ,19.7748   ,.137290   /
      DATA TBX(148,5),TBX(148,6),TBX(148,7),TBX(148,8),TBX(148,9)
     *    / 12.1506   ,8.36225   ,5.17379   ,17.8974   ,1.17613   /
      DATA TBX(148,10)/68./
      DATA TBXC(149),TBX(149,1),TBX(149,2),TBX(149,3),TBX(149,4)
     *    /'TM  '   ,28.1819   ,2.02859   ,15.8851   ,.238849   /
      DATA TBX(149,5),TBX(149,6),TBX(149,7),TBX(149,8),TBX(149,9)
     *    / 15.1542   ,10.9975   ,2.98706   ,102.961   ,6.75621   /
      DATA TBX(149,10)/69./
      DATA TBXC(150),TBX(150,1),TBX(150,2),TBX(150,3),TBX(150,4)
     *    /'TM+3'   ,27.3083   ,1.78711   ,19.3320   ,.136974   /
      DATA TBX(150,5),TBX(150,6),TBX(150,7),TBX(150,8),TBX(150,9)
     *    / 12.3339   ,7.96778   ,5.38348   ,17.2922   ,1.63929   /
      DATA TBX(150,10)/69./
      DATA TBXC(151),TBX(151,1),TBX(151,2),TBX(151,3),TBX(151,4)
     *    /'YB  '   ,28.6641   ,1.98890   ,15.4345   ,.257119   /
      DATA TBX(151,5),TBX(151,6),TBX(151,7),TBX(151,8),TBX(151,9)
     *    / 15.3087   ,10.6647   ,2.98963   ,100.417   ,7.56672   /
      DATA TBX(151,10)/70./
      DATA TBXC(152),TBX(152,1),TBX(152,2),TBX(152,3),TBX(152,4)
     *    /'YB+2'   ,28.1209   ,1.78503   ,17.6817   ,.159970   /
      DATA TBX(152,5),TBX(152,6),TBX(152,7),TBX(152,8),TBX(152,9)
     *    / 13.3335   ,8.18304   ,5.14657   ,20.3900   ,3.70983   /
      DATA TBX(152,10)/70./
      DATA TBXC(153),TBX(153,1),TBX(153,2),TBX(153,3),TBX(153,4)
     *    /'YB+3'   ,27.8917   ,1.73272   ,18.7614   ,.138790   /
      DATA TBX(153,5),TBX(153,6),TBX(153,7),TBX(153,8),TBX(153,9)
     *    / 12.6072   ,7.64412   ,5.47647   ,16.8153   ,2.26001   /
      DATA TBX(153,10)/70./
      DATA TBXC(154),TBX(154,1),TBX(154,2),TBX(154,3),TBX(154,4)
     *    /'LU  '   ,28.9476   ,1.90182   ,15.2208   ,9.98519   /
      DATA TBX(154,5),TBX(154,6),TBX(154,7),TBX(154,8),TBX(154,9)
     *    / 15.1000   ,.261033   ,3.71601   ,84.3298   ,7.97628   /
      DATA TBX(154,10)/71./
      DATA TBXC(155),TBX(155,1),TBX(155,2),TBX(155,3),TBX(155,4)
     *    /'LU+3'   ,28.4628   ,1.68216   ,18.1210   ,.142292   /
      DATA TBX(155,5),TBX(155,6),TBX(155,7),TBX(155,8),TBX(155,9)
     *    / 12.8429   ,7.33727   ,5.59415   ,16.3535   ,2.97573   /
      DATA TBX(155,10)/71./
      DATA TBXC(156),TBX(156,1),TBX(156,2),TBX(156,3),TBX(156,4)
     *    /'HF  '   ,29.1440   ,1.83262   ,15.1726   ,9.59990   /
      DATA TBX(156,5),TBX(156,6),TBX(156,7),TBX(156,8),TBX(156,9)
     *    / 14.7586   ,.275116   ,4.30013   ,72.0290   ,8.58154   /
      DATA TBX(156,10)/72./
      DATA TBXC(157),TBX(157,1),TBX(157,2),TBX(157,3),TBX(157,4)
     *    /'HF+4'   ,28.8131   ,1.59136   ,18.4601   ,.128903   /
      DATA TBX(157,5),TBX(157,6),TBX(157,7),TBX(157,8),TBX(157,9)
     *    / 12.7285   ,6.76232   ,5.59927   ,14.0366   ,2.39699   /
      DATA TBX(157,10)/72./
      DATA TBXC(158),TBX(158,1),TBX(158,2),TBX(158,3),TBX(158,4)
     *    /'TA  '   ,29.2024   ,1.77333   ,15.2293   ,9.37046   /
      DATA TBX(158,5),TBX(158,6),TBX(158,7),TBX(158,8),TBX(158,9)
     *    / 14.5135   ,.295977   ,4.76492   ,63.3644   ,9.24354   /
      DATA TBX(158,10)/73./
      DATA TBXC(159),TBX(159,1),TBX(159,2),TBX(159,3),TBX(159,4)
     *    /'TA+5'   ,29.1587   ,1.50711   ,18.8407   ,.116741   /
      DATA TBX(159,5),TBX(159,6),TBX(159,7),TBX(159,8),TBX(159,9)
     *    / 12.8268   ,6.31524   ,5.38695   ,12.4244   ,1.78555   /
      DATA TBX(159,10)/73./
      DATA TBXC(160),TBX(160,1),TBX(160,2),TBX(160,3),TBX(160,4)
     *    /'W   '   ,29.0818   ,1.72029   ,15.4300   ,9.22590   /
      DATA TBX(160,5),TBX(160,6),TBX(160,7),TBX(160,8),TBX(160,9)
     *    / 14.4327   ,.321703   ,5.11982   ,57.0560   ,9.88750   /
      DATA TBX(160,10)/74./
      DATA TBXC(161),TBX(161,1),TBX(161,2),TBX(161,3),TBX(161,4)
     *    /'W+6 '   ,29.4936   ,1.42755   ,19.3763   ,.104621   /
      DATA TBX(161,5),TBX(161,6),TBX(161,7),TBX(161,8),TBX(161,9)
     *    / 13.0544   ,5.93667   ,5.06412   ,11.1972   ,1.01074   /
      DATA TBX(161,10)/74./
      DATA TBXC(162),TBX(162,1),TBX(162,2),TBX(162,3),TBX(162,4)
     *    /'RE  '   ,28.7621   ,1.67191   ,15.7189   ,9.09227   /
      DATA TBX(162,5),TBX(162,6),TBX(162,7),TBX(162,8),TBX(162,9)
     *    / 14.5564   ,.350500   ,5.44174   ,52.0861   ,10.4720   /
      DATA TBX(162,10)/75./
      DATA TBXC(163),TBX(163,1),TBX(163,2),TBX(163,3),TBX(163,4)
     *    /'OS  '   ,28.1894   ,1.62903   ,16.1550   ,8.97948   /
      DATA TBX(163,5),TBX(163,6),TBX(163,7),TBX(163,8),TBX(163,9)
     *    / 14.9305   ,.382661   ,5.67589   ,48.1647   ,11.0005   /
      DATA TBX(163,10)/76./
      DATA TBXC(164),TBX(164,1),TBX(164,2),TBX(164,3),TBX(164,4)
     *    /'OS+4'   ,30.4190   ,1.37113   ,15.2637   ,6.84706   /
      DATA TBX(164,5),TBX(164,6),TBX(164,7),TBX(164,8),TBX(164,9)
     *    / 14.7458   ,.165191   ,5.06795   ,18.0030   ,6.49804   /
      DATA TBX(164,10)/76./
      DATA TBXC(165),TBX(165,1),TBX(165,2),TBX(165,3),TBX(165,4)
     *    /'IR  '   ,27.3049   ,1.59279   ,16.7296   ,8.86553   /
      DATA TBX(165,5),TBX(165,6),TBX(165,7),TBX(165,8),TBX(165,9)
     *    / 15.6115   ,.417916   ,5.83377   ,45.0011   ,11.4722   /
      DATA TBX(165,10)/77./
      DATA TBXC(166),TBX(166,1),TBX(166,2),TBX(166,3),TBX(166,4)
     *    /'IR+3'   ,30.4156   ,1.34323   ,15.8620   ,7.10909   /
      DATA TBX(166,5),TBX(166,6),TBX(166,7),TBX(166,8),TBX(166,9)
     *    / 13.6145   ,.204633   ,5.82008   ,20.3254   ,8.27903   /
      DATA TBX(166,10)/77./
      DATA TBXC(167),TBX(167,1),TBX(167,2),TBX(167,3),TBX(167,4)
     *    /'IR+4'   ,30.7058   ,1.30923   ,15.5512   ,6.71983   /
      DATA TBX(167,5),TBX(167,6),TBX(167,7),TBX(167,8),TBX(167,9)
     *    / 14.2326   ,.167252   ,5.53672   ,17.4911   ,6.96824   /
      DATA TBX(167,10)/77./
      DATA TBXC(168),TBX(168,1),TBX(168,2),TBX(168,3),TBX(168,4)
     *    /'PT  '   ,27.0059   ,1.51293   ,17.7639   ,8.81174   /
      DATA TBX(168,5),TBX(168,6),TBX(168,7),TBX(168,8),TBX(168,9)
     *    / 15.7131   ,.424593   ,5.78370   ,38.6103   ,11.6883   /
      DATA TBX(168,10)/78./
      DATA TBXC(169),TBX(169,1),TBX(169,2),TBX(169,3),TBX(169,4)
     *    /'PT+2'   ,29.8429   ,1.32927   ,16.7224   ,7.38979   /
      DATA TBX(169,5),TBX(169,6),TBX(169,7),TBX(169,8),TBX(169,9)
     *    / 13.2153   ,.263297   ,6.35234   ,22.9426   ,9.85329   /
      DATA TBX(169,10)/78./
      DATA TBXC(170),TBX(170,1),TBX(170,2),TBX(170,3),TBX(170,4)
     *    /'PT+4'   ,30.9612   ,1.24813   ,15.9829   ,6.60834   /
      DATA TBX(170,5),TBX(170,6),TBX(170,7),TBX(170,8),TBX(170,9)
     *    / 13.7348   ,.168640   ,5.92034   ,16.9392   ,7.39534   /
      DATA TBX(170,10)/78./
      DATA TBXC(171),TBX(171,1),TBX(171,2),TBX(171,3),TBX(171,4)
     *    /'AU  '   ,16.8819   ,  .4611   ,18.5913   , 8.6216   /
      DATA TBX(171,5),TBX(171,6),TBX(171,7),TBX(171,8),TBX(171,9)
     *    / 25.5582   , 1.4826   , 5.8600   ,36.3956   ,12.0658   /
      DATA TBX(171,10)/79./
      DATA TBXC(172),TBX(172,1),TBX(172,2),TBX(172,3),TBX(172,4)
     *    /'AU+1'   ,28.0109   ,1.35321   ,17.8204   ,7.73950   /
      DATA TBX(172,5),TBX(172,6),TBX(172,7),TBX(172,8),TBX(172,9)
     *    / 14.3359   ,.356752   ,6.58077   ,26.4043   ,11.2299   /
      DATA TBX(172,10)/79./
      DATA TBXC(173),TBX(173,1),TBX(173,2),TBX(173,3),TBX(173,4)
     *    /'AU+3'   ,30.6886   ,1.21990   ,16.9029   ,6.82872   /
      DATA TBX(173,5),TBX(173,6),TBX(173,7),TBX(173,8),TBX(173,9)
     *    / 12.7801   ,.212867   ,6.52354   ,18.6590   ,9.09680   /
      DATA TBX(173,10)/79./
      DATA TBXC(174),TBX(174,1),TBX(174,2),TBX(174,3),TBX(174,4)
     *    /'HG  '   ,20.6809   ,  .5450   ,19.0417   , 8.4484   /
      DATA TBX(174,5),TBX(174,6),TBX(174,7),TBX(174,8),TBX(174,9)
     *    / 21.6575   , 1.5729   , 5.9676   ,38.3246   ,12.6089   /
      DATA TBX(174,10)/80./
      DATA TBXC(175),TBX(175,1),TBX(175,2),TBX(175,3),TBX(175,4)
     *    /'HG+1'   ,25.0853   ,1.39507   ,18.4973   ,7.65105   /
      DATA TBX(175,5),TBX(175,6),TBX(175,7),TBX(175,8),TBX(175,9)
     *    / 16.8883   ,.443378   ,6.48216   ,28.2262   ,12.0205   /
      DATA TBX(175,10)/80./
      DATA TBXC(176),TBX(176,1),TBX(176,2),TBX(176,3),TBX(176,4)
     *    /'HG+2'   ,29.5641   ,1.21152   ,18.0600   ,7.05639   /
      DATA TBX(176,5),TBX(176,6),TBX(176,7),TBX(176,8),TBX(176,9)
     *    / 12.8374   ,.284738   ,6.89912   ,20.7482   ,10.6268   /
      DATA TBX(176,10)/80./
      DATA TBXC(177),TBX(177,1),TBX(177,2),TBX(177,3),TBX(177,4)
     *    /'TL  '   ,27.5446   ,.655150   ,19.1584   ,8.70751   /
      DATA TBX(177,5),TBX(177,6),TBX(177,7),TBX(177,8),TBX(177,9)
     *    / 15.5380   ,1.96347   ,5.52593   ,45.8149   ,13.1746   /
      DATA TBX(177,10)/81./
      DATA TBXC(178),TBX(178,1),TBX(178,2),TBX(178,3),TBX(178,4)
     *    /'TL+1'   ,21.3985   ,1.47110   ,20.4723   ,.517394   /
      DATA TBX(178,5),TBX(178,6),TBX(178,7),TBX(178,8),TBX(178,9)
     *    / 18.7478   ,7.43463   ,6.82847   ,28.8482   ,12.5258   /
      DATA TBX(178,10)/81./
      DATA TBXC(179),TBX(179,1),TBX(179,2),TBX(179,3),TBX(179,4)
     *    /'TL+3'   ,30.8695   ,1.10080   ,18.3841   ,6.53852   /
      DATA TBX(179,5),TBX(179,6),TBX(179,7),TBX(179,8),TBX(179,9)
     *    / 11.9328   ,.219074   ,7.00574   ,17.2114   ,9.80270   /
      DATA TBX(179,10)/81./
      DATA TBXC(180),TBX(180,1),TBX(180,2),TBX(180,3),TBX(180,4)
     *    /'PB  '   ,31.0617   ,  .6902   ,13.0637   , 2.3576   /
      DATA TBX(180,5),TBX(180,6),TBX(180,7),TBX(180,8),TBX(180,9)
     *    / 18.4420   , 8.6180   , 5.9696   ,47.2579   ,13.4118   /
      DATA TBX(180,10)/82./
      DATA TBXC(181),TBX(181,1),TBX(181,2),TBX(181,3),TBX(181,4)
     *    /'PB+2'   ,21.7886   ,1.33660   ,19.5682   ,.488383   /
      DATA TBX(181,5),TBX(181,6),TBX(181,7),TBX(181,8),TBX(181,9)
     *    / 19.1406   ,6.77270   ,7.01107   ,23.8132   ,12.4734   /
      DATA TBX(181,10)/82./
      DATA TBXC(182),TBX(182,1),TBX(182,2),TBX(182,3),TBX(182,4)
     *    /'PB+4'   ,32.1244   ,1.00566   ,18.8003   ,6.10926   /
      DATA TBX(182,5),TBX(182,6),TBX(182,7),TBX(182,8),TBX(182,9)
     *    / 12.0175   ,.147041   ,6.96886   ,14.7140   ,8.08428   /
      DATA TBX(182,10)/82./
      DATA TBXC(183),TBX(183,1),TBX(183,2),TBX(183,3),TBX(183,4)
     *    /'BI  '   ,33.3689   ,  .7040   ,12.9510   , 2.9238   /
      DATA TBX(183,5),TBX(183,6),TBX(183,7),TBX(183,8),TBX(183,9)
     *    / 16.5877   , 8.7937   , 6.4692   ,48.0093   ,13.5782   /
      DATA TBX(183,10)/83./
      DATA TBXC(184),TBX(184,1),TBX(184,2),TBX(184,3),TBX(184,4)
     *    /'BI+3'   ,21.8053   ,1.23560   ,19.5026   ,6.24149   /
      DATA TBX(184,5),TBX(184,6),TBX(184,7),TBX(184,8),TBX(184,9)
     *    / 19.1053   ,.469999   ,7.10295   ,20.3185   ,12.4711   /
      DATA TBX(184,10)/83./
      DATA TBXC(185),TBX(185,1),TBX(185,2),TBX(185,3),TBX(185,4)
     *    /'BI+5'   ,33.5364   ,.916540   ,25.0946   ,.039042   /
      DATA TBX(185,5),TBX(185,6),TBX(185,7),TBX(185,8),TBX(185,9)
     *    / 19.2497   ,5.71414   ,6.91555   ,12.8285   ,-6.7994   /
      DATA TBX(185,10)/83./
      DATA TBXC(186),TBX(186,1),TBX(186,2),TBX(186,3),TBX(186,4)
     *    /'PO  '   ,34.6726   ,.700999   ,15.4733   ,3.55078   /
      DATA TBX(186,5),TBX(186,6),TBX(186,7),TBX(186,8),TBX(186,9)
     *    / 13.1138   ,9.55642   ,7.02588   ,47.0045   ,13.6770   /
      DATA TBX(186,10)/84./
      DATA TBXC(187),TBX(187,1),TBX(187,2),TBX(187,3),TBX(187,4)
     *    /'AT  '   ,35.3163   ,.685870   ,19.0211   ,3.97458   /
      DATA TBX(187,5),TBX(187,6),TBX(187,7),TBX(187,8),TBX(187,9)
     *    / 9.49887   ,11.3824   ,7.42518   ,45.4715   ,13.7108   /
      DATA TBX(187,10)/85./
      DATA TBXC(188),TBX(188,1),TBX(188,2),TBX(188,3),TBX(188,4)
     *    /'RN  '   ,35.5631   ,  .6631   ,21.2816   , 4.0691   /
      DATA TBX(188,5),TBX(188,6),TBX(188,7),TBX(188,8),TBX(188,9)
     *    /  8.0037   ,14.0422   , 7.4433   ,44.2473   ,13.6905   /
      DATA TBX(188,10)/86./
      DATA TBXC(189),TBX(189,1),TBX(189,2),TBX(189,3),TBX(189,4)
     *    /'FR  '   ,35.9299   ,.646453   ,23.0547   ,4.17619   /
      DATA TBX(189,5),TBX(189,6),TBX(189,7),TBX(189,8),TBX(189,9)
     *    / 12.1439   ,23.1052   ,2.11253   ,150.645   ,13.7247   /
      DATA TBX(189,10)/87./
      DATA TBXC(190),TBX(190,1),TBX(190,2),TBX(190,3),TBX(190,4)
     *    /'RA  '   ,35.7630   ,.616341   ,22.9064   ,3.87135   /
      DATA TBX(190,5),TBX(190,6),TBX(190,7),TBX(190,8),TBX(190,9)
     *    / 12.4739   ,19.9887   ,3.21097   ,142.325   ,13.6211   /
      DATA TBX(190,10)/88./
      DATA TBXC(191),TBX(191,1),TBX(191,2),TBX(191,3),TBX(191,4)
     *    /'RA+2'   ,35.2150   ,.604909   ,21.6700   ,3.57670   /
      DATA TBX(191,5),TBX(191,6),TBX(191,7),TBX(191,8),TBX(191,9)
     *    / 7.91342   ,12.6010   ,7.65078   ,29.8436   ,13.5431   /
      DATA TBX(191,10)/88./
      DATA TBXC(192),TBX(192,1),TBX(192,2),TBX(192,3),TBX(192,4)
     *    /'AC  '   ,35.6597   ,.589092   ,23.1032   ,3.65155   /
      DATA TBX(192,5),TBX(192,6),TBX(192,7),TBX(192,8),TBX(192,9)
     *    / 12.5977   ,18.5990   ,4.08655   ,117.020   ,13.5266   /
      DATA TBX(192,10)/89./
      DATA TBXC(193),TBX(193,1),TBX(193,2),TBX(193,3),TBX(193,4)
     *    /'AC+3'   ,35.1736   ,.579689   ,22.1112   ,3.41437   /
      DATA TBX(193,5),TBX(193,6),TBX(193,7),TBX(193,8),TBX(193,9)
     *    / 8.19216   ,12.9187   ,7.05545   ,25.9443   ,13.4637   /
      DATA TBX(193,10)/89./
      DATA TBXC(194),TBX(194,1),TBX(194,2),TBX(194,3),TBX(194,4)
     *    /'TH  '   ,35.5645   ,.563359   ,23.4219   ,3.46204   /
      DATA TBX(194,5),TBX(194,6),TBX(194,7),TBX(194,8),TBX(194,9)
     *    / 12.7473   ,17.8309   ,4.80703   ,99.1722   ,13.4314   /
      DATA TBX(194,10)/90./
      DATA TBXC(195),TBX(195,1),TBX(195,2),TBX(195,3),TBX(195,4)
     *    /'TH+4'   ,35.1007   ,.555054   ,22.4418   ,3.24498   /
      DATA TBX(195,5),TBX(195,6),TBX(195,7),TBX(195,8),TBX(195,9)
     *    / 9.78554   ,13.4661   ,5.29444   ,23.9533   ,13.3760   /
      DATA TBX(195,10)/90./
      DATA TBXC(196),TBX(196,1),TBX(196,2),TBX(196,3),TBX(196,4)
     *    /'PA  '   ,35.8847   ,.547751   ,23.2948   ,3.41519   /
      DATA TBX(196,5),TBX(196,6),TBX(196,7),TBX(196,8),TBX(196,9)
     *    / 14.1891   ,16.9235   ,4.17287   ,105.251   ,13.4287   /
      DATA TBX(196,10)/91./
      DATA TBXC(197),TBX(197,1),TBX(197,2),TBX(197,3),TBX(197,4)
     *    /'U   '   ,36.0228   ,  .5293   ,23.4128   , 3.3253   /
      DATA TBX(197,5),TBX(197,6),TBX(197,7),TBX(197,8),TBX(197,9)
     *    / 14.9491   ,16.0927   , 4.1880   ,100.613   ,13.3966   /
      DATA TBX(197,10)/92./
      DATA TBXC(198),TBX(198,1),TBX(198,2),TBX(198,3),TBX(198,4)
     *    /'U+3 '  ,35.5747   ,.520480   ,22.5259   ,3.12293   /
      DATA TBX(198,5),TBX(198,6),TBX(198,7),TBX(198,8),TBX(198,9)
     *    / 12.2165   ,12.7148   ,5.37073   ,26.3394   ,13.3092   /
      DATA TBX(198,10)/92./
      DATA TBXC(199),TBX(199,1),TBX(199,2),TBX(199,3),TBX(199,4)
     *    /'U+4 '  ,35.3715   ,.516598   ,22.5326   ,3.05053   /
      DATA TBX(199,5),TBX(199,6),TBX(199,7),TBX(199,8),TBX(199,9)
     *    / 12.0291   ,12.5723   ,4.79840   ,23.4582   ,13.2671   /
      DATA TBX(199,10)/92./
      DATA TBXC(200),TBX(200,1),TBX(200,2),TBX(200,3),TBX(200,4)
     *    /'U+6 '  ,34.8509   ,.507079   ,22.7584   ,2.89030   /
      DATA TBX(200,5),TBX(200,6),TBX(200,7),TBX(200,8),TBX(200,9)
     *    / 14.0099   ,13.1767   ,1.21457   ,25.2017   ,13.1665   /
      DATA TBX(200,10)/92./
      DATA TBXC(201),TBX(201,1),TBX(201,2),TBX(201,3),TBX(201,4)
     *    /'NP  '   ,36.1874   ,.511929   ,23.5964   ,3.25396   /
      DATA TBX(201,5),TBX(201,6),TBX(201,7),TBX(201,8),TBX(201,9)
     *    / 15.6402   ,15.3622   ,4.18550   ,97.4908   ,13.3573   /
      DATA TBX(201,10)/93./
      DATA TBXC(202),TBX(202,1),TBX(202,2),TBX(202,3),TBX(202,4)
     *    /'NP+3'   ,35.7074   ,.502322   ,22.6130   ,3.03807   /
      DATA TBX(202,5),TBX(202,6),TBX(202,7),TBX(202,8),TBX(202,9)
     *    / 12.9898   ,12.1449   ,5.43227   ,25.4928   ,13.2544   /
      DATA TBX(202,10)/93./
      DATA TBXC(203),TBX(203,1),TBX(203,2),TBX(203,3),TBX(203,4)
     *    /'NP+4'   ,35.5103   ,.498626   ,22.5787   ,2.96627   /
      DATA TBX(203,5),TBX(203,6),TBX(203,7),TBX(203,8),TBX(203,9)
     *    / 12.7766   ,11.9484   ,4.92159   ,22.7502   ,13.2116   /
      DATA TBX(203,10)/93./
      DATA TBXC(204),TBX(204,1),TBX(204,2),TBX(204,3),TBX(204,4)
     *    /'NP+6'   ,35.0136   ,.489810   ,22.7286   ,2.81099   /
      DATA TBX(204,5),TBX(204,6),TBX(204,7),TBX(204,8),TBX(204,9)
     *    / 14.3884   ,12.3300   ,1.75669   ,22.6581   ,13.1130   /
      DATA TBX(204,10)/93./
      DATA TBXC(205),TBX(205,1),TBX(205,2),TBX(205,3),TBX(205,4)
     *    /'PU  '   ,36.5254   ,.499384   ,23.8083   ,3.26371   /
      DATA TBX(205,5),TBX(205,6),TBX(205,7),TBX(205,8),TBX(205,9)
     *    / 16.7707   ,14.9455   ,3.47947   ,105.980   ,13.3812   /
      DATA TBX(205,10)/94./
      DATA TBXC(206),TBX(206,1),TBX(206,2),TBX(206,3),TBX(206,4)
     *    /'PU+3'   ,35.8400   ,.484938   ,22.7169   ,2.96118   /
      DATA TBX(206,5),TBX(206,6),TBX(206,7),TBX(206,8),TBX(206,9)
     *    / 13.5807   ,11.5331   ,5.66016   ,24.3992   ,13.1991   /
      DATA TBX(206,10)/94./
      DATA TBXC(207),TBX(207,1),TBX(207,2),TBX(207,3),TBX(207,4)
     *    /'PU+4'   ,35.6493   ,.481422   ,22.6460   ,2.89020   /
      DATA TBX(207,5),TBX(207,6),TBX(207,7),TBX(207,8),TBX(207,9)
     *    / 13.3595   ,11.3160   ,5.18831   ,21.8301   ,13.1555   /
      DATA TBX(207,10)/94./
      DATA TBXC(208),TBX(208,1),TBX(208,2),TBX(208,3),TBX(208,4)
     *    /'PU+6'   ,35.1736   ,.473204   ,22.7181   ,2.73848   /
      DATA TBX(208,5),TBX(208,6),TBX(208,7),TBX(208,8),TBX(208,9)
     *    / 14.7635   ,11.5530   ,2.28678   ,20.9303   ,13.0582   /
      DATA TBX(208,10)/94./
      DATA TBXC(209),TBX(209,1),TBX(209,2),TBX(209,3),TBX(209,4)
     *    /'AM  '   ,36.6706   ,.483629   ,24.0992   ,3.20647   /
      DATA TBX(209,5),TBX(209,6),TBX(209,7),TBX(209,8),TBX(209,9)
     *    / 17.3415   ,14.3136   ,3.49331   ,102.273   ,13.3592   /
      DATA TBX(209,10)/95./
      DATA TBXC(210),TBX(210,1),TBX(210,2),TBX(210,3),TBX(210,4)
     *    /'CM  '   ,36.6488   ,.465154   ,24.4096   ,3.08997   /
      DATA TBX(210,5),TBX(210,6),TBX(210,7),TBX(210,8),TBX(210,9)
     *    / 17.3990   ,13.4346   ,4.21665   ,88.4834   ,13.2887   /
      DATA TBX(210,10)/96./
      DATA TBXC(211),TBX(211,1),TBX(211,2),TBX(211,3),TBX(211,4)
     *    /'BK  '   ,36.7881   ,.451018   ,24.7736   ,3.04619   /
      DATA TBX(211,5),TBX(211,6),TBX(211,7),TBX(211,8),TBX(211,9)
     *    / 17.8919   ,12.8946   ,4.23284   ,86.0030   ,13.2754   /
      DATA TBX(211,10)/97./
      DATA TBXC(212),TBX(212,1),TBX(212,2),TBX(212,3),TBX(212,4)
     *    /'CF  '   ,36.9185   ,.437533   ,25.1995   ,3.00775   /
      DATA TBX(212,5),TBX(212,6),TBX(212,7),TBX(212,8),TBX(212,9)
     *    / 18.3317   ,12.4044   ,4.24391   ,83.7881   ,13.2674   /
      DATA TBX(212,10)/98./
      DATA TBM(  1),TBM(  2),TBM(  3),TBM(  4),TBM(  5),TBM(  6)
     *    /1.0079  ,4.0026  ,6.94    ,9.01218 ,10.81   ,12.011   /
      DATA TBM(  7),TBM(  8),TBM(  9),TBM( 10),TBM( 11),TBM( 12)
     *    /14.0067 ,15.999  ,18.9984 ,20.17   ,22.98977,24.305   /
      DATA TBM( 13),TBM( 14),TBM( 15),TBM( 16),TBM( 17),TBM( 18)
     *    /26.98154,28.08   ,30.97376,32.06   ,35.453  ,39.94    /
      DATA TBM( 19),TBM( 20),TBM( 21),TBM( 22),TBM( 23),TBM( 24)
     *    /39.09   ,40.08   ,44.9559 ,47.9    ,50.941  ,51.996   /
      DATA TBM( 25),TBM( 26),TBM( 27),TBM( 28),TBM( 29),TBM( 30)
     *    /54.9380 ,55.84   ,58.9332 ,58.7    ,63.54   ,65.377   /
      DATA TBM( 31),TBM( 32),TBM( 33),TBM( 34),TBM( 35),TBM( 36)
     *    /69.72   ,72.5    ,74.9216 ,78.9    ,79.904  ,83.80    /
      DATA TBM( 37),TBM( 38),TBM( 39),TBM( 40),TBM( 41),TBM( 42)
     *    /85.467  ,87.62   ,88.9059 ,91.22   ,92.9064 ,95.9     /
      DATA TBM( 43),TBM( 44),TBM( 45),TBM( 46),TBM( 47),TBM( 48)
     *    /98.9062 ,101.0   ,102.9055,106.4   ,107.868 ,112.40   /
      DATA TBM( 49),TBM( 50),TBM( 51),TBM( 52),TBM( 53),TBM( 54)
     *    /114.82  ,118.6   ,121.70  ,127.60  ,126.9045,131.30   /
      DATA TBM( 55),TBM( 56),TBM( 57),TBM( 58),TBM( 59),TBM( 60)
     *    /132.9054,137.30  ,138.905 ,140.12  ,140.9077,144.20   /
      DATA TBM( 61),TBM( 62),TBM( 63),TBM( 64),TBM( 65),TBM( 66)
     *    /147.    ,150.40  ,151.96  ,157.20  ,158.9254,162.50   /
      DATA TBM( 67),TBM( 68),TBM( 69),TBM( 70),TBM( 71),TBM( 72)
     *    /164.9304,167.20  ,168.9342,173.0   ,174.97  ,178.4    /
      DATA TBM( 73),TBM( 74),TBM( 75),TBM( 76),TBM( 77),TBM( 78)
     *    /180.947 ,183.8   ,186.2   ,190.2   ,192.2   ,195.0    /
      DATA TBM( 79),TBM( 80),TBM( 81),TBM( 82),TBM( 83),TBM( 84)
     *    /196.9665,200.5   ,204.3   ,207.2   ,208.9804,210.     /
      DATA TBM( 85),TBM( 86),TBM( 87),TBM( 88),TBM( 89),TBM( 90)
     *    /210.    ,222.    ,223.    ,226.0254,227.    ,232.0381 /
      DATA TBM( 91),TBM( 92),TBM( 93),TBM( 94),TBM( 95),TBM( 96)
     *    /231.0359,238.029 ,237.0482,242.    ,243.    ,247.     /
      DATA TBM( 97),TBM( 98)
     *    /247.    ,249.    /
      END 
      subroutine gotoer
      write (6,10)
10    format(1x,'subroutine gotoer was called')
      stop
      end
      SUBROUTINE SORT(IPHASE)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      DIMENSION L(IRS),TEMP(IRS),ITEMP(IRS)
      EQUIVALENCE (ITEMP(1),TEMP(1))
      IOF=0
      IF(IPHASE.GE.2.AND.IPHASE.NE.1000) THEN
      DO 4100 IIPHAS=2,IPHASE
4100      IOF = IOF + ICR(IIPHAS-1)
      END IF
      IF (IPHASE.EQ.1000) THEN
        IC=0
C        DO 4101 IIPHAS=1,8                !cp 03aug97
        DO 4101 IIPHAS=1,99
4101      IC  = IC  + ICR(IIPHAS)
      ELSE
        IC=ICR(IPHASE)
      END IF
      DO 2 I=1,IC
2       L(I)=I
      M=IC-1
      DO 15 LL=1,5000
      K=0
      DO 10 J=1,2 
        DO 10 I=J,M,2
          INP=L(I)+IOF
          INP1=L(I+1)+IOF
          R=REFS(INP1,2)-REFS(INP,2)
          IF(R) 20,10,10
20          IZ=L(I) 
          L(I)=L(I+1)
          L(I+1)=IZ
          K=K+1
10          CONTINUE
      IF(K.EQ.0) GO TO 5
15      CONTINUE
5     DO 30 I=1,IC
30      ITEMP(I)=IREFS(I+IOF) 
      DO 31 I=1,IC
      LX=L(I)
      IREFS(I+IOF)=ITEMP(LX)
31      CONTINUE
      DO 33 J=1,3
        DO 32 I=1,IC
32        TEMP(I)=REFS(I+IOF,J)
        DO 33 I=1,IC
          LX=L(I)
33        REFS(I+IOF,J)=TEMP(LX)
      DO 34 I=1,IC
34      ITEMP(I)=FMGNTD(I+IOF)
      DO 35 I=1,IC
        LX=L(I)
35      FMGNTD(I+IOF)=ITEMP(LX)
      IF (NPROF.EQ.8) THEN
      DO 40 I=1,IC
40        TEMP(I)=HALFL(I+IOF)
      DO 41 I=1,IC
        LX=L(I)
41        HALFL(I+IOF)=TEMP(LX)
      DO 50 I=1,IC
50        TEMP(I)=HALFG(I+IOF)
      DO 51 I=1,IC
        LX=L(I)
51        HALFG(I+IOF)=TEMP(LX)
      DO 60 I=1,IC
60        TEMP(I)=GAM(I+IOF)
      DO 61 I=1,IC
        LX=L(I)
61        GAM(I+IOF)=TEMP(LX) 
      END IF
      RETURN
      END
      SUBROUTINE EXPUT
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      REAL LAMDA
      INTEGER PTR
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
C      COMMON/ALLP/FINAL(8*MSZ,2),ILOC
      COMMON/ALLP/FINAL(nfinal,2),ILOC
      CHARACTER SYMB(99,20)*1,ISPACE*1,IDOT*1,ISTAR*1,IPLUS*1,IMINUS*1 
     * ,IBAR*1
      CHARACTER TITLE*70,PHSNM(99)*50
      CHARACTER*1 IOUT(150)
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),XMAS(16)
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ),
     * BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/PRFX/DELT,TL,GAM1,GAM2,PRFDER,iph,delta
      COMMON /PVII/ TF1,TF2,TF4,TF6,TF8,TF9,C4
      common/spvii/rl,da1l,da2l,da3l,da4l,da5l,da6l,da7l
      common/spvii/rh,da1h,da2h,da3h,da4h,da5h,da6h,da7h
C
      COMMON/FONDI/BKCOM(IDSZ),BKDIS(IDSZ),BKAM(IDSZ),BKPOL(IDSZ)
      common/sizestrain/sizeG(15),strainG(15),sizeL(15),strainL(15)
     *,siz(15),strain(15),NsizeStrain
      CHARACTER*8 TITLE1,TITLE2,TITLE3,TITLE4,TITLE6,TITLE7
      CHARACTER*8 TITLE8
C     CHARACTER*8 TITLE8,TITLE9
C
      DIMENSION DEL(10),KI(10),DUMMY(2*MSZ+1)
C     DIMENSION RJAC(MSZ,MSZ)
C     EQUIVALENCE (RJAC(1,1),KR(1))
C     DIMENSION LABEL(12),KS(50)
      DIMENSION LABEL(12)
      COMMON/STRUPHASE/aphase(IRS),tavix(irs),srix(irs)
      LOGICAL VERT    
      CHARACTER*80 OUTFILE
      dimension phs(idsz)     
      equivalence(phs(1),kr(1))
      common/maxint/xmaxint     
C      DATA ISPACE,IDOT,ISTAR,IPLUS,IMINUS,IBAR/' ','.','*','+','-','0'/
      DATA ISPACE,IDOT,ISTAR,IPLUS,IMINUS,IBAR
     1       /' ', '.',  '*',  '+',   '0', '-'/
      TITLE1 = 'OSS.COR"'
      TITLE2 = 'TOTALE "'
      TITLE3 = 'COM.TOT"'
      TITLE4 = 'DIS.TOT"'
      TITLE6 = 'AM.TOT "'
      TITLE7 = 'POL.TOT"'
      TITLE8 = 'ALL.TOT"'
C      TITLE9 = 'AIR.TOT"'
      IF(NXT.EQ.0.OR.JOBTYP.GE.3) NPROF=NPROF-1
      RAD = 45./ATAN(1.)  
      zer = glb(1)
      IF(IPL.NE.0)READ(5,19,END=10000)ISCALE,IDIF 
19    FORMAT(BZ,2I8)
c10000 IF(IPL2.NE.0)READ(5,17,END=10001)IFY,IFYC,IFM,IFD,IFB,EXPAND
c10000 IF(JOBTYP.GT.2)IFY=0
C      IF(JOBTYP.GT.2)IFD=0
c17    FORMAT(BZ,5I1,3X,F8.0)
10000   IF(NXT.NE.0.AND.JOBTYP.LT.3)then
C        CALL REWRIT(ISCALE,IDIF,IFY,IFYC,IFM,IFD,IFB,EXPAND)
        CALL REWRIT(ISCALE,IDIF)
        endif
C          if(i2d94.NE.0)THEN
C            call write94(iscale,idif)
C          END IF
C  
      if(ipl2.ne.0)then
      open(69,file='bragg.dat',status='unknown')
      open(690,file='xy-int.dat',status='unknown')
      write(690,960)    
960   format(2x,'2theta_i',6x,'y_o',7x,'y_c',8x,'yo-yc',2x,'w(yo-yc)^2'
     1)
      endif
      IF(JOBTYP.GT.2)GOTO 37
      WRITE(6,10)
10    FORMAT(/42H0AVERAGE INTENSITY DIFFERENCE FOR PATTERN,,/,
     * 37H0GIVEN FOR BLOCKS OF 20 OBSERVATIONS.)
ccccc    APAGUE APENAS 2 vezes o C para efetivar as linhas abaixo
cccc      if(ipl2.ne.0)then
cccc      open(69,file='bragg.dat',status='unknown')
ccccc      write(69,96)       
ccccc96    format(1x,'2theta_B2',1x,'y-axis_position')      
cccc      open(690,file='xy-int.dat',status='unknown')
cccc      write(690,960)    
cccc960   format(2x,'2theta_i',6x,'y_o',7x,'y_c',8x,'yo-yc',2x,'w(yo-yc)^2'
cccc     1)
cccc      endif
      DO 15 I=1,NPTS,200
      DO 14 J=1,10
        J20=20*(J-1)
        DEL(J)=0
        DO 13 K=1,20
          J20IK=J20+I+K-1
13          DEL(J)=DEL(J)+Y(J20IK)-BK(J20IK)-YC(J20IK)
        KI(J)=J+I/20
        DEL(J)=DEL(J)/20.
        IF(KI(J)*20.GE.NPTS)GOTO 15
14        CONTINUE
      J=10
15      WRITE(6,20)(KI(JJ),DEL(JJ),JJ=1,J)
20    FORMAT(1H0,6X,10(I4,2X,F5.1))
37    DO 47 I=1,NPTS
      YC(I)=YC(I)+BK(I)
47      IF(JOBTYP.GT.2)Y(I)=YC(I)
      IF(IPC.EQ.0)GOTO 36
       INQUIRE(UNIT=6,NAME=OUTFILE)
       LBKSL = 0
       DO I=1,80
cvv         IF ( OUTFILE(I:I).EQ.'\' ) LBKSL=I
         IF ( OUTFILE(I:I).EQ.'A' ) LBKSL=I
       END DO
!      WRITE(6,*) 'Outfile:',OUTFILE
       OUTFILE(LBKSL+1:) = 'plotinfo.bin'
!      WRITE(6,*) 'Outfile:',OUTFILE
       OPEN(9,FILE=OUTFILE(1:LBKSL+8),status='unknown')
       open(10,file=OUTFILE(1:LBKSL+12),form='binary',
     1    status='unknown')
!      OPEN(9,FILE='PLOTINFO',status='unknown')
!      open(10,file='plotinfo.bin',form='UNFORMATTED',status='unknown')
       write(10) (bk(i),i=1,npts)
       WRITE(9,'(A70)') TITLE
       WRITE(9,8730) NPHASE,(ICR(IIPHAS),IIPHAS=1,NPHASE)
8730   FORMAT (1X, 'NO. OF PHASESZ', 2X,I4,/,1X,
     *  ' NO. OF REFLECTIONS IN EACH PHASEQ',2X,15I4)
       WRITE(9,8731)
8731   FORMAT (1X,'BRAGG POSITIONSZ')
      DO 400 K=1,NPHASE
      do 401 i=1,npts
401     phs(i) = 0.
      IF (MOD(K,2).EQ.1) ILOC = ILOC + 1
C       GAM1=PAR(K,17)
      WRITE (6,4289) K,PHSNM(K)
4289    FORMAT (/'1PHASE NO. = ',I3,5X,'PHASE NAME ',A50)
      WRITE (9,4298) PHSNM(K)
4298    FORMAT (/1X,A50)
      T2OBS=0.0
      TDOBS=0.0
      RFNMR=0.0
      RFDNR=0.0
      ICZ=0
      DO 4300 IIPHAS=1,NPHASE
4300      ICZ = ICZ + ICR(IIPHAS)
      IXX=0
      ixxx=0
      DO 481 IX=1,ICZ
        IF( K.NE.IREFS(IX)/(256*256*256*8) )GOTO 481
        SHIFT = GLB(10) * COS(REFS(IX,2)/2./57.2958)
     *               + GLB(11) * SIN(REFS(IX,2)/57.2958)
        IXX=IXX+1
        DO 485 J=1,NEXCRG
C-----CHECK FOR SPLIT PEARSON VII PROFILE
        if (nprof+1.eq.5) then
          IF ((REFS(IX,2)+SHIFT+glb(1)).GE.(ALOW(J)-WDT*fwhm(IX,1))
     *                                   .AND.
     *          (REFS(IX,2)+SHIFT+glb(1)).LE.(AHIGH(J)+WDT*fwhm(IX,2)))
     *                                  GOTO 481
C-----FOR ALL OTHER PROFILES
        else
          IF ((REFS(IX,2)+SHIFT+glb(1)).GE.(ALOW(J)-WDT*REFS(IX,1))
     *                                   .AND.
     *          (REFS(IX,2)+SHIFT+glb(1)).LE.(AHIGH(J)+WDT*REFS(IX,1)))
     *                                  GOTO 481
        end if
485         CONTINUE
      IXXx=IXXx+1                       !test !cp 29 jun 98
        IRL=MOD(IREFS(IX),256)-128
        IRK=MOD(IREFS(IX)/256,256)-128
        IRH=MOD(IREFS(IX)/(256*256),256)-128
        IRC=MOD(IREFS(IX)/(256*256*256),8)
C-----CHECK FOR THE SPLIT PEARSON VII PROFILE
C-----IF SO CHANGE THE PROFILE LIMITS
        if (nprof+1.eq.5) then
          RMIN=REFS(IX,2)+glb(1)+SHIFT-WDT*fwhm(IX,1)
          RMAX=REFS(IX,2)+glb(1)+SHIFT+WDT*fwhm(IX,2)
C-----FOR ALL OTHER PROFILES
        else
          RMIN=REFS(IX,2)+glb(1)+SHIFT-WDT*REFS(IX,1)
          RMAX=REFS(IX,2)+glb(1)+SHIFT+WDT*REFS(IX,1)
        end if
        MIN=(RMIN-THMIN)/STEP+1.5
        MAX=(RMAX-THMIN)/STEP+1.5
        MIN=MAX0(MIN,1)
        MIN=MIN0(MIN,NPTS)
        MAX=MIN0(MAX,NPTS)
C-----PATCH TO CALCULATE R-BRAGG
        TL=REFS(IX,1)
        VERT=REFS(IX,2).LE.RLIM
        IF (NPROF+1.EQ.6) THEN
           GAM1=PAR(K,17) + PAR(K,18) * REFS(IX,2)
        ELSE IF (NPROF+1.EQ.7) THEN
           GAM1=PAR(K,17)+(PAR(K,18)+PAR(K,19)/REFS(IX,2))/REFS(IX,2)
           CALL PRSVII(GAM1)
        ELSE IF (NPROF+1.EQ.8) THEN 
           GAM1=GAM(IX)
           TLG=HALFG(IX)
           TLL=HALFL(IX)
        else if (nprof+1.eq.5) then
           rl=par(k,17)+(par(k,18)+par(k,19)/refs(ix,2))/refs(ix,2)
           rh=par(k,24)+(par(k,25)+par(k,26)/refs(ix,2))/refs(ix,2)
           call mspvii(par(k,27),tl)
        END IF
        BB=TL*TL
        TIOBS =0.0
        TIC=0.0
        FMGNTD(IX) = FMGNTD(IX)*REFS(IX,3)*PAR(K,1)
        DO 410 IS=MIN,MAX
          TH=THMIN+FLOAT(IS-1)*STEP
          IF(NEXCRG.GT.0) THEN
            DO 1450 IEXC=1,NEXCRG
            IF (TH.GE.ALOW(IEXC).AND.TH.LE.AHIGH(IEXC)) GO TO 410 
1450            CONTINUE
          END IF
          TH2 = 0.5*TH/RAD
          DELTA=TH-REFS(IX,2)-glb(1)-SHIFT
          DELT=DELTA*DELTA
C     NEXT LINE IS NECESSEARY FOR 2 PHASES WITH VERY DIFFERENT FWHM.
          IF (DELT/BB.GT.WDT*WDT) GO TO 410
          IF(.NOT.VERT)GOTO 4
            if (iasym.eq.0) then
               YX=DELT*SIGN(1.,DELTA)
               Z=1.-PAR(K,14)*YX/TAN(TH2)
             else
               YX=SIGN(1.0,DELTA)*DELTA/(2*TL)
               IF(TH2.GT.(45.0/RAD)) THEN
                   TH2 = TH2-(90.0/RAD)
               END IF
               Z=(PAR(K,14)/TAN(TH2)) * (2.0*(DELTA/(2*TL))*EXP(-YX))
               Z=1+Z
             end if
          IF(Z.LE.0.0)Z=0.0001
          GOTO 5
4           Z=1.
5           OMEGA= Z*PROFIL(NPROF+1,DELT/BB)
            if (nprof+1.eq.5) then
              phs(is) = phs(is) + omega*FMGNTD(ix)
            else 
              phs(is) = phs(is) + omega*FMGNTD(ix)/tl
            end if
          TIC   = TIC + OMEGA
          TIOBS=TIOBS+OMEGA*(Y(IS)-BK(IS))/AMAX1((YC(IS)-BK(IS)),1.)
410       CONTINUE
        TIC   = TIC   * FMGNTD(IX)/TL 
        TIOBS = TIOBS * FMGNTD(IX)/TL
C ! FROM ITALIAN CODE  !cp ap 97  
C        ***************************************************************
C        NEXT LINE IS FOR NOT EVALUATING R-BRAGG FOR REFLECTIONS
C        WHICH ARE OUTSIDE THE MEASUREMENTS RANGE BUT WHOSE TAILS ARE
C        IN THE PATTERN
C        ***************************************************************
      IF((REFS(IX,2)+SHIFT+glb(1)).LE.THMAX) THEN 
        T2OBS = T2OBS+TIOBS 
        TDOBS = TDOBS + ABS(TIOBS-TIC)
C  below is to correct the Fobs and Fcal for alpha_2 !cp 24.May.98
C          if(irc.eq.1) then
C                  xlor=1.0
C                   else
C                  xlor=ratio(2)
C          end if
        IF (IPC.EQ.2.OR.IPC.eq.3) THEN
      TFCAL= SQRT(ABS(TIC/REFS(IX,3)/PAR(K,1)/tavix(ix)/srix(ix)))
      TFOBS= SQRT(ABS(TIOBS/REFS(IX,3)/PAR(K,1)/tavix(ix)/srix(ix)))
      Afcal=tfcal*cos(aphase(ix))
      Bfcal=tfcal*sin(aphase(ix))
      Afobs=tfobs*cos(aphase(ix))
      Bfobs=tfobs*sin(aphase(ix))
C        Line below also from italian code   !cp ap 97 
            IF(REFS(IX,2).LE.THMAX)THEN 
                RFNMR = RFNMR + ABS(TFOBS-TFCAL)
                RFDNR = RFDNR + ABS(TFOBS)
            END IF       
        end if
      end if            
        IF (NPROF+1.EQ.8) GO TO 9222
        if (nprof+1.eq.5) go to 9221
        IF(MOD(IXXx-1,60).EQ.0) THEN
c alterei aqui para criar arquivo com dados para determinacao de estruturas e graficos ! cp 22abril01
        IF(IPC.EQ.1) then
                     WRITE(6,482)
                     if(ipl2.ne.0)write(69,4821)
        end if             
        IF(IPC.EQ.2) then
                     WRITE(6,8482)
                     if(ipl2.ne.0)write(69,84821)
        end if
        IF(IPC.EQ.3) then
                     WRITE(6,8483)
                     if(ipl2.ne.0)write(69,84831)
        end if
        END IF
482     FORMAT(32H0NO.  CODE    H   K   L   HW     ,
     *     16HPOSN      ICALC ,4X,5HIOBS )
4821    FORMAT('   H   K   L   HW     POSN      ICALC       IOBS   Bragg
     1_Pos')
8482    FORMAT('0NO.  CODE     H   K   L    FWHM  ',
     *     '  POSN    ICALC ',5X,'IOBS',6X,'FCALC',5X,'FOBS   PHASE_C')
84821 FORMAT('   H   K   L     FWHM    POSN    ICALC    IOBS      FCALC 
     1     FOBS   PHASE_C  Bragg_Pos')
8483    FORMAT('0NO.  CODE     H   K   L     FWHM  ',
     *     '  POSN    ICALC ',5X,'IOBS',5X,'A_CALC    B_CALC',
     *        5X,'A_OBS     B_OBS')
84831 FORMAT('   H   K   L     FWHM    POSN    ICALC       IOBS     A_CA
     1L  B_CALC      A_OBS     B_OBS  Bragg_Pos')
        HW=REFS(IX,1)
      IF (IPC.EQ.1) THEN
        WRITE(6,470)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *          TIC,TIOBS               
          if(irc.eq.1)then
          if (ipl2.ne.0)  write(69,691) irh,irk,irl,hw,
     1         refs(ix,2)+glb(1)+shift,tic,tiobs,(-k*xmaxint/10) 
          end if
      else if (ipc.eq.2) THEN
              IF (IRC.EQ.1)THEN                        
        WRITE(6,470)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *       TIC,TIOBS,TFCAL,TFOBS,rad*aphase(ix)
          if (ipl2.ne.0)  write(69,692) irh,irk,irl,hw,
     1         refs(ix,2)+glb(1)+shift,tic,tiobs,
     1         tfcal,tfobs,rad*aphase(ix),(-k*xmaxint/10) 
              else
        WRITE(6,470)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *       TIC,TIOBS
              end if
      else if (ipc.eq.3) THEN   
              IF (IRC.EQ.1)THEN                        
        WRITE(6,4701)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *         TIC,TIOBS,Afcal,Bfcal,Afobs,Bfobs 
          if (ipl2.ne.0)  write(69,691) irh,irk,irl,hw,
     1         refs(ix,2)+glb(1)+shift,tic,tiobs,
     1         Afcal,Bfcal,Afobs,Bfobs,(-k*xmaxint/10) 
              else
        WRITE(6,4701)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *         TIC,TIOBS 
              end if
      END IF   
C      WRITE (9,8147) REFS(IX,2)+glb(1)+SHIFT
      WRITE (9,8147) REFS(IX,2)+glb(1)+SHIFT,irc,irh,irk,irl,FMGNTD(ix)  
c      if (ipl2.ne.0)  write(69,691) REFS(IX,2)+glb(1)+SHIFT,k*(-100) 
c                                                      !cp 22abril01
691   format(1x,3(i3,1x),2(f8.4,1x),2(f8.0,1x),f9.0)      
692   format(1x,3(i3,1x),2(f8.4,1x),4(f8.1,1x),f9.4,1x,f9.0)      
693   format(1x,3(i3,1x),2(f8.4,1x),2(f8.0,1x),4(f8.2,1x),f9.0)      
694   format(1x,3(i3,1x),3(f8.4,1x),2(f8.0,1x),f9.0)      
695   format(1x,3(i3,1x),3(f8.4,1x),2(f8.0,1x),2(f8.4,1x),f9.4,f9.0)      
696   format(1x,3(i3,1x),3(f8.4,1x),2(f8.0,1x),4(f8.2,1x),f9.0)      
697   format(1x,3(i3,1x),2(f8.4,1x),2(f8.0,1x),3(f8.4,1x),f9.4,f9.0) 
698   format(1x,3(i3,1x),2(f8.4,1x),2(f8.0,1x),3(f8.4,1x),2(f8.4,1x),
     1f9.4,f9.0)
699   format(1x,3(i3,1x),2(f8.4,1x),2(f8.0,1x),3(f8.4,1x),4(f8.4,1x),
     1f9.0)
       GO TO 481
9221      if(mod(ixxx-1,60).eq.0) then
          if(ipc.eq.1) then
                 write(6,9481)          
                 if(ipl2.ne.0)write(69,94811)
          end if
          if(ipc.eq.2) then
                 write(6,9483)          
                 if(ipl2.ne.0)write(69,94831)
                 end if
          if(ipc.eq.3) then
                 write(6,9485)
                 if(ipl2.ne.0)write(69,94851) 
          end if
        end if
9481    format(41H0NO.  CODE    H   K   L      HWL     HWH  ,
     *     17HPOSN       ICALC ,4X,5HIOBS ' PHASE')
c             x123x123*123*12345678*12345678*12345678*123456778*12345678*12345678*
94811 format('   H   K   L      HWL       HWH      POSN     ICALC      I
     1OBS    PHASE  Bragg_Pos')
9483    format(42H0NO.  CODE    H   K   L      HWL     HWH  ,
     *     16HPOSN      ICALC ,4X,5HIOBS ,7X,'FCALC',4X,'FOBS',7x,
     *     'PHASE_C')
c             x123x123*123*12345678*12345678*12345678*123456778*12345678*12345678*12345678*12345678*
94831 format('   H   K   L      HWL       HWH      POSN     ICALC      I
     1OBS    FCALC     FOBS     PHASE_C  Bragg_Pos')
9485    format(42H0NO.  CODE    H   K   L      HWL     HWH  ,
     *     16HPOSN      ICALC ,4X,5HIOBS ,7X,'A_CALC',5X,'B_CALC',5x,
     *     'A_OBS',5X,'B_OBS')
c             x123x123*123*12345678*12345678*12345678*123456778*12345678*12345678*12345678*12345678*
94851 format( '   H   K   L      HWL      HWH      POSN     ICALC      I
     1OBS   A_CALC   B_CALC   A_OBS    B_OBS  Bragg_Pos')

      if(ipc.EQ.1) then
        write(6,471)ixx,irc,irh,irk,irl,fwhm(ix,1),fwhm(ix,2),
     1         refs(ix,2)+glb(1)+shift,tic,tiobs 
          if(irc.eq.1)then
          if (ipl2.ne.0)  write(69,694) irh,irk,irl,fwhm(ix,1),
     1         fwhm(ix,2),refs(ix,2)+glb(1)+shift,tic,tiobs,
     1         (-k*xmaxint/10) 
          end if
      else if(ipc.eq.2) THEN
              IF (IRC.EQ.1)THEN                        
        write(6,471)ixx,irc,irh,irk,irl,fwhm(ix,1),fwhm(ix,2),
     1      refs(ix,2)+glb(1)+shift,tic,tiobs,tfcal,tfobs,rad*aphase(ix)
          if (ipl2.ne.0)  write(69,695) irh,irk,irl,fwhm(ix,1),
     1         fwhm(ix,2),refs(ix,2)+glb(1)+shift,tic,tiobs,
     1         tfcal,tfobs,rad*aphase(ix),(-k*xmaxint/10) 
              else
        write(6,471)ixx,irc,irh,irk,irl,fwhm(ix,1),fwhm(ix,2),
     1      refs(ix,2)+glb(1)+shift,tic
              end if
      else if(ipc.eq.3) THEN
              IF (IRC.EQ.1)THEN                        
       write(6,4711)ixx,irc,irh,irk,irl,fwhm(ix,1),fwhm(ix,2),
     1 refs(ix,2)+glb(1)+shift,tic,tiobs,Afcal,Bfcal,Afobs,Bfobs
          if (ipl2.ne.0)  write(69,696) irh,irk,irl,fwhm(ix,1),
     1         fwhm(ix,2),refs(ix,2)+glb(1)+shift,tic,tiobs,
     1         Afcal,Bfcal,Afobs,Bfobs,(-k*xmaxint/10) 
              else
       write(6,4711)ixx,irc,irh,irk,irl,fwhm(ix,1),fwhm(ix,2),
     1 refs(ix,2)+glb(1)+shift,tic,tiobs
              end if
      end if
471    format(1x,2i4,3x,3i4,3f8.3,2f10.0,2f10.3,1X,F8.2)
4711   format(1x,2i4,3x,3i4,3f8.3,2f10.0,4f10.3)
C      write (9,8147) refs(ix,2)+glb(1)+shift
      write (9,8147) refs(ix,2)+glb(1)+shift,irc,irh,irk,irl,FMGNTD(ix)
c      if (ipl2.ne.0)  write(69,691) REFS(IX,2)+glb(1)+SHIFT,k*(-100)  !cp 22abril01
        go to 481
9222      CONTINUE
        IF(MOD(IXXx-1,60).EQ.0) THEN
        IF (IPC.EQ.1) then
                      WRITE(6,9482)
                     if(ipl2.ne.0)write(69,94821)
        end if
        IF (IPC.EQ.2) then
                      WRITE(6,19482)
                     if(ipl2.ne.0)write(69,94822)
        end if
        if (ipc.eq.3) then
                      write(6,29482)
                     if(ipl2.ne.0)write(69,94823)
        end if
        END IF    
C                 123456789*123456789*123456789*12345    123456789*1234
9482    FORMAT(35H0  NO. CODE     H   K   L     HW     ,
     *         14HPOSN     ICALC,4X,5HIOBS ,6X,
     *         'HG      HL      ETA   PHASE_C')
94821 FORMAT('   H   K   L     HW      POSN     ICALC    IOBS    HG     
     1HL      ETA   PHASE_C')
C                 123456789*123456789*123456789*12345    123456789*1234
19482   FORMAT('0 NO. CODE     H   K   L     HW     ',
     *     'POSN    ICALC ',5X,5HIOBS ,4X,' HG      HL      ETA   ',
     *     '  FCALC     FOBS   PHASE_C')
94822 FORMAT('   H   K   L     HW       POSN    ICALC    IOBS    HG     
     1HL      ETA   FCALC     FOBS   PHASE_C')
c
29482   FORMAT('0 NO. CODE     H   K   L     HW     ',
     *     'POSN    ICALC ',5X,5HIOBS ,4X,' HG      HL      ETA   ',
     *     ' A_CALC    B_CALC     A_OBS     B_OBS')
94823 FORMAT('   H   K   L     HW       POSN    ICALC    IOBS    HG     
     1HL      ETA   A_CALC    B_CALC     A_OBS     B_OBS')
c
        HW=REFS(IX,1)
      IF (IPC.EQ.1) THEN
      WRITE(6,9470)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,TIC
     *         ,TIOBS,TLG,TLL,GAM1,rad*aphase(ix)  
      if(ipl2.ne.0)write(69,697)IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *         tic,TIOBS,TLG,TLL,GAM1,rad*aphase(ix),(-k*xmaxint/10)
      ELSE IF(IPC.EQ.2) THEN
              IF (IRC.EQ.1)THEN                        
        WRITE(6,9470)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *        TIC,TIOBS,TLG,TLL,GAM1,TFCAL,TFOBS,aphase(ix) 
      if(ipl2.ne.0)WRITE(69,698)IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *        TIC,TIOBS,TLG,TLL,GAM1,TFCAL,TFOBS,aphase(ix),
     *        (-k*xmaxint/10)
              else
        WRITE(6,9470)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     *        TIC,TIOBS,TLG,TLL,GAM1
              end if              
      ELSE IF(IPC.EQ.3) THEN
              IF (IRC.EQ.1)THEN                        
        WRITE(6,94701)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     1     TIC,TIOBS,TLG,TLL,GAM1,Afcal,Bfcal,Afobs,Bfobs
      if(ipl2.ne.0)WRITE(69,699)IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     1     TIC,TIOBS,TLG,TLL,GAM1,Afcal,Bfcal,Afobs,Bfobs,
     1     (-k*xmaxint/10)
              else
        WRITE(6,94701)IXX,IRC,IRH,IRK,IRL,HW,REFS(IX,2)+glb(1)+SHIFT,
     1     TIC,TIOBS,TLG,TLL,GAM1
              end if
      END IF
9470    FORMAT(1X,2I4,3X,3I4,2F8.3,2F10.0,3F8.3,2F10.3,1x,f8.2)
94701   FORMAT(1X,2I4,3X,3I4,2F8.3,2F10.0,3F8.3,4F10.3)
C      WRITE (9,8147) REFS(IX,2)+glb(1)+SHIFT
4810   WRITE (9,8147) REFS(IX,2)+glb(1)+SHIFT,irc,irh,irk,irl,FMGNTD(ix)
c      if (ipl2.ne.0)  write(69,691) REFS(IX,2)+glb(1)+SHIFT,k*(-100)       !cp 22abril01
c8147  FORMAT (1X,F8.3)
C                                                version I - 
c8147  FORMAT (1X,F8.3,2x,i1,3i4,g10.4)
C version II - format for compatibility with Sakthevil's PLOT program
8147  FORMAT (1X,F8.3,2h K,i1,sp,3i4.3,g10.4,s) 
c
481       CONTINUE
      write(10) (phs(i),i=1,npts)
470     FORMAT(1X,2I4,3X,3I4,2F8.3,2F10.0,2F10.3,1X,2F8.2) 
4701    FORMAT(1X,2I4,3X,3I4,2F8.3,2F10.0,5F10.3) 
      TDOBS=100.*TDOBS/T2OBS
      WRITE(6,483)TDOBS
      IF (IPC.NE.1)  WRITE(6,10483) 100.*RFNMR/RFDNR
483     FORMAT(26H0 DERIVED BRAGG R-FACTOR = ,F8.2)
10483   FORMAT(  '0 DERIVED R-F            = ',F8.2)
      FINAL(ILOC,2-MOD(K,2)) = TDOBS
            if(NsizeStrain .eq. 9) call size(k)
400     CONTINUE
36    CONTINUE
      WRITE (9,8476) NPTS, THMIN, STEP
8476  FORMAT (1X,'NPTSZ',I5,/, 1X, 'THMINZ',F8.4,/,1X,'STEPZ',F8.4,/,
     *1X,4HYOBS,4X,6HYCALCZ)
      DO 31 I=1,NPTS
      BK(I)=THMIN+FLOAT(I-1)*STEP
c31    WRITE(9,8477)Y(I),YC(I)
c     !incluidas as linhas open(69) e open(690) e tudo o mais
c       com write e etc.. essa linha com"31" tambem teve o31 mudado de lugar. 02set00
      WRITE(9,8477)Y(I),YC(I)
8477  FORMAT (1X,2F8.0)
31    if(ipl2.ne.0)write(690,6692)THMIN+FLOAT(I-1)*STEP,y(i),yc(i),
     !y(i)-yc(i),((y(i)-yc(i))**2)/y(i)
c31    if(ipl2.ne.0)write(690,6692)THMIN+FLOAT(I-1)*STEP,y(i),yc(i),
c     !y(i)-yc(i)
6692   format(1x,f9.4,4(1x,f9.0))
      CLOSE (9)
      close (10)        
      if(ipl2.ne.0)then
      close(69)
      close(690)
      endif
C     IT BUILDS THE OBSERVED DATA FILE CORRECTED FOR ABSORPTION
      IF(IPLOSS.EQ.1.AND.IPBIG.EQ.0) THEN
      OPEN(31,FILE='PLOTOSS.COR')
      WRITE(31,12002) TITLE,TITLE1
12002 FORMAT(1X,'"',A70,A8)
      WRITE(31,12004) NPTS,STEP,THMIN,THMAX,1
12004 FORMAT(I6,3F15.5,I5)
      DO 12031 I=1,NPTS
      WRITE(31,12008) Y(I)
12008 FORMAT(F15.5)
12031 CONTINUE
      CLOSE (31)
      ENDIF
C     IT BUILDS THE CALCULATED DATA FILE
C     (BRAGG+COMPTON+DISORDINE+AMORPHOUS)
      IF(IPLCAL.EQ.1.AND.IPBIG.EQ.0)THEN
      OPEN(32,FILE='PLOTCAL.TOT')
      WRITE(32,12002) TITLE,TITLE2
      WRITE(32,12004) NPTS,STEP,THMIN,THMAX,1
      DO 12032 I=1,NPTS
      WRITE(32,12008) YC(I)
12032 CONTINUE
      CLOSE (32)
      ENDIF
C     IT BUILDS THE TOTAL COMPTON SCATTERING FILE
C     FOR ALL PHASES
      IF(IPLCOM.EQ.1.AND.IPBIG.EQ.0)THEN
      OPEN(33,FILE='PLOTCOM.TOT')
      WRITE(33,12002) TITLE,TITLE3
      WRITE(33,12004) NPTS,STEP,THMIN,THMAX,1
      DO 12033 I=1,NPTS
      WRITE(33,12008) BKCOM(I)
12033 CONTINUE
      CLOSE (33)
      ENDIF
C     IT BUILDS THE TOTAL DISORDER SCATTERING FILE
C     FOR ALL PHASES
      IF(IPLDIS.EQ.1.AND.IPBIG.EQ.0)THEN
      OPEN(34,FILE='PLOTDIS.TOT')
      WRITE(34,12002) TITLE,TITLE4
      WRITE(34,12004) NPTS,STEP,THMIN,THMAX,1
      DO 12034 I=1,NPTS
      WRITE(34,12008) BKDIS(I)
12034 CONTINUE
      CLOSE (34)
      ENDIF
C     IT BUILDS THE POLYNOMIAL BACKGROUND FILE
      IF(IPLPOL.EQ.1.AND.IPBIG.EQ.0)THEN
      OPEN(37,FILE='PLOTPOL.TOT')
      WRITE(37,12002) TITLE,TITLE7
      WRITE(37,12004) NPTS,STEP,THMIN,THMAX,1
      DO 12037 I=1,NPTS
      WRITE(37,12008) BKPOL(I)
12037 CONTINUE
      CLOSE (37)
      ENDIF
C     IT BUILDS THE AMORPHOUS FILE
      IF(IPLAM.EQ.1.AND.IPBIG.EQ.0)THEN
        OPEN(36,FILE='PLOTAM.TOT')
        WRITE(36,12002) TITLE,TITLE6
        WRITE(36,12004) NPTS,STEP,THMIN,THMAX,1
        DO 12036 I=1,NPTS
        WRITE(36,12008) BKAM(I)
12036     CONTINUE
        CLOSE (36)
      ENDIF
C     IT BUILDS THE TOTAL PLOT FILE
      IF(IPBIG.EQ.1)THEN
      OPEN(38,FILE='PLOTBIG.DAT')
      WRITE(38,12012) TITLE,TITLE8
      WRITE(38,12014) 'ANG','OSS','CAL','+AM','+POL','+DIS',
     * '+COM','RES'
      DO 12038 I=1,NPTS
      XXXX=THMIN+FLOAT(I-1)*STEP
      WRITE(38,12018) XXXX,Y(I),YC(I),BKAM(I),
     *BKAM(I)+BKPOL(I),BKAM(I)+BKPOL(I)+BKDIS(I),
     *BKAM(I)+BKPOL(I)+BKDIS(I)+BKCOM(I),
     *((Y(I)-YC(I))**2)/Y(I)
12038  CONTINUE
       CLOSE (38)
      ENDIF
12012 FORMAT(1X,'"',A70,A8,'"')
12014 FORMAT(A10,',',A10,',',A10,',',A10,',',A10,',',A10,
     *',',A10,',',A10,',',A10)
12018 FORMAT(F10.3,',',F10.3,',',F10.3,',',F10.3,',',F10.3,
     *',',F10.3,',',F10.3,',',F10.3,',',F10.3)
12045 IF(IOT.EQ.0)GOTO 30
      NP=NPTS/240
      WRITE(6,101)(((BK(240*N-300+60*I+J),Y(240*N-300+60*I+J),YC(
     * 240*N-300+60*I+J),VAR(240*N-300+60*I+J),I=1,4),J=1,60),N=1,NP) 
101   FORMAT((2H1 ,4(6H2THETA,3X,4HYOBS,4X,5HYCALC,1X,8HVARIANCE,1X),/,
     *  60(4(1X,F7.3,3F8.0),/)))
      NPTS2=NPTS-NP*240
      IF(NPTS2.EQ.0)GOTO 30
      NCOL=NPTS2/60
      WRITE(6,103)
103   FORMAT(2H1 ,4(6H2THETA,3X,4HYOBS,4X,5HYCALC,1X,8HVARIANCE,1X))
      NLINES=NPTS2-NCOL*60
      IF(NLINES.EQ.0)GOTO 33
      NCOL1=NCOL+1
      NP=240*NP-60
      DO 34 J=1,NLINES
34      WRITE(6,102)(BK(NP+60*I+J),Y(NP+60*I+J),YC(NP+60*I+J),
     *      VAR(NP+60*I+J),I=1,NCOL1)
102   FORMAT(4(1X,F7.3,3F8.0))
      NP=NP+NLINES
      NLINES=60-NLINES
33    DO 35 J=1,NLINES
35      WRITE(6,102)(BK(NP+60*I+J),Y(NP+60*I+J),YC(NP+60*I+J),
     *      VAR(NP+60*I+J),I=1,NCOL)
30    IF(IPL.EQ.0)GOTO 45
      DO 210 I=1,12
210     LABEL(I)=10*I*ISCALE
      WRITE(6,100) TITLE,LABEL
100   FORMAT(1H1,30X,A70//1X,12I10)
      DO 11 I=1,12
11      LABEL(I)=-60*IDIF+10*I*IDIF
      WRITE(6,203) LABEL
203   FORMAT(1H ,12I10)
      DO 60 I=1,NPTS
1       DO 3 J=1,120
3         IOUT(J)=ISPACE
      IOUT(1)=IDOT
      IOUT(61)=IDOT
      IOUT(120)=IDOT
200     IY=INT(Y(I))/ISCALE+2 
      IY=MAX0(MIN0(IY,120),2)
      IOUT(IY-1)=IBAR
      IOUT(IY)=IPLUS
      IF(IY.LE.119)IOUT(IY+1)=IBAR
      IY=INT(YC(I))/ISCALE+2
      IY=MAX0(MIN0(IY,120),2)
      IOUT(IY)=IMINUS
      IY=INT(Y(I)-YC(I))/IDIF+61
      IY=MAX0(MIN0(IY,120),2)
      IOUT(IY)=ISTAR
      IF(MOD(I-1,10).NE.0) GOTO 55
      TLABEL=THMIN+STEP*FLOAT(I-1)
      WRITE(6,104)(IOUT(J),J=1,113),TLABEL,IOUT(120)
104     FORMAT(1X,113A1,F6.2,A1)
      GO TO 60
55      WRITE(6,201) (IOUT(J),J=1,120)
201     FORMAT(1X,120A1)
60      CONTINUE
45    IF(JOBTYP.LT.3)GOTO 72
      REWIND 4
      WRITE(4,70)THMIN,STEP,THMAX
70    FORMAT(3F8.4) 
      WRITE(4,71)(YC(I),I=1,NPTS)
71    FORMAT(8(F7.0,1X))
72      continue
C 72    IF(IPL2.NE.0)CALL CALPLT(IFY,IFYC,IFM,IFD,IFB,EXPAND,glb(1))
      IF (IPLST.NE.0.AND.MAXS.NE.0) THEN
      INUMB = ICYRUN + 1
      NPAGES = MAXS/12
      IF (MOD(MAXS,12).NE.0) NPAGES = NPAGES+1
C     LIST PARAMETERS IN EACH CYCLE
      WRITE (6,2002)
2002    FORMAT ('1PARAMETERS IN EACH CYCLE')
      DO 7924 J=1,NPAGES
      REWIND 8
        ISTART = 1 + (J-1)*12
        IFINIS  = MIN0(MAXS,12 + (J-1)*12)
        WRITE (6,2001) (L,L=ISTART,IFINIS)
2001      FORMAT ('0CYCLE',12(3X,I2,5X))
        DO 7924 K=1,INUMB
          READ  (8,ERR=99992) (DUMMY(I),I=1,2*MAXS+4)
7924        WRITE (6,2000,ERR=99990) K-1,(DUMMY(I),I=ISTART,IFINIS)
2000    FORMAT (1X,I2,')',2X,12E10.4)
C     LIST R-VALUES  IN EACH CYCLE
      REWIND 8
      WRITE (6,2012)
2012    FORMAT ('1R-VALUE VARIATION WITH CYCLE')
      WRITE (6,2011)
2011    FORMAT ('0CYCLE',4X,'R-P',5X,'R-WP',2X,'    S   ',3X,'D-W D') 
      DO 7944 K=1,INUMB
        READ  (8,ERR=99992) (DUMMY(I),I=1,2*MAXS+4)
7944      WRITE (6,2010,ERR=99990) K-1,(DUMMY(I),I=2*MAXS+1,2*MAXS+4) 
2010    FORMAT (1X,I2,')',2X, 4F8.2)
C TESTE PARA GERAR SAIDA PARA 9411
          if(i2d94.NE.0)THEN
            call write94(iscale,idif)
          END IF
C     LIST PARAMETER SHIFTS IN EACH CYCLE
      IF (INUMB.EQ.1) GO TO 88880
      WRITE (6,2022)
2022    FORMAT ('1APPLIED PARAMETER SHIFT IN EACH CYCLE')
      DO 7964 J=1,NPAGES
      REWIND 8
        ISTART = 1 + MAXS+ (J-1)*12
        IFINIS = MIN0(2*MAXS,12 + MAXS+ (J-1)*12)
        WRITE (6,2021) (L,L=ISTART-MAXS,IFINIS-MAXS)
2021      FORMAT ('0CYCLE',12(3X,I2,5X))
        DO 7964 K=1,INUMB-1 
          READ  (8,ERR=99992) (DUMMY(I),I=1,2*MAXS+4)
          IF (K.EQ.1) GO TO 7964
7964        WRITE (6,2020,ERR=99990) K,(DUMMY(I),I=ISTART,IFINIS)
2020    FORMAT (1X,I2,')',2X,12E10.4)
      END IF
88880 CONTINUE
C     CODE FOR PRINTING PARAMETERS AND STD. DEV. IN THE FINAL CYCLE
      IF (IPLST.EQ.2.and.maxs.ne.0) THEN
      WRITE (6,5000) TITLE
5000  FORMAT ('1PARAMETERS AND STANDARD DEVIATIONS IN THE FINAL CYCLE'
     1 ,' FOR DATA BASE ',/,' " ',A70,' " ')
      WRITE (6,7445,ERR=99991) ((FINAL(I,J),J=1,2),I=1,ILOC)
7445  FORMAT (12(1X,E10.4))
      END IF
      RETURN
99990 STOP 'ERROR IN WRITING TO FILE 6 IN SUBROUTINE EXPUT'
99991 STOP 'ERROR WRITING PARAMS,ST.DEV. IN SUBROUTINE EXPUT'
99992 STOP 'ERROR IN READING FROM UNIT 8 IN SUBROUTINE EXPUT'
      END 
C      SUBROUTINE REWRIT(ISCALE,IDIF,IFY,IFYC,IFM,IFD,IFB,EXPAND)
      SUBROUTINE REWRIT(ISCALE,IDIF)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      REAL LAMDA
      INTEGER PTR,fondo
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94 
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),XMAS(16)
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/MULTIP/TMASSA(99),MLTPHASE,xmltp(99),MURT(NATS),SAQF(99) 
     *,wtis(99)
      common/sizestrain/sizeG(15),strainG(15),sizeL(15),strainL(15)
     *,siz(15),strain(15),NsizeStrain
      common/codebck/ibckcode
       REWIND 5
C line 1
      WRITE(5,1)TITLE
1     FORMAT(A70)
      JOBTYP=JOBTYP-1
      NPROF=NPROF-1 
      nplof = nprof  
      IF (NsizeStrain.eq.9) nplof = 9
      INSTRM = INSTRM-1
c line 2.1
c line 2.1 changed due to size-strain calculation (NsizeStrain)
      WRITE(5,2)JOBTYP,nplof,NPHASE,ibckcode,NEXCRG,NSCAT,INSTRM,IPREF,
     1          iasym,iabsr,idata,isphase,I2D94
2     FORMAT(13I4, 5x,'LINE 2.1')
C line 2.2
      if(ibckcode.eq.-1)write(5,202)IAS,fondo
202   format(2i4,49x,'LINE 2.2')
C line 3
      WRITE(5,14)IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2,IPLST,
     *IPLOSS,IPLCAL,IPLPOL,IPLCOM,IPLDIS,IPLAM,IPBIG
Cc 14    FORMAT(19I1,38x,'LINE 3')
14    FORMAT(3(5i1,1x),3i1,36x,'LINE 3')
C line 4
      WRITE(5,23)LAMDA,RATIO(2),BKPOS,WDT,CTHM,TMV,RLIM,sw
23    FORMAT(3f8.5,6F8.4)
C line 5
      WRITE(5,30)MCYCLE,EPS,RELAX
30    FORMAT(1I4,5F4.2,33x,'CYCLS EPS RELAX P_CALC')
      IF(NBCKGD.LT.2)GOTO 120 
C line 6(*)
      WRITE(5,34)(POS(I),BCK(I),I=1,NBCKGD)
34    FORMAT(2F8.2)
120   IF(NEXCRG.LE.0)GOTO 122
C line 7(*)
      WRITE(5,341)(ALOW(I),AHIGH(I),I=1,NEXCRG)
341     format(2f8.2,41x,'EXCLUDED REGION')
122   IF(NSCAT.LE.0)GOTO 124
      DO 125 I=1,NSCAT
           IF (JOBTYP .EQ. 1.OR.JOBTYP.EQ.3) GOTO 1228
C line 8.1 XRD (*)
        WRITE(5,38)NAM(I),DFP(I),DFPP(I),XMAS(I),NSCAT
38        FORMAT(A4,3F8.4,29X,'SCATTERING SET ',I2)
        GOTO 126
C line 8.1 ND(*)
1228        WRITE(5,3838)NAM(I),DFP(I),XMAS(I),NSCAT
3838        FORMAT(A4,2F8.4,37X,'SCATTERING SET ',I2)
C line 8.2 XRD(*)
126       IF(JOBTYP.EQ.0.OR.JOBTYP.EQ.2)WRITE(5,39)(AC(J,I),J=1,9) 
39          FORMAT(9F8.5)
125     CONTINUE
124   CONTINUE
C line 9
      WRITE(5,43)MAXS
43    FORMAT(I8,49X,'PARAMS REFINED')
      N=0
      DO 4310 IIPHAS=1,NPHASE 
4310    N=N+NATOM(IIPHAS)
      DO 84 I=1,N
      DO 84 J=1,11
84        A(I,J)=SIGN(1.,A(I,J))*(FLOAT(10*LP(I,J))+ABS(A(I,J)))
      DO 87 I=1,NPHASE
      DO 88 J=1,6 
88        PAR(I,J+5)=SAVE(I,J)
      DO 87 J=1,27
87        APAR(I,J)=SIGN(1.,APAR(I,J))*(FLOAT(10*LPAR(I,J))+ABS(APAR(I,J
     *     )))
      DO 89 J=1,20
89      AGLB(J)=SIGN(1.,AGLB(J))*(FLOAT(10*LGLB(J))+ABS(AGLB(J)))
C line 10.1
      WRITE(5,45)GLB(1),GLB(10),GLB(11),glb(8),glb(9),glb(12),glb(13)
C line 10.11
      WRITE(5,451)AGLB(1),AGLB(10),AGLB(11),aglb(8),aglb(9),aglb(12),
     *            aglb(13)
45    FORMAT(7F8.4,1x,'ZER DISP TRANS p q r t')
451   format(7f8.4,1X,'CODEWORDS')
C new lines in the ICF !cp may 10 97 and (ibgd test) !cp jun 97
      if (ibgd.eq.1) goto 4600 
C      if (iax.eq.0) then
C line 10.2
         write(5,634)glb(20),glb(18),glb(19)
634      format(3f8.4,33X,'AM MON1 MON2')
C line 10.21
         write(5,6342)aglb(20),aglb(18),aglb(19)
6342     format(3f8.4,33x,'CODEWORS')
C         write(5,634)glb(17),glb(20),glb(18),glb(19),CTIME
c634      format(5f8.4,17X,'AIR AM MON1 MON2 CTIME')
C         write(5,6342)aglb(17),aglb(20),aglb(18),aglb(19)
c6342     format(4f8.4,25x,'CODEWORS')
C      end if
C      write(5,635)AI1,AI2,AI3,AI4,AI5,AI6,AI7,AI8,AI9,AI10
c635   format(f8.5,f8.6,f8.3,2f8.4,2f8.5,f8.3,f8.5,f8.6)
4600   IF(NBCKGD.EQ.0)WRITE(5,46)(GLB(J),J=2,7),(AGLB(J),J=2,7)
46    FORMAT(6f9.2,3x,'BACKGROUND',/,6F9.4,3X,'CODEWORDS') 
477   DO 81 K=1,NPHASE
      IOF=0
      IF(K.GT.1) THEN
        DO 4311 IIPHAS=2,K
4311        IOF = IOF + NATOM (IIPHAS-1)
      END IF
C line 11.1
      WRITE(5,465)PHSNM(K),K
465     FORMAT(A50,7X,'PHASE NUMBER ',I2)
C line 11.2
      WRITE(5,48)NATOM(K),NMOL(K), saqf(k),(PREF(K,I),I=1,3),wtis(k)
48    FORMAT(2I4,f7.4,1x,3F4.1,f7.2,22X,'#ATMS #FU AFQPA PREFDIR ISWT')
      N=NATOM(K)
C line 11.3
      WRITE(5,460)(SYMB(K,I),I=1,20)
460     FORMAT(20A1,37X,'SPACE GROUP')
C Changing N to 'so' !cp oct 96 
           do 3947 isof=1,n
           xl(isof+iof,5)=xl(isof+iof,5)*xmltp(k)/murt(isof+iof)
3947         continue
C !cp oct 96 #6 murt parametrs included below. FORMAT modified...
      WRITE(5,65)(ATEXT(I+IOF),murt(i+iof),NTYP(I+IOF),
     *     (XL(I+IOF,J),J=1, 5),(A(I+IOF,J),J=1, 5),
     *        (XL(I+IOF,J),J=6,11),(A(I+IOF,J),J=6,11),
     *            I=1,N)
c65      FORMAT(A4,1x,i4,1x,a4,2X,5F8.5,2x,'At M At x y z B So'/,
C     *      16X,5F8.2,2x,'CODEWORDS'/,6F8.5,10X,'BETAS',/,
C     *      6F8.2,10x,'CODEWORDS')
65      FORMAT(A4,1x,i4,1x,a4,2X,5F8.5,2x,'LBL M NTYP x y z B So'/,
     *      16X,5F8.2,2x,'CODEWORDS'/,6F8.5,10X,'BETAS',/,
     *      6F8.2,10x,'CODEWORDS')
      WRITE(5,50) PAR(K,1),PAR(K,2),
     *                 APAR(K,1),APAR(K,2),
     *PAR(K,3),PAR(K,4),PAR(K,5),PAR(K,21),par(k,20),
     *         PAR(K,15),PAR(K,16),
     *APAR(K,3),APAR(K,4),APAR(K,5),APAR(K,21),apar(k,20),
     *         APAR(K,15),APAR(K,16),
     *              (PAR(K,I),I=6,11), (APAR(K,I),I=6,11), 
     *              PAR(K,12),PAR(K,13),PAR(K,14),
     *                 APAR(K,12),APAR(K,13),APAR(K,14),
     *              PAR(K,17),PAR(K,18),PAR(K,19),
     *                 APAR(K,17),APAR(K,18),APAR(K,19),
     1              par(k,24),par(k,25),par(k,26),
     1                 apar(k,24),apar(k,25),apar(k,26),
     1              par(k,27),
     1                 apar(k,27)
C !cp oct 96. End of FORMAT MODIFICATION
50      FORMAT(G8.3,F8.4,41X,'SCALE Bo(OVERALL)',/,2F8.2,/,
     *        7F8.5,1X,'U V W CT Z X Y',/,7F8.2,/,
     *        6F8.4,9X,'CELL PARAMETERS',/,6F8.2,/,
     *        3F8.5,33X,'PREF1 PREF2 R/RCF_ASYM',/,3F8.2,/,
     *        3F8.4,33X,'NA NB NC (MIX_PARAMS)',/,3F8.2,/,
     1        3F8.4,33X,'NA NB NC (HIGH SIDE)',/,3F8.2,/,
     1         f8.4,49x,'PEARSON ASYM.FACTOR',/,F8.2)
81      CONTINUE
      IF(IPL.NE.0)WRITE(5,191)ISCALE,IDIF
191   FORMAT(2I8,41X,'LINE PRINTER INFO')
C lines commented to avoid LINE 13 in the ICF    !cp jul 97
C      IF(IPL2.NE.0)WRITE(5,171)IFY,IFYC,IFM,IFD,IFB,EXPAND
C      IF(IPL2.EQ.0)WRITE(5,171)  2,   1,  1,  1,0,0.95
c171   FORMAT(5I1,3X,F8.4,41x,'CALCOMP INFO')
C                       JOBTYP  must be reproduced !!!
151   JOBTYP=JOBTYP+1
      RETURN
      END
      SUBROUTINE SUMMAT(IPM,CSK,DISK,DYCDD,ISODER,TOTCS)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      REAL LAMDA, LAMDAM
      INTEGER PTR,FONDO
      COMMON/F1/RJAC(MSZ,MSZ),V1(MSZ)
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),XMAS(16)
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
C     COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ),BK(IDSZ),NPTS
      COMMON/G3/COND,IORD1,IORD2,TH,NUM 
      COMMON/G4/TANN(NOV),DERSTO(NOV,MSZ)
      COMMON /PVII/ TF1,TF2,TF4,TF6,TF8,TF9,C4
      common/spvii/rl,da1l,da2l,da3l,da4l,da5l,da6l,da7l
      common/spvii/rh,da1h,da2h,da3h,da4h,da5h,da6h,da7h
      COMMON/PRFX/DELT,TL,GAM1,GAM2,PRFDER,iph,delta
C      common/prfx/iph,delta
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      COMMON/FONDI/BKCOM(IDSZ),BKDIS(IDSZ),BKAM(IDSZ),BKPOL(IDSZ)
      DIMENSION DERIV(MSZ)
      DIMENSION CSK(99),DISK(99),DYCDD(99),ISODER(NATS)
      LOGICAL VERT
      REAL ISODER
      shift = 0.
      DO 2 J=1,MSZ
2       DERIV(J)=0.0 
      YCALC=0.0
      IL=0
C-----CALCULATE THE CONTRIBUTION OF THE REFLECTIONS ORD1 TO ORD2 TO THE
C-----DERIVATIVES W.R.T. THE PROFILE INTENSITY YOBS
      IF(IORD1.EQ.0)GOTO 12
      DO 33 I=IORD1,IORD2
      IPH=IREFS(I)/(256*256*256*8)
      IL=IL+1
      j=mod(i,nov)+1
C test for asymmetry function !cp ap 12 97
      if (iasym.eq.0) then
        VERT = REFS(I,2).LE.RLIM
      else  
        VERT = ABS(REFS(I,2)-90.0).GE.RLIM
      end if
      SHIFT = GLB(10) * COS(REFS(I,2)/2./57.2958)
     *           + GLB(11) * SIN(REFS(I,2)/57.2958)
      DELTA=TH-REFS(I,2)-GLB(1)-SHIFT 
      tanth=tan(th*3.14159265359/360.)
      DELT=DELTA*DELTA
      TL=REFS(I,1)
      IF (NPROF.EQ.6) GAM1 = PAR(IPH,17) + PAR(IPH,18) * REFS(I,2)
      IF (NPROF.EQ.7) GAM1 = PAR(IPH,17) + PAR(IPH,18) / REFS(I,2)
     *                              + PAR(IPH,19)/REFS(I,2)/REFS(I,2) 
      IF (NPROF.EQ.7) CALL PRSVII(gam1)
      if (nprof.eq.5) then
        rl=par(iph,17)+(par(iph,18)+par(iph,19)/refs(i,2))/refs(i,2)
        rh=par(iph,24)+(par(iph,25)+par(iph,26)/refs(i,2))/refs(i,2)
        call mspvii(par(iph,27),tl)
      end if
      IF(NPROF.EQ.8) THEN
        TLL = HALFL(I)
        GAM1 = GAM(I)
        TLR = TLL/TL
      END IF
      BB=TL*TL
c-----NEXT LINE IS NECESSEARY FOR 2 PHASES WITH VERY DIFFERENT FWHM.
          IF (DELT/BB.GT.WDT*WDT) GO TO 33
          IF(VERT)then
C       test for asymmetry model               !cp may 01 97
             if (iasym.eq.0)then
                YX=DELT*SIGN(1.0,DELTA)
                Z=1.0-PAR(IPH,14)*YX/TANN(J)
                IF ( Z.LE.0.0 ) Z=0.0001
             else
                 YX=SIGN(1.0,DELTA)*DELTA/(2*TL)
                 TANNJ=TANN(J)
                 IF (TANNJ.GE.1.0) THEN
                   TANNJ=TAN(ATAN(TANNJ)-3.14159265359/2)
                 END IF
                 Z=(PAR(IPH,14)/TANNJ) * (2.0*(DELTA/(2*TL))*EXP(-YX))
                 Z=1+Z
                 IF ( Z.LE.0.0 ) Z=.0001
             end if
          else
C                          GOTO 5
c4                       Z=1.0
             Z=1.0
          end if
5     PRTEMP = PROFIL(NPROF,DELT/BB)
      if (nprof.eq.5) then
        OMEGA = REFS(I,3)*Z*PRTEMP*PAR(IPH,1)
      else
        OMEGA = REFS(I,3)*Z*PRTEMP*PAR(IPH,1)/TL
      end if
      YCALC = YCALC+OMEGA*FMGNTD(I)
      IF ( JOBTYP.GT.2 ) GOTO 33
      X = PRFDER*2.*DELT/BB-1.0
      DO 3 K=1,MAXS
        DER=1.0
c-----Broadening Coeficients Derivatives
      if(nprof.ne.5) then
        DO 6 M=3,5
          IF(LPAR(IPH,M).EQ.K) DER=X/TL/2.
6         CONTINUE
      end if
        IF (LPAR(IPH,20).EQ.K) DER=X/TL/2.
        X1=0.0
c-----Asymmetry Derivative
        IF (VERT) then
           if (iasym.eq.0)then
                X1=PAR(IPH,14)*SIGN(1.,DELTA)*BB/TANN(J)/Z 
           else 
                X1=-PAR(IPH,14)*EXP(-YX)*(TL/(2*DELTA)-
     &              SIGN(1.0,DELTA)*1.0/4)/TANNJ/Z
           end if
        end if
c-----Zero, Displacement, and Transparancy Derivative
        IF ( LGLB(1).EQ.K ) DER=DELTA*(PRFDER+X1)
        IF ( LGLB(10).EQ.K ) DER=DELTA*(PRFDER+X1)
        IF ( LGLB(11).EQ.K ) DER=DELTA*(PRFDER+X1)
        IF ( (LPAR(IPH,14).EQ.K) .AND. VERT ) THEN
          IF ( IASYM.EQ.0 ) THEN
            DER = YX/Z
          ELSE
            DER = -2.0*(DELTA/(2.0*TL))*EXP(-YX)/Z
          END IF
        END IF
        IF ( NPROF.EQ.8 ) GO TO 8
c-----Pseudo-Voigt Shape Derivatives
        IF(NPROF.EQ.6) THEN
          KRP1=LPAR(IPH,17)
          IF(K.EQ.KRP1)DER=(.636619772/(1.+4.*DELT/BB)-
     *             .939437279*EXP(-2.772588722*DELT/BB))/PRTEMP
          KRP1=LPAR(IPH,18) 
          IF(K.EQ.KRP1)DER=(.636619772/(1.+4.*DELT/BB)-
     *             .939437279*EXP(-2.772588722*DELT/BB))/PRTEMP
        END IF
c-----Pearson VII Shape Derivatives
      IF (NPROF.EQ.7) THEN
          KRP1=LPAR(IPH,17) 
          IF(K.EQ.KRP1)
     *      DER=-ALOG(1.+TF9*DELT/BB)+TF4*(DELT/BB)/(1.+TF9*DELT/BB)
     1            +TF8
          KRP1=LPAR(IPH,18)
          IF(K.EQ.KRP1)
     *      DER=-ALOG(1.+TF9*DELT/BB)+TF4*(DELT/BB)/(1.+TF9*DELT/BB)
     1            +TF8
          KRP1=LPAR(IPH,19)
          IF(K.EQ.KRP1)
     *      DER=-ALOG(1.+TF9*DELT/BB)+TF4*(DELT/BB)/(1.+TF9*DELT/BB)
     1            +TF8
      END IF
c-----Lattice Parameter Derivatives
8         DO 7 M=6,11
          IF(LPAR(IPH,M).EQ.K) DER=(PRFDER+X1)*DELTA
7         CONTINUE
        DERIV(K)=DERSTO(J,K)*DER*OMEGA+DERIV(K)
3        continue
c----TCHZ Profile Derivatives
      IF(NPROF.EQ.8) THEN
        OMEGA8 = Z*PAR(IPH,1)*REFS(I,3)/TL
         DO 1000 K = 1,MAXS
          DO 1001 M=3,5
            IF (LPAR(IPH,M).NE.K) GO TO 1001
            DERIV(K) = DERIV(K)+ OMEGA8*DERSTO(J,K)/2.*
     *                  (.939437279*EXP(-2.772588722*DELT/BB) -
     *                            .636619772/(1.+4.*DELT/BB)) *
     *                                  (1.36603*TLR/TL-.95438*TLR*TLR/T
     *         L+.33348*TLR**3./TL)
1001          CONTINUE
            IF (LPAR(IPH,20).EQ.K)
     *        DERIV(K) = DERIV(K)+ OMEGA8*DERSTO(J,K)/2.*
     *                  (.939437279*EXP(-2.772588722*DELT/BB) -
     *                            .636619772/(1.+4.*DELT/BB)) *
     *                                  (1.36603*TLR/TL-.95438*TLR*TLR/T
     *         L+.33348*TLR**3./TL)
          IF (LPAR(IPH,15).NE.K) GO TO 1002
          DERIV(K) = DERIV(K)+ OMEGA8*
     *              (.939437279*EXP(-2.772588722*DELT/BB) - 
     *                      .636619772/(1.+4.*DELT/BB)) *
     *                          ((1.36603*TLR/TL-.95438*TLR*TLR/TL+.3334
     *       8*TLR**3./TL)*     DERSTO(J,K)/2. - FMGNTD(I)*TANN(J)*
     *                                        (1.36603/TL-.95438*TLR/TL+
     *       .33348*TLR*TLR/TL))
1002        CONTINUE
          IF (LPAR(IPH,16).NE.K) GO TO 1003
          DERIV(K) = DERIV(K)+ OMEGA8*
     *              (.939437279*EXP(-2.772588722*DELT/BB) - 
     *                      .636619772/(1.+4.*DELT/BB)) *
     *                          ((1.36603*TLR/TL-.95438*TLR*TLR/TL+.3334
     *       8*TLR**3./TL)* DERSTO(J,K)/2. - FMGNTD(I)*SQRT(1+TAN
     *       N(J)*TANN(J))*                   (1.36603/TL-.95438*TLR/TL+
     *       .33348*TLR*TLR/TL))
1003        CONTINUE
1000    CONTINUE
      END IF
33      CONTINUE
c-----FORM SUMS
12    IF(NBCKGD.NE.0)GOTO 11
      DO 10 II=2,7
      IF(LGLB(II).EQ.0.)GOTO 10
      KM=LGLB(II)
      IF(II.EQ.2)DERIV(KM)=DERIV(KM)+1.
      IF(II.EQ.2)GOTO 10
      DERIV(KM)=DERIV(KM)+((THMIN+FLOAT(IPM-1)*STEP)/BKPOS-1.)**(II-2)
10      CONTINUE
11    YC(IPM)=YCALC 
      IF (JOBTYP.GT.2) GOTO 20
c
c
      DO 100 K = 1,NPHASE
c
c
c
      LK1  = LPAR(K,1)
      LK2  = LPAR(K,2)
c
c-----UPDATING GLOBAL SCALE DERIVATE FOR BKG CONTRIBUTE
c
      IF(FONDO.EQ.1.OR.FONDO.EQ.2) THEN
      IF(LK1.NE.0) DERIV(LK1)=DERIV(LK1)+GCOM(K)*CSK(K)+GCOM(K)*DISK(K)
      END IF
c
c
c-----UPDATING DERIVATE OF Q OVERALL FOR BKG CONTRIBUTE
c
      IF(FONDO.EQ.2) THEN
        IF(LK2.NE.0) DERIV(LK2) = DERIV(LK2) +  DYCDD(K)
      END IF
c
c-----UPDATING DERIVATE OF ISOTROPIC THERMAL PARAMETERS FOR BKG CONTRIBUTE
c
       IF(FONDO.EQ.1) THEN
        IOF = 0
        IF(K.GT.1) THEN
          DO 13 I = 2,K
13          IOF = IOF + NATOM(I-1)
        ENDIF
c
        DO 14 I = 1,NATOM(K)
         IISO=LP(I+IOF,4)
         IF(IISO.NE.0) DERIV(IISO) = DERIV(IISO) + ISODER(I+IOF)
         ISODER(I+IOF)=0.0
14        CONTINUE
       ENDIF
100   CONTINUE
c-----DYC RESPECT TO AMORPHOUS SCALE FACTOR
      LK = LGLB(20)
      IF(LK.NE.0) DERIV(LK) = DERIV(LK) + AMORPHOUS(IPM)
c-----MONOCHROMATOR PARAMETERS DERIVATIVES
      LAMDAM=(LAMDA(1)*RATIO(1)+LAMDA(2)*RATIO(2))/(RATIO(1)+RATIO(2))
      ESSE=2*SIN((TH-GLB(1)-SHIFT)*0.00872665)/LAMDAM
      LK=LGLB(18)
      IF ( LK.NE.0 ) THEN
c-------NEXT LINE FOR A LORENTZIAN MONOCHROMATOR BASS-BAND  FUNCTION
        ASS5=1/(1+GLB(18)*ESSE**GLB(19))
        DERMON=-(ESSE**(GLB(19)))/ ((1+GLB(18)*ESSE**GLB(19))**2 )
c--------------------------------------------------------------------
        DERIV(LK)=DERIV(LK)+TOTCS/ASS5*DERMON
      ENDIF
C   !cp ap 20 97
      LK=LGLB(19)
      ASS5=1/(1+GLB(18)*ESSE**GLB(19))                 
      IF ( LK.NE.0 ) THEN
        DERMON=-GLB(18)*ALOG(ESSE)*(ESSE**GLB(19))       
     *       /(1+GLB(18)*(ESSE**GLB(19)))**2                   
        DERIV(LK)=DERIV(LK)+TOTCS/ASS5*DERMON
      ENDIF
c-----FORM THE UPPER TRIANGULAR OF "RJAC" = MATRIX OF NORMAL EQUATIONS
      DELTA = Y(IPM)-BK(IPM)-YCALC
      DO 9 J=1,MAXS 
        X = DERIV(J)/(VAR(IPM)) 
        V1(J) = V1(J)+DELTA*X
        DO 9 KK=J,MAXS
9         RJAC(J,KK) = RJAC(J,KK)+X*DERIV(KK)
20    RETURN
      END 
      SUBROUTINE REFGEN(IPHASE,ZERO,DIS,TRANS,PREF,PREFOR)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      COMMON /HKLCTL/ IHKL(3,48),AZ(48),NC1(10,7,99),ICHKL(99),
     1    N1HKL(99),IER
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      REAL LAMDA
      INTEGER H1,H2,H3,fondo
      LOGICAL ORH1,ORH2,ORH3
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      COMMON/CELLX/A,B,C,ALPHA,BETA,GAMMA,AL(3,3) 
      DIMENSION PREF(99,3)
      EQUIVALENCE (AL(1,1),ASTAR),(AL(2,2),BSTAR),(AL(3,3),CSTAR),
     * (AL(1,2),BACOSC),(AL(1,3),CACOSB),(AL(2,3),BCCOSA)
      CHARACTER*8 LAU(14)
      DATA LAU(1),LAU(2),LAU(3),LAU(4),LAU(5),LAU(6),LAU(7),LAU(8),
     * LAU(9),LAU(10),LAU(11),LAU(12),LAU(13),LAU(14)
     *  /'1BAR','2/M','MMM','4/M','4/MMM','3BAR   R','3BAR M R','3BAR',
     *   '3BAR M 1','3BAR 1 M','6/M','6/MMM','M3','M3M'/  
     
      SIN1(A) = SIN(A*6.28318531)
      
!     RAD = ATAN(1.0)/45.0
      RAD = 3.14159265359/360.0 
      IF ( NAXIS.GT.3 ) STOP 5001 
      KXIS = NAXIS
      NAXIS = 1 
1002  N = 0
      J = 0
      I2D = 0
      I3D = 0
      I1D = 1
      I12D = 1
      I13D = 1
      I123D = 1
      I23D = 1
      I2DEL = 1
      I3DEL = 1
      DO 8000 I=1,8 
        IF ( NCONT(I).EQ.0 ) GO TO 9000
        IF ( NCONT(I)-8192 ) 8002,8001,8000 
8001    I3DEL = 3
        I1D = 2
        I123D = 3
        GO TO 9000
8002    IF ( NCONT(I)-576 ) 8004,8003,8004
8003    I3DEL = 2
        I23D = 2
        GO TO 8000
8004    IF ( NCONT(I)-516 ) 8008,8009,8008
8009    I3DEL = 2
        I13D = 2
        GO TO 8000
8008    IF ( NCONT(I)-68 ) 8012,8013,8012
8012    IF ( NCONT(I)-580 ) 8015,8016,8015
8013    I2DEL = 2
        I12D = 2
        GO TO 8000
8016    I123D = 2
        I3DEL = 2
8015    CONTINUE
8000    CONTINUE
9000  N = J-1
97    TMIN = (SIN1(THMIN/720.0)/LAMDA(1))**2
      IC = 0
      IF ( IPHASE.GE.2 ) THEN
        DO 4100 IIPHAS=2,IPHASE
4100      IC = IC+ICR(IIPHAS-1)
      END IF
C     IF(IPHASE.EQ.2)IC=ICR(1)
      SMAX = (SIN1(THMAX/720.0)/LAMDA(1))**2
c
C     **************************************************************
C     SMAX IS CHANGED TO ACCOUNT FOR THE LEFT TAILS OF THE REFLECTIONS
C     THAT ARE PRESENT AT ANGLES GREATER THAN THMAX
C     **************************************************************
C     PRINT *,'THMAX=',THMAX
C      THMAX1=(U*(TAN(THMAX*RAD))**2+V*TAN(THMAX*RAD)+W+ZZZ*(1+
C     *(TAN(THMAX*RAD))**2))
C also incorporating the cotg^2 term  !cp may 01 97
      THMAX1=(U*(TAN(THMAX*RAD))**2+V*TAN(THMAX*RAD)+W+ZZZ*(1+
     *    (TAN(THMAX*RAD))**2) + uc/(tan(thmax*rad)))
      IF ( THMAX1.GT.0. ) THEN
        THMAX1= WDT*SQRT(THMAX1)/2
      ELSE
        WRITE (6,3456) POS,IPHASE
        WRITE (7,3456) POS,IPHASE
        STOP 'SQUARE OF FWHM IS NEGATIVE'
      END IF
c
      ANGTTETA=THMAX+THMAX1
      IF ( (THMAX+THMAX1).GE.180.0 ) ANGTTETA=180.0
      SMAX1 = (SIN1((ANGTTETA)/720.0)/LAMDA(1))**2
      SMAX = SMAX1
c
      THMAXX=ANGTTETA
      IF ( LAMDA(2).GT.LAMDA(1) ) TMIN=TMIN*(LAMDA(1)/LAMDA(2))**2
      IF ( LAMDA(2).GT.0. .AND. LAMDA(2).LT.LAMDA(1) ) 
     *     SMAX=SMAX*(LAMDA(1)/LAMDA(2))**2 
      CALL OP1(IPHASE)
      WRITE(6,36)LAU(NSPGRP)
36    FORMAT(15H LAUE SYMMETRY ,A8,33H WILL BE USED TO GENERATE INDICES,
     */,'        ------------------------------------------')
      I1MAX=A*2.0*SQRT(SMAX)
      IF ( KXIS.LT.3 ) THEN
        I2MAX=B*2.0*SQRT(SMAX)
        I3MAX=C*2.0*SQRT(SMAX)
      ELSE
        I2MAX = 2.0*C*SQRT(SMAX)
        I3MAX = 2.0*B*SQRT(SMAX)
      END IF
      L1 = NSPGRP
      IF ( L1.GE.13.AND.NCONT(1).EQ.580 ) L1=L1+2
      IF ( I12D+I13D+I23D.NE.5 ) GO TO 2131
      I12D = 2
      I13D = 1
      I23D = 2
      I2DEL = 2
      I3DEL = 2
2131  IF ( L1.EQ.6 ) I2D=1
      IF ( L1.EQ.7 ) I2D=2
      I3D = I2D
      I1 = -1
2200  I1 = 1+I1
      IF ( I1-I1MAX ) 2201,2201,2303
2201  H1 = I1
      IF ( I2D.GT.0 ) I2MAX=I1
      GO TO (2203,2204,2204,2205,2206,2208,2209,2207,2204,2206,2204,
     *  2206,2206,2206,2206,2206),L1
      CALL GOTOER
2203  I2= -MIN0(I2MAX,I1*I2MAX)
      GO TO 2210
2204  I2 = 0
      GO TO 2210
2205  I2= MIN0(I1,1)
      GO TO 2210
2206  I2 = I1
      GO TO 2210
2207  I2= MIN0(-I1+1,0)
      GO TO 2210
2208  I2 = -2*I1
      GO TO 2210
2209  I2 = -I1/2
2210  I2= I2DEL*(I2/I2DEL)+MOD(I1,I12D) 
      GO TO 2221
2220  I2 = I2DEL+I2 
2221  H2 = I2
      IF ( I2-I2MAX ) 2222,2222,2200
2222  GO TO (2223,2223,2224,2224,2224,2225,2225,2224,2224,2223,2224,
     *  2224,2226,2227,2236,2235),L1
      CALL GOTOER
2223  I3= -MIN0(I3MAX,(I1+IABS(I2))*I3MAX)
c2223 I3= -MIN0(I3MAX,I1*I3MAX)  !(ALLEN MADE THIS CHANGE,REVERTED LATER)
      GO TO 2228
2224  I3 = 0
      GO TO 2228
2225  I3 = -I1-I2
      GO TO 2228
2226  I3= MIN0(I2,I1+I3DEL)
      GO TO 2228
2227  I3 = I2
2228  I3= I3DEL*(I3/I3DEL)+MOD(MOD(I1+I1D*I2,I123D)+I123D,I123D) +MOD(I1
     * ,I13D)+MOD(I2,I23D)
      IF ( I3D-1 ) 2232,2229,2231
2229  I3MAX = I1
      IF ( I2-I1 ) 2230,2232,2230
2230  I3MAX = I3MAX-1
      GO TO 2232
2231  I3MAX = I2
      GO TO 2232
2235  I3= I2+MOD(I1,2)
      GO TO 2232
2236  I3 = I1+2-MOD(I2,2)
      IF ( I1.EQ.I2.AND.MOD(I1,2).EQ.0 ) I3=I1
      GO TO 2232
2233  I3 = I3DEL+I3 
2232  H3 = I3
      IF ( KXIS.NE.3 ) GO TO 2240
      H3 = I2
      H2 = I3
2240  CONTINUE
      IF ( I3-I3MAX ) 113,113,2220
113   SQ=FLOAT(H1*H1)*ASTAR+FLOAT(H2*H2)*BSTAR+FLOAT(H3*H3)*CSTAR
     * +FLOAT(H1*H2)*BACOSC+FLOAT(H1*H3)*CACOSB+FLOAT(H2*H3)*BCCOSA
      SQ=SQ/4.
2234  IF ( SQ-SMAX ) 3000,3000,2233
2303  ICR(IPHASE)=IC
      IF (IPHASE.GE.2) THEN
      DO 4101 IIPHAS=2,IPHASE
4101      ICR(IPHASE) = ICR(IPHASE)-ICR(IIPHAS-1) 
      END IF
C     IF(IPHASE.EQ.2)ICR(IPHASE)=IC-ICR(1)
      CALL SORT(IPHASE)
      RETURN
3000  IF ( SQ - TMIN ) 2233,3117,3117
3117  IF ( IC.GT.IRS-2 ) GO TO 6001
C     NEXT IF BLOCK FOR PHASE WITH DISTINCT ORIENTATION ALONG PREF(NAXIS,I)
      IF (PREFOR.GT.99.) THEN 
c
      ORH1 = IFIX(PREF(IPHASE,1)).EQ.H1
      IF (.NOT.(ORH1).AND.H1.NE.0.AND.IFIX(PREF(IPHASE,1)).NE.0)
     *        ORH1 = MOD(H1,IFIX(PREF(IPHASE,1))).EQ.0
      IF (.NOT.ORH1) GO TO 2233
c
      ORH2 = IFIX(PREF(IPHASE,2)).EQ.H2
      IF (.NOT.(ORH2).AND.H2.NE.0.AND.IFIX(PREF(IPHASE,2)).NE.0)
     *        ORH2 = MOD(H2,IFIX(PREF(IPHASE,2))).EQ.0
      IF (.NOT.ORH2) GO TO 2233
c
      ORH3 = IFIX(PREF(IPHASE,3)).EQ.H3
      IF (.NOT.(ORH3).AND.H3.NE.0.AND.IFIX(PREF(IPHASE,3)).NE.0)
     *        ORH3 = MOD(H3,IFIX(PREF(IPHASE,3))).EQ.0
      IF (.NOT.ORH3) GO TO 2233
      IORDR1=0
      IORDR2=0
      IORDR3=0
      IF (H1.NE.0) IORDR1 = INT(H1/IFIX(PREF(IPHASE,1)))
      IF (H2.NE.0) IORDR2 = INT(H2/IFIX(PREF(IPHASE,2)))
      IF (H3.NE.0) IORDR3 = INT(H3/IFIX(PREF(IPHASE,3)))
      IF (IORDR1.EQ.IORDR2.AND.IORDR2.EQ.IORDR3) GO TO 9257
      IF (IORDR1.EQ.0.AND.(IORDR2.EQ.IORDR3)) GO TO 9257
      IF (IORDR2.EQ.0.AND.(IORDR3.EQ.IORDR1)) GO TO 9257
      IF (IORDR3.EQ.0.AND.(IORDR1.EQ.IORDR2)) GO TO 9257
      IF (IORDR1.EQ.0.AND.IORDR2.EQ.0) GO TO 9257
      IF (IORDR2.EQ.0.AND.IORDR3.EQ.0) GO TO 9257
      IF (IORDR3.EQ.0.AND.IORDR1.EQ.0) GO TO 9257
      END IF
9257  CONTINUE
C     SEPARATE PHASE FOR ORIENTATION COMPLETE
      IHKL(1,1)=H1
      IHKL(2,1)=H2
      IHKL(3,1)=H3
      AZ(1)=0.0
      IER=0
      CALL SMTRY2(IPHASE)
      IF(IER.NE.0)GOTO 2233
      LXN=2
      IF(LAMDA(2).EQ.0..OR.LAMDA(2).EQ.LAMDA(1))LXN=1
      DO 3118 LX1=1,LXN
      SQH=SQ*LAMDA(LX1)*LAMDA(LX1)
      TAN2=SQH/(1.-SQH)
      IF(SQH.GE.1.)GOTO 3118
      TANX=SQRT(TAN2)
C     SHIFT DUE TO SAMPLE DISPLACEMENT AND TRANSPARENCY  
      SHIFT =  DIS*SQRT(1-SQH)+TRANS*
     *                    SQRT(1.-(1.-2.*SQH)*(1.-2.*SQH))
      POS=ATAN(TANX)/RAD+ ZERO + SHIFT
C     IF(POS.GT.THMAX.OR.POS.LT.THMIN)GOTO 3118 
      IF ( POS.GT.THMAXX .OR. POS.LT.THMIN ) GOTO 3118 
      IC=IC+1
      XABS=1.0
      IF(TMV.LE.0.000001)XABS=1.
      PLOR=1./(2.*SQH*SQRT(1.-SQH))*XABS
      IF (INSTRM.EQ.2) THEN 
        PLOR = PLOR * (0.95+0.05*(1.-2.*SQH)*(1.-2.*SQH)) 
        GO TO 4000
      END IF
      IF(JOBTYP.EQ.3)PLOR=PLOR*(1.+(1.-2.*SQH)*(1.-2.*SQH)*CTHM)
      IF(JOBTYP.EQ.1)PLOR=PLOR*(1.+(1.-2.*SQH)*(1.-2.*SQH)*CTHM)
4000    CONTINUE
      IREFS(IC)=256*(256*(256*(8*IPHASE+LX1)+128+H1)+128+H2)+128+H3 
      FMGNTD(IC)=MULT(H1,H2,H3,KXIS)
c
c------CALCULATE FWHM FOR PSEUDOVOIGT WITH GAUSS AND LORENTZ
      IF (NPROF.EQ.8) THEN
        HALFG(IC) =    (U*TAN2+V*TANX+W+ZZZ*(1+TAN2))
        IF (HALFG(IC).GT.0.) THEN
          HALFG(IC) = SQRT(HALFG(IC)) 
        ELSE
          WRITE (6,3456) POS,IPHASE
          WRITE (7,3456) POS,IPHASE
          STOP 'SQUARE OF FWHM IS NEGATIVE'
3456        FORMAT (3X,'SQUARE OF FWHM NEGATIVE AT TWO-THETA',F8.3,
     *             'FOR PHASE NO. ',I4) 
        END IF
        HALFL(IC) = ULOR*TANX+VLOR/SQRT(1.-SQH) 
        REFS(IC,1) = (HALFG(IC)**5.+2.69269*HALFG(IC)**4.*HALFL(IC)+
     *              2.42843*HALFG(IC)**3.*HALFL(IC)**2.+
     *         4.47163*HALFG(IC)**2.*HALFL(IC)**3.+
     *              0.07842*HALFG(IC)*HALFL(IC)**4.+HALFL(IC)**5.)**0.2
        TLR = HALFL(IC)/REFS(IC,1)
        GAM(IC) = 1.36603*TLR-0.47719*TLR*TLR+0.11116*TLR**3.
      else if (nprof.eq.5) then
        REFS(IC,1) = (U*TAN2+V*TANX+W)
      else   
C incorporating cotg^2  !cp may 01 97
        REFS(IC,1)=(U*TAN2+V*TANX+W+ZZZ*(1+TAN2)+uc/tan2)
      end if
        IF (REFS(IC,1).GT.0.0) THEN
          REFS(IC,1)= SQRT(REFS(IC,1))
        ELSE
          WRITE (6,3456) POS,IPHASE
          WRITE (7,3456) POS,IPHASE
          STOP 'SQUARE OF FWHM IS NEGATIVE'
        END IF
7000    REFS(IC,2)=POS
      REFS(IC,3)=PLOR
3118    CONTINUE
      GOTO 2233
6001  WRITE(6,24) IRS
      WRITE(7,24) IRS
      WRITE(*,24) IRS
24    FORMAT (//' TOO MANY REFLECTIONS (',I6,'). INCREASE *IRS* IN THE',
     *          ' PARAMETER STATEMENT IN THE SOURCE CODE ') 
      STOP
C     GO TO 2303
      END
C  SUBROUTINE FINDC and SUBROUTINE COMPTON subroutine DISORDER: by Canton et all.
C Added by CPS between March-May 1997
      SUBROUTINE FINDC(K,NSCAT)
      include 'param.inc'
C     PARAMETER(IDSZ=8000,IRS=4096,NATS=16,MSZ=63 ,NOV=2048,NFINAL=1024)
C ----THIS SUBROUTINE ASSIGNS AUTOMATICALLY TO ATOMS
C     THE RIGHT CONSTANTS TO CALCULATE THE COMPTON SCATTERING
      INTEGER PTC
      CHARACTER*4 TABNC(128),TBXC(256)
      COMMON/TABCHR/TABNC,TBXC
      COMMON/TABLES/TBX(256,10),TBD(128,20),TABN(128),tbm(128)
      COMMON/INC/TCS(95,4)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     $NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCATx
      COMMON/COMP/CC(4,16),ZEFF(16),PTC(NATS)
      CHARACTER*4 NAM(16),NOME
      COMMON/COEFC/NAM
      CHARACTER*1 NOM(4),PIU,MENO
      DATA PIU,MENO/'+','-'/
      DO 7 II=1,16
7     CONTINUE
      IOF = 0
      NS  = 0
c-----K = NUMBER OF PHASE UNDER CONSIDERATION
      IF(K.GT.1) THEN
c-----CALCULATE IOF = ALL ATOMS OF THE K-1 PHASES
      DO 5 I = 2,K
5     IOF = IOF + NATOM(I-1)
      NS = NSAVE
      ENDIF
c-----DEFINE NK = NUMBER OF ATOMS OF THE K-TH PHASE
      NK = NATOM(K)
      DO 9999 I = 1,NK
C                      IF(NS.EQ.0) GO TO 25
      DO 20 J = 1,NSCAT
      IF(NTYP(I+IOF).EQ.NAM(J)) GO TO 30
20    CONTINUE
c-----212 = ALL THE POSSIBLE NAMES OF ATOMS AND IONS
25    DO 40 J = 1,212
      IF(NTYP(I+IOF).EQ.TBXC(J)) GO TO 50
40    CONTINUE
      WRITE(6,41) NTYP(I+IOF)
41    FORMAT(' COMPTON SCATTERING COEFFICIENT NOT FOUND FOR ',A4)
      STOP ' COMPTON SCATTERING DATA MISSING'
30    PTC(I+IOF) = J
      GO TO 9999
50    NOME = TBXC(J)
c-----FIND NA = THE ATOMIC NUMBER OF I-TH ATOM
      IF(J.EQ.1.OR.J.EQ.2.OR.J.EQ.3) THEN
      NA = TBX(J,10)
      ELSE
      NA = TBX(J,10) + 1.0
      ENDIF
      NS = NS + 1
      PTC(I+IOF) = NS
c-----PUT IN CC THE 4 COMPTON COEFFICIENTS
      DO 60 L = 1,4
      CC(L,NS) = TCS(NA,L)
60    CONTINUE
      NAM(NS) = NOME
c
      READ(NOME,65) (NOM(L),L=1,4)
65    FORMAT(4A1)
C     FIND IN WHAT COLUMN THERE IS + OR -
      DO 80 L= 1,4
      IF(NOM(L).EQ.PIU)  GO TO 90
      IF(NOM(L).EQ.MENO) GO TO 95
80    CONTINUE
C CASE WITH NOR + NOR -
      ZEFF(NS) = NA
      GO TO 9999
c-----CASE WITH PLUS
   90 L = L + 1
      IF(L.GT.4) STOP '  SOMETHING IS WRONG IN ATOMIC NAME '
      READ(NOM(L),1000) NN
1000  FORMAT(I1)
      ZEFF(NS) = NA - NN
      GO TO 9999
c-----CASE WITH MINUS
   95 L = L + 1
      IF(L.GT.4) STOP '  SOMETHING IS WRONG IN ATOMIC NAME '
      READ(NOM(L),1000) NN
      ZEFF(NS) = NA + NN
 9999 CONTINUE
      NSAVE = NS
c
      RETURN
      END
      SUBROUTINE COMPTON(K,STH,CISK)
      include 'param.inc'
C     PARAMETER(IDSZ=8000,IRS=4096,NATS=16,MSZ=63 ,NOV=2048,NFINAL=1024)
c-----THIS SUBROUTINE IS USED TO COMPUTE:
C     CISK = COMPTON INTENSITY SCATTERED BY THE K-TH CRYSTALLINE PHASE
C            AT THE IPM-TH POINT OF THE X-RAY SPECTRUM.
C     THE METHOD FOLLOWED IS THAT REPORTED IN :
c*****A NEW ANALYTIC APPROXIMATION TO ATOMIC INCOHERENT X-RAY SCATTERING
C     INTENSITIES.
C     BY VEDENE H. SMITH JR. AJIT J. THAKKAR AND DOUGLAS C. CHAPMAN
C     ACTA CRYST. (1975). A31, 391-392.
c-----COMPTON SCATTERING OF THE AMORPHOUS PHASE IS CONTAINED INTO AM(S)
C     MOREOVER CISK IS MULTIPLIED BY THE FOLLOWING GEOMETRICAL FACTOR:
C     SEE:
C     RULAND, ACTA CRYST. 1961, 14, 1180.
C     ASS1 = CORRECTION FOR ABSORPTION OF COMPTON SCATTERING BY THE SAMPLE
C     BDF  = BREIT-DIRAC RECOIL FACTOR
C     ASS2 = CORRECTION FOR ABSORPTION BY AIR CONTAINED BETWEEN THE SAMPLE
C            AND THE COUNTER
C     ASS3 = CORRECTION FOR THE CUTTING UP OF COMPTON SCATTERING BY THE
C            MONOCHROMATOR ( ALSO COMPTON SCATTERING IS DIFFRACTING )
      REAL LAMDA
      INTEGER PTC
      DIMENSION CSP(2)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     *NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      COMMON /RTMTX/ IVEC(99,192),MLTPHS(99),ICNTPHS(99)
      COMMON/COMP/CC(4,16),ZEFF(16),PTC(NATS)
c-----HMC = H/MC
      DATA HMC/0.024263935/
c-----NK = ATOMS IN THE ASSIMETRIC UNIT OF K-TH PHASE
      NK = NATOM(K)
c-----IRL = NUMBER OF EQUIVALENT POSITION - ALSO THE IDENTITY
C           AND IF IT IS PRESENT THE SIMMETRY CENTRE
      IRL = MLTPHS(K)
c-----CALCULATE IOF = ALL ATOMS OF THE K-1 PHASES
      IOF = 0
      IF(K.GT.1) THEN
      DO 110 I = 2,K
  110 IOF = IOF + NATOM(I-1)
      ENDIF
      DO 120 ICX = 1,2
      S  =  STH / LAMDA(ICX)
      S2 =  S * S
      S4 = S2 * S2
      CSP(ICX) = 0.0
      DO 130 I = 1,NK
      J = PTC(I+IOF)
      FNUM    =   1.0 + CC(1,J) * S2 + CC(2,J) * S4
      FDEN    = ( 1.0 + CC(3,J) * S2 + CC(4,J) * S4 )**2
      CSP(ICX)=CSP(ICX) + ZEFF(J) * XL(I+IOF,5) * (1.0-(FNUM/FDEN))
  130 CONTINUE
      CSP(ICX) = CSP(ICX) * FLOAT(IRL)
c-----UPDATE S
      S  = 2.0 * S
      S2 =   S * S
c-----COMPUTE ASS1
      ASS1 =   1.0 + 0.75 * HMC * LAMDA(ICX) * S2
c-----COMPUTE BDF
      BDF  = ( 1.0 + 0.50 * HMC * LAMDA(ICX) * S2 ) ** 2
c-----COMPUTE ASS2 = EXP (-DELTA(MU)*D)
C      WHERE:
C      MU       = ABSORPTION COEFFICIENT OF AIR
C      DELTA(MU)= VARIATION OF MU WITH LAMDA, WAS CALCULATED BY RIELLO
C                 USING REGRESSION ANALYSIS ON ESTIMATED MU BY ASSUMING
C                 AIR COMPOSITION 20% O2 AND 80% N2 AT 300¯K AND 1 ATM.:
C                 MU = (E ** 3.40089)*(1.0E-04)*( LAMDA ** 2.79287)
C      D = 17.3 CM. DISTANCE SPECIMENT-COUNTER OR RADIUS OF THE CAMERA
      ASS2 = EXP(-1.5*HMC*17.3*(29.99E-04)*(LAMDA(ICX)**3.79287)*S2)
c-----COMPUTE ASS3
C     FORMULA MUST BE CHANGED FOR OTHER RADIATION AND MONOCHROMATOR
c-----  ASS3 IS A LORENTZIAN FUNCTION
        ASS3=1/(1+GLB(18)*S**GLB(19))
c----- ASS3 IS A GAUSSIAN FUNCTION
C       ASS3=EXP(-GLB(18)*S2)
C     ASS3=0.38*EXP(-4.0*S2)+0.62*EXP(-0.15*S2)
      CSP(ICX) = CSP(ICX) * ASS2 * ASS3  / ( ASS1 * BDF )
  120 CONTINUE
C     COMPUTE CISK = COMPTON INTENSITY SCATTERED BY THE K-TH PHASE
      CISK = RATIO(1) * CSP(1) + RATIO(2) * CSP(2)
      RETURN
      END
      SUBROUTINE DISORDER(K,STH,IDERIV,SDK,DYC,FONDO,DERISO)
      include 'param.inc'
C     PARAMETER(IDSZ=8000,IRS=4096,NATS=16,MSZ=63 ,NOV=2048,NFINAL=1024)
c
c-----THIS SUBROUTINE COMPUTES SDK = THE SCATTERING DUE TO THE THERMAL
C     OR LATTICE DISORDER IN THE SAMPLE FOR K-TH PHASE AT THE IPM POINT,
C     (THEN ALSO IN STEP POINTS THAT DO NOT CONTAIN BRAGG REFLECTIONS),
C     P.RIELLO, G. FAGHERAZZI, D. CLEMENTE AND P.CANTON IN
C     J. APPL. CRYST. (1995) 28,115-120
      REAL    LAMDA
      INTEGER PTR,FONDO
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),xmas(16)
      COMMON /RTMTX/ IVEC(99,192),MLTPHS(99),ICNTPHS(99)
      COMMON/ATFAT/FI2(2)
      DIMENSION DIS(2),STHL2(2),DER(2),DERDIS(NATS,2)
      DIMENSION DERISO(NATS)
c-----COMPUTE IOF = ATOM NUMBER IN THE K-1 PHASE
      IOF = 0
      IF(K.GT.1) THEN
      DO 110 I = 2,K
110   IOF = IOF + NATOM(I-1)
      ENDIF
      DO 111 I=1,NATS
      DERDIS(I,1)=0.0
      DERDIS(I,2)=0.0
111   CONTINUE
c-----IRL = NUMBER OF EQUIVALENT POSITIONS- ALSO THE IDENTITY POSITION
C           AND THE SIMMETRY CENTRE IF PRESENT
      IRL = MLTPHS(K)
c-----NK = NUMBER OF ATOMS IN THE K-TH PHASE
      NK  = NATOM(K)
      DO 120 ICX = 1,2
      STHL2(ICX) = ( STH / LAMDA(ICX) ) ** 2
      IF(FONDO.EQ.1) THEN
c-----COMPUTE FI2 = SUM OF SQUARE SCATTERING FACTORS DUE TO ALL ATOMS
C                   IN THE CELL AT LAMDA(ICX) CORRECTED FOR THE ISOTROPIC
C                    THERMAL FACTORS.
      FI2(ICX) = 0.0
      DO 131 I = 1,NK
      NI = PTR(I+IOF)
      FI = 0.0
      AC(10,NI) = 0.0
      DO 141 II = 1,9,2
141   FI = FI + AC(II,NI)*EXP(-AC(II+1,NI)*STHL2(ICX))
      FI = FI + DFP(NI)
      FI2(ICX) = FI2(ICX) + XL(I+IOF,5) * FLOAT(IRL)*
     $ ( 1.0 - EXP(-XL(I+IOF,4)*2.0*STHL2(ICX)) )* (FI**2+DFPP(NI)**2)
c-----NEXT LINES EVALUATE THE DERIVATES
      IF(LP(I+IOF,4).NE.0.AND.IDERIV.EQ.2) THEN
      DERDIS(I+IOF,ICX) = XL(I+IOF,5) * 2.0 * STHL2(ICX) *
     $ EXP(-XL(I+IOF,4)*2.0*STHL2(ICX))*FLOAT(IRL)*(FI**2+DFPP(NI)**2)
      ELSE
      DERDIS(I+IOF,ICX)=0.0
      ENDIF
131   CONTINUE
      DIS(ICX) = FI2(ICX)
      END IF
      IF(FONDO.EQ.2) THEN
c-----COMPUTE FI2 = SUM OF SQUARE SCATTERING FACTORS DUE TO ALL ATOMS
C                   IN THE CELL AT LAMDA(ICX) CORRECTED FOR THE OVERALL
C                    THERMAL FACTORS.
      FI2(ICX) = 0.0
      DO 130 I = 1,NK
      NI = PTR(I+IOF)
      FI = 0.0
      AC(10,NI) = 0.0
      DO 140 II = 1,9,2
140   FI = FI + AC(II,NI)*EXP(-AC(II+1,NI)*STHL2(ICX))
      FI = FI + DFP(NI)
      FI2(ICX) = FI2(ICX) + XL(I+IOF,5) * (FI**2+DFPP(NI)**2)
130   CONTINUE
      FI2(ICX) = FI2(ICX) * FLOAT(IRL)
      DIS(ICX) = ( 1.0 - EXP(-PAR(K,2)*2.0*STHL2(ICX)) ) * FI2(ICX)
      LK =LPAR(K,2)
      IF(LK.NE.0.AND.IDERIV.EQ.2) THEN
      DER(ICX) = 2.0 * STHL2(ICX) * FI2(ICX) *
     $                   EXP(-PAR(K,2)*2.0*STHL2(ICX))
      ELSE
      DER(ICX) = 0.0
      ENDIF
c
c
      END IF
  120 CONTINUE
c-----COMPUTE SDK = SCATTERING DISORDER DUE TO THE K-TH PHASE
      SDK = RATIO(1) * DIS(1) + RATIO(2) * DIS(2)
c-----COMPUTE  DERIVATIVE OF YC RESPECT TO ISOTROPIC THERMAL FACTORS
C                   XL(I+IOF,4) IN THE K-TH PHASE
      IF(FONDO.EQ.1) THEN
      DO 151 I=1,NK
      IF(LP(I+IOF,4).NE.0.AND.IDERIV.EQ.2) THEN
        DERISO(I+IOF)=RATIO(1) * DERDIS(I+IOF,1) +
     $                        RATIO(2) * DERDIS(I+IOF,2)
       ELSE
        DERISO(I+IOF)=0.0
      END IF
  151 CONTINUE
      END IF
c-----COMPUTE  DERIVATIVE OF YC RESPECT TO OVERALL THERMAL FACTOR
C                   PAR(K,2) IN THE K-TH PHASE
      IF(FONDO.EQ.2) THEN
        IF(LK.NE.0.AND.IDERIV.EQ.2) THEN
        DYC = RATIO(1) * DER(1) + RATIO(2) * DER(2)
        ELSE
        DYC = 0.0
        ENDIF
      ENDIF
      RETURN
      END
      BLOCK DATA COMTAB
C                          *** COMPTON TABLE ***
c-----THE VALUES OF THE COSTANTS FOR  COMPTON SCATTERING ARE TAKEN
C     FROM:
C     V.H. SMITH JR. A.J. THAKKAR AND C. CHAPMAN
C     ACTA CRYST. (1975), A31, 391.
      COMMON/INC/TCS(95,4)
C                    -----H -----
      DATA TCS(1,1),TCS(1,2),TCS(1,3),TCS(1,4)
     $/0.65390,-0.18456,8.2356,12.438/
C                    -----HE-----
      DATA TCS(2,1),TCS(2,2),TCS(2,3),TCS(2,4)
     $/0.72391,-0.21464,9.1019,15.566/
C                    -----LI-----
      DATA TCS(3,1),TCS(3,2),TCS(3,3),TCS(3,4)
     $/26.076,6.8817,26.299,0.88062/
C                    -----BE-----
      DATA TCS(4,1),TCS(4,2),TCS(4,3),TCS(4,4)
     $/16.966,1117.7,40.948,103.99/
C                    -----B -----
      DATA TCS(5,1),TCS(5,2),TCS(5,3),TCS(5,4)
     $/7.2997,272.81,22.693,34.403/
C                    -----C -----
      DATA TCS(6,1),TCS(6,2),TCS(6,3),TCS(6,4)
     $/4.3353,93.125,14.656,14.628/
C                    -----N -----
      DATA TCS(7,1),TCS(7,2),TCS(7,3),TCS(7,4)
     $/4.5051,40.965,11.061,7.3225/
C                    -----O -----
      DATA TCS(8,1),TCS(8,2),TCS(8,3),TCS(8,4)
     $/3.2434,19.377,8.2735,4.0087/
C                    -----F -----
      DATA TCS(9,1),TCS(9,2),TCS(9,3),TCS(9,4)
     $/2.7771,10.031,6.5952,2.3142/
C                     -----NE-----
      DATA TCS(10,1),TCS(10,2),TCS(10,3),TCS(10,4)
     $/3.1880,5.6531,5.7748,1.3790/
C                     -----NA-----
      DATA TCS(11,1),TCS(11,2),TCS(11,3),TCS(11,4)
     $/16.285,45.169,13.167,5.4569/
C                     -----MG-----
      DATA TCS(12,1),TCS(12,2),TCS(12,3),TCS(12,4)
     $/115.98,3227.0,70.762,77.161/
C                     -----AL-----
      DATA TCS(13,1),TCS(13,2),TCS(13,3),TCS(13,4)
     $/107.62,2740.9,67.151,59.872/
C                     -----SI-----
      DATA TCS(14,1),TCS(14,2),TCS(14,3),TCS(14,4)
     $/80.719,1508.4,52.054,36.898/
C                     -----P -----
      DATA TCS(15,1),TCS(15,2),TCS(15,3),TCS(15,4)
     $/55.933,703.48,37.402,20.872/
C                     -----S -----
      DATA TCS(16,1),TCS(16,2),TCS(16,3),TCS(16,4)
     $/43.695,416.47,30.100,13.430/
C                     -----CL-----
      DATA TCS(17,1),TCS(17,2),TCS(17,3),TCS(17,4)
     $/34.592,252.66,24.481,8.8480/
C                     -----AR-----
      DATA TCS(18,1),TCS(18,2),TCS(18,3),TCS(18,4)
     $/27.153,150.48,19.697,5.8774/
C                     -----K -----
      DATA TCS(19,1),TCS(19,2),TCS(19,3),TCS(19,4)
     $/30.545,193.71,21.865,6.5615/
C                     -----CA-----
      DATA TCS(20,1),TCS(20,2),TCS(20,3),TCS(20,4)
     $/36.945,289.98,25.950,8.1578/
C                     -----SC-----
      DATA TCS(21,1),TCS(21,2),TCS(21,3),TCS(21,4)
     $/33.212,232.59,23.472,6.8604/
C                     -----TI-----
      DATA TCS(22,1),TCS(22,2),TCS(22,3),TCS(22,4)
     $/29.486,182.05,20.962,5.7117/
C                     -----V -----
      DATA TCS(23,1),TCS(23,2),TCS(23,3),TCS(23,4)
     $/26.193,143.17,18.683,4.8164/
C                     -----CR-----
      DATA TCS(24,1),TCS(24,2),TCS(24,3),TCS(24,4)
     $/18.929,72.508,13.819,3.0415/
C                     -----MN-----
      DATA TCS(25,1),TCS(25,2),TCS(25,3),TCS(25,4)
     $/21.292,94.848,15.188,3.6204/
C                     -----FE-----
      DATA TCS(26,1),TCS(26,2),TCS(26,3),TCS(26,4)
     $/19.276,77.511,13.788,3.1076/
C                     -----CO-----
      DATA TCS(27,1),TCS(27,2),TCS(27,3),TCS(27,4)
     $/17.616,64.657,12.621,2.7056/
C                     -----NI-----
      DATA TCS(28,1),TCS(28,2),TCS(28,3),TCS(28,4)
     $/16.346,55.896,11.687,2.4272/
C                     -----CU-----
      DATA TCS(29,1),TCS(29,2),TCS(29,3),TCS(29,4)
     $/12.851,33.834,9.3436,1.7194/
C                     -----ZN-----
      DATA TCS(30,1),TCS(30,2),TCS(30,3),TCS(30,4)
     $/14.619,45.462,10.338,2.0820/
C                     -----GA-----
      DATA TCS(31,1),TCS(31,2),TCS(31,3),TCS(31,4)
     $/17.176,64.6920,11.820,2.5707/
C                     -----GE-----
      DATA TCS(32,1),TCS(32,2),TCS(32,3),TCS(32,4)
     $/19.343,83.421,13.080,2.9211/
C                     -----AS-----
      DATA TCS(33,1),TCS(33,2),TCS(33,3),TCS(33,4)
     $/19.471,84.720,13.126,2.8165/
C                     -----SE-----
      DATA TCS(34,1),TCS(34,2),TCS(34,3),TCS(34,4)
     $/20.359,92.997,13.650,2.8488/
C                     -----BR-----
      DATA TCS(35,1),TCS(35,2),TCS(35,3),TCS(35,4)
     $/20.316,92.534,13.623,2.7028/
C                     -----KR-----
      DATA TCS(36,1),TCS(36,2),TCS(36,3),TCS(36,4)
     $/18.986,80.312,12.812,2.3584/
C                     -----RB-----
      DATA TCS(37,1),TCS(37,2),TCS(37,3),TCS(37,4)
     $/21.269,101.61,14.209,2.6226/
C                     -----SR-----
      DATA TCS(38,1),TCS(38,2),TCS(38,3),TCS(38,4)
     $/24.497,135.93,16.192,3.0174/
C                     -----Y -----
      DATA TCS(39,1),TCS(39,2),TCS(39,3),TCS(39,4)
     $/24.539,136.20,16.237,2.9014/
C                     -----ZR-----
      DATA TCS(40,1),TCS(40,2),TCS(40,3),TCS(40,4)
     $/23.692,126.53,15.732,2.6710/
C                     -----NB-----
      DATA TCS(41,1),TCS(41,2),TCS(41,3),TCS(41,4)
     $/19.178,81.654,12.952,19.678/
C                     -----MO-----
      DATA TCS(42,1),TCS(42,2),TCS(42,3),TCS(42,4)
     $/17.443,67.055,11.867,1.6907/
C                     -----TC-----
      DATA TCS(43,1),TCS(43,2),TCS(43,3),TCS(43,4)
     $/18.256,73.766,12.363,1.7500/
C                     -----RU-----
      DATA TCS(44,1),TCS(44,2),TCS(44,3),TCS(44,4)
     $/15.431,51.934,10.601,1.3610/
C                     -----RH-----
      DATA TCS(45,1),TCS(45,2),TCS(45,3),TCS(45,4)
     $/14.236,43.914,9.8385,1.1954/
C                     -----PD-----
      DATA TCS(46,1),TCS(46,2),TCS(46,3),TCS(46,4)
     $/12.125,31.377,8.4946,0.94341/
C                     -----AG-----
      DATA TCS(47,1),TCS(47,2),TCS(47,3),TCS(47,4)
     $/11.910,30.295,8.3388,0.90612/
C                     -----CD-----
      DATA TCS(48,1),TCS(48,2),TCS(48,3),TCS(48,4)
     $/12.283,32.428,8.5514,0.93212/
C                     -----IN-----
      DATA TCS(49,1),TCS(49,2),TCS(49,3),TCS(49,4)
     $/12.893,36.009,8.9130,0.98158/
C                     -----SN-----
      DATA TCS(50,1),TCS(50,2),TCS(50,3),TCS(50,4)
     $/13.547,40.031,9.3037,1.0326/
C                     -----SB-----
      DATA TCS(51,1),TCS(51,2),TCS(51,3),TCS(51,4)
     $/13.486,39.657,9.2631,1.0057/
C                     -----TE-----
      DATA TCS(52,1),TCS(52,2),TCS(52,3),TCS(52,4)
     $/14.490,46.190,9.8671,1.0869/
C                     -----I -----
      DATA TCS(53,1),TCS(53,2),TCS(53,3),TCS(53,4)
     $/14.975,49.498,10.161,1.1116/
C                     -----XE-----
      DATA TCS(54,1),TCS(54,2),TCS(54,3),TCS(54,4)
     $/14.885,48.883,10.102,1.0775/
C                     -----CS-----
      DATA TCS(55,1),TCS(55,2),TCS(55,3),TCS(55,4)
     $/16.159,57.996,10.886,1.1773/
C                     -----BA-----
      DATA TCS(56,1),TCS(56,2),TCS(56,3),TCS(56,4)
     $/17.835,71.158,11.921,13.128/
C                     -----LA-----
      DATA TCS(57,1),TCS(57,2),TCS(57,3),TCS(57,4)
     $/18.225,74.363,12.167,1.3191/
C                     -----CE-----
      DATA TCS(58,1),TCS(58,2),TCS(58,3),TCS(58,4)
     $/17.590,69.295,11.741,1.2501/
C                     -----PR-----
      DATA TCS(59,1),TCS(59,2),TCS(59,3),TCS(59,4)
     $/16.613,61.934,11.071,1.1678/
C                     -----ND-----
      DATA TCS(60,1),TCS(60,2),TCS(60,3),TCS(60,4)
     $/16.148,58.598,10.747,1.1195/
C                     -----PM-----
      DATA TCS(61,1),TCS(61,2),TCS(61,3),TCS(61,4)
     $/15.741,55.777,10.458,10.774/
C                     -----SM-----
      DATA TCS(62,1),TCS(62,2),TCS(62,3),TCS(62,4)
     $/15.359,53.212,10.182,1.0396/
C                     -----EU-----
      DATA TCS(63,1),TCS(63,2),TCS(63,3),TCS(63,4)
     $/15.027,51.053,9.9393,1.0059/
C                     -----GD-----
      DATA TCS(64,1),TCS(64,2),TCS(64,3),TCS(64,4)
     $/14.888,50.133,9.8426,0.97728/
C                     -----TB-----
      DATA TCS(65,1),TCS(65,2),TCS(65,3),TCS(65,4)
     $/14.449,47.248,9.5464,0.93253/
C                     -----DY-----
      DATA TCS(66,1),TCS(66,2),TCS(66,3),TCS(66,4)
     $/13.795,43.137,9.1027,0.87906/
C                     -----HO-----
      DATA TCS(67,1),TCS(67,2),TCS(67,3),TCS(67,4)
     $/13.447,41.045,8.8604,0.84486/
C                     -----ER-----
      DATA TCS(68,1),TCS(68,2),TCS(68,3),TCS(68,4)
     $/13.146,39.287,8.6467,0.81464/
C                     -----TM-----
      DATA TCS(69,1),TCS(69,2),TCS(69,3),TCS(69,4)
     $/12.882,37.803,8.4552,0.78869/
C                     -----YB-----
      DATA TCS(70,1),TCS(70,2),TCS(70,3),TCS(70,4)
     $/12.670,36.643,8.2968,0.76627/
C                     -----LU-----
      DATA TCS(71,1),TCS(71,2),TCS(71,3),TCS(71,4)
     $/12.576,36.129,8.2264,0.74720/
C                     -----HF-----
      DATA TCS(72,1),TCS(72,2),TCS(72,3),TCS(72,4)
     $/12.392,35.077,8.1062,0.72092/
C                     -----TA-----
      DATA TCS(73,1),TCS(73,2),TCS(73,3),TCS(73,4)
     $/11.966,32.657,7.8420,0.67752/
C                     -----W -----
      DATA TCS(74,1),TCS(74,2),TCS(74,3),TCS(74,4)
     $/11.429,29.720,7.5123,0.62773/
C                     -----RE-----
      DATA TCS(75,1),TCS(75,2),TCS(75,3),TCS(75,4)
     $/10.841,26.660,7.1515,0.57682/
C                     -----OS-----
      DATA TCS(76,1),TCS(76,2),TCS(76,3),TCS(76,4)
     $/10.553,25.229,6.9734,0.54767/
C                     -----IR-----
      DATA TCS(77,1),TCS(77,2),TCS(77,3),TCS(77,4)
     $/10.255,23.781,6.7889,0.51891/
C                     -----PT-----
      DATA TCS(78,1),TCS(78,2),TCS(78,3),TCS(78,4)
     $/9.1444,18.754,6.1116,0.43980/
C                     -----AU-----
      DATA TCS(79,1),TCS(79,2),TCS(79,3),TCS(79,4)
     $/8.6692,16.793,5.8176,0.40417/
C                     -----HG-----
      DATA TCS(80,1),TCS(80,2),TCS(80,3),TCS(80,4)
     $/8.8187,17.414,5.9035,0.40768/
C                     -----TL-----
      DATA TCS(81,1),TCS(81,2),TCS(81,3),TCS(81,4)
     $/9.0707,18.476,6.0519,0.41721/
C                     -----PB-----
      DATA TCS(82,1),TCS(82,2),TCS(82,3),TCS(82,4)
     $/9.3622,19.739,6.2252,0.42863/
C                     -----BI-----
      DATA TCS(83,1),TCS(83,2),TCS(83,3),TCS(83,4)
     $/9.5574,20.605,6.3414,0.43371/
C                     -----PO-----
      DATA TCS(84,1),TCS(84,2),TCS(84,3),TCS(84,4)
     $/9.9256,22.287,6.5624,0.44878/
C                     -----AT-----
      DATA TCS(85,1),TCS(85,2),TCS(85,3),TCS(85,4)
     $/10.267,23.902,6.7681,0.46167/
C                     -----RN-----
      DATA TCS(86,1),TCS(86,2),TCS(86,3),TCS(86,4)
     $/10.425,24.663,6.8643,0.46342/
C                     -----FR-----
      DATA TCS(87,1),TCS(87,2),TCS(87,3),TCS(87,4)
     $/11.077,27.947,7.2583,0.49341/
C                     -----RA-----
      DATA TCS(88,1),TCS(88,2),TCS(88,3),TCS(88,4)
     $/11.906,32.420,7.7600,0.53275/
C                     -----AC-----
      DATA TCS(89,1),TCS(89,2),TCS(89,3),TCS(89,4)
     $/12.372,35.070,8.0440,0.55039/
C                     -----TH-----
      DATA TCS(90,1),TCS(90,2),TCS(90,3),TCS(90,4)
     $/12.708,37.029,8.2501,0.56015/
C                     -----PA-----
      DATA TCS(91,1),TCS(91,2),TCS(91,3),TCS(91,4)
     $/11.931,32.573,7.7684,0.51436/
C                     -----U------
      DATA TCS(92,1),TCS(92,2),TCS(92,3),TCS(92,4)
     $/11.596,30.746,7.5590,0.49185/
C                     -----NP-----
      DATA TCS(93,1),TCS(93,2),TCS(93,3),TCS(93,4)
     $/11.277,29.055,7.3579,0.47097/
C                     -----PU-----
      DATA TCS(94,1),TCS(94,2),TCS(94,3),TCS(94,4)
     $/10.569,25.479,6.9123,0.43325/
C                     -----AM-----
      DATA TCS(95,1),TCS(95,2),TCS(95,3),TCS(95,4)
     $/10.243,23.918,6.7044,0.41402/
c
      END
      SUBROUTINE INPTR
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 )
      REAL LAMDA, lamdam
      INTEGER PTR,fondo
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
C      COMMON/AIRPAR/AI1,AI2,AI3,AI4,AI5,AI6,AI7,AI8,AI9,AI10
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      CHARACTER*72 DATAID
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      CHARACTER SYMB(99,20)*1,SPG(20)*1
C      CHARACTER TITLE*70,PHSNM(8)*50
      CHARACTER TITLE*70,PHSNM(99)*50
C      CHARACTER*70 TITLE
C      character*50 PHSNM(99)
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      COMMON/BLNK1/NJNK(22),FILL2(751)
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),xmas(16)
      CHARACTER*5 DATE
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/CELLX/AA,B,C,ALPHA,BETA,GAMMA,AL(3,3)
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      COMMON/MULTIP/TMASSA(99),MLTPHASE,xmltp(99),MURT(NATS),SAQF(99)  
     *,wtis(99)
      common/codebck/ibckcode     
      common/simoper/isimop
      DIMENSION XRYZ(10),BACK(6),FBACK(6)
      EQUIVALENCE (GLB(2),BACK(1))
     1    ,(AGLB(2),FBACK(1))
      common/sizestrain/sizeG(15),strainG(15),sizeL(15),strainL(15)
     *,siz(15),strain(15),NsizeStrain
      common/convert/icnvt   
      common/maxint/xmaxint     
      CHARACTER*26 LOWER,UPPER       
      DATA UPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LOWER/'abcdefghijklmnopqrstuvwxyz'/
C 2 lines below were commented to include the code from Ian madsen
C      DATA XRYZ(1),XRYZ(2),XRYZ(3),XRYZ(4),XRYZ(5)
C     * /2.28962,1.93597,1.54051,0.70926,0.556363/ 
C Adding the code from Ian Madsen (10 feb 98)  (Au,Cr,Fe,Co,Cu,
C                                               Mo,Ag,Ta,W ,Au)
      DATA XRYZ(1), XRYZ(2),  XRYZ(3),  XRYZ(4),  XRYZ(5)
     *   /2.748510,2.289620, 1.935970, 1.788965, 1.540520/
      DATA XRYZ(6), XRYZ(7),  XRYZ(8),  XRYZ(9), XRYZ(10)
     * /  0.709260,0.559360, 0.215947, 0.209010, 0.180195/
cc      icnvt=0
C line 1
      READ(5,1,END=99999)TITLE
1001  FORMAT(1X,A70)
c1     FORMAT(BZ,A70)
1     FORMAT(A70)
C line 2
      READ(5,2,END=99999)JOBTYP,NPROF,NPHASE,NBCKGD,NEXCRG,NSCAT,INSTRM,
     *    IPREF,iasym,IABSR,idata,isphase,i2d94  
2     FORMAT(BZ,13I4)
              if(nprof.lt.0)then
                  icnvt=1
                  write(6,8000)
                  write(*,8000)
                  call conv94
8000  format(/' $ INPUT FILE CONVERSION FROM DBWS9411 TO DBWS9807',
     */,' $ CHECK ALL PARAMETERS',/,' $ SPECIAL ATTENTION REQUIRED TO'
     *' ATOM MULTIPLICITIES AND OCCUPANCIES',/,
     *33H $ See USER's GUIDE,  LINE 11-41 ,/,
     *'                  GOOD LUCK!',//)
C        1234567890123456789012345678901234567890
                  GOTO 88088
              else
                  icnvt = 0
                    IF(ISPHASE.GT.NPHASE) THEN
                       WRITE(6,7334)
7334  FORMAT(' Internal Standard Phase does not exist.',
     *       ' Check its number in line 2, column 12.',/
     *       ' No Internal Standard will be used in the QPA')
                       isphase = 0
                    ENDIF
              endif
C     open file to write +/- dbws9006 & dbws9411 format (readable by ATOMS and ZORTEP)
      IF(I2D94.EQ.1)open(53,file='ICF94.ICF',status='unknown')
C
       if (nprof.eq.9) then
            nprof = 7
            NsizeStrain = 9
       end if
        if(NBCKGD.eq.-1)then
            ibckcode=nbckgd
            nbckgd = 0  
            ibdg   = 0
        else
            fondo  = 0
            ibgd   = 1
C            iax    = 0
            ibckcode=nbckgd
        endif     
      IF(NPHASE.EQ.0)NPHASE=1 
      INSTRM=INSTRM+1
      JOBTYP=JOBTYP+1
      NPROF=NPROF+1 
      WRITE(6,3)
3     FORMAT('1RIETVELD ANALYSIS PROGRAM DBWS-9807, RELEASE ',
     *' 18.12.99,',/,'   COPYRIGHT 1998 BY R.A.YOUNG,',/,
     *'   SEE APPENDIX B ON COPYRIGHT-AND-FAIR-USE IN USERS GUIDE',
     *' FOR CONDITIONS,',/,
     *'   THIS PROGRAM SHOULD BE CITED IN ALL WORK TO WHICH IT',
     *' HAS CONTRIBUTED AS AN UPGRADE OF:',/,
     *'  "DBWS-9411 - an upgrade of the DBWS programs for Rietveld',/, 
     *'   Refinement with PC and mainframe computers",',/,
     *'   J. Appl. Cryst., Vol. 28, 366-7; 1995')       
      write(6,311) idsz,irs,nats,msz,nov
311   format(' PROGRAM PARAMETERS:',/,' IDSZ=',i5,4x,'IRS=',i5,4x,
     +          'NATS=',i4,5x,'MSZ=',i3,5x,'NOV=',i5)
      WRITE(6,1001)TITLE
      WRITE(*,3)
      write(*,311) idsz,irs,nats,msz,nov
      WRITE(*,1001)TITLE
      IF(JOBTYP.EQ.1)WRITE(6,4)
4     FORMAT(15H FOR X-RAY DATA )
      IF(JOBTYP.EQ.1.AND.INSTRM.EQ.2) WRITE (6,9490)
9490  FORMAT (40H COLLECTED IN SYNCHROTRON AT NSLS OR SRS)
      IF(JOBTYP.EQ.2)WRITE(6,5)
      IF(JOBTYP.EQ.2.AND.INSTRM.EQ.2) WRITE (6,9493)
9493  FORMAT(1X,37H VARYING NO. OF COUNTERS AT EACH STEP )
      IF(JOBTYP.EQ.3)WRITE(6,6)
6     FORMAT(25H PATTERN CALCULATION,XRAY  )
      IF(JOBTYP.EQ.4)WRITE(6,7)
7     FORMAT(28H PATTERN CALCULATION,NEUTRON  )
      IF(JOBTYP.LT.1.OR.JOBTYP.GT.4)STOP 7777
5     FORMAT(43H FOR NEUTRON DATA, NUCLEAR INTENSITIES ONLY )
      if(idata.eq.0)write(6,2120)
      if(idata.eq.1)write(6,2121)
      if(idata.eq.2)write(6,2122)
      if(idata.eq.3)write(6,2123)
      if(idata.eq.4)write(6,2124)
      if(idata.eq.5)write(6,2125)
      if(idata.eq.6)write(6,2126)
2120  format(' READ DATA IN TRADITIONAL DBWS FORMAT')      
2121  format(' READ DATA IN FREE FORMAT')      
2122  format(' READ DATA IN GSAS STD FORMAT')      
2123  format(' READ DATA IN PHILIPS UDF FORMAT')
2124  format(' READ DATA IN SCINTAG TXT FORMAT')
2125  format(' READ DATA IN SIEMENS UXD FORMAT')
2126  format(' READ DATA IN RIGAKU ASC FORMAT')
      WRITE(6,9)NPHASE,NEXCRG,NSCAT
9     FORMAT(19H NUMBER OF PHASES=  ,I4,/,
     * 29H NUMBER OF EXCLUDED REGIONS=  ,I4,/,
     *  28H NUMBER OF SCATTERING SETS=  ,I4)
      IF(NBCKGD-1)110,111,112 
110   WRITE(6,10)
10    FORMAT(27H BACKGROUND TO BE REFINED   )
      GOTO 113
111   WRITE(6,11)
11    FORMAT(37H BACKGROUND DATA TO BE READ FROM FILE )
      GOTO 113
112   WRITE(6,12)NBCKGD
12    FORMAT(40H BACKGROUND CORRECTION BY INTERPOLATION  ,
     * 12HBETWEEN THE  ,I4,13H POINTS GIVEN )
c113   WRITE(6,13) NPROF-1
c13    FORMAT(5H THE  ,I4,32HTH PROFILE FUNCTION WAS SELECTED )
113   CONTINUE
      IF(NPROF-1.EQ.0) THEN
       WRITE(6,14001) NPROF-1
14001  FORMAT (' GAUSSIAN PROFILE, NPROF = ',I4)
       ELSE IF(NPROF-1.EQ.1) THEN
       WRITE(6,14002)NPROF-1
14002  FORMAT (' LORENTZIAN (CAUCHY) PROFILE, NPROF = ',I4)
       ELSE IF(NPROF-1.EQ.2) THEN
       WRITE(6,14003) NPROF-1
14003  FORMAT (' MOD 1 LORENTZIAN PROFILE, NPROF = ',I4)
       ELSE IF(NPROF-1.EQ.3) THEN
       WRITE(6,14004) NPROF-1
14004  FORMAT (' MOD 2 LORENTZIAN PROFILE, NPROF = ',I4)
       ELSE IF(NPROF-1.EQ.4) THEN
       WRITE(6,14005) NPROF-1
14005  FORMAT (' SPLIT PEARSON VII PROFILE, NPROF =  ',I4)
       ELSE IF(NPROF-1.EQ.5) THEN
       WRITE(6,14006) NPROF-1
14006  FORMAT (' PSEUDO-VOIGT (PV) PROFILE, NPROF = ',I4)
       ELSE IF(NPROF-1.EQ.6) THEN
       WRITE(6,14007) NPROF-1
14007  FORMAT (' PEARSON VII PROFILE, NPROF = ',I4)
       ELSE IF(NPROF-1.EQ.7) THEN
       WRITE(6,14008)NPROF-1
14008  FORMAT (' THOMPSON-COX-HASTINGS (PV) PROFILE, NPROF = ',I4)
      END IF
      IF (IPREF.EQ.0) THEN
       WRITE(6,7445)
      ELSE
       WRITE(6,7446)
      ENDIF
7445  FORMAT(1X,'IPREF=0, UDA-RIETVELD PREFERRED ORIENTATION FUNCTION')
7446  FORMAT(1X,'IPREF=1, MARCH-DOLLASE PREFERRED ORIENTATION FUNCTION')
C  SURFACE ROUGHNESS
C !cp jun 97 7448 introduced line below
      if (iabsr.eq.1) write(6,9004)
9004  format(1x,'IABSR=1, CORRECTION OF SURFACE ROUGHNESS BY YOUNG')
      if (iabsr.eq.2) write(6,9002)
9002  format(1x,'IABSR=2, CORRECTION OF SURFACE ROUGHNESS BY, SPARKS, ET
     * AL.')
      if (iabsr.eq.3) write(6,9003)
9003  format(1x,'IABSR=3, CORRECTION OF SURFACE ROUGHNESS BY SUORTTI')
      if (iabsr.eq.4) write(6,9001)
9001  format(1x,'IABSR=4, CORRECTION OF SURFACE ROUGHNESS BY PITSCHKE, H
     *ERMANN, AND MATTERN')
C  asymmetry correction  Test for asymmetry model included !cp ap 97
      if (iasym.eq.0) write(6,8049)
8049  format(1x,'IASYM=0, Usual Rietveld Asymmetry')      
      if (iasym.eq.1) write(6,8048)
8048  format(1x,'IASYM=1, Asymmetry by Riello et al.:',
     *' Powder Diffraction,10,204-206,1995')
c-----FONDO=0 :BKG EVALUATED USING STANDARD METHODS
c-----FONDO=1 :BKG EVALUATED USING THE ISOTROPIC THERMAL FACTORS
c-----FONDO=2 :BKG EVALUATED USING THE OVERAL THERMAL FACTORS
c-----AIR SCATTERING,IAXX= 1 AIR SCATTERING ADDED TO THE BACKGROUND IF SCAIR<>0
C                    IAXX= -1 AIR SCATT. EVALUATED AND SUBTRACTED FROM DATA
C                    IAXX= 0 AIR SCATT NOT EVALUATED and ibgd = 1
c-----LINEAR ABSORP.CORR., IAS= 1  DATA CORRECTED FOR LINEAR ABSORPTION
c
c-----FI = AIR FRACTION INSIDE THE SAMPLE   
        if(ibckcode.eq.-1) then
C line 2.2
          read(5,202,END=99999)ias, fondo
202     FORMAT(BZ,2I4)
        endif
C          read(5,202,END=99999)fondo,iaxX   
C              IF(IAXX)5330,5331,5332
C      !cp jul 97 start
c5330          IAX = 1      
C              ibgd = 0
C              GOTO 5333
c5331          IBGD = 1
C              iax = 0               !not used in calculations !cp jul 97
C              GOTO 5333
c5332          IAX = 0
C              ibgd = 0    
5333      if (ibgd.eq.1) then
          write(6,7447)
7447      format(' NO AMOPHOUS AND COMPTON CORRECTION TO THE BGD')
          goto 7448    
        end if
C      IF(IAx.EQ.1) WRITE(6,7456)
c7456  FORMAT(1X,'IAX = 1, AIR SCATTERING SUBTRACTED FROM THE DATA')
      IF(IAS.EQ.0) WRITE(6,7458)
7458  FORMAT(1X,'IAS = 0, NO LINEAR ABSORPTION CORRECTION')
      IF(IAS.EQ.1) WRITE(6,7460)
7460  FORMAT(1X,'IAS = 1, LINEAR ABSORPTION CORRECTION IS APPLIED')
      IF(FONDO.EQ.1) WRITE(6,7461)
7461  FORMAT(1X,'FONDO = 1,ISOTROPIC B FACTOR USED FOR BKG EVALUATION')
      IF(FONDO.EQ.2) WRITE(6,7462)
7462  FORMAT(1X,'FONDO = 2,OVERALL Q USED FOR BKG EVALUATION')
C !cp jul 97 stop
C line 3
c7448  READ(5,14,END=99999)IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2
C     * ,IPLST,IPLOSS,IPLCAL,IPLPOL,IPLCOM,IPLDIS,IPLAM,IPLAIR,IPBIG
c14    FORMAT(BZ,19I1)
C line 3
7448  READ(5,14,END=99999)IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2
     * ,IPLST,IPLOSS,IPLCAL,IPLPOL,IPLCOM,IPLDIS,IPLAM,IPBIG
        if(ibckcode.ne.-1) then
            IPLOSS=0
            IPLCAL=0
            IPLPOL=0
            IPLCOM=0
            IPLDIS=0
            IPLAM =0
C            IPLAIR=0
            IPBIG =0
        endif
          isimop=0
14    FORMAT(BZ,4(5I1,1x) )
      IF(IOT.NE.0)WRITE(6,15) 
15    FORMAT(46H OUTPUT OBSERVED AND CALCULATED INTENSITIES ON,
     * 11H LAST CYCLE )
      IF(IPL.NE.0)WRITE(6,16) 
16    FORMAT(27H GENERATE LINE PRINTER PLOT )
C ipl2 is UNUSED in this release
C      IF(IPL2.NE.0)WRITE(6,17)
c17    FORMAT(22H GENERATE CALCOMP PLOT )
      IF(IPL2.NE.0)write(6,167)
167   format(' Generate file with obs, calc, dif intensiies, weighted di
     !ference, and Bragg peak positions')
      IF(IPLST.EQ.1)WRITE(6,717)
717   FORMAT(24H GENERATE PARAMETER LIST )
      IF(IPLST.EQ.2)WRITE(6,718)
718   FORMAT(24H GENERATE PARAMETER LIST,/,
     *' GENERATE PARAMETERS AND STD.DEV. IN FINAL CYCLE FOR DATA BASE')
      IF(IPC.EQ.1)WRITE(6,181) 
      IF(IPC.EQ.2)WRITE(6,182) 
      IF(IPC.EQ.3)WRITE(6,183) 
181    FORMAT(20H OUTPUT INTENSITIES )
C              123456789*123456789*123456789*123456789*123456789*
182    FORMAT(' OUTPUT ABSOLUTE VALUES OF STRUCTURE FACTORS +',
     *       ' PHASE ANGLE' )
183    FORMAT(47H OUTPUT A AND B (CALC + OBS) STRUCTURE FACTORS )
      IF(MAT.NE.0.AND.JOBTYP.LT.3)WRITE(6,19)
19    FORMAT(26H OUTPUT CORRELATION MATRIX )
      IF(NXT.NE.0.AND.JOBTYP.LT.3)WRITE(6,20)
20    FORMAT(24H GENERATE NEW INPUT FILE )
      IF(LST1.NE.0)WRITE(6,21)
21    FORMAT(22H PRINT REFLECTION LIST )
      IF(NPHASE.EQ.1)LST3=0
      IF(LST3.NE.0)WRITE(6,8) 
8     FORMAT(' Print merged reflection list')
      IF(LST2.NE.0)WRITE(6,22)
22    FORMAT(' Print corrected data') 
C line 4
      READ(5,23,END=99999)LAMDA,RATIO(2),BKPOS,WDT,CTHM,TMV,RLIM,sw
23    FORMAT(BZ,9F8.0)
C        !cp ap 21 97
      RATIO(1) = 1.0
      LAMDAM=(RATIO(1)*LAMDA(1)+RATIO(2)*LAMDA(2))/(RATIO(1)+RATIO(2))
      WRITE(6,24)LAMDA(1),LAMDA(2),LAMDAM
24    FORMAT(15H WAVELENGTHS=  ,2F9.5,16H  LAMDA MEAN =  ,F9.5)
C  stop !cp ap 21 97
      WRITE(6,223)RATIO(2)
223   FORMAT(23H ALPHA2:ALPHA1 RATIO =  ,F8.4)
C      DO 109 IX=1,5     (changed to 10 wavelengths)
      DO 109 IX=1,10
C      IF(LAMDA(1)+0.1.GT.XRYZ(IX))GOTO 119
      IF(1.03*LAMDA(1).GT.XRYZ(IX))GOTO 119
109     CONTINUE
119   IXRAY=IX
      WRITE(6,26) WDT
26    FORMAT(23H BASE OF PEAK = 2.0*HW*,F8.2)
      WRITE(6,27)CTHM
27    FORMAT(27H MONOCHROMATOR CORRECTION =,F8.4) 
C      WRITE(6,28)TMV
c28    FORMAT(28H ABSORPTION CORRECTION=-0.1*,F8.4,11H*SIN 2THETA)
      WRITE(6,28)TMV,sw
28    FORMAT(1X,'ABSORPTION CORRECTION COEFFICIENT = ',F8.4,' CM-1',/,
     $       1X,'SLAB-WIDTH = ',F8.4,' CM.')
      if (iasym.eq.0)then
        write(6,290)rlim
      else      
        WRITE(6,291)90.0 - RLIM,90.0 + rlim
290     FORMAT(' RIETVELD ASYMMETRY CORRECTION FOR ANGLES LESS THAN ',
     *    F8.3, 8H DEGREES )   
291    FORMAT(43H ASYMMETRY CORRECTION FOR ANGLES LESS THAN ,F8.3, 
     *  8H DEGREES,
     * /25x,18H AND GREATER THAN ,F8.3,
     * 8H DEGREES)  
C     * /,1X,'AIR FRACTION INSIDE THE SAMPLE= ',F8.4)  
      end if
C line 5
      READ(5,30,END=99999)MCYCLE,EPS,RELAX,THMIN,STEP,THMAX 
30    FORMAT(BZ,1I4,5F4.0,3F8.0)
      IF(JOBTYP.GT.2)MCYCLE=1 
      ICYRUN = MCYCLE
      WRITE(6,31)MCYCLE
31    FORMAT(20H NUMBER OF CYCLES =  ,I4)
      WRITE(6,32)RELAX
32    FORMAT(19H RELAXATION FACTORS,/,17H FOR COORDINATES=,F5.2,
     * /,36H FOR ANISOTRPIC TEMPERATURE FACTORS= ,F5.2,/,
     *  22H FOR FWHM PARAMETERS= ,F5.2,/,
     *   23H FOR LATTICE CONSTANTS= ,F5.2)
      WRITE(6,33)EPS
33    FORMAT(11H EPS-VALUE=,F6.1)
      IF(NBCKGD.LT.2)GOTO 120 
C line 6(*)
      READ(5,34,END=99999)(POS(I),BCK(I),I=1,NBCKGD)
34    FORMAT(BZ,2F8.2)
      WRITE(6,35)(POS(I),BCK(I),I=1,NBCKGD)
C     WRITE(*,35)(POS(I),BCK(I),I=1,NBCKGD)
35    FORMAT(12H BACKGROUND ,/,9H POSITION,4X,9HINTENSITY,/,
     *  100(2F9.4,/)) 
120   IF(NEXCRG.LE.0)GOTO 122 
C line 7(*)
      READ(5,34,END=99999)(ALOW(I),AHIGH(I),I=1,NEXCRG)
      WRITE(6,37)(ALOW(I),AHIGH(I),I=1,NEXCRG)
37    FORMAT(18H EXCLUDED REGIONS /5H FROM,5X,2HTO/,20(2F9.4,/))
122   IF(NSCAT.LE.0)GOTO 124
      DO 125 I=1,NSCAT
       IF(JOBTYP.EQ.2.OR.JOBTYP.EQ.4) THEN
C line 8.1 XRD(*)
         READ(5,3838,END=99999)NAM(I),DFP(I),XMAS(I) 
3838       FORMAT(BZ,A4,2F8.0)
         GOTO 125
       ENDIF
C line 8.1 ND(*)
      READ(5,38,END=99999)NAM(I),DFP(I),DFPP(I),XMAS(I) 
38      FORMAT(BZ,A4,3F8.0)
      K=0
C line 8.2 XRD(*)
126     READ(5,39,END=99999)(AC(J,I),J=1,9)
39      FORMAT(BZ,9F8.0)
      IF(AC(1,I).EQ.-100.)CALL COEF(I,K)
      IF(AC(3,I).NE.0.)GOTO 125
      K=K+1
      POSI(K)=AC(1,I)
      SCAT(K)=AC(2,I)
      IF(K.LE.29)GOTO126
      WRITE(6,40) 
40      FORMAT(34H TOO MANY SCATTERING TABLE ENTRIES )
      STOP 7700
125     CONTINUE
      IF(JOBTYP.EQ.1.OR.JOBTYP.EQ.3)GOTO 129
      WRITE(6,422)
422   FORMAT(19H SCATTERING LENGTHS)
      WRITE(6,423)(NAM(I),DFP(I),I=1,NSCAT)
423   FORMAT(5H FOR ,A4,5H     ,F10.6)
      DO 424 I=1,NSCAT
      DO 425 J=1,9
425       AC(J,I)=0.0
424     DFPP(I)=0.0
      GOTO 124
129   WRITE(6,42)
42    FORMAT(12H FORMFACTORS )
      WRITE(6,41)(NAM(I),DFP(I),DFPP(I),(AC(J,I),J=1,9),I=1,NSCAT)
41    FORMAT(5H FOR ,A4,5H DFP=,F10.6,6H DFPP=,F10.6,/,
     * 14H COEFFICIENTS= ,9F10.6)
124   CONTINUE
C line 9
      READ(5,43,END=99999)MAXS
43    FORMAT(BZ,I8) 
      if(maxs.gt.msz) then
      write(6,7733)
      write(*,7733)
7733  format(
     * /, ' * YOU HAVE DECLARED MORE CODEWORDS THAN WILL FIT INTO *'
     * /, ' * THE -MSZ- ARRAY.  EITHER DECREASE THE # OF CODEWORD *'
     * /, ' * OR  INCREASE  THE  -MSZ-  ARRAY  SIZE AND RECOMPILE *')      
      goto 99995
      END IF
      write(*,431)mcycle,maxs
431   format(1x,'INPUT:    CYCLES =',I4,5X,'REFINABLE PARAMETERS =',I4)
C     CHECK DIMENSIONING FOR SOME EQUIVALENCED ARRAYS
      IF (IDSZ.LT.MSZ*MAXS) THEN
      WRITE (6,7989) IDSZ, MSZ,MAXS
      WRITE (7,7989) IDSZ, MSZ,MAXS
7989    FORMAT (' CHANGE IDSZ OR MSZ SO THAT IDSZ IS EQUAL TO OR'
     *   , ' GREATER THAN MSZ*MAXS',/,
     *        ' IDSZ, MAXIMUM NO. OF DATA POINTS = ',I6,/,
     *           ' MSZ,  MATRIX SIZE                = ',I6,/,
     *              ' MAXS, NO. OF PARAMETERS VARIED   = ',I6,/,
     *                 ' *** JUST A DIMENSIONING ERROR ***'/)
      STOP 'IDSZ IS LESS THAN MSZ*MAXS'
      END IF
c
      IF(JOBTYP.GT.2)MAXS=0
      WRITE(6,44)MAXS
44    FORMAT(29H NUMBER OF PARAMETERS VARIED= ,I5)
C line 10.1
      READ(5,45,END=99999)glb(1),glb(10),glb(11),glb(8),glb(9),glb(12),
     *                    glb(13) 
C line 10.11
      READ(5,45,END=99999)aglb(1),aglb(10),aglb(11),aglb(8),aglb(9),
     1                    aglb(12),aglb(13)
45    FORMAT(BZ,10F8.0)
      WRITE(6,46)glb(1),aglb(1)
46    FORMAT(32H GLOBAL PARAMETERS AND CODEWORDS ,/,
     * 11H ZEROPOINT= ,2F8.2) 
c-----PMON1,PMON2=PARAMETER OF THE MONOCHROMATOR
c-----IT READS PARAMETER OF THE MONOCHROMATOR AND AIR SCALE
c-----IF THE MONOCHROMATOR WORKS ON THE INCIDENT BEAM PUT :
c-----PMON1=1;PMON2=0;        
C next 2 'READ' are for amorphous bkg codes
C SCAIR=glb(17), FLSCAIR=aglb(17), SCAM =glb(20), FLSCAM=aglb(20)
C PMON1=glb(18), FLMON1 =aglb(18), PMON2=glb(19), FLMON2=aglb(19)
C      READ(5,455,END=99999)SCAIR,FLSCAIR,SCAM,FLSCAM,
C     *  PMON1,FLMON1,PMON2,FLMON2  
c
C      READ(5,455,END=99999)glb(17),aglb(17),glb(20), aglb(20),
C     *  glb(18),aglb(18),glb(19), aglb(19)
c455     FORMAT(BZ,8F8.0)
C  !cp jun 95 start
      if (ibgd.eq.1) goto 4555     
C      if (iax.eq.0) then 
C line 10.2(*) and 10.21(*)
C         READ(5,455,END=99999)glb(17),glb(20),glb(18),glb(19),ctime,
C     *                       aglb(17),aglb(20),aglb(18),aglb(19)
c455     FORMAT(BZ,5F8.0/4f8.0)
C line 10.2 and line 10.21
         READ(5,455,END=99999)glb(20), glb(18), glb(19),
     *                       aglb(20),aglb(18),aglb(19)
455     FORMAT(BZ,3F8.0/3f8.0)
C      end if
C       IF(IAx.EQ.1.AND.(glb(17).NE.0.0.OR.aglb(17).NE.0.0)) goto 88887
C      WRITE(6,466)glb(17),aglb(17),glb(20),aglb(20) 
C      write(*,466)glb(17),aglb(17),glb(20),aglb(20) 
c466   FORMAT(26H AIR BACKGROUND CODEWORDS ,/,
C     * ' AIRSCALE= ',2F8.2,/,
C     * ' AMORPHOUS SCALE= ',2F8.2)
      WRITE(6,466)glb(20),aglb(20) 
C               123456789012345678901234567890
466   FORMAT(31H AMORPHOUS SCALE and CODEWORD= ,2F14.4)
      WRITE(6,467)glb(18),aglb(18),glb(19),aglb(19)
467   FORMAT(48H MONOCROMATOR BANDPASS PARAMETERS AND CODEWORDS,/,
     *' PARAMETERS MONOC=',2F8.4,/,
     *'                  ',2F8.4)
c-----IT READS THE AIR SCATTERING FITTING PARAMETER
C line 10.3(*)
C       READ(5,456,END=99999)AI1,AI2,AI3,AI4,AI5,AI6,AI7,AI8,AI9,AI10
c456    FORMAT(BZ,10F8.0)
C       WRITE(6,457)AI1,AI2,AI3,AI4,AI5,AI6,AI7,AI8,AI9,AI10
c457   FORMAT(26H AIR SCATTERING PARAMETERS ,/,
C     * '   AI1 = ',F8.4,/,
C     * '   AI2 = ',F8.4,/,
C     * '   AI3 = ',F8.4,/,
C     * '   AI4 = ',F8.4,/,
C     * '   AI5 = ',F8.4,/,
C     * '   AI6 = ',F8.4,/,
C     * '   AI7 = ',F8.4,/,
C     * '   AI8 = ',F8.4,/,
C     * '   AI9 = ',F8.4,/,
C     * '   AI10= ',F8.4)
c
4555  IF(NBCKGD.NE.0)GOTO 49
C line 10.4(*)
C line 10.3 and line 10.31
      READ(5,51,END=99999)BACK,FBACK
51    FORMAT(BZ,6F9.2)
      WRITE(6,52)BKPOS,BACK,FBACK
52    FORMAT(36H BACKGROUND PARAMETERS AND CODEWORDS,/
     * ' ORIGIN OF BACKGROUND POLYNOMIAL AT TWO-THETA = ',F8.3,'DEGREES'
     *   ,/,6G12.6,/,6(F8.3,4x))
49    WRITE (6,7100) glb(10),aglb(10),glb(11),aglb(11)
7100  FORMAT (1X,'DISPLACEMENT PEAKSHIFT PARAMETER AND CODEWORD',2F8.2
     * ,/,' TRANSPARENCY PEAKSHIFT PARAMETER AND CODEWORD',2F8.2)
      write (6,7101) glb(8),aglb(8),glb(9),aglb(9),glb(12),aglb(12),
     *               glb(13),aglb(13) 
7101  FORMAT (' SURFACE ROUGHNESS P PARAMETER AND CODEWORD',2F9.4,
     * /,     ' SURFACE ROUGHNESS Q PARAMETER AND CODEWORD',2F9.4,
     * /,     ' SURFACE ROUGHNESS R PARAMETER AND CODEWORD',2F9.4,
     * /,     ' SURFACE ROUGHNESS T PARAMETER AND CODEWORD',2F9.4,/)
      IF(JOBTYP.GT.2)GOTO 483 
c-----IF PATTERN CALCULATION ONLY FOR SYNCHROTRON X-RAY DATA GO TO 501
      IF(JOBTYP.EQ.1.AND.INSTRM.EQ.2) GOTO 501
c-----IF PATTERN CALCULATION ONLY FOR MULTIPLE NEUTRON DATA GO TO 601
      IF(JOBTYP.EQ.2.AND.INSTRM.EQ.2) GOTO 601
C 471      FORMAT(BZ,3F8.2,A56)     !Some step size require the format below !cp 30jun99
471    format(BZ,F8.5,f8.7,f8.5,A56)
C read DBWS formated
      if(idata.eq.0.or.idata.eq.1) then
         READ(4,471,END=99998)THMIN,STEP,THMAX,DATAID
         write(*,4711)thmin,thmax,step
4711  format(1x,'DATA RANGE (2THETA):  START =',F8.3,', STOP =',
     *      F8.3,', STEP =',F8.3,/)
      WRITE (6,7967) DATAID
7967  FORMAT (5X,'DATA ID ',A56)
      write(6,4711)thmin,thmax,step
      NPTS=(THMAX-THMIN)/STEP+1.5
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      WRITE (*,1000) IDSZ,NPTS
1000    FORMAT(1X,'PROGRAM CAN HANDLE ',I5,'POINTS',/,
     *      1X,I5,' POINTS WERE INPUT',/,
     *        ' INCREASE IDSZ IN PARAMETER STATEMENT')
      STOP ' TOO MANY DATA POINTS '
      ENDIF
C read GSAS data file (read start,stop,step and data)
      else if(idata.eq.2) then
        CALL GSASREAD
        GOTO 3333
        END IF
C read Philips data file (read start,stop,step and data)
      IF(IDATA.EQ.3) THEN
        CALL PHILIPSREAD
        GOTO 3333
        END IF
C read SCINTAG data file (read start,stop,step and data)
      IF(IDATA.EQ.4) THEN
        CALL scintag
        GOTO 3333
        END IF
C read  SIEMENS UXD data file (read start,stepsize,stepcount and data)
      IF(IDATA.EQ.5)THEN
        CALL SIEMENSREAD
        GOTO 3333
        END IF
C read  RIGAKU data file (read start,stepsize,stop,stepcount and data)
      IF(IDATA.EQ.6)THEN
        CALL rigakuread
        GOTO 3333
        END IF
C END SPECIAL FORMATS DATA FILE        
      if(idata.eq.0)READ(4,472,END=99998)(Y(I),I=1,NPTS)
472   FORMAT(BZ,8(F7.0,1X))                              
      if(idata.eq.1.)read(4,*,END=99998)(Y(I),I=1,NPTS)
C     CHECK FOR NON-ZERO COUNTS AT ALL POINTS
      xmaxint = 0.
3333  do 444 i=1,npts
      if(y(i).le.1.0E-6)y(i)=1.
      if(y(i).gt.xmaxint)xmaxint=y(i)
444   continue
      DO 4484 I=1,NPTS
C     BUILD UP AMORPHOUS-VECTOR
C     IF REQUIRED MAKE ABSORPTION CORRECTION
c-----COMPUTE TWO-THETA       
         TH = THMIN + FLOAT(I-1) * STEP
         VAR(I)=Y(I)
C      end if
c-----COMPUTE SAMPLE ABSORPTION CORRECTION
      IF(IAS.EQ.1) THEN
        CALL ABSORP(TMV,SW,TH,ABC)
        Y(I) = Y(I) / ABC
        VAR(I)=VAR(I)/ABC
      ENDIF
C  4484    VAR(I) = Y(I)
4484  continue
      GOTO 484
C       READ DATA FROM SYNCHROTRON SOURCE AND CORRECT DATA FOR DEAD TIME,
C      CALCULATE VARIANCE FOR EACH OF THE DATA POINTS
C     DATE IS OCT85,FEB86,AUG86,SRS83,SRS91
C     NRANGE IS THE NO OF BLOCKS IN WHICH THE DATA ARE GIVEN
C     OFSTI0 - DARK BEAM CURRENT, READING WITH NO ELECTRONS IN THE CHAMBER
C     OFSTI1 - DETECTOR DARK BEAM CURRENT
C     CHMBRI - ALL CURRENTS ARE  NORMALISED TO CHMBRI
501   READ (4,796) DATE,DATAID
796   FORMAT(A5,A56)
      WRITE (6,7967) DATAID
      READ (4,502,END=99998) NRANGE,CHMBRI,TAUK
502   FORMAT (BZ,I8, F8.0,G10.0)
C     NPTS IS THE COUNTER FOR TOTAL NO OF POINTS IN ALL THE RANGES
      NPTS = 0
      DO 530 J=1,NRANGE
C     READ INFORMATION ABOUT EACH RANGE 
      READ (4,500,END=99998) ANGMIN,STEP,ANGMAX,STPTIM,OFSTI0,OFSTI1
500     FORMAT (BZ,6F8.0)
c
      IPTS= (ANGMAX-ANGMIN)/STEP +1.5 
C     FIND MAXIMUM AND MINIMUM TWO THETA IN ALL RANGES
      IF (J.EQ.1) THMIN = ANGMIN
      IF (J.EQ.NRANGE) THMAX = ANGMAX 
      IF (J.GT.1) THEN
        READ (4,510)
        IPTS = IPTS - 1
      END IF
      IF (NPTS+IPTS.GT.IDSZ) THEN
        WRITE(6,1000) IDSZ,NPTS+IPTS
        WRITE (7,1000) IDSZ,NPTS+IPTS 
        STOP ' TOO MANY DATA POINTS ' 
      ENDIF
C     VAR(I) IS USED FOR TWO QUANTITIES, JUST TO SAVE SOME SPACE.
      IF (DATE.EQ.'OCT85') THEN
        READ (4,510,END=99998) (Y(I+NPTS),VAR(I+NPTS),I=1,IPTS)
      ELSE IF (DATE.EQ.'FEB86') THEN
        READ (4,511,END=99998) (Y(I+NPTS),VAR(I+NPTS),I=1,IPTS)
      ELSE IF (DATE.EQ.'AUG86') THEN
        READ (4,512,END=99998) (VAR(I+NPTS),Y(I+NPTS),I=1,IPTS)
      ELSE IF (DATE.EQ.'SRS83') THEN
        READ (4,513,END=99998) (VAR(I+NPTS),Y(I+NPTS),I=1,IPTS)
      ELSE IF (DATE.EQ.'SRS91') THEN
        READ (4,514,END=99998) (VAR(I+NPTS),Y(I+NPTS),I=1,IPTS)
      ELSE
        WRITE (6,513)
      END IF
510     FORMAT (BZ,24X,F7.0,1X,F7.0)
C     USE THIS FOR DATA TAKEN IN FEB 86 
511     FORMAT (BZ,24X,F7.0,1X,F7.0)
C     USE THIS FOR DATA TAKEN IN AUG 86 
512     FORMAT (BZ,51X,F10.0,1X,F10.0)
C     USE THIS FOR DATA TAKEN ON THE SRS AT 8.3 IN FALL 1989
513     FORMAT (BZ,37X,F9.0,1X,F9.0)
C     USE THIS FOR DATA TAKEN ON THE SRS AT 9.1 IN FALL 1989
514     FORMAT (BZ,29X,F9.0,10X,F9.0)
515     FORMAT (5X,'WHEN AND WHERE WERE THESE DATA TAKEN. OCT85,FEB86,
     *  AUG86,SRS83,SRS91')
C     TAUK = 2.4621E-6
      DO 520 I=NPTS+1,IPTS+NPTS
        Y(I) = Y(I)/STPTIM
        FLEET1 = VAR(I)-OFSTI0
        FLEET2= FLEET1/CHMBRI
        FLEET2= (FLEET2*(1-TAUK*Y(I)))**2.
        VAR(I) = (STPTIM*Y(I))/FLEET2 
        Y(I) = stptim*(Y(I)/(1-TAUK*Y(I))-OFSTI1)*CHMBRI/FLEET1
        IF (ABS(Y(I)).LT.0.01) Y(I)=1.
520       CONTINUE
530     NPTS=NPTS+IPTS
      write(*,5301)thmin,thmax,step
      write(6,5301)thmin,thmax,step
5301  format(1x,'DATA RANGE (2THETA):  START =',F8.2,', STOP =',
     *      F8.2,', STEP =',F8.2,/)
      GOTO 484
C       ENDS SYNCHROTRON DATA MODIFICATION
C     BEGIN VARIANCE CALCULATION FOR VARYING NO. OF COUNTERS AT EACH STEP
601   READ (4,471,END=99998) THMIN,STEP,THMAX,DATAID
      WRITE (6,7967) DATAID
      write(*,6011)thmin,thmax,step
      write(6,6011)thmin,thmax,step
6011  format(1x,'DATA RANGE (2THETA):  START =',F8.2,', STOP =',
     *      F8.2,', STEP =',F8.2,/)
      NPTS = (THMAX-THMIN)/STEP+1.5
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      STOP ' TOO MANY DATA POINTS '
      ENDIF
      READ (4,602,END=99998) (VAR(I),Y(I),I=1,NPTS)
602   FORMAT (10(F2.0,F6.0))
      DO 603 I=1,NPTS
603     VAR(I)=Y(I)/VAR(I)
      GO TO 484
c
483   NPTS=(THMAX-THMIN)/STEP+1.5
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      STOP ' TOO MANY DATA POINTS '
      ENDIF
484   CONTINUE
      IF(NBCKGD-1)474,475,476 
474   TH=THMIN-STEP 
      DO 473 I=1,NPTS
      TH=TH+STEP
      THX=TH/BKPOS-1.0
      BK(I)=BACK(1)
      DO 473 J=2,6
473       BK(I)=BK(I)+BACK(J)*THX**(J-1)
      GOTO 477
475   open (3,file=' ',status='unknown')
      READ(3,472,END=99997)(BK(I),I=1,NPTS)
      GOTO 477
476   DIFB=POS(1)-THMIN
      IF(DIFB)4760,4761,4761
4761  NBX=DIFB/STEP+2.5
      DO 4762 I=1,NBX
4762    BK(I)=BCK(1)
      NXX=1
      GOTO 4764
4760  NBX=(POS(2)-THMIN)/STEP+1.5
      BK(1)=BCK(1)-DIFB/(POS(2)-POS(1))*(BCK(2)-BCK(1))
      BSTEP=STEP/(POS(2)-POS(1))*(BCK(2)-BCK(1))
      DO 4763 I=2,NBX
4763    BK(I)=BK(I-1)+BSTEP
      NXX=2
4764  NBC=NBCKGD-1
      IF(POS(NBCKGD).GE.THMAX)GOTO 4767 
      POS(NBCKGD+1)=THMAX
      BCK(NBCKGD+1)=BCK(NBCKGD)
      NBC=NBC+1
4767  DO 4765 J=NXX,NBC
      BSTEP=STEP*(BCK(J+1)-BCK(J))/(POS(J+1)-POS(J))
      NINC=(POS(J+1)-POS(J))/STEP+1.5 
      N2X=MIN0(NPTS,NBX+NINC)
      DO 4766 I=NBX,N2X
4766      BK(I)=BK(I-1)+BSTEP 
      NBX=N2X
      IF(NBX.EQ.NPTS)GOTO 477
4765    CONTINUE
477   CONTINUE
      DO 4100 K=1,NPHASE
4100    ICR(K)=0
C start loop on phase
      DO 81 K=1,NPHASE
C line 11.1
      READ(5,465,END=99999)PHSNM(K)
465     FORMAT(BZ,A50)
      WRITE(6,47)K,PHSNM(K) 
47      FORMAT(7H1PHASE  ,I2,/,1X,A50)
C line 11.2
      READ(5,48,END=99999)NATOM(K),NMOL(K),saqf(k),(PREF(K,I),I=1,3),
     *wtis(k)
      N=NATOM(K)
         IF (saqf(k).LE.0.)THEN
               saqf(k)=1.0
               write(*,478)k
         ENDIF
      ISTEST = ISPHASE
         IF (K.EQ.isphase) THEN
                IF (wtis(K).EQ.0.) then
                   isphase = 0
                   write(*,578) istest
                   write(6,5781)
                END IF
         END IF
478   format(' Particle absorption factor (phase ',i2, ') is now 1.')
578   format(' WARNING: Internal Standard WT% = ZERO. Check ISWT in'
     *' line 11.2 for phase ',i2,'.',/
     *' ISPHASE (line 2) turned to ZERO. Amourphous content will not',
     *' be calculated.')
5781  format(15x,' WARNING: Wt% for internal Standard is ZERO.',/,15x,
     *           ' Check ISWT in line 11.2 for this phase.',/,15x,
     *           ' ISPHASE (line 2) turned to ZERO.',/,15x,
     *           ' AMORPHOUS CONTENT WILL NOT BE CALCULATED.')
      if (k.eq.isphase)then
            WRITE(6,4801)N,NMOL(K),saqf(k),(PREF(K,I),I=1,3),wtis(k)
4801     FORMAT(18H NUMBER OF ATOMS= ,I4,/,
     *         40H NUMBER OF FORMULA UNITS PER UNIT CELL= ,I4,/, 
     *         30H PARTICLE ABSORPTION FACTOR = , F8.4,/,
     *         31H PREFERRED ORIENTATION VECTOR=  ,3F8.4,/
     *           ' MASS% IN THE SAMPLE= ',F7.2)
      else
        WRITE(6,480)N,NMOL(K),saqf(k),(PREF(K,I),I=1,3)
480       FORMAT(18H NUMBER OF ATOMS= ,I4,/,
     *         40H NUMBER OF FORMULA UNITS PER UNIT CELL= ,I4,/, 
     *         30H PARTICLE ABSORPTION FACTOR = , F8.4,/,
     *         31H PREFERRED ORIENTATION VECTOR=  ,3F8.4)  
      end if
48      FORMAT(BZ,2I4,F7.0,1x,3F4.0,f7.2)
c line 11.3
      READ(5,460,END=99999)(SYMB(K,I),I=1,20)
460     FORMAT(BZ,20A1)
      DO 135 I=1,20   
          SPG(i)=SYMB(K,I)
C convert lower case to upper case
          DO 22371 ik=1,26
          IF (spg(i).EQ.LOWER(ik:ik)) spg(i)=UPPER(ik:ik)
22371     CONTINUE
c finish conversion
135       continue
      CALL SPGP(SPG)
C getting multiplicity of each phase !cp jun 96)
      isimop=1
       call rtmt(ipl1,k)
       xmltp(k)=MLTPHASE
       WRITE(6,3393)MLTPHASE
3393  format(' The multiplicity of the general site is ', i3)
c-----READ AND PRINT FOR EACH ATOM
      IOF=0
      IF(K.GT.1) THEN
        DO 4325 IIPHAS=2,K
4325        IOF = IOF + NATOM(IIPHAS-1) 
      END IF
      WRITE(6,57) 
C line 11-4i
C READ and WRITE, and respectives FORMAT command lines below were
C changed to incorporate the parameter MURT(I+IOF) * !cp jun 96
      READ(5,65,END=99999)(ATEXT(I+IOF),murt(I+IOF),NTYP(I+IOF),
     *      (XL(I+IOF,J),J=1, 5), (A(I+IOF,J),J=1, 5),
     *      (XL(I+IOF,J),J=6,11), (A(I+IOF,J),J=6,11),I=1,N)
65      FORMAT(BZ,A4,1X,I4,1X,A4,2X,5F8.5,/,16X,5F8.2,/,6F8.5,/,6F8.2) 
C convert lower case to upper case
        DO ITN = 1,N
          DO ITIPO = 1,2
            DO 22378 ik=1,26
             IF (NTYP(ITN+IOF)(ITIPO:ITIPO).EQ.LOWER(ik:ik)) 
     1           NTYP(ITN+IOF)(ITIPO:ITIPO)=UPPER(ik:ik)
22378       CONTINUE
          END DO
        END DO 
c finish conversion
      WRITE(6,170)(ATEXT(I+IOF),murt(I+IOF),NTYP(I+IOF),(XL(I+IOF,J),
     *   J=1,11),I=1,N)
170     FORMAT(1X,A4,1X,i4,1x,A4,7x,5F10.5/22X,6F10.5) 
C  !cp jun 96 ... (CONVERT sof MULTIPLICITY)(also changed in OUTPTR)
      do 3331 i=1,n
        if(int(a(i+iof,5)/10).ne.0.and.xl(i+iof,5).eq.0)then
             xl(i+iof,5)=1e-6
        end if
            xl(i+iof,5) = xl(i+iof,5) * murt(i+iof) / xmltp(k)
3331    continue
C par(k,21) introduced below. It is for the term cot**2 in the pv-5 FWHM !cp Aug 95
C line 11-5, line 11-6, line 11-7, line 11-8 and line 11-9
      READ(5,50,END=99999) PAR(K,1),PAR(K,2),                   !S  O_B (line 11-5)
     *                          APAR(K,1),APAR(K,2),
     * PAR(K,3),PAR(K,4),PAR(K,5),PAR(K,21),par(k,20),          !FWHM (line 11-6)
     * PAR(K,15),PAR(K,16),
     *     APAR(K,3),APAR(K,4),APAR(K,5),APAR(K,21),apar(k,20),
     *     APAR(K,15),APAR(K,16),
     *                  (PAR(K,I),I=6,11), (APAR(K,I),I=6,11),  !Unit cell (line 11-7)
     *                  PAR(K,12),PAR(K,13),PAR(K,14),          !G1 G2 P (line 11-8)
     *                     APAR(K,12),APAR(K,13),APAR(K,14),
     *                  PAR(K,17),PAR(K,18),PAR(K,19),          !NA NB NC (line 11-91)
     *                     APAR(K,17),APAR(K,18),APAR(K,19),
     *                  PAR(K,24),PAR(K,25),PAR(K,26),          !NA NB NC HS (line 11-93)
     *                     APAR(K,24),APAR(K,25),APAR(K,26),
     *                  par(k,27),apar(k,27)                    !s-PVII (line 11-95)
50      FORMAT(BZ,G8.2,F8.0,/,2F8.2,/7F8.0,/,7F8.2,/,
     *         6F8.0,/,6F8.2,/,3F8.0,/,3F8.2,/,3F8.4,/,3F8.2,/,
     *         3F8.4,/,3f8.2,/,f8.4,/,f8.2)    
C checking for zeros if TCH-PV is being used
      if(nprof.eq.8)then
         if(par(k,15).eq.0)par(k,15)=1e-8
         if(par(k,16).eq.0)par(k,16)=1e-8
         if(par(k,3).eq.par(k,4).and.par(k,4).eq.par(k,5)
     *.and.par(k,5).eq.par(k,20).and.par(k,3).eq.0) par(k,5)=1e-8
      end if
C checking for zeros in PV #5
      if(nprof.eq.6)then
        if(par(k,17).eq.0)par(k,17)=1e-6
      end if
C      if(int(apar(k,20)/10).ne.0.and.par(k,20).eq.0)par(k,20)=1e-9
C !cp Aug 95 introducing par(k,21)=ct
C CHECKING FOR NON-REFINABLE PARAMETERS
cCCC                        FOR  CT  IN TCHZ AND SPVII FUNCTIONS
      if(nprof.eq.8.or.nprof.eq.5)then       
            par(k,21)=0.0
          if(apar(k,21).ne.0.)then
              write(*,9477)
              write(6,9477)
9477  format(1x,'NON-REFINABLE PARAMETER TURNED ON (CT) WITH' 
     *          ' SPLIT PEARSON VII OR TCHZ PROFILE FUNCTION')
              stop
           end if
      end if               
cCCC                             FOR  X,Y,Z IN NON-TCHZ FUNCTION
      IF(NPROF.NE.8) THEN
         if(par(k,20).ne.0.or.par(k,16).ne.0.or.apar(k,15).ne.0)then
            write(6,94781)
            par(k,20)=0.0
            par(k,16)=0.0
            par(k,15)=0.0
94781 format(6x,'NON-REFINABLE PARAMETER RESET TO ZERO (Z,X,Y) FOR'/ 
     *      ,6X,'NON TCHZ PROFILE FUNCTION')  
         end if
         if(apar(k,20).ne.0.or.apar(k,16).ne.0.or.apar(k,15).ne.0)then
            write(*,94782)
            write(6,94782)
94782 format(1x,'NON-REFINABLE PARAMETER TURNED ON (Z,X,Y) WITH', 
     *' NOT TCHZ PROFILE FUNCTION')  
              stop     
            end if         
      END IF
cCCCC                            FOR  RIET_ASYM,X,Y,Z,CT IN SPVII
      IF (NPROF.EQ.5) THEN
        do 9772, kkS=1,3
          if(par(k,13+kkS).ne.0.0) then
            par(k,13+kkS)=0.0
            write(6,9773)
9773       format(1X,' RIET_ASYM X Y RESET TO ZERO FOR SPVII FUNCTION')
          end if
          if(apar(k,13+kkS).ne.0.0) then
            write(*,9774)
            write(6,9774)
            stop
          end if
9774  format(1x,'NON-REFINEABLE PARAMETER TURNED ON (X,Y,R/RCF_ASYM)',
     *  /,1X,', WITH SPLIT PEARSON VII PROFILE FUNCTION')
9772    continue
c
        do 9775, kkS=1,2
          if ( par(k,19+kkS).ne.0.0 ) then
            par(k,19+kkS)=0.0
            write(6,9776)
9776       format(1X,' CT,Z  RESET TO ZERO FOR SPVII FUNCTION')
          end if
          if ( apar(k,19+kkS).ne.0.0 ) then
            write(*,9777)
            write(6,9777)
            stop
          end if
9777  format(1x,'NON-REFINEABLE PARAMETER TURNED ON (CT,Z) ',
     *   /,1X,', WITH SPLIT PEARSON VII PROFILE FUNCTION')
9775    continue
      END IF
cCCCCC             FOR HIGH SIDE PARAMETERS IN NON-PVII PROFILE FUNCTION
      if (nprof.ne.5) then
        do 772, kk=1,4
         if ( par(k,23+kk).ne.0.0 ) then
          par(k,23+kk)=0.0
          write(6,9459)
9459     format(' NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY RESET TO ',
     *     'ZERO FOR A,'/',            NON-SPLIT PEARSON VII FUNCTION')
        end if
        if(apar(k,23+kk).ne.0.0) then
          write(*,9460)
          write(6,9460)
          stop
        end if
9460     format (' NA NB NC (HIGH SIDE) AND PEARSON ASYMMETRY ONLY ',
     *    'REFINABLE FOR A,'/',            SPLIT PEARSON VII FUNCTION')
772        continue
      end if
cCCCCC END OF CHECKING NON-REFINABLE PARAMETERS
      WRITE(6,66) PAR(K,1),PAR(K,2),
     *             (PAR(K,I),I=6,11),
     *              PAR(K,12),PAR(K,13),PAR(K,14)
66      FORMAT(' OVERALL SCALE FACTOR=',G12.6,/,
     *           ' OVERALL TEMP. FACTOR=',F12.5,/,
     *            ' DIRECT CELL PARAMETERS=',6F9.4,/,
     *             ' PREFERRED ORIENTATION PARAMETERS=',2F7.3,/,
     *              ' ASYMMETRY PARAMETER=',F8.4)
c-----CHECK FOR SPLIT PEARSON PROFILE
       if (nprof.eq.5) then
       write(6,774) par(k,17),par(k,18),par(k,19),
     1                    par(k,24),par(k,25),par(k,26),
     1                        par(k,27)
774       format('  LOW SIDE EXPONENT COEFFICIENTS=',3F12.4,/,
     1            '  HIGH SIDE EXPONENT COEFFICIENTS=',3F12.4,/,
     1             '  SPLIT PEARSON VII ASYMMETRY PARAMETER=',3F8.4,/)
c-----IF NOT THE SPLIT PEARSON VII PROFILE
       else
       write(6,779) PAR(K,17),PAR(K,18),PAR(K,19)
779      format(' MIXING PARAMETERS = ',3(G10.3,1x)/)
       end if
       write(6,773) PAR(K,3),PAR(K,4),PAR(K,5),
     *                  PAR(K,21),par(k,20),PAR(K,15),PAR(K,16)
773      format(' FWHM PARAMETERS (U,V,W,CT,Z,X,Y)=',7F9.4,/)
c
      AA=PAR(K,6) 
      B=PAR(K,7)
      C=PAR(K,8)
      ALPHA=PAR(K,9)
      BETA=PAR(K,10)
      GAMMA=PAR(K,11)
C     CALL CELL2
        CALL CELL2(K,LAMDAM)
C ************************************** !cp ap 97 (from It code)
      IF(FONDO.EQ.1.OR.FONDO.EQ.2) THEN
      SCABKG(K) = GCOM(K) * PAR(K,1)
      ENDIF
C ************************************************************
C        WRITE(6,68) K,VOLi(K),K,GCOM(K)
        WRITE(6,68) K,VOLi(K)
c68      FORMAT(' CELL VOLUME PHASE(',I2,' ) = ',F12.4,/,
C     &         ' GCOM(',I2,' ) = ',F12.4)
68      FORMAT(' CELL VOLUME PHASE(',I2,' ) = ',F12.4)
      DO 72 I=1,6 
72        SAVE(K,I)=PAR(K,I+5)
      DO 73 I=1,3 
73        PAR(K,I+5)=AL(I,I)
      PAR(K,9)=AL(2,3)
      PAR(K,10)=AL(1,3)
      PAR(K,11)=AL(1,2)
      WRITE(6,82) 
c-----PRINT CODEWORDS FOR ATOMIC PARAMETERS
      WRITE(6,85)(ATEXT(I+IOF),(A(I+IOF,J),J=1,11),I=1,N) 
85      FORMAT(1X,A4,17X,5F10.2/22X,6F10.2)
      DO 83 I=1,N 
        DO 84 J=1,11
          X=A(I+IOF,J)
          IYY=INT(ABS(X)/10.)
          IF(IYY.GT.MSZ) GO TO 99996
          LP(I+IOF,J)=IYY
84          A(I+IOF,J)=(ABS(X)-10.*FLOAT(IYY))*SIGN(1.,X)
83        CONTINUE
c-----PRINT CODEWORDS FOR PROFILE PARAMETERS
      WRITE(6,89) APAR(K,1),APAR(K,2),
     *             (APAR(K,I),I=6,11),
     *              APAR(K,12),APAR(K,13),APAR(K,14)
89      FORMAT(' OVERALL SCALE FACTOR=',f8.2,/,
     *           ' OVERALL TEMP. FACTOR=',F8.2,/,
     *            ' DIRECT CELL PARAMETERS=',6F8.2,/,
     *             ' PREFERRED ORIENTATION PARAMETERS=',2F8.2,/,
     *              ' ASYMMETRY PARAMETER=',F8.4)
C !cp ap 97 (from It code)
      IF(FONDO.EQ.1.AND.(PAR(K,2).NE.0.0.OR.APAR(K,2).NE.0.0)) THEN
      GOTO   88888
      END IF
      IF(FONDO.EQ.2.AND.PAR(K,2).EQ.0.0.AND.APAR(K,2).EQ.0.0) THEN
      GOTO   88889
      END IF
c-----CHECK FOR SPLIT PEARSON PROFILE
      if (nprof.eq.5) then
        write(6,775) apar(k,17),apar(k,18),apar(k,19),
     1                  apar(k,24),apar(k,25),apar(k,26),
     1                    apar(k,27)
775     format(' LOW SIDE EXPONENT COEFFICIENTS=',3f8.2,/,
     1            ' HIGH SIDE EXPONENT COEFFICIENTS=',3f8.2,/,
     1             ' SPLIT PEARSON VII ASSYMETRY PARAMETER=',3F8.2)
c-----IF NOT THE SPLIT PEARSON VII PROFILE
      else
        write(6,781) APAR(K,17),APAR(K,18),APAR(K,19)
781     format(' MIXING PARAMETERS = ',3(g10.3,1x),/)
      end if
      write(6,776) APAR(K,3),APAR(K,4),APAR(K,5),
     *             APAR(K,21),apar(k,20),APAR(K,15),APAR(K,16)
776   format(' FWHM PARAMETERS (U,V,W,CT,Z,X,Y)=',7F8.2,/)
      DO 87 I=1,27
        X=APAR(K,I)
        IYY=INT(ABS(X)/10.) 
        LPAR(K,I)=IYY
87      APAR(K,I)=(ABS(X)-10.*FLOAT(IYY))*SIGN(1.,X)
151   CALL LOOKUP(K,NATOM(K),NSCAT,IXRAY,JOBTYP)
      IF(FONDO.EQ.1.OR.FONDO.EQ.2) CALL FINDC(K,NSCAT)
      U=PAR(K,3)
      V=PAR(K,4)
      W=PAR(K,5)
      ZZZ = PAR(K,20)
      uc = par(k,21)
      IF (NPROF.EQ.8) THEN
        ULOR=PAR(K,15)
        VLOR=PAR(K,16)
      END IF
      CALL REFGEN(K,glb(1),glb(10),glb(11),PREF,PAR(K,12))
      isimop=0
      CALL RTMT(IPL1,K)
C         mltp(k)=MLTPHASE
      ICY=1
      ICZ=0
      DO 4101 IIPHAS=1,K
4101      ICZ = ICZ + ICR(IIPHAS)
      IF(K.GE.2) ICY=1+ICZ-ICR(K)
      IF(LST1.NE.1)GOTO 479 
      IXDEL=0
      DO 481 IXX=ICY,ICZ
        IX=IXX-IXDEL
        IRL=MOD(IREFS(IXX),256)-128
        IRK=MOD(IREFS(IXX)/256,256)-128
        IRH=MOD(IREFS(IXX)/(256*256),256)-128
        IRC=MOD(IREFS(IXX)/(256*256*256),8)
        MLTT = NINT(FMGNTD(IXX))
        IF (NPROF.EQ.8) THEN
          IF(MOD(IX-1,60).EQ.0)WRITE(6,7482)
          WRITE(6,7470)IX,IRC,IRH,IRK,IRL,MLTT,
     1       (REFS(IXX,JX),JX=1,3),HALFL(IXX),HALFG(IXX),GAM(IXX)
7482      FORMAT(38H1NO.  CODE    H   K   L  MULT   HW     ,
     *            16HPOSN      FACTOR,'       HWL       HWG     ETA '/)
7470      FORMAT(1X,2I4,3X,3I4,I6,2F8.3,F10.6,3F8.3)
        else if (nprof.eq.5) then
          fwhm(ixx,1)=2.0*(REFS(IXX,1))*par(k,27)/(1.0+par(k,27))
          fwhm(ixx,2)=2.0*(REFS(IXX,1))/(1.0+par(k,27))
          IF(MOD(IX-1,60).EQ.0)WRITE(6,777)
          WRITE(6,778)IX,IRC,IRH,IRK,IRL,MLTT,
     1            fwhm(ixx,1),fwhm(ixx,2),(REFS(IXX,JX),JX=1,3)
777       FORMAT(49H1NO.  CODE    H   K   L  MULT      HWL    HWH    ,
     *            22H FWHM   POSN    FACTOR/)
778       FORMAT(1X,2I4,3X,3I4,I6,4F8.3,F10.6)
        else
          IF ( MOD(IX-1,60).EQ.0 ) WRITE(6,482)
          WRITE(6,470)IX,IRC,IRH,IRK,IRL,MLTT,
     1                            (REFS(IXX,JX),JX=1,3)
482       FORMAT(40H1NO.  CODE    H   K   L  MULT     HW     ,
     *         16HPOSN    FACTOR  /)
470       FORMAT(1X,2I4,3X,3I4,I6,2F8.3,F10.6)
        end if
481     CONTINUE
479     DO 499 IX=ICY,ICZ
499       REFS(IX,3)=REFS(IX,3)*FMGNTD(IX) ! FLOAT(MLTT(IX))
81    CONTINUE       
C end of great loop on phases
C FIND CODEWORDS FOR GLOBAL PARAMETERS: LGLB(I) & AGLB(I)
      DO 88 I=1,20
      X=AGLB(I)
      IYY=INT(ABS(X)/10.)
      LGLB(I)=IYY 
88      AGLB(I)=(ABS(X)-10.*FLOAT(IYY))*SIGN(1.,X)
57    FORMAT(/25H ***INITIAL PARAMETERS***/' ATOM    M NTYP        ',6X
     * ,1HX,9X,1HY,9X,1HZ,9X,1HB,8X,2HSo,/,27X,3HB11
     *  ,7X,3HB22,7X,3HB33,7X,3HB12,7X,3HB13,7X,3HB23)
67    FORMAT(BZ,4F8.2)
76    FORMAT(' Cell constants=',6F12.8)
77    FORMAT(BZ,2F8.4)
82    FORMAT('0***Coding of variables***'/' ATOM',24X,'X',9X,'Y',9X,'Z'
     * ,9X,'B',9X,'So',/27X,'B11',7X,'B22',7X,
     *   'B33',7X,'B12',7X,'B13',7X,'B23')
86    FORMAT(BZ,8F8.2)
        call qpainit
      IF (JOBTYP.GT.2) RETURN 
      NATOMS = 0
      DO 7150 K=1,NPHASE
7150    NATOMS = NATOMS + NATOM(K)
C     COUNT NO. OF USES OF SAME LOCATION IN NORMAL MATRIX
      DO 17200 I=1,MAXS
        LCOUNT = 0
        DO 17201 J =1,20
          IF ( I.EQ.LGLB(J) ) LCOUNT = LCOUNT + 1
17201     CONTINUE 
        DO 17202 K=1,NPHASE
          DO 17202 J=1,27
            IF ( I.EQ.LPAR(K,J) )  LCOUNT = LCOUNT + 1
17202       CONTINUE
        DO 17203 K=1,NATOMS
          DO 17203 J=1,11
            IF(I.EQ.LP(K,J)) LCOUNT = LCOUNT + 1
17203       CONTINUE
        IF (LCOUNT.GE.2)   WRITE (6,17204) I, LCOUNT
17204   FORMAT ('0******* CODEWORD ',I3,' is used ',I3,' times **')
17200   CONTINUE
C     END COUNT NO. OF USES OF SAME LOCATION IN NORMAL MATRIX
C     CHECK FOR HOLES IN THE NORMAL MATRIX
      ISTOP=0
      DO 7200 I=1,MAXS
c
      DO 7201 J =1,20
        IF (I.EQ.LGLB(J)) GO TO 7200
7201      CONTINUE
C  HOLES IN PHASE PARAMETER
      DO 7202 K=1,NPHASE
        DO 7202 J=1,27
          IF(I.EQ.LPAR(K,J)) GO TO 7200
7202        CONTINUE
C HOLES IN ATOMS PARAMETERS
      DO 7203 K=1,NATOMS
        DO 7203 J=1,11
          IF(I.EQ.LP(K,J)) GO TO 7200 
7203        CONTINUE
      WRITE (7,7204) I
      WRITE (6,7204) I
      WRITE (*,7204) I
7204    FORMAT (5X,/,'***** HOLE IN THE MATRIX ******. ELEMENT ',I3,1X, 
     *     'IN THE NORMAL MATRIX IS MISSING')
      ISTOP = 1
7200    CONTINUE
      IF (ISTOP.EQ.1) STOP
C      IF (ISTOP.EQ.1) STOP 'HOLE IN THE MATRIX'
C     END CHECKING FOR HOLE IN THE NORMAL MATRIX
C     CHECK FOR NON-ZERO COUNTS AT ALL POINTS
C                                             not necessary anymore  !nov98
c
C      DO 5678 I=1,NPTS
C      IF (VAR(I).LE.1.0E-6) THEN
C        ANGLEP = THMIN + FLOAT(I-1) * STEP
C        WRITE (6,5679) I,ANGLEP
C        WRITE (7,5679) I,ANGLEP
c5679   FORMAT (5X,'ZERO COUNTS AT STEP NO.',I7,'  AT TWO THETA ',F8.3)
C        STOP '*** ZERO COUNTS IN DATA FILE AT A STEP ***' 
C      ENDIF
c5678    CONTINUE
      RETURN
99999 STOP 'END OF FILE TAPE5'
99998 STOP 'END OF FILE TAPE4'
99997 STOP 'END OF FILE TAPE3'
99996 STOP 'MATRIX SIZE IS TOO SMALL'
99995 stop 'MAXS > MSZ'
88888 STOP 'WITH FONDO=1 YOU MUST USE ISOTROPIC THERMAL FACTORS'
88889 STOP 'IF FONDO=2 YOU MUST USE THE OVERALL THERMAL FACTOR'
88088 close (4)
      close (5)
      close (6)
      close (7)
      close (8,status='DELETE') 
      stop
C      stop 'CHECK MULTIPLICITIES IN NEWINP.INP * CONVERSION DONE'
c88887 STOP 'IF IAX =1 YOU MUST SET AIR AND FLAIR EQUAL ZERO'
      END 
      SUBROUTINE CALCUL(NN)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      REAL LAMDA
      INTEGER PTR,fondo
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT 
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON /HKLCTL/ IHKL(3,48),AZ(48),NC1(10,7,99),ICHKL(99),
     1    N1HKL(99),IER
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),xmas(16)
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ),
     * BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON /RTMTX/ IVEC(99,192),MLTPHS(99),ICNTPHS(99)
      COMMON/G4/TANN(NOV),DERSTO(NOV,MSZ)
      COMMON/PRFX/DELT,TL,GAM1,GAM2,PRFDER,iph,delta
C      common/prfx/iph,delta
      common/spvii/rl,da1l,da2l,da3l,da4l,da5l,da6l,da7l
      common/spvii/rh,da1h,da2h,da3h,da4h,da5h,da6h,da7h
      DIMENSION HNN(3),XI(14) 
      DIMENSION DERIV(MSZ),SNX(NATS),SINTL(NOV),SA(NATS),SB(NATS)
      DIMENSION TEMP(NATS),SNEX(NATS),T(3),SUMAX(NATS,9),SUMBX(NATS,9)
      DIMENSION H(3),AL(3,3),SM(3,3)       
      COMMON/STRUPHASE/aphase(IRS),tavix(irs),srix(irs)
      LOGICAL VERT,PAC
      PI = 90.0/ATAN(1.0)              ! 360./3.14159265359          
      IPH = IREFS(NN)/(256*256*256*8)  ! 256/256/256/8
      IOF = 0
      IF ( IPH.GT.1 ) THEN
      DO 4340 IIPHAS=2,IPH
4340  IOF = IOF + NATOM(IIPHAS-1)
      ENDIF
      IRL = MLTPHS(IPH)
      N = NATOM(IPH)
      ICENT = ICNTPHS(IPH)
c-----ZEROIZE THE DERIVATIVES OF THIS REFLECTION W.R.T. TO PARAMETERS 
      NX=0
      DO 4310 IIPHAS=1,NPHASE 
4310    NX = NX+NATOM(IIPHAS)
      DO 12 I=1,MSZ
12      DERIV(I)=0.0 
      DO 21 I=1,NX
        DO 21 J=1,9
          SUMBX(I,J) = 0.
21        SUMAX(I,J) = 0.
      CV=0.0
      DV=0.0
      AV=0.0
      BV=0.0
      PAC = PAR(IPH,12).NE.0.0 .OR. LPAR(IPH,12).NE.0
     * .OR. PAR(IPH,13).NE.0.0 .OR. LPAR(IPH,13).NE.0
      DO 303 I=1,3
303     AL(I,I)=PAR(IPH,I+5)
      AL(3,2)=PAR(IPH,9)
      AL(2,3)=AL(3,2)
      AL(3,1)=PAR(IPH,10)
      AL(1,3)=AL(3,1)
      AL(2,1)=PAR(IPH,11)
      AL(1,2)=AL(2,1)
      AH=PAR(IPH,3)
      BH=PAR(IPH,4)
      CH=PAR(IPH,5)
      DH=PAR(IPH,20)
      if (nprof.eq.6) then
        eh=par(iph,21)
      end if
      IF (NPROF.EQ.8) THEN
        AH2=PAR(IPH,15)
        BH2=PAR(IPH,16)
      END IF
      B1=PAR(IPH,12)
      B2=PAR(IPH,13)
      NM=MOD(NN,NOV)+1
      ICX=MOD(IREFS(NN)/(256*256*256),8)
      SLABDA=LAMDA(ICX)*LAMDA(ICX)/4.
      N=NATOM(IPH)
      FMGNTD(NN)=0.0
      PAKNN=0.0
      PRECOR=1.0
      TR=0.0
      HNN(3)=MOD(IREFS(NN),256)-128
      HNN(2)=MOD(IREFS(NN)/256,256)-128 
      HNN(1)=MOD(IREFS(NN)/(256*256),256)-128
      DO 1 I=1,3
1       H(I)=HNN(I) 
c-----CALCULATION OF TEMP.FACTOR,POSITION AND FWHM
      SS=0.
      DO 2 I=1,3
        IHKL(I,1) = HNN(I)
        DO 2 J=I,3
2         SS = HNN(I)*AL(I,J)*HNN(J)+SS                         ! GSAS DH2
      IF ( PAC ) THEN
        TT = 0.
        DO 16 I=1,3
         DO 16 J=I,3 
16          TT = TT+PREF(IPH,I)*AL(I,J)*PREF(IPH,J)             ! GSAS DP2
        CALL SMTRY2(IPH) 
        PRECOR = 0.0
        PREXP = 0.0
        DPRECOR = 0.0
        DO IJ=1,ICHKL(IPH)
          PAK = 0.
          DO 33 I=1,3
            DO 33 J=I,3 
33            PAK = PREF(IPH,I)*AL(I,J)*FLOAT(IHKL(J,IJ))+PAK     ! GSAS CA
          PAK = PAK*PAK/(TT*SS)
          PAKNN = (PI/2.0)**2
          IF (PAK.NE.0) PAKNN=ATAN(SQRT(ABS((1.0-PAK)/PAK)))**2
          IF (IPREF.EQ.0) THEN
            PREXPX = EXP(B1*PAKNN)
            PRECORX = B2+(1.0-B2)*PREXPX
            DPRECORX = PAKNN*(1.0-B2)*PREXPX/PRECORX
          ELSE
            PREXPX = B1*B1*PAK+(1.0-PAK)/B1
            PRECORX = 1.0/PREXPX**1.5
            DPRECORX = -1.5*(2.*B1*PAK-(1.-PAK)/B1/B1)/
     1        (PREXPX**2.5*PRECORX)
          END IF
          PREXP = PREXP+PREXPX
          PRECOR = PRECOR+PRECORX
          DPRECOR = DPRECOR+DPRECORX
C         print '(5x,3f10.5,i3)',prexpx,precorx,dprecorx,ij     ! ******
        END DO
!       write(6,'(a,4i3,a,i2)') 'hkl',(ihkl(i,1),i=1,3),nprof,' iph',iph  
!ACL debug 980520
!       write (6,*) 'tt=',tt,' SS=',SS,' Prexp=',prexp             !ACL debug 980520
!0000000011111111112222222222333333333344444444445555555555666666666677777777778
!2345678901234567890123456789012345678901234567890123456789012345678901234567890
        PREXP = PREXP/FLOAT(ICHKL(IPH))
        PRECOR = PRECOR/FLOAT(ICHKL(IPH))
        DPRECOR = DPRECOR/FLOAT(ICHKL(IPH))
!       write (6,*) 'Precor=',precor,' DPrecor=',dprecor        !ACL debug 980520
      END IF
13    SINTL(NM) = SQRT(SS)
      SSNN = 0.25*SS
      TAV = EXP(-2.0*PAR(IPH,2)*SSNN)*PRECOR
      SINTH = SLABDA*SS
      COSTH = 1.0-SINTH
      TANTH = SQRT(SINTH/COSTH) 
      TANN(NM) = TANTH
C   Correction of microabsorption
      if (iabsr.eq.1) then
         isith=(1.0)/(sqrt(sinth))
       sr = glb(12)*(1.0-glb(8)*exp(-glb(9))+glb(8)*exp(-glb(9)/
     *     sqrt(sinth)))+(1.0-glb(12))*(1+glb(13)*(asin(sqrt(sinth)))-
     *     1.5707963268)
      else if (iabsr.eq.2) then
         sith = sqrt(sinth)
      sr = 1.0-glb(13)*(asin(sith)-1.5707963268)
      else if (iabsr.eq.3) then
         isith=(1.0)/(sqrt(sinth))
      sr = 1.0-glb(8)*exp(-glb(9))+glb(8)*exp(-glb(9)*isith)
      else if (iabsr.eq.4) then
         isith=(1.0)/(sqrt(sinth))
      sr = 1.0-glb(8)*glb(9)*(1.0-glb(9))-isith*glb(8)*glb(9)*
     *          (1.0-glb(9)*isith)
      end if
c
      REFS(NN,2)=ATAN(TANTH)*PI
      HALFG(NN)=(AH*TANTH*TANTH+BH*TANTH+CH+DH/COSTH+eh/(tanth*tanth))
      IF (HALFG(NN).GT.0.0) THEN
      HALFG(NN) = SQRT(HALFG(NN))
      ELSE
      WRITE (6,3456) REFS(NN,2),IPH
      WRITE (7,3456) REFS(NN,2),IPH
      WRITE (*,3456) REFS(NN,2),IPH
C      STOP 'SQUARE OF FWHM IS NEGATIVE'
      stop
3456    FORMAT (/,3X,'SQUARE OF FWHM NEGATIVE AT TWO-THETA=',F8.3,
     *     ' FOR PHASE NO. ',I4)
      END IF
      IF (NPROF.EQ.8) THEN
        HALFL(NN) = AH2*TANTH+BH2*SQRT(1.+TANTH*TANTH)
        BB = (HALFG(NN)**5.+2.69269*HALFG(NN)**4.*HALFL(NN)+
     *     2.42843*HALFG(NN)**3.*HALFL(NN)**2.+
     *        4.47163*HALFG(NN)**2.*HALFL(NN)**3. +
     *           0.07842*HALFG(NN)*HALFL(NN)**4. + HALFL(NN)**5.)**0.2
        TLR = HALFL(NN)/BB
        GAM(NN) = 1.36603*TLR-0.47719*TLR*TLR+0.11116*TLR**3.
      ELSE
        BB = HALFG(NN)
      END IF
      TL=BB
      REFS(NN,1)=BB 
      BB=BB*BB
c-----VERT=.TRUE. IF ASYMMETRY CORRECTION IS TO BE CALCULATED
      VERT=.FALSE. 
      if(iasym.eq.0) then
         IF(REFS(NN,2).LE.RLIM.AND.NPROF.NE.5)VERT=.TRUE. 
      else
         IF (ABS(REFS(NN,2)-90.0).GE.RLIM)VERT=.TRUE.
      end if
c-----CALCULATION OF COS(H.X),SIN(H.X) AND TEMP. FACTOR FOR EACH ATOM 
3     DO 8 I=1,N
      SNXI=0.0
      SAI=0.0
      SBI=0.0
      DO 9 J=1,11
9         XI(J)=XL(I+IOF,J)
      DO 10 IR=1,IRL
        DO 301 J=1,3
          IV=IVEC(IPH,IR)/32768/32**(3-J)
          IV=MOD(IV,32)
          SM(J,1)=IV/9-1
          SM(J,2)=MOD(IV/3,3)-1
          SM(J,3)=MOD(IV,3)-1
301         T(J)=FLOAT(MOD(IVEC(IPH,IR)/32**(3-J),32)-16)/12.0
        X=0.0
        DO 6 II=1,3
          YY=0.0
          X=T(II)*HNN(II)+X 
          DO 7 J=1,3
7             YY=SM(J,II)*HNN(J)+YY
6           H(II)=YY
        TR=X
        ARG=TR
        DO 11 J=1,3
11          ARG=H(J)*XI(J)+ARG
        ARG=6.28318530718*ARG
        ARG2=H(1)*H(1)*XI(6)+H(2)*H(2)*XI(7)+
     *        H(3)*H(3)*XI(8)+2.*H(1)*H(2)*XI(9)+2.*
     *              H(1)*H(3)*XI(10)+2.*H(2)*H(3)*XI(11)
        EXPARG=EXP(-ARG2)
        COSA=COS(ARG)*EXPARG
        SINA=SIN(ARG)*EXPARG
        SAI=SAI+COSA
        IF(ICENT.EQ.1)SBI=SINA+SBI
        DO 14 JJ=1,3
          SUMAX(I,JJ)=SUMAX(I,JJ)+H(JJ)*SINA
14          IF(ICENT.EQ.1)SUMBX(I,JJ)=SUMBX(I,JJ)+H(JJ)*COSA
        SUMAX(I,4)=SUMAX(I,4)+H(1)*H(1)*COSA
        SUMAX(I,5)=SUMAX(I,5)+H(2)*H(2)*COSA
        SUMAX(I,6)=SUMAX(I,6)+H(3)*H(3)*COSA
        SUMAX(I,7)=SUMAX(I,7)+H(1)*H(2)*COSA
        SUMAX(I,8)=SUMAX(I,8)+H(1)*H(3)*COSA
        SUMAX(I,9)=SUMAX(I,9)+H(2)*H(3)*COSA
        IF(ICENT.NE.1)GOTO 10
        SUMBX(I,4)=SUMBX(I,4)+H(1)*H(1)*SINA
        SUMBX(I,5)=SUMBX(I,5)+H(2)*H(2)*SINA
        SUMBX(I,6)=SUMBX(I,6)+H(3)*H(3)*SINA
        SUMBX(I,7)=SUMBX(I,7)+H(1)*H(2)*SINA
        SUMBX(I,8)=SUMBX(I,8)+H(1)*H(3)*SINA
        SUMBX(I,9)=SUMBX(I,9)+H(2)*H(3)*SINA
10      CONTINUE
      TEMP(I)=EXP(-XL(I+IOF,4)*SSNN)
      SA(I)=SAI
      SB(I)=SBI
      NI=PTR(I+IOF)
      BNI=DFPP(NI)
      FFX=DFP(NI) 
      AC(10,NI)=0.0
      DO 302 II=1,9,2
302       FFX=FFX+AC(II,NI)*EXP(-AC(II+1,NI)*SSNN)
      SNEXI=FFX*XL(I+IOF,5)*TEMP(I)
      SNXI=BNI*XL(I+IOF,5)*TEMP(I)
c-----CALCULATE A AND B OF F
        AV=SNEXI*SAI+AV
        BV=SNEXI*SBI+BV 
        SNEX(I)=2.*SNEXI*TAV*RATIO(ICX) 
        SNX(I)=2.*SNXI*TAV*RATIO(ICX)
        CV=CV+SNXI*SAI
        DV=DV+SNXI*SBI
8       CONTINUE
      FNN=RATIO(ICX)*(CV*CV+AV*AV+DV*DV+BV*BV)*TAV*sr
C PREPARING PHASE and struc fact TO BE PRINTED  !cp june 98
      TAVix(nn)=tav
      SRix(nn)=sr
      if(av.eq.0)av=1e-6
      aphase(NN)=atan(BV/AV)
C      write(*,2233)av,bv
c2233  format('AV = ',f12.4,' BV = ', f12.4)      
C      if(AV.lt.0.and.BV.GT.0)aphase(nn)=aphase(nn)+3.1415927359
C      if(AV.lt.0.and.BV.lt.0)aphase(nn)=aphase(nn)+3.1415927359
C      if(AV.lt.0.and.BV.eq.0)aphase(nn) =-3.1415927359
      if(AV.lt.0.)aphase(nn)=aphase(nn)+3.1415927359
C      if(AV.gt.0.and.BV.eq.0)aphase(nn) =0.00
      if(BV.lt.0.and.AV.eq.0)aphase(nn) =1.5*3.1415927359
      if(BV.gt.0.and.AV.eq.0)aphase(nn) =0.5*3.1415927359
120   FMGNTD(NN)=FNN
!      write(6,*) 'FMGNTD(nn)=',FMGNTD(nn)                   !ACL debug 980520
      IF(MAXS.EQ.0)RETURN
c-----CALCULATE DERIVATIVES
      DO 17 I=1,N
        SNEXI=SNEX(I)
        SNXI=SNX(I) 
        SAI=SA(I)
        SBI=SB(I)
        DO 22 J=1,11
          K=LP(I+IOF,J)
          IF(K.EQ.0) GOTO22
          IF(J.GT.5) GOTO 221 
          IF(J.GT.3)GOTO 29
          SUMA=SUMAX(I,J)
          SUMB=SUMBX(I,J)
          DER=-((AV*SUMA-BV*SUMB)*SNEXI+(CV*SUMA-DV*SUMB)*SNXI)*6.283185
     *     3071
          GOTO 26
29        IF(J.GT.4)GOTO 31
          DER=-((SAI*AV+SBI*BV)*SNEXI+(SAI*CV+SBI*DV)*SNXI)*SSNN
          GOTO 26
31        DER=((SAI*AV+SBI*BV)*SNEXI+(SAI*CV+SBI*DV)*SNXI)/XL(I+IOF,5)
          GOTO 26
221       SUMA=SUMAX(I,J-2)
          SUMB=SUMBX(I,J-2)
          DER=-((AV*SUMA+BV*SUMB)*SNEXI+(CV*SUMA+DV*SUMB)*SNXI)
          IF(J.GE.9)DER=2.*DER
26        DERIV(K) = SIGN(1.,A(I+IOF,J))*DER+DERIV(K)
22        CONTINUE
17      CONTINUE
c-----CALCULATE DERIVATIVES
c-----Preferred Orientation Derivatives
      K = LPAR(IPH,12)
      IF ( K.NE.0 ) THEN
          DERIV(K) = DERIV(K)+FNN*DPRECOR
C         print '(3x,i3,3f10.5)',k,fnn,dprecor,deriv(k)         ! ********
      END IF
      K = LPAR(IPH,13)
      IF (IPREF.EQ.0) THEN
        IF(K.NE.0)DERIV(K)=DERIV(K)+(1.-PREXP)*FNN/PRECOR
      ELSE
          IF (K.NE.0) THEN
             WRITE (6,7607)
             STOP
          ENDIF
      ENDIF
7607    FORMAT (1X,'G2 IS NOT A REFINABLE PARAMETER FOR IPREF = 1')
c-----Derivatives for microabsorption parameter
      if (iabsr.eq.1) then
        k = lglb(13)
         srd = (1.0-glb(12))*(asin(sqrt(sinth))-1.5707963268)
         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
        k = lglb(12)
         srd = 1.0 -glb(8)*exp(-glb(9))+glb(8)*exp(-glb(9)/
     *         sqrt(sinth))-1.0-glb(13)*(asin(sqrt(sinth))-1.5707963268)
         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
        k = lglb(9)
         srd = glb(8)*glb(12)*(exp(-glb(9))-exp(-glb(9)/sqrt(sinth))/
     *          sqrt(sinth))
         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
        k = lglb(8)
         srd = -glb(12)*(exp(-glb(9))+exp(-glb(9)/sqrt(sinth)))
         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
      else if (iabsr.eq.2) then
        kkl = lglb(12)
        k  = lglb(9)
        kl = lglb(8)
         if (k.ne.0.or.kl.ne.0.or.kkl.ne.0) then
          write (6,9598)
          write (*,9598)
          stop
         end if
9598   format (' P AND/OR Q AND/OR R ARE NOT REFINABLE PARAMETERS FOR ',
     *         'IABSR=2')
9599   format (' R AND/OR T IS NOT A REFINABLE PARAMETER FOR THE IABSR',
     *         ' CHOICE')
        k = lglb(13)
         srd = 1.5707963268 - asin(sqrt(sinth))
         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
      else if (iabsr.eq.3) then
        kkl = lglb(13)
        k = lglb(12)
         if (k.ne.0.or.kkl.ne.0) then
          write (6,9599)
          write (*,9599)
          stop
         end if
        k = lglb(9)
         srd = glb(8)*exp(-glb(9))-glb(8)*isith*exp(-glb(9)*isith)
C         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn
        k = lglb(8)
         srd = -exp(-glb(9))+exp(-glb(9)*isith)
C         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
         if (k.ne.0) deriv(k) = deriv(k) + srd*fnn
      else if (iabsr.eq.4) then
       kkl = lglb(13)
       k = lglb(12)
        if (k.ne.0.or.kkl.ne.0) then
          write (6,9599)
          write (*,9599)
          stop
        end if
       k = lglb(9)
        srd = glb(8)*(2*glb(9)-1)+glb(8)*isith*(2*glb(9)*isith-1)
        if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
       k = lglb(8)
        srd = glb(9)*(glb(9)-1)-glb(9)*isith*(1-glb(9)*isith)
        if (k.ne.0) deriv(k) = deriv(k) + srd*fnn/sr
      end if
c----Overall Temperature and Scale Factor
      K=LPAR(IPH,2)
      IF(K.NE.0)DERIV(K)=DERIV(K)-2.*SSNN*FNN
      K=LPAR(IPH,1)
      IF(K.NE.0)DERIV(K)=DERIV(K)+FNN/PAR(IPH,1)
      SINTH=FNN*PI*SLABDA/(SQRT(SINTH*COSTH)*BB)
      SS=FNN/TL
      X=TANTH*TANTH
c-----Broadening Derivatives
      if (nprof.eq.5.or.nprof.eq.8) goto 9212
      DO 78 J=3,5
        K=LPAR(IPH,J)
        IF(K.EQ.0) GOTO 78
        DERIV(K)=X*SS+DERIV(K)
78    X=X/TANTH   
         k=lpar(iph,21)
         if(k.eq.0) goto 9212
C for cot^2 case                                    !cp nov 29 96
         deriv(k)=ss/(tanth*tanth)+deriv(k)
c-----Split Pearson VII Broadening Derivatives
9212  if (nprof.eq.5) then
      if (delta.lt.0.0) then
         da3 = da3l
         da1 = da1l
      else
         da3 = da3h
         da1 = da1h
      end if
      DO 780 J=3,5
          K=LPAR(IPH,J)
          IF(K.EQ.0) GOTO 780
          DERIV(K)=DERIV(K)+x*ss*((da3*delt/(1.0+da1*delt))-(1.0/tl))
780     X=X/TANTH
      end if
9211    CONTINUE
c-----TCHZ Broadening Derivatives
      IF (NPROF.EQ.8) THEN
        TL = REFS(NN,1)
        TLG = HALFG(NN)
        TLL = HALFL(NN)
        DHDHG = 0.2/TL**4.*(5.*TLG**4.+10.77076*TLG**3.*TLL+
     *        7.28529*TLG*TLG*TLL*TLL+8.94326*TLG*TLL**3.
     *         + 0.07842*TLL**4.)
        DHDHL = 0.2/TL**4.*(2.69269*TLG**4.+ 4.85686*TLG**3.*TLL
     *      +13.41489*TLG*TLG*TLL*TLL + 0.31368*TLG*TLL**3.+5.*TLL**4.)
        DO 9078 J=3,5
          K=LPAR(IPH,J)
          IF(K.EQ.0)GOTO 9078
          DERIV(K)=DHDHG*X*SS+DERIV(K)
9078      X=X/TANTH
        K=LPAR(IPH,20)
        IF  (K.NE.0) DERIV(K)=DHDHG*SS/COSTH+DERIV(K)
        K = LPAR(IPH,15)
        IF (K.EQ.0) GO TO 9213
        DERIV(K) = 2.*FNN*DHDHL*TANTH+DERIV(K)
9213    CONTINUE
        K = LPAR(IPH,16)
        IF (K.EQ.0) GO TO 9214
        DERIV(K) = 2.*FNN*DHDHL/SQRT(COSTH) + DERIV(K)
9214    CONTINUE
      END IF
c-----Profile Shape Derivatives
      K = LPAR(IPH,17)
      IF(K.NE.0.AND.(NPROF.EQ.6.OR.NPROF.EQ.7)) DERIV(K)=DERIV(K)+ FNN
      K = LPAR(IPH,18)
      IF(K.NE.0.AND.NPROF.EQ.6) DERIV(K)=DERIV(K)+ FNN * REFS(NN,2)
      IF(K.NE.0.AND.NPROF.EQ.7) DERIV(K)=DERIV(K)+ FNN / REFS(NN,2)
      K = LPAR(IPH,19)
      IF(K.NE.0.AND.NPROF.EQ.7)DERIV(K)=DERIV(K)+FNN/REFS(NN,2)
     *               /REFS(NN,2)
c-----Split Pearson VII Shape Derivative
      if (nprof.eq.5) then
        k = lpar(iph,17)
          if (k.ne.0.0) then
             if (delta.lt.0.0) then
                deriv(k) = deriv(k)+fnn*(-alog(1.0+da1l*delt/bb)+
     *                da7l*delt/bb/(1.0+da1l*delt/bb))
             else
                deriv(k) = deriv(k)+fnn*da6l
          end if
      end if
        k = lpar(iph,18)
        if (k.ne.0.0) then
           if (delta.lt.0.0) then
               deriv(k) = deriv(k)+fnn*(-alog(1.0+da1l*delt/bb)+
     *                da7l*delt/bb/(1.0+da1l*delt/bb))/refs(nn,2)
           else
               deriv(k) = deriv(k)+fnn*da6l/refs(nn,2)
           end if
        end if
      k = lpar(iph,19)
        if (k.ne.0.0) then
           if (delta.lt.0.0) then
              deriv(k) = deriv(k)+fnn*(-alog(1.0+da1l*delt/bb)+
     *        da7l*delt/bb/(1.0+da1l*delt/bb))/refs(nn,2)/refs(nn,2)
           else
              deriv(k) = deriv(k)+fnn*da6l/refs(nn,2)/refs(nn,2)
           end if
        end if
      k = lpar(iph,24)
        if (k.ne.0.0) then
           if (delta.lt.0.0) then
                 deriv(k) = deriv(k)+fnn*da6h
           else
                 deriv(k) = deriv(k)+fnn*(-alog(1.0+da1h*delt/bb)+
     *                da7h*delt/bb/(1.0+da1h*delt/bb))
           end if
        end if
      k = lpar(iph,25)
        if (k.ne.0.0) then
           if (delta.lt.0.0) then
            deriv(k) = deriv(k)+fnn*da6h/refs(nn,2)
           else
            deriv(k) = deriv(k)+fnn*(-alog(1.0+da1h*delt/bb)+
     *                da7h*delt/bb/(1.0+da1h*delt/bb))/refs(nn,2)
           end if
        end if
      k = lpar(iph,26)
        if (k.ne.0.0) then
           if (delta.lt.0.0) then
              deriv(k) = deriv(k)+fnn*da6h/refs(nn,2)/refs(nn,2)
           else
              deriv(k) = deriv(k)+fnn*(-alog(1.0+da1h*delt/bb)+
     *          da7h*delt/bb/(1.0+da1h*delt/bb))/refs(nn,2)/refs(nn,2)
           end if
        end if
      end if
c-----Split Pearson VII Asymmetry Derivative
      k = lpar(iph,27)
      if (delta.lt.0.0) then
          da1 = da1l
          da4 = da4l
          da5 = da5l
      else
          da1 = da1h
          da4 = da4h
          da5 = da5h
      end if
      if (k.ne.0.and.nprof.eq.5) then
         deriv(k)=deriv(k)+fnn*(da4+da5*delt/bb/(1.0+da1*delt/bb))
      end if
c-----Zero, Displacement, and Transparancy Derivatives
      K=LGLB(1)
      IF(K.NE.0)DERIV(K)=DERIV(K)+2.0*FNN/BB
      K=LGLB(10)
      IF (K.NE.0) DERIV(K)=DERIV(K)+2.*FNN/BB*SQRT(COSTH) 
      K=LGLB(11)
      IF (K.NE.0) DERIV(K)=DERIV(K)+2.*FNN/BB*SIN(REFS(NN,2)/57.2958)
c-----Lattice Parameter Derivatives
      DO 79 J=1,6
      K=LPAR(IPH,J+5)
      IF(K.EQ.0) GOTO 79
      IF(J.LT.4) X=HNN(J)*HNN(J)
      IF(J.EQ.4) X=HNN(2)*HNN(3)
      IF(J.EQ.5) X=HNN(1)*HNN(3)
      IF(J.EQ.6) X=HNN(1)*HNN(2)
      DERIV(K)= X*SINTH + DERIV(K)
79    CONTINUE
c-Asymmetry Derivative.  Test for asymmetry model included !cp may 97  
          K=LPAR(IPH,14)
      IF((K.NE.0).AND.VERT)then
             if (iasym.eq.0) then
               DERIV(K)=-FNN/TANTH+DERIV(K)
             else
               TANTHE=TANTH
               IF (TANTHE.GE.1.)TANTHE=TAN(ATAN(TANTHE-3.14159265359/2))
               DERIV(K)=-FNN/TANTH+DERIV(K)
             end if
      END IF
c-----STORE DERIVATIVES FOR LIMO REFLECTIONS AT A TIME
      DO 80 I=1,MAXS
80      DERSTO(NM,I)=DERIV(I) 
48    RETURN
      END 
      FUNCTION PROFIL(N,X)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
C      COMMON/PRFX/DELT,TL,GAM1,GAM2,PRFDER,IPH
      COMMON/PRFX/DELT,TL,GAM1,GAM2,PRFDER,iph,delta
C      common/prfx/iph,delta
      COMMON /PVII/ TF1,TF2,TF4,TF6,TF8,TF9,C4
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      common/spvii/rl,da1l,da2l,da3l,da4l,da5l,da6l,da7l
      common/spvii/rh,da1h,da2h,da3h,da4h,da5h,da6h,da7h
      GOTO(10,20,30,40,50,60,70,80),N
      STOP 'ILLEGAL PROFILE FUNCTION REQUEST'
10    PROFIL=.939437279*EXP(-2.772588722*X)
      PRFDER=2.772588722
      RETURN
20    PROFIL=.636619772/(1.+4.*X)
      PRFDER=4./(1.+4.*X)
      RETURN
30    PROFIL=.819449653/(1.+1.656854248*X)**2.
      PRFDER=3.313708496/(1.+1.656854248*X)
      RETURN
40    PROFIL=.766420937/(1.+2.349604208*X)**1.5
      PRFDER=3.5244063/(1.+2.349604208*X)
      RETURN
50    if (delta.lt.0.0) then
        profil = (1.0+da1l*x)**(-rl)*da2l
        prfder = da1l*rl/(1.0+da1l*x)
      else
        profil = (1.0+da1h*x)**(-rh)*da2h
        prfder = da1h*rh/(1.0+da1h*x)
      end if
      RETURN
60    PROFIL=GAM1*.636619772/(1.+4.*X)+ 
     * (1-GAM1)*.939437279*EXP(-2.772588722*X)
      PRFDER=(GAM1*2.546479088/(1.+4.*X)**2+(1-GAM1)*2.6046732048*
     * EXP(-2.772588722*X))/PROFIL
      RETURN
70    PROFIL=C4/(1.+4.*(2**(1./GAM1)-1)*X)**GAM1
      PRFDER=TF6*GAM1/(1.+4.*(2**(1./GAM1)-1)*X)
      RETURN
80    PROFIL=GAM1*.636619772/(1.+4.*X)+ 
     * (1-GAM1)*.939437279*EXP(-2.772588722*X)
      PRFDER=(GAM1*2.546479088/(1.+4.*X)**2+(1-GAM1)*2.6046732048*
     * EXP(-2.772588722*X))/PROFIL
      RETURN
      END 
      SUBROUTINE  PRSVII(T)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
C     CALCULATES THE COEFFICIENT C4 AND ITS DERIVATIVE WRT T FOR PHASE K
C     THIS SUBROUTINE OBTAINED FROM IMMIRZI'S CODE
C     RGD = GAMMA(T)/GAMMA(T-0.5)
C     DGD = D(RGD)/D(T)
C     COMMON/PRFX/DELT,TL,GAM1,GAM2,PRFDER,IPH
      COMMON/PRFX/DELT,TL,GAM1,GAM2,PRFDER,iph,delta
C      common/prfx/iph,delta
      COMMON /PVII/ TF1,TF2,TF4,TF6,TF8,TF9,C4
      DIMENSION  RGD(52),DGD(52)
      DATA  RGD/   .156535,  .183771,  .209926,  .235089,  .259339,
     1   .282749,  .305380,  .327289,  .348527,  .369141,  .389170,
     2   .408654,  .427625,  .446115,  .464153,  .481764,  .498972,
     3   .515799,  .532265,  .548390,  .564190,  .579680,  .594877,
     4   .609794,  .624443,  .638837,  .652986,  .666902,  .680594,
     5   .694071,  .707342,  .720415,  .733297,  .754996,  .758519,
     6   .770871,  .783059,  .795089,  .806966,  .818696,  .830282,
     7   .841731,  .853045,  .864230,  .875290,  .886227,  .897046,
     8   .907750,  .918344,  .928828,  .939207,  .949484/
      DATA  DGD/  1.390576, 1.334036, 1.282282, 1.234747, 1.190975,
     1  1.150537, 1.113079, 1.078304, 1.045950, 1.015831,  .987556,
     2   .961144,  .936314,  .913013,  .891052,  .870228,  .850726,
     3   .832193,  .814684,  .798032,  .782274,  .767112,  .752807,
     4   .739098,  .726022,  .713542,  .701621,  .690110,  .679232,
     5   .668690,  .658669,  .648871,  .639521,  .630580,  .621825,
     6   .613555,  .605434,  .597760,  .590086,  .582896,  .575967,
     7   .569075,  .562482,  .555962,  .549890,  .543743,  .538044,
     8   .532381,  .526905,  .521578,  .516400,  .511296/
      DATA  PIG/3.1415926/, AL2/0.6931472/
      IF (T.LT.0.6) THEN
      WRITE (6,2) T
2       FORMAT (//,' UNACCEPTABLE M VALUE IN CALCULATING PEARSON VII'
     *,' FUNCTION ',F6.2,//)
      STOP
      END IF
      IT=IFIX(T-0.6)
      FT=T-FLOAT(IT)
      N1=(FT-0.6)/0.02+1.0001 
      DG=DGD(N1)+(DGD(N1+1)-DGD(N1))*(FT-0.58-0.02*FLOAT(N1))/0.02
      RG=RGD(N1)+(RGD(N1+1)-RGD(N1))*(FT-0.58-0.02*FLOAT(N1))/0.02
      IF(IT.GT.0) THEN
      DO 4 I=1,IT
      FL1=FT+FLOAT(I-1)
      DG=DG*(FL1)/(FL1-0.5)-0.5*RG/(FL1-0.5)**2
4     RG=RG*FL1/(FL1-0.5)
      ENDIF
      TF1=SQRT(2.0**(1.0/T)-1.0)
      C4=2.0*RG*TF1/SQRT(PIG) 
      TF2=AL2*2.0**(1.0/T)/(2.0*TF1**2) 
C     DC4DT= 2./SQRT(PIG)*(TF1*DG - RG*TF2/T**2)
      TF4=4.0*AL2*2.0**(1.0/T)/T
      TF6=4.0*T*TF1**2
      TF8=DG/RG-TF2/T**2
      TF9=4.0*TF1**2
C     TF3=360.0*TF6/PIG
      RETURN
      END
C    Split Pearson VII coding for function and derivatives
C    Based on Toraya's Code from Profit, Ver 1.22N
      subroutine mspvii(a,w)
      common/spvii/rl,da1l,da2l,da3l,da4l,da5l,da6l,da7l
      common/spvii/rh,da1h,da2h,da3h,da4h,da5h,da6h,da7h
      dl = gamma(rl-0.5)/(sqrt(2.0**(1.0/rl)-1.0)*gamma(rl))
      dh = gamma(rh-0.5)/(sqrt(2.0**(1.0/rh)-1.0)*gamma(rh))
      dld = gamma(rl-0.5+0.001)/(sqrt(2.0**(1.0/(rl+0.001))-1.0)*
     *           gamma(rl+0.001))
      dhd = gamma(rh-0.5+0.001)/(sqrt(2.0**(1.0/(rh+0.001))-1.0)*
     *           gamma(rh+0.001))
      da1l = (((1.0+a)/a)**2.0)*(2.0**(1.0/rl)-1.0)
      da1h = ((1.0+a)**2.0)*(2.0**(1/rh)-1.0)
      da2l = 1.128379167*(1.0+a)*(1.0/(a*dl+dh))/w
      da2h = da2l
      da3l = 2.0*rl*da1l/w
      da3h = 2.0*rh*da1h/w
      da4l = (1.0/(1.0+a))-(dl/(a*dl+dh))
      da4h = da4l
      da5l = 2.0*rl*((1+a)/a)**2.0*(1.0+a)/(a**3)
      da5h = -2.0*rh*(1+a)**2.0*(1+a)
      da6l = -a*1000.0*(dld-dl)/(a*dl+dh)
      da6h = -1000.0*(dhd-dh)/(a*dl+dh)
      da7l = alog(2.0)*((1.0+a)/a)**2.0*(2.0**(1.0/rl)-1.0+1.0)/rl
      da7h = alog(2.0)*(1.0+a)**2.0*(2.0**(1.0/rh)-1.0+1.0)/rh
      return
      end
c
      SUBROUTINE outscr(I,R2,r3,x) 
      print '(a,i2)',' CYCLE NUMBER - ',i
      print '(4(a,f8.2,a))', ' R-P =',r2,'%','    R-WP =',r3,'%',
     +  '       R-EXPECTED =',x,'%','  S = ',r3/x
      RETURN
      END 
C 
      SUBROUTINE OUTPTR(ICYCLE)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      INTEGER PTR,fondo
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
C      COMMON/ALLP/FINAL(8*MSZ,2),ILOC  
C       size of FINAL is too small in case of many atoms !!!
      COMMON/ALLP/FINAL(nfinal,2),ILOC
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/G3/COND,IORD1,IORD2,TH,NUM 
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,PLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/G2/S1,S2,ss2,S3,ss4,D1,D2,D4,R1,R2,R3,r2nobk,r3nobk
      COMMON/F1/SMM(MSZ,MSZ),V(MSZ)
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),xmas(16)
      COMMON/MULTIP/TMASSA(99),MLTPHASE,xmltp(99),MURT(NATS),SAQF(99)
     *,wtis(99)
      common/dirCV/ DCSM(6,6),DCV(6)      
      DIMENSION SY(30),SZ(30),DUMMY(2*MSZ+4)
      DIMENSION VOL(99), DVOL(99), W(99), DW(99), XMASS(99), DMASS(99)
      DIMENSION FR(99), FRP(99), dfr(99)  
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      dimension stmassa(99),smass(99)
      common/sizestrain/sizeG(15),strainG(15),sizeL(15),strainL(15)
     *,siz(15),strain(15),NsizeStrain
      DO 8123 I=1,2*MSZ+4
8123    DUMMY(I) = 0.
      WRITE(6,57)R2,R3,r3nobk
      WRITE(7,57)R2,R3,r3nobk
c      WRITE(6,57)R2,r2nobk,R3,r3nobk                      !removed 08jan2001
c      WRITE(7,57)R2,r2nobk,R3,r3nobk
57    FORMAT(/,
     *  ' R-P        = ',F8.2,'%',/
     *  ' R-WP       = ',F8.2,'%',5x,
     *  ' R-WP(Background Removed) = ',F8.2,'%')
c57    FORMAT(/,                                            !removed 08jan2001
c     *  ' R-P        = ',F8.2,'%',5x,
c     *  ' R-P(Background Removed) = ',F8.2,'%',/,
c     *  ' R-WP       = ',F8.2,'%',5x,
c     *  ' R-WP(Background Removed) = ',F8.2,'%')
      if ( maxs.eq.0 .and. maxsx.ne.0 ) then
        X=100.*SQRT((FLOAT(NUM)-FLOAT(MAXSX))*1./D2)
      else
        X=100.*SQRT((FLOAT(NUM)-FLOAT(MAXS))*1./D2)
      end if
      IF ( MAXS.EQ.0 .AND. MCYCLE.EQ.1 ) THEN
        CALL outscr(MCYCLX,R2,R3,X)
C        CALL outscr(MCYCLX,R3)
      ELSE
        CALL outscr(ICYCLE-1,R2,R3,X)
C        CALL outscr(ICYCLE-1,R3)
      ENDIF
C     CALL DISP(ICYCLE,R3)
      WRITE(6,58)X
58    FORMAT(' R-EXPECTED = ',F8.2,'%')
      WRITE (6,66) SQRT(S2/FLOAT(NUM-MAXS)), ss4/ss2
66    FORMAT (
     *' S          = ',F8.2,5X,'SQRT(RESIDUAL/N-P)','GOODNESS OF FIT',/,
     *' D - W D    = ',F8.2,5X,'UNWEIGHTED DURBIN-WATSON STATISTIC D')
      I=NUM-MAXS
      WRITE(6,62)I
      WRITE(6,63)S1,D1,S3,D2,S2,COND
62    FORMAT( 13H0N-P        =,1X,I8)
63    FORMAT(/,7X,'SUMYDIF', 7X,'SUMYOBS', 6X,'SUMYCALC', 6X,'SUMWYOBS'
     * ,'SQ',6X,'RESIDUAL',6X,'CONDITION',/,1X,6E14.4/)
      DUMMY(2*MSZ+1) = R2
      DUMMY(2*MSZ+2) = R3
      DUMMY(2*MSZ+3) = SQRT(S2/FLOAT(NUM))
      DUMMY(2*MSZ+4) = ss4/ss2
C     FINAL PARAMETERS AND R-FACTORS
      IF (IPLST.NE.0.AND.MAXS.EQ.0) THEN
        DO 410 IP=1,NPHASE
          IOF=0
          IF(IP.GT.1) THEN
            DO 44340 IIPHAS=2,IP
44340         IOF = IOF + NATOM(IIPHAS-1)
          END IF
          N=NATOM(IP)
          DO 604 I=1,N
            DO 604 J=1,5
              KM=LP(I+IOF,J)
              IF(KM.NE.0) THEN
C  !cp jun 96 start   
                 if(j.eq. 5) then
                   dummy(km) = xl(i+iof,j)*xmltp(ip)/murt(i+iof)
                 else  
C !cp jun 96 stop
                   DUMMY(KM)  =XL(I+IOF,J)         
C !cp sept 96 start
                 endif                                    
C !cp sept 96 stop
              ENDIF
604           CONTINUE
        DO  607 I=1,N
          DO  607 J=6,11
            KM=LP(I+IOF,J)
C              IF(KM.NE.0) DUMMY(KM)  =XL(I+IOF,J) 
            IF(KM.NE.0) then                            
C  !cp jun 96 start   
               if(j.eq. 5) then
                 dummy(km) = xl(i+iof,j)*xmltp(ip)/murt(i+iof)
               endif  
C !cp jun 96 stop              
               DUMMY(KM)  =XL(I+IOF,J)
            endif                                        
607         CONTINUE
        DO 6017 J=1,27
          KM=LPAR(IP,J)
          IF(KM.NE.0) DUMMY(KM) = PAR(IP,J)
6017        CONTINUE
      CALL DIRECT(DCSM,DCV,IP)
      DO 30030 I=1,6
        KM = LPAR(IP,I+5)
        MATCH = 0
        IF (KM.NE.0) THEN
          IF (I.GE.4) THEN
            DO 30031 MMM=1,3
            IF (KM.EQ.LPAR(IP,MMM+5)) MATCH=1
30031         CONTINUE
          ENDIF
          IF (MATCH.EQ.0) DUMMY(KM) = SAVE(IP,I)
        END IF
30030     continue
410       CONTINUE
C      DO 6019 J=1,13   !cp may 97
      DO 6019 J=1,20
C        IF(NPROF.EQ.8.AND.(J.EQ.8.OR.J.EQ.9)) GO TO 6019
C        IF(NPROF.EQ.7.AND.J.EQ.8) GO TO 6019
C        IF(NPROF.EQ.6.AND.J.EQ.8) GO TO 6019
C        if(nprof.eq.5.and.(j.eq.8.or.j.eq.9)) go to 6019
C glb(8) and gkb(9) are now used to surface roughness
C and glb(14), glb(15), glb(16) aree unused  !cp may 97
        IF(NPROF.EQ.8.AND.(J.EQ.14.OR.J.EQ.15.or.j.eq.16)) GO TO 6019
        IF(NPROF.EQ.7.AND.(J.EQ.14.OR.J.EQ.15.or.j.eq.16)) GO TO 6019
        IF(NPROF.EQ.6.AND.(J.EQ.14.OR.J.EQ.15.or.j.eq.16)) GO TO 6019
        if(nprof.eq.5.and.(J.EQ.14.OR.J.EQ.15.or.j.eq.16)) go to 6019
        KM=LGLB(J)
        IF(KM.NE.0) DUMMY(KM) = GLB(J)
6019      CONTINUE
      WRITE (8,ERR=99990) (DUMMY(I),I=1,MAXSX)
     *       ,(DUMMY(I),I=MSZ+1,MSZ+MAXSX),(DUMMY(I),I=2*MSZ+1,2*MSZ+4)
      DO 9487 I=1,4
      IF (MOD(I,2).EQ.1) ILOC = ILOC + 1
9487    FINAL(ILOC,2-MOD(I,2)) = DUMMY(2*MSZ+I)
      END IF
      IF (MAXS.EQ.0.AND.MCYCLE.EQ.1) RETURN
      ILOC = 0
C      DO 8127 I=1,8*MSZ
      do 8127 i=1,nfinal
      DO 8127 J=1,2
8127      FINAL(I,J) = 0.
      WRITE(6,1)ICYCLE
      WRITE(7,111)ICYCLE
111   FORMAT(/14H CYCLE NUMBER=,I4)
1     FORMAT(1H1,120(1H+)//14H CYCLE NUMBER=,I4)
      DO 10 IP=1,NPHASE
      IOF=0
      IF(IP.GT.1) THEN
        DO 4340 IIPHAS=2,IP 
4340        IOF = IOF + NATOM(IIPHAS-1) 
      END IF
      WRITE(6,3)IP,PHSNM(IP)
3       FORMAT(7H0PHASE ,I2,2H: ,A50,/' NEW PARAMETERS, SHIFTS, AND ',
     *     'STANDARD DEVIATIONS='//5H ATOM,5X,1HX,8X,2HDX,7X,2HSX,8X,1HY
     *   ,8X, 2HDY,7X,2HSY,8X,1HZ,8X,2HDZ,7X,2HSZ,8X,1HB,6X,2HDB,5X,2HSB
     *   ,       6X,2HSo,5X,3HDSo,4X,3HSSo)
      N=NATOM(IP) 
      tmassa(ip)=0.   
      stmassa(ip)=0.
      DO  4 I=1,N 
        DO  5 J=1,5
          ILOC = ILOC + 1
          KM=LP(I+IOF,J)
          IF(KM.NE.0)GOTO  6
          SY(J)=0.
          SZ(J)=0.
C  !cp jun 97 start
      if(j.eq.5) then
        final(iloc,1) = xl(i+iof,j)*xmltp(ip)/murt(i+iof)
      else
          FINAL(ILOC,1) = XL(I+IOF,J) 
      endif
C  !cp jun 97 stop
          GOTO  5 
6           SZ(J)=SQRT(ABS(SMM(KM,KM))) 
          SY(J)=V(KM)*A(I+IOF,J)*RELAX(1)
C  !cp jun 96 start   
           if(j.eq. 5) then
              dummy(km) = xl(i+iof,j)*xmltp(ip)/murt(i+iof)
              XL(I+IOF,J)=XL(I+IOF,J)+SY(J)
              DUMMY(KM+MSZ)  = SY(J)*xmltp(ip)/murt(i+iof)
              FINAL(ILOC,1) = xl(i+iof,j)*xmltp(ip)/murt(i+iof)
              final(iloc,2) = sz(j)*xmltp(ip)/murt(i+iof)
           else
C !cp jun 96 stop
              DUMMY(KM)  =XL(I+IOF,J)
              XL(I+IOF,J)=XL(I+IOF,J)+SY(J)
              DUMMY(KM+MSZ)  = SY(J)
              FINAL(ILOC,1) = XL(I+IOF,J)                       
              FINAL(ILOC,2) = SZ(J)
                      smass(ip)=sz(5)
           endif
5       CONTINUE

C !cp jun 96 start
4     WRITE(6,100)ATEXT(I+IOF),(XL(I+IOF,J),SY(J),SZ(J),j=1,4),
     *          xl(i+iof,5)*xmltp(ip)/murt(i+iof),SY(5)*xmltp(ip)/
     *          murt(i+iof),SZ(5)*xmltp(ip)/murt(i+iof)
C !cp jun 96 stop
c100     FORMAT(1X,A4,2X,5(3F7.4,1X))
c !CP 10 jun 99. format below changed from 5((3F7.4,1X)) to that
100     FORMAT(1X,A4,2X,3(3F9.6,1X),2(3F7.4,1X) )
      WRITE(6,101)
101     FORMAT(5H0ATOM,7X,3HB11,6X,4HDB11,6X,4HSB11,7X,3HB22,6X,4HDB22,
     *     6X,4HSB22,7X,3HB33,6X,4HDB33,6X,4HSB33,/,12X,3HB12,6X,4HDB12,
     *   6X,  4HSB12,7X,3HB13,6X,4HDB13,6X,4HSB13,7X,3HB23,6X,4HDB23,6X,
     *   4HSB23)
      DO  7 I=1,N 
        DO  8 J=6,11
          ILOC = ILOC + 1
          KM=LP(I+IOF,J)
          IF(KM.NE.0)GOTO  9
          SY(J)=0.0
          SZ(J)=0.0
          FINAL(ILOC,1) = XL(I+IOF,J) 
          GOTO  8 
9         SZ(J)=SQRT(ABS(SMM(KM,KM))) 
          SY(J)=V(KM)*A(I+IOF,J)*RELAX(2)
C  !cp jun 96 start   
           if(j.eq. 5) then
              dummy(km) = xl(i+iof,j)*xmltp(ip)/murt(i+iof)
              XL(I+IOF,J)=XL(I+IOF,J)+SY(J)
              DUMMY(KM+MSZ)  = SY(J)*xmltp(ip)/murt(i+iof)
              FINAL(ILOC,1) = xl(i+iof,j)*xmltp(ip)/murt(i+iof)
              final(iloc,2)=sz(j)*xmltp(ip)/murt(i+iof)    
             else
              DUMMY(KM)  =XL(I+IOF,J)
              XL(I+IOF,J)=XL(I+IOF,J)+SY(J)
              DUMMY(KM+MSZ)  = SY(J)
              FINAL(ILOC,1) = XL(I+IOF,J) 
              FINAL(ILOC,2) = SZ(J)
             endif
C !cp jun 96 stop
8           CONTINUE
7         WRITE(6,102)ATEXT(I+IOF),(XL(I+IOF,J),SY(J),SZ(J),J=6,11)
102     FORMAT(1X,A4,9F10.6,/,5X,9F10.6)
      DO 17 J=1,27
        IF ( J.EQ.12 ) THEN
          ILOC = ILOC + 1
          FINAL(ILOC,1) = PREF(IP,1)
          FINAL(ILOC,2) = PREF(IP,2)
          ILOC = ILOC + 1
          FINAL(ILOC,1) = PREF(IP,3)
        END IF
        ILOC = ILOC + 1
        IF ( J.EQ.6 ) ILOC1 = ILOC
C         IF (J.EQ.11) ILOC2 = ILOC
        KM = LPAR(IP,J)
        IF ( KM.EQ.0 ) THEN
          SY(J) = 0.
          SZ(J) = 0.
          FINAL(ILOC,1) = PAR(IP,J)
        ELSE
          SZ(J) = SQRT(ABS(SMM(KM,KM)))
          SY(J) = V(KM)*APAR(IP,J)*RELAX(3)
          DUMMY(KM) = PAR(IP,J)
          PAR(IP,J) = PAR(IP,J)+SY(J)
          DUMMY(KM+MSZ)  = SY(J)
          FINAL(ILOC,1) = PAR(IP,J)
          FINAL(ILOC,2) = SZ(J)
        END IF
17      CONTINUE
      CALL DIRECT(DCSM,DCV,IP)
      ILOC2 = ILOC
      ILOC = ILOC1-1
      DO 30 I=1,6 
        ILOC = ILOC + 1
        KM = LPAR(IP,I+5)
        MATCH=0
        SY(I+5)=DCV(I)-SAVE(IP,I)
        IF (KM.NE.0) THEN
         IF (I.GE.4) THEN
            DO 30032 MMM=1,3
            IF (KM.EQ.LPAR(IP,MMM+5)) MATCH=1
30032         CONTINUE
          ENDIF
          IF (MATCH.EQ.0) THEN
            DUMMY(KM) = SAVE(IP,I)
            DUMMY(KM+MSZ) = SY(I+5)
          END IF
        END IF
        FINAL(ILOC,1) = DCV(I)
        FINAL(ILOC,2) = DCSM(I,I)
30        SAVE(IP,I)=DCV(I)
      ILOC = ILOC2
c==============================================      
      if(dcv(1).eq.dcv(2).and.dcv(2).ne.dcv(3))then
c         if(dcv(4).eq.dcv(5).and.dcv(5).eq.90.) then
         if(dcv(4).eq.dcv(5).and.dcv(5).ne.90.0) then
            if(dcv(6).eq.90.0.or.dcv(6).eq.120.0) then
                ratiodcv=dcv(3)/dcv(1)   
                sratiodcv=(dcv(1)*dcsm(3,3)+dcv(3)*dcsm(1,1))/(dcv(3)
     1                    *dcv(3))
                WRITE(6,1031) PAR(IP,1),SY(1),SZ(1),
     *                 PAR(IP,2),SY(2),SZ(2),
     *                   (DCV(I),SY(I+5),DCSM(I,I),I=1,6),
     *                        ratiodcv,sratiodcv,
     *                      (PAR(IP,I),SY(I),SZ(I),I=12,16)
1031     FORMAT( '0OVERALL SCALE FACTOR=',3G9.3,/,
     *         '0OVERALL TEMP. FACTOR=',3F9.4,/,
     *        '0CELL PARAMETERS=',/,3(3F11.6,/),3(3F11.4,/),/,
     *        '             c/a= ' ,1f11.6, 2x,'+/-',1f11.6,//,
     *      '0PREFERRED ORIENTATION PARAMETERS=',/,2(1X,3F8.5,/),
     *    '0ASYMMETRY PARAMETER=',/,1X,3F8.4,/,
     *  '0LORENTZIAN HALF WIDTH PARAMS (X AND Y) ',/,2(1X,3F8.5,/))     
                goto 1032
            endif  
         endif 
      endif
c      endif
c cpaiva 15set2000      
c=============================================
      WRITE(6,103) PAR(IP,1),SY(1),SZ(1),
     *                 PAR(IP,2),SY(2),SZ(2),
     *                   (DCV(I),SY(I+5),DCSM(I,I),I=1,6),
     *                      (PAR(IP,I),SY(I),SZ(I),I=12,16)
103     FORMAT( '0OVERALL SCALE FACTOR=',3G9.3,/,
     *         '0OVERALL TEMP. FACTOR=',3F9.4,/,
     *        '0CELL PARAMETERS=',/,3(3F11.6,/),3(3F11.4,/),
     *      '0PREFERRED ORIENTATION PARAMETERS=',/,2(1X,3F8.5,/),
     *    '0ASYMMETRY PARAMETER=',/,1X,3F8.4,/,
     *  '0LORENTZIAN HALF WIDTH PARAMS (X AND Y) ',/,2(1X,3F8.5,/))
c-----CHECK FOR THE SPLIT PEARSON VII PROFILE
1032   if(nprof.eq.5) then
         write(6,106) (par(ip,i),sy(i),sz(i),i=17,19),
     1               (par(ip,i),sy(i),sz(i),i=24,26),
     1             (par(ip,i),sy(i),sz(i),i=27,27)
106      format( '0LOW SIDE EXPONENT PARAMETERS (NA, NB, NC)=',/,
     1    3(1x,3g10.4,/),'0HIGH SIDE EXPONENT PARAMETERS (NA, NB, NC)='
     1    ,/,3(1x,3g10.4,/),'0SPLIT PEARSON VII ASSYMETRY PARAMETER=',
     1    1x,3f10.6,/)
       else
c-----IF NOT A SPLIT PEARSON VII PROFILE
         write(6,105) (par(ip,i),sy(i),sz(i),i=17,19)
c105        format('0MIXING PARAMETERS (NA, NB, NC)=',/,3(1X,3G10.4,/))  !cp may 97
105        format('0MIXING PARAMETERS ',/,' NA= ',3g10.3,/
     *                                    ' NB= ',3g10.3,/,
     *                                    ' NC= ',3g10.3,/)
       end if
         write(6,107) (par(ip,i),sy(i),sz(i),i=3,5),
     1                    par(ip,21),sy(21),sz(21),
     1                    par(ip,20),sy(20),sz(20)
c107      format( '0FWHM PARAMETERS (U,V,W,Z,CT)=',/,5(1X,3F10.6,/))

107      format( '0FWHM PARAMETERS=',/,' U = ',3F10.6,/,
     *                                 ' V = ',3f10.6,/,
     *                                 ' W = ',3f10.6,/,
     *                                 ' CT= ',3f10.6,/,
     *                                 ' Z = ',3f10.6,/)

c-----Modification introduced by Carlos O. Paiva-Santos to perform
c-----Quantitative phase analysis, 03/94. Added 05/94, T.S. Moss   
c-----CHANGES TO INCORPORATE THE REFINED OCCUPANCY. Paiva-Santos (Feb-Mar/95)
          do  428 i=1,n 
          icoco=ptr(i+iof)
          tmassa(IP) = tmassa(IP) + xl(i+iof,5)*xmas(icoco)*xmltp(ip)
          stmassa(ip) = stmassa(ip) + smass(ip)*xmas(icoco)*xmltp(ip)   
428         continue
C      print '(a,f12.5, I2)', ' stmassa = ', stmassa(ip), ip  
      XFAC = 3.141592654 / 180.000000
      DCV(4) = XFAC * DCV(4)
      DCSM(4,4) = DCSM(4,4) * XFAC
      DCV(5) = XFAC * DCV(5)
      DCSM(5,5) = DCSM(5,5) * XFAC
      DCV(6) = DCV(6) * XFAC
      DCSM(6,6) = DCSM(6,6) * XFAC
c-----Calculations of VOLUME and SVZM (=W) for each phase
c-----and respectives standard deviations
c-----New standard deviation code introduced in nov 96 !cp
      ARGCOS= 1-(COS(DCV(4)))**2-(COS(DCV(5)))**2-(COS(DCV(6)))**2
     1 + 2 * (COS(DCV(4))) * (COS(DCV(5))) * (COS(DCV(6)))
      V0 = DCV(1) * DCV(2) * DCV(3)
      VOL(IP) = V0 * SQRT(ARGCOS)
      VOSQ = 0.5*vol(ip)/argcos
      ARG1 = VOSQ*(2 * COS(DCV(4)) * SIN(DCV(4)) -
     1   2*SIN(DCV(4)) *COS(DCV(5)) *COS(DCV(6))) *DCSM(4,4)
      ARG2 = VOSQ*(2 * COS(DCV(5)) * SIN(DCV(5)) -
     1   2*SIN(DCV(5)) *COS(DCV(4)) *COS(DCV(6))) *DCSM(5,5)
      ARG3 = VOSQ*(2 * COS(DCV(6)) * SIN(DCV(6)) -
     1   2*SIN(DCV(6)) *COS(DCV(4)) *COS(DCV(5))) *DCSM(6,6)
      DVOL(IP) = SQRT((VOL(IP) * DCSM(1,1) / DCV(1))**2 
     1 + (VOL(IP) * DCSM(2,2) / DCV(2))**2 +
     1 (VOL(IP) * DCSM(3,3) / DCV(3))**2 + ARG1**2 + ARG2**2 + ARG3**2)
C standard deviations are calculed below                      !cp nov 96
      W(IP) = PAR(IP,1) * tmassa(IP) * VOL(IP)/SAQF(IP)
                          dw(ip) =  (sz(1)/par(ip,1))     +
     *                              (dvol(ip)/vol(ip))    +
     *                              (stmassa(ip)/tmassa(ip))/saqf(ip)    
C   end of std
      WRITE(6,2705)VOL(IP), DVOL(IP), tmassa(IP),
     * 1.66113*tmassa(ip)/vol(IP)
2705  FORMAT('Volume= ', F9.3, '(+/-)', F7.3,' UCW= ',F7.2,
     *' U.C.Density = 'f7.3,' gr/cm^3',/,
     *'                     _____________________________')
10      CONTINUE
C ****** QUANTITATIVE ANALYSIS ***************
      WTOTAL = 0.000000
      DWTOTAL = 0.000000
      DO 2707 I = 1, NPHASE
      WTOTAL = WTOTAL + W(I)
C      DWTOTAL = DWTOTAL + DW(I) 
2707  CONTINUE
      DO 2708 I = 1, NPHASE
C      !cp nov 10 96
          if (nphase.eq.1)then
              XMASS(I) = 100. * W(I) / WTOTAL
              dmass(1)= 0.
          else
              XMASS(I) = 100. * W(I) / WTOTAL
              dmass(i) = 100*dw(i)
          end if
C Print comand below is for testing the codes. !cp 13june98
C      print '(a,f12.5,2x,a,f12.5)', ' xmass = ',xmass(i),'+/-',
C     ! dmass(i) 
2708  CONTINUE
      IINNMOL = 0
      DO 2713 I = 1, NPHASE
      IF (IINNMOL.EQ.1) GOTO 2713
      IF (NMOL(I).EQ.0) IINNMOL=1
2713    CONTINUE
      IF (IINNMOL.EQ.1) THEN
         DO 2714 I = 1, NPHASE
C ** printing results
         WRITE(6,2741) I, XMASS(I), DMASS(I)
2741     FORMAT(' PHASE = ',I2,' => %MASS = ',F6.2,'(+/-)',F6.2,2X,
     1     '%MOLAR = NOT COMPUTED')
2714     CONTINUE
      ELSE
C ****    CALCULATION OF MOLAR FRACTION  ****
       FT = 0.0000000
       DO 2709 I = 1, NPHASE
       FRP(I) = XMASS(I) * NMOL(I) / tmassa(I)
       FT = FT + FRP(I)
2709     CONTINUE
       DO 2710 I = 1, NPHASE
       FR(I) = 100. * FRP(I) / FT
C      dfr(I) = dmass(i)*nmol(i)/tmassa(i)
       dfr(I) = fr(I) * dmass(i)/xmass(i)
C ** printing results
       WRITE(6,2711) I, XMASS(I), DMASS(I), FR(I), dfr(I)  
2711     FORMAT(' PHASE = ',I2,' => %MASS = ',F6.2,'(+/-)',F6.2,2X,
     1     '%MOLAR = ',F6.2,'(+/-)',F6.2)
2710   CONTINUE
      ENDIF
      if(isphase.ne.0)then
       write(6,2717)
2717    format(/,' Considering Amorphous Content:')    
         sfIS=wtis(isphase)/xmass(isphase)
          if(sfIS.gt.1.)then
            write(6,2718)ISPHASE, wtis(isphase)
2718        format(' PROBLEM:Amount of Internal Standard (Phase #',i2,
     *             ') is less than the specified ' f6.2,'%.'/,
     *             ' Amorphous content not computed. Check ISWT in'
     *             ' line 11.2 for this phase')
            go to 2720
          end if
         do 2715 i=1,nphase
           write(6,2716)i,xmass(i)*sfIS  
2716       format(' PHASE = ',I2,' => %MASS = ',F6.2)
2715     continue   
       write(6,2719)100*(1.-sfIS)
2719       format(' AMORPHOUS  => %MASS = ',F6.2)       
      end if
2720  write(6,2712)
2712  format(1x,/)
c291   DO 19 J=1,13    !cp may 10 97
291   DO 19 J=1,20
        ILOC = ILOC +1
        KM=LGLB(J)
        IF(KM.NE.0)GOTO 20
        SY(J)=0.
        SZ(J)=0.
        FINAL(ILOC,1) = GLB(J)
        GOTO 19
20      SZ(J)=SQRT(ABS(SMM(KM,KM)))
        SY(J)=V(KM)*AGLB(J)*RELAX(4)
        DUMMY(KM) = GLB(J)
        GLB(J)=GLB(J)+SY(J)
        DUMMY(KM+MSZ)  = SY(J)
        FINAL(ILOC,1) = GLB(J)
        FINAL(ILOC,2) = SZ(J) 
19      CONTINUE
      if (ILOC .gt. nfinal) Stop
     +          ' Parameter NFINAL in PARAM.INC file too small'
      WRITE(6,109)GLB(1),SY(1),SZ(1)
109   FORMAT('0GLOBAL PARAMETERS',/,
     *       ' ZEROPOINT (ZER)             :',3F8.4)
      WRITE(6,181) (GLB(JJ),SY(JJ),SZ(JJ),JJ=10,11),
     *              (glb(kk),sy(kk),sz(kk),kk= 8, 9),
     *              (glb(ll),sy(ll),sz(ll),ll=12,12),
     *              (glb(mm),sy(mm),sz(mm),mm=13,13)
181   FORMAT(' SAMPLE DISPLACEMENT (DISP)  :',3F8.4,/,
     *       ' SAMPLE TRANSPARENCY (TRANSP):',3F8.4,/,
     *       ' ROUGHNESS PARAMETERS        :'      ,/,
     *       '              P              :',3F8.4,/,
     *       '              Q              :',3F8.4,/,
     *       '              R              :',3F8.4,/,
     *       '              T              :',3F8.4,/)
C !cp ap 20 97  !from It. codes
      WRITE(6,1884)GLB(20),SY(20),SZ(20)
1884  FORMAT( ' AMORPHOUS SCALE (SCAM):',3F11.4)
C !cp ap 20 97   !from It. codes
      WRITE(6,1882) (GLB(J),SY(J),SZ(J),J=18,19)
1882  FORMAT (/' MONOCROMATOR BANDPASS PARAMETERS (PMONI) ',/,6F8.4)
      IF(NBCKGD.NE.0)GOTO 90
      WRITE(6,84)(GLB(J),SY(J),SZ(J),J=2,7)
84    FORMAT(/,'0BACKGROUND PARAMETERS',/,(1X,3G12.6))
90    CONTINUE
      WRITE (8,ERR=99990) (DUMMY(I),I=1,MAXS)
     *   ,(DUMMY(I),I=MSZ+1,MSZ+MAXS),(DUMMY(I),I=2*MSZ+1,2*MSZ+4)
      RETURN
99990 STOP 'ERROR IN WRITING TO UNIT 8 IN OUTPTR'
      END
      FUNCTION GAMMA(X)
      DOUBLE PRECISION DU1,DU2,DU4
  201 FORMAT(/' X IN GAMMA(X) IS LESS THAN 0.1   X =',F8.4)
  202 FORMAT(/' X IN GAMMA(X) IS GREATER THAN 8.0   X =',F8.4)
      IF(X-1.0) 11,13,13
   11 IF(X-0.1) 97,12,12
   12 GX=X+1.0
      GO TO 21
   13 IF(X-2.0) 14,14,15
   14 GX=X
      GO TO 21
   15 IF(X-8.0) 16,16,98
   16 DO 17 I=1,8
      IN=I
      GX=X-FLOAT(I)
      IF(GX.GE.1.0.AND.GX.LE.2.0) GO TO 21
   17 CONTINUE
   21 DU1=GX-1.0
      DU2=DU1*DU1
      DU4=DU2*DU2
      GG=1.0D0-0.5771917D0*DU1+0.9882059D0*DU2-0.8970569D0*DU1*DU2
     1  +0.9182069D0*DU4-0.7567041D0*DU4*DU1+0.4821994D0*DU4*DU2
     2  -0.1935278D0*DU1*DU2*DU4+0.0358683D0*DU4*DU4
      IF(X-1.0) 31,32,32
   31 GAMMA=GG/X
      GO TO 99
   32 IF(X-2.0) 33,33,34
   33 GAMMA=GG
      GO TO 99
   34 GAMMA=X-1.0
      IF(IN-1) 35,35,36
   35 GAMMA=GAMMA*GG
      GO TO 99
   36 DO 37 I=2,IN
   37 GAMMA=GAMMA*(X-FLOAT(I))
      GAMMA=GAMMA*GG
      GO TO 99
   97 WRITE(*,201) X
      GAMMA=9.513508
      GO TO 99
   98 WRITE(6,202) X
      GAMMA=5040.0
   99 RETURN
      END 
C  SUBROUTINE ABSORP and SUBROUTINE ARIA: by Canton et all. Added by cps between
C  march-may 1997
      SUBROUTINE ABSORP(MU,SW,TH,ABC)
c-----THIS SUBROUTINE CORRECTS THE EXPERIMENTAL INTENSITIES FOR THE
C                ABSORPTION EFFECTS AS REPORTED BY:
C     1) H. P. KLUG & L. E. ALEXANDER, X-RAY DIFFRACTION PROCEDURES,
C        1970, PAG.487.
C     2) A. IMMIRZI, ACTA CRYST. , 1980, B36, 2378-2385.
c
C     IT IS WRITTEN TAKING INTO ACCOUNT THE SYMMETRIC REFLECTION ARRANGEMENT
C     ( KLUG & ALEXANDER, 1970, FIG. 5-52 PAG. 390) AND  THE  FINITE
C     THICKNESS OR WIDTH OF THE SLAB, IN THIS SITUATION THE ABSORPTION
C     IS GENERALLY LOW AND INCREASES SLIGHTLY WITH 2 THETA.
C     MU = LINEAR ABSORPTION COEFFICIENT IN CM-1.
C     SW = SAMPLE THICKNESS IN CM.
      REAL MU
      EX  = ( 2.0 * MU * SW ) / SIN (TH * 0.008726646)
      ABC = 1.0 - EXP (-EX)
      RETURN
      END
cccccc subroutine inpam: To read the amorphous data file !cp may 10 97
      subroutine inpam
      include 'param.inc'
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
C unit 11 = amorphous file
        OPEN(11, FILE=' ', status='unknown')
c-----CONTROL IF THERE IS THE AMORPHOUS FILE .
c-----AND IF IT IS ON THE SAME POINTS OF THE DATA FILE
      READ(11,8500,END=99999)THMIN1,STEP1,THMAX1,TMV1,SW1,DATAID1
8500    FORMAT(BZ,5F8.3,A16)
      WRITE (6,79677) DATAID1
79677  FORMAT (5X,'DATA AMORPHOUS ',A16)
       IF(THMIN1.NE.THMIN)THEN
       WRITE (6,150)THMIN1,THMIN
150    FORMAT(1X,'AMORPHOUS THMIN=',F8.2,
     * 'IS DIFFERENT FROM THE DATA THMIN=',F8.2)
       END IF
       IF(STEP1.NE.STEP)THEN
       WRITE (6,1514)STEP1,STEP
1514    FORMAT(1X,'AMORPHOUS STEP=',F8.2,
     * 'IS DIFFERENT FROM THE DATA STEP=',F8.2)
       END IF
       IF(THMAX1.NE.THMAX)THEN
       WRITE (6,152)THMAX1,THMAX
152    FORMAT(1X,'AMORPHOUS THMAX=',F8.2,
     * 'IS DIFFERENT FROM THE DATA THMAX=',F8.2)
       END IF
C      READ(11,4481,END=99999)(AMORPHOUS(I),I=1,NPTS)
c4481   FORMAT(BZ,8(F7.0,1X))
C      read the rest of the file in free format
      READ(11,*,END=99999)(AMORPHOUS(I),I=1,NPTS)
       DO 448 I = 1,NPTS
       TH=THMIN+(I-1)*STEP
c------AMORPHOUS CORRECTION FOR AIR SCATTERING
C       CALL ARIA(TMV1,SW1,FI1,TH,SCA)
C       AMORPHOUS(I)=0.0
      IF(IAS.EQ.1) THEN
      CALL ABSORP(TMV1,SW1,TH,ABC)
      AMORPHOUS(I) = AMORPHOUS(I) / ABC
      ENDIF
504    FORMAT(F8.0)
448    CONTINUE
      return               
99999 write(6,806) 11
  806 format(1x,' END OF FILE ENCOUNTERED AT START ON TAPE =',I2,'OF
     $PHASE AMORPHUS')
      STOP   '   *** END OF FILE ENCOUNTERED FOR AMORPHOUS FILE *** '
      end    
C subroutine convert from dbws9411 to dbws9807
      SUBROUTINE conv94
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 )
      INTEGER PTR,fondo
      REAL LAMDA
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/REFLS/IREFS(IRS),REFS(IRS,3),FMGNTD(IRS),ICR(99)
     *  ,HALFL(IRS),HALFG(IRS),GAM(IRS),fwhm(irs,2)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      CHARACTER SYMB(99,20)*1,SPG(20)*1
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      CHARACTER TITLE*70,PHSNM(99)*50
C      CHARACTER*70 TITLE
C      character*50 PHSNM(99)
C      CHARACTER SYMB(99,20)*1
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/BLNK1/NJNK(22),FILL2(751)
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),xmas(16)
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/CELLX/AA,B,C,ALPHA,BETA,GAMMA,AL(3,3)
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      COMMON/MULTIP/TMASSA(99),MLTPHASE,xmltp(99),MURT(NATS),SAQF(99)  
     *,wtis(99)
      COMMON /SPGCOM/ NSPGRP,NAXIS,NCONT(10),NC,NPOL,MULTP
      common/codebck/ibckcode     
      common/simoper/isimop      
      DIMENSION XRYZ(10),BACK(6),FBACK(6),MLTT(IRS)
      EQUIVALENCE (MLTT(1),FMGNTD(1)),
     1            (GLB(2),BACK(1)),
     1            (AGLB(2),FBACK(1))
      common/convert/icnvt
      CHARACTER*26 LOWER,UPPER
      DATA UPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LOWER/'abcdefghijklmnopqrstuvwxyz'/
C      CHARACTER*80 OUTFILE
      DATA XRYZ(1), XRYZ(2),  XRYZ(3),  XRYZ(4),  XRYZ(5)
     *   /2.748510,2.289620, 1.935970, 1.788965, 1.540520/
      DATA XRYZ(6), XRYZ(7),  XRYZ(8),  XRYZ(9), XRYZ(10)
     * /  0.709260,0.559360, 0.215947, 0.209010, 0.180195/
C      DATA XRYZ(1),XRYZ(2),XRYZ(3),XRYZ(4),XRYZ(5)
C     * /2.28962,1.93597,1.54051,0.70926,0.556363/ 
      xnothing = xryz(1)
C      icnvt=0        
C      INQUIRE(UNIT=6,NAME=OUTFILE)
C      LBKSL = 0
C      DO I=1,80
C        IF ( OUTFILE(I:I).EQ.'\' ) LBKSL=I
C      END DO
c!     WRITE(6,*) 'Outfile:',OUTFILE
C      OUTFILE(LBKSL+1:) = 'newicf.inp'
c!     WRITE(6,*) 'Outfile:',OUTFILE
C      open(51, file=OUTFILE(1:LBKSL+10), status='unknown')
c!     open(51, file='newicf.inp', status='unknown')
      rewind 5
      write(6,642)
642   format('*** New ICF starts here with next line ',
     !       '*** Delete everything above it')
C line 1
      READ(5,1,END=99999)TITLE
1     FORMAT(BZ,A70)
      write(6,801)TITLE
801   FORMAT(A70)
      WRITE(*,1001)TITLE
1001  FORMAT(///,1X,A70)
C line 2.1
      READ(5,2,END=99999)JOBTYP,NPROF,NPHASE,NBCKGD,NEXCRG,NSCAT,INSTRM
     *                   ,IPREF,IABSR
2     FORMAT(BZ,12I4)
      iasym = 0
      nprof = abs(nprof)
      idata=0     
      isphase=0
      i2d94 =0
      write(6,802)JOBTYP,NPROF,NPHASE,nbckgd,NEXCRG,NSCAT,INSTRM,IPREF,
     1          iasym,iabsr,idata,isphase,i2d94
802    FORMAT(13I4,5x,'LINE 2.1')
C line 3
      READ(5,14,END=99999)IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2
     * ,IPLST
14    FORMAT(BZ,11I1)
      iploss = 0
      iplcl  = 0
      iplcal = 0
      iplpol = 0
      iplcom = 0
      ipldis = 0
      iplam  = 0
      ipbig  = 0
      write(6,803)IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2,IPLST,
     *IPLOSS,IPLCAL,IPLPOL,IPLCOM,IPLDIS,IPLAM,IPBIG
803    FORMAT(3(5i1,1x),3i1,36x,'LINE 3')
c
      READ(5,23,END=99999)LAMDA,RATIO(2),BKPOS,WDT,CTHM,TMV,RLIM
23    FORMAT(BZ,8F8.0)
      sw = 0.
      write(6,804)LAMDA,RATIO(2),BKPOS,WDT,CTHM,TMV,RLIM,sw
804   FORMAT(3f8.5,6F8.4)
C  line 5
      READ(5,30,END=99999)MCYCLE,EPS,RELAX,THMIN,STEP,THMAX 
30    FORMAT(BZ,1I4,5F4.0,3F8.0)
C 
      write(6,805)MCYCLE,EPS,RELAX,THMIN,STEP,THMAX
805   FORMAT(1I4,5F4.2,3f8.3,9x,'CYCLS EPS RELAX P_CALC')
c
      IF(NBCKGD.LT.2)GOTO 120 
C line 6(*)
      READ(5,34,END=99999)(POS(I),BCK(I),I=1,NBCKGD)
34    FORMAT(BZ,2F8.2)
      write(6,806)(POS(I),BCK(I),I=1,NBCKGD)
806   FORMAT(2F8.2)   
c
C line 7(*)
120   IF(NEXCRG.LE.0)GOTO 122 
      READ(5,34,END=99999)(ALOW(I),AHIGH(I),I=1,NEXCRG)
      write(6,807)(ALOW(I),AHIGH(I),I=1,NEXCRG)
807     format(2f8.2,41x,'EXCLUDED REGION')
C  line 8(*) read
122   IF(NSCAT.LE.0)GOTO 124
      DO 125 I=1,NSCAT
         IF(JOBTYP.EQ.1.OR.JOBTYP.EQ.3) THEN
           READ(5,3838,END=99999)NAM(I),DFP(I),XMAS(I) 
3838       FORMAT(BZ,A4,2F8.0)
           GOTO 125
         ENDIF
        READ(5,38,END=99999)NAM(I),DFP(I),DFPP(I),XMAS(I) 
38      FORMAT(BZ,A4,3F8.0)
        K=0
126     READ(5,39,END=99999)(AC(J,I),J=1,9)
39      FORMAT(BZ,9F8.0)
        IF(AC(1,I).EQ.-100.)CALL COEF(I,K)
        IF(AC(3,I).NE.0.)GOTO 125
        K=K+1
        POSI(K)=AC(1,I)
        SCAT(K)=AC(2,I)
        IF(K.LE.29)GOTO126
        WRITE(6,40) 
40      FORMAT(34H TOO MANY SCATTERING TABLE ENTRIES )
        STOP 7700
125     CONTINUE
C     DO 424 I=1,NSCAT
C       DO 425 J=1,9
c425       AC(J,I)=0.0
c424     DFPP(I)=0.0
      GOTO 124
124   CONTINUE
C  line 8 (write)
      IF(NSCAT.LE.0)GOTO 813
      DO 808 I=1,NSCAT
           IF (JOBTYP .EQ. 1.OR.JOBTYP.EQ.3) GOTO 809
C line 8.1 XRD (*)
        write(6,810)NAM(I),DFP(I),DFPP(I),XMAS(I),NSCAT
810        FORMAT(A4,3F8.4,29X,'SCATTERING SET ',I2)
        GOTO 814
C line 8.1 ND(*)
809        write(6,811)NAM(I),DFP(I),XMAS(I),NSCAT
811        FORMAT(A4,2F8.4,37X,'SCATTERING SET ',I2)
C line 8.2 XRD(*)
814       IF(JOBTYP.EQ.0.OR.JOBTYP.EQ.2)write(6,812)(AC(J,I),J=1,9) 
812          FORMAT(9F8.5)
808   CONTINUE
813   CONTINUE
C line 9
      READ(5,43,END=99999)MAXS
43    FORMAT(BZ,I8) 
      write(6,815)MAXS
815   FORMAT(I8,49X,'PARAMS REFINED')
C  line 10 (read)
      READ(5,45,END=99999)glb(1),glb(10),glb(11),glb(8),glb(9),glb(12),
     *                    glb(13) 
      READ(5,45,END=99999)aglb(1),aglb(10),aglb(11),aglb(8),aglb(9),
     1                    aglb(12),aglb(13)
45    FORMAT(BZ,10F8.0)
C line 10.1 (write)
      write(6,816)GLB(1),GLB(10),GLB(11),glb(8),glb(9),glb(12),glb(13)
C line 10.11 (write)
      write(6,817)AGLB(1),AGLB(10),AGLB(11),aglb(8),aglb(9),aglb(12),
     *            aglb(13)
816   FORMAT(7F8.4,1x,'ZER DISP TRANS p q r t')
817   format(7f8.4,1X,'CODEWORDS')
C line 10.3 (read & write)
      IF(NBCKGD.NE.0)GOTO 49
      READ(5,51,END=99999)BACK,FBACK
51    FORMAT(BZ,6F9.4)
      write(6,818)(GLB(J),J=2,7),(AGLB(J),J=2,7)
818   FORMAT(6f9.2,3x,'BACKGROUND',/,6F9.4,3X,'CODEWORDS') 
C read AND rewrite for phases
c
49      DO 81 K=1,NPHASE
C line 11.1
        READ(5,465,END=99999)PHSNM(K)
465     FORMAT(BZ,A50)
c465     FORMAT(BZ,A57)
        write(6,819)PHSNM(K),K
819     FORMAT(A50,7X,'PHASE NUMBER ',I2)
c
C line 11.2
        READ(5,48,END=99999)NATOM(K),NMOL(K),
     1            (PREF(K,I),I=1,3)
48      FORMAT(BZ,2I4,8x,3F4.0)
      saqf(k) = 1.
      write(6,820)NATOM(K),NMOL(K), saqf(k),(PREF(K,I),I=1,3),wtis(K)
820   FORMAT(2I4,f7.4,1x,3F4.1,F7.2,22X,'#ATMS #FU AFQPA PREFDIR ISWT')
c
      N=NATOM(K)
C line 11.3
      READ(5,460,END=99999)(SYMB(K,I),I=1,20)
460   FORMAT(BZ,20A1)
      write(6,821)(SYMB(K,I),I=1,20)
821   FORMAT(20A1,37X,'SPACE GROUP')
      DO 135 I=1,20
          SPG(I)=SYMB(K,I)
C convert lower case to upper case
          DO 22371 ik=1,26
          IF (spg(i).EQ.LOWER(ik:ik)) spg(i)=UPPER(ik:ik)
22371     CONTINUE
c finish conversion
135      continue
      CALL SPGP(SPG)
C getting multiplicity of each phase !cp jun 96)
      isimop=1
       call rtmt(ipl1,k)
       xmltp(k)=MLTPHASE 
       WRITE(*,33)MLTPHASE, K
C       WRITE(6,33)MLTPHASE, K
33    format(' General position multiplicity is ', i4,' for phase ',i2)
C line 11.4
c-----READ  FOR EACH ATOM
        IOF=0
        IF(K.GT.1) THEN
          DO 4325 IIPHAS=2,K
4325        IOF = IOF + NATOM(IIPHAS-1) 
        END IF
C line 11.41 and 11.42
        READ(5,65,END=99999)(ATEXT(I+IOF),NTYP(I+IOF),
     *      (XL(I+IOF,J),J=1, 5), (A(I+IOF,J),J=1, 5),
     *         (XL(I+IOF,J),J=6,11), (A(I+IOF,J),J=6,11),
     *               I=1,N)
65      FORMAT(BZ,2A4,8X,5F8.5,/,16X,5F8.2,/,6F8.5,/,6F8.2) 
C Creating MURT and 'So' 
      do 3947 isof=1,n
C         murt(isof+iof) = ifix(xl(isof+iof,5)*xmltp(k))
         xl(isof+iof,5)=1.
3947  continue
C      write(6,822)(ATEXT(I+IOF),murt(i+iof),NTYP(I+IOF),
      write(6,822)(ATEXT(I+IOF),NTYP(I+IOF),
     *     (XL(I+IOF,J),J=1, 5),(A(I+IOF,J),J=1, 5),
     *        (XL(I+IOF,J),J=6,11),(A(I+IOF,J),J=6,11),
     *            I=1,N)
822   FORMAT(A4,1x,'   #',1x,a4,2X,5F8.5,2x,'LBL M NTYP x y z B So'/,
     *      16X,5F8.2,2x,'CODEWORDS'/,6F8.5,10X,'BETAS',/,
     *      6F8.2,10x,'CODEWORDS')
c
        READ(5,50,END=99999) PAR(K,1),PAR(K,2),
     *                          APAR(K,1),APAR(K,2),
     * PAR(K,3),PAR(K,4),PAR(K,5),PAR(K,20),PAR(K,15),PAR(K,16),
     * APAR(K,3),APAR(K,4),APAR(K,5),APAR(K,20),APAR(K,15),APAR(K,16),
     *                  (PAR(K,I),I=6,11), (APAR(K,I),I=6,11), 
     *                  PAR(K,12),PAR(K,13),PAR(K,14),
     *                     APAR(K,12),APAR(K,13),APAR(K,14),
     *                  PAR(K,17),PAR(K,18),PAR(K,19),
     *                     APAR(K,17),APAR(K,18),APAR(K,19),
     *                  PAR(K,24),PAR(K,25),PAR(K,26),
     *                     APAR(K,24),APAR(K,25),APAR(K,26),
     *                  par(k,27),
     *                     apar(k,27)
50      FORMAT(BZ,G8.2,F8.0,/,2F8.2,/6F8.0,/,6F8.2,/,
     *         6F8.0,/,6F8.2,/,3F8.0,/,3F8.2,/,3F8.4,/,3F8.2,/,
     *         3F8.4,/,3f8.2,/,f8.4,/,f8.2)
      par(k,21)  = 0.0
      apar(k,21) = 0.0
      write(6,823) PAR(K,1),PAR(K,2),
     *                 APAR(K,1),APAR(K,2),
     *PAR(K,3),PAR(K,4),PAR(K,5),PAR(K,21),par(k,20),
     *         PAR(K,15),PAR(K,16),
     *APAR(K,3),APAR(K,4),APAR(K,5),APAR(K,21),apar(k,20),
     *         APAR(K,15),APAR(K,16),
     *              (PAR(K,I),I=6,11), (APAR(K,I),I=6,11), 
     *              PAR(K,12),PAR(K,13),PAR(K,14),
     *                 APAR(K,12),APAR(K,13),APAR(K,14),
     *              PAR(K,17),PAR(K,18),PAR(K,19),
     *                 APAR(K,17),APAR(K,18),APAR(K,19),
     1              par(k,24),par(k,25),par(k,26),
     1                 apar(k,24),apar(k,25),apar(k,26),
     1              par(k,27),
     1                 apar(k,27)
823     FORMAT(G8.3,F8.4,41X,'SCALE Bo(OVERALL)',/,2F8.2,/,
     *        7F8.5,1X,'U V W CT Z X Y',/,7F8.2,/,
     *        6F8.4,9X,'CELL PARAMETERS',/,6F8.2,/,
     *        3F8.5,33X,'PREF1 PREF2 R/RCF_ASYM',/,3F8.2,/,
     *        3F8.4,33X,'NA NB NC (MIX_PARAMS)',/,3F8.2,/,
     1        3F8.4,33X,'NA NB NC (HIGH SIDE)',/,3F8.2,/,
     1         f8.4,49x,'PEARSON ASYM.FACTOR',/,F8.2)
81      CONTINUE
      RETURN
99999 STOP 'END OF FILE TAPE5'
      END
C SUBROUTINE TO READ PHILIPS UDF DATA FILE
      SUBROUTINE PHILIPSREAD
      include 'param.inc'
      character*56 test,dataid
      character*10 athmin,athmax
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT 
     * ,ULOR,VLOR,ZZZ,UC
      open (96,file='trash',status='unknown')
      read(4,'(a56)')dataid
10      read(4,1)test
1     format(a56)
      do 3 i=1 ,56
        if(test(I:I).eq.',')then
             lb1=i
             goto 31
        end if
3      continue
31      if(test(1:lb1-1).ne.'DataAngleRange')goto 10
C                            DataAngleRange
      do 4 i=lb1+1,56
        if(test(I:I).eq.',')then
          lb2=i
          goto 41
        end if
4     continue
c
41    do 5 i=lb2+1,56
        if(test(I:I).eq.',')then
          lb3=i
          goto 51
        end if
5     continue
C      
51    athmin=test(lb1+1:lb2-1)
      athmax=test(lb2+1:lb3-1)
C      write(*,11)athmin,athmax
c11    format(/' thmax= ',a10,'  thmin= ',a10)
      read(4,2)test,step
2     format(A13,f8.3)
C      write(*,12)step
c12    format(\ ' Step= ',f8.3)      
      write(96,110)athmin,athmax,step
110   format(2(a10,1x),f8.3)
30     read(4,1)test
      do i=1 ,56
      if(test(I:I).eq.',')lb1=i
      end do
C         
      if(test(1:(lb1-1)).ne.'RawScan')goto 30
C                            RawScan
      rewind 96
      read(96,*)thmin,thmax,step

      write(*,4711)thmin,thmax,step
4711  format(1x,'DATA RANGE (2THETA):  START =',F8.3,', STOP =',
     *      F8.3,', STEP =',F8.3,/)
      WRITE (6,7967) DATAID
7967  FORMAT (5X,'DATA ID ',A56)
      write(6,4711)thmin,thmax,step
      NPTS=(THMAX-THMIN)/STEP+1.5
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      WRITE (*,1000) IDSZ,NPTS
1000    FORMAT(1X,'PROGRAM CAN HANDLE ',I5,'POINTS',/,
     *      1X,I5,' POINTS WERE INPUT',/,
     *        ' INCREASE IDSZ IN PARAMETER STATEMENT')
      STOP ' TOO MANY DATA POINTS '
      ENDIF

c
      read(4,*,END=99998)(Y(I),I=1,NPTS)
C        
      close (96,status='DELETE') 
      return
99998 STOP 'END OF FILE unit=4'
C      close (96,status='DELETE') 
      end
C SUBROUTINE GSASREAD * READ GSAS FORMATTED DATA FILE
      subroutine gsasread     
      include 'param.inc'
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT 
     * ,ULOR,VLOR,ZZZ,UC
      character*72 test,dataid
      dimension lab(20)
      dimension nctr(IDSZ)
      open (96,file='trash',status='unknown') 
C      open (12,file='fromGsas.rit')
      read(4,'(a72)')dataid
C      write(*,'(3x,a72)')dataid
      j=0
      k=1
10    read(4,1)test
1     format(a72)
      if(test(1:4).ne.'BANK')goto 10
31    do 2 i=k ,72
        if(test(i:i).eq.' ')then
              j=j+1
              lab(j)=i-1
              goto 21
        end if
2      continue
c
21    do 3 k=i,66
        if(test(k:k).ne.' ')then
             j=j+1
             lab(j)=k
             goto 31
        end if
3      continue
      l4=lab(4)
      l5=lab(5)
      l10=lab(10)
      l11=lab(11)
      l12=lab(12)
      l13=lab(13)
C      write(*,110)test(l4:l5),test(l10:l11),test(l12:l13)
      write(96,110)test(l4:l5),test(l10:l11),test(l12:l13)
C                   npts        start*100     step*100
110   format(3(a10,1x))
      REWIND 96
      read(96,*)npts,thmin,step
      thmin=thmin/100
      step=step/100
      thmax =thmin+npts*step
      write(*,4711)thmin,thmax,step
4711  format(1x,'DATA RANGE (2THETA):  START =',F8.3,', STOP =',
     *      F8.3,', STEP =',F8.3,/)
      WRITE (6,7967) DATAID
7967  FORMAT (5X,'DATA ID ',A56)
      write(6,4711)thmin,thmax,step
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      WRITE (*,1000) IDSZ,NPTS
1000    FORMAT(1X,'PROGRAM CAN HANDLE ',I5,'POINTS',/,
     *      1X,I5,' POINTS WERE INPUT',/,
     *        ' INCREASE IDSZ IN PARAMETER STATEMENT')
      STOP ' TOO MANY DATA POINTS '
      ENDIF
c
      read(4,411,END=99998)(nctr(I),Y(I),I=1,NPTS)
411   format(10(i2,f6.0))
      do 412 i=1,npts
      if (nctr(i).eq.0)nctr(i)=1
      if(y(i).eq.0.)y(i)=0.
C      y(i)=nctr(i)*y(i)
412   continue      
C         write(12,121)thmin,step,thmax,dataid,(y(i),i=1,npts)
c121      format(3f8.3,1x,a56,/,8(f7.0,1x))
      close (96,status='DELETE') 
C      close (96) 
      RETURN             
99998 STOP 'END OF FILE unit=4'
      END            
C subroutine to read SCINTAG TXT file      
      subroutine scintag 
      include 'param.inc'
      character*56 test,dataid
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT 
     * ,ULOR,VLOR,ZZZ,UC
      character*15 athmin,athmax,astep
      common/labels/lb1,lb2,lb3,test
      open (96,file='trash',status='unknown')
C     open (97,file='trash.DAT',status='unknown')
      read(4,'(a56)')dataid
C      write(*,'(3x,a72)')dataid
1     format(a56)  
96    format(1x,a15)
10    read(4,1)test
        if(test(1:12).ne.'Start Angle:') goto 10
C                         Start Angle:      
            call readasc
      athmin=test(lb2:lb3)
20    read(4,1)test
        if(test(1:11).ne.'Stop Angle:')goto 20
C                         Stop Angle:      
            call readasc
      athmax=test(lb2:lb3)
30    read(4,1)test
        if(test(1:10).ne.'Step Size:')goto 30
C                         Step Size:      
            call readasc
      astep=test(lb2:lb3)
C      write(*,11)athmin,athmax,astep
c11    format(/' thmin= ',a15,' thmax= ',a15,' step= ',a15)
      write(96,110)athmin,athmax,astep
110   format(3(a15,1x))
C         
70    read(4,1)test
      if(test(1:5).ne.'Range')goto 70
      rewind 96
      read(96,*)thmin,thmax,step
c
      NPTS=(THMAX-THMIN)/STEP+1.5
C      write(*,51)thmin,thmax,step,npts
c51    format(3(f8.3,2x),i5)
c
      write(*,4711)thmin,thmax,step
4711  format(1x,'DATA RANGE (2THETA):  START =',F8.3,', STOP =',
     *      F8.3,', STEP =',F8.3,/)
      WRITE (6,7967) DATAID
7967  FORMAT (5X,'DATA ID ',A56)
      write(6,4711)thmin,thmax,step
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      WRITE (*,1000) IDSZ,NPTS
1000    FORMAT(1X,'PROGRAM CAN HANDLE ',I5,'POINTS',/,
     *      1X,I5,' POINTS WERE INPUT',/,
     *        ' INCREASE IDSZ IN PARAMETER STATEMENT')
      STOP ' TOO MANY DATA POINTS '
      ENDIF
c
      read(4,*,END=99998)(xxx,Y(I),xxx,xxx,I=1,NPTS)
C     
C     CHECK FOR NON-ZERO COUNTS AT ALL POINTS
c6     FORMAT(8(F8.0,1X))
C      write(97,62)thmin,step,thmax,dataid
c62    format(3f8.3,a56)      
C     WRITE(97,61)(Y(I),I=1,NPTS)
c1    format(8(f7.0,1x))            
C        
      close (96,status='DELETE') 
C     close (96) 
      RETURN
99998 STOP 'END OF FILE unit=4'
      end        
      subroutine readasc
      character*56 test
      common/labels/lb1,lb2,lb3,test
      do 3 i=1 ,56
        if(test(I:I).eq.':')then
             lb1=i
             goto 31
        end if
3      continue
c
31    do 4 i=lb1+1,56
        if(test(I:I).ne.' ')then
          lb2=i
          goto 41
        end if
4     continue
c
41    do 5 i=lb2,56
        if(test(I:I).eq.' ')then
          lb3=i-1
          goto 51
        end if
5     continue
51      return
!      stop
      end                                                                 
C SUBROUTINE TO READ SIEMENS UXD DATA FILE
      SUBROUTINE SIEMENSREAD
      include 'param.inc'
      character*72 test,dataid
      character*15 athmin,astep,anpts
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT 
     * ,ULOR,VLOR,ZZZ,UC
      open (96,file='trash',status='unknown')
      read(4,'(a72)')dataid
10    read(4,1,end=99999)test
1     format(a72)
      if(test(2:9).eq.'STEPSIZE') then
           astep=test(12:20)
C           write(*,'(1x,a15)') astep
           goto 10  
      else if(test(2:7).eq.'2THETA') then
           athmin=test(10:16) 
C           write(*,'(1x,a15)') athmin
           goto 10
      else if(test(2:10).eq.'STEPCOUNT') then
           anpts=test(13:19)
C           write(*,'(1x,a15)') anpts
           goto 10
      else if(test(2:7).eq.'COUNTS') then
          goto 20
      else
          goto 10
      end if
C      write(*,11)athmin,athmax
c11    format(/' thmax= ',a10,'  thmin= ',a10)
20    write(96,110)athmin,astep,anpts
110   format(1x,3(a15,1x))
      rewind 96
      read(96,*)thmin,step,npts
      thmax=thmin + step*npts
      write(*,4711)thmin,thmax,step
4711  format(1x,'DATA RANGE (2THETA):  START =',F8.3,', STOP =',
     *      F8.3,', STEP =',F8.3,/)
      WRITE (6,7967) DATAID
7967  FORMAT (5X,'DATA ID ',A72)
      write(6,4711)thmin,thmax,step
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      WRITE (*,1000) IDSZ,NPTS
1000    FORMAT(1X,'PROGRAM CAN HANDLE ',I5,'POINTS',/,
     *      1X,I5,' POINTS WERE INPUT',/,
     *        ' INCREASE IDSZ IN PARAMETER STATEMENT')
      STOP ' TOO MANY DATA POINTS '
      ENDIF

c
      read(4,*,END=99998)(Y(I),I=1,NPTS)
      close (96,status='DELETE') 
C      close (96)
      return
99998 STOP 'END OF FILE unit=4'
99999 stop 'IS THE FILE NOT UXD SIEMENS FORMAT?'
      end
C SUBROUTINE TO READ RIGAKU DATA FILE
      SUBROUTINE rigakuread
      include 'param.inc'
      character*72 test
      character*15 athmin,astep,athmax,anpts
      COMMON/DATAX/Y(IDSZ),VAR(IDSZ),YC(IDSZ),KR(IDSZ)
     * ,BK(IDSZ),NPTS,AMORPHOUS(IDSZ)
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT 
     * ,ULOR,VLOR,ZZZ,UC
      open (96,file='trash',status='unknown')
C      read(4,'(a72)')dataid
10    read(4,1,end=99999)test
1     format(a)
      if(test(1:6).eq.'*START') then
        do 30 i=1 ,56
          if(test(I:I).eq.'=') then
             athmin=test(i+1:i+11)
C             write(*,'(1x,a15)') athmin
             goto 10  
           end if
30       continue
      end if
      if(test(1:5).eq.'*STOP') then
        do 31 i=1 ,56
          if(test(I:I).eq.'=') then
             athmax=test(i+1:I+11) 
C             write(*,'(1x,a15)') athmax
             goto 10         
           end if
31        continue
      end if
      if(test(1:5).eq.'*STEP') then
        do 32 i=1 ,56
          if(test(I:I).eq.'=') then
             astep=test(i+1:I+8) 
C             write(*,'(1x,a15)') astep
             goto 10         
           end if
32        continue
      end if
      if(test(1:6).eq.'*COUNT'.and.test(1:7).ne.'*COUNTE') then 
        do 33 i=1 ,56
          if(test(I:I).eq.'=') then
             anpts=test(i+1:I+8) 
C             write(*,'(1x,a15)') anpts
             goto 20         
           end if
33        continue
      end if
          goto 10
c20      write(*,11)athmin,astep,athmax,anpts
c11    format(/' thmax= ',a10,' step= ',a10,'  thmin= ',a10,'npts= ',a10)
20    write(96,110)athmin,astep,athmax,anpts
110   format(1x,4(a15,1x))
      rewind 96
      read(96,*)thmin,step,thmax,npts
      write(*,4711)thmin,thmax,step
4711  format(1x,'DATA RANGE (2THETA):  START =',F8.3,', STOP =',
     *      F8.3,', STEP =',F8.3,/)
      WRITE (6,7967) 
7967  FORMAT (5X,'DATA ID: RIGAKU DATA FILE ')
      write(6,4711)thmin,thmax,step
      IF (NPTS.GT.IDSZ) THEN
      WRITE(6,1000) IDSZ,NPTS
      WRITE (7,1000) IDSZ,NPTS
      WRITE (*,1000) IDSZ,NPTS
1000    FORMAT(1X,'PROGRAM CAN HANDLE ',I5,'POINTS',/,
     *      1X,I5,' POINTS WERE INPUT',/,
     *        ' INCREASE IDSZ IN PARAMETER STATEMENT')
      STOP ' TOO MANY DATA POINTS '
      ENDIF
C read the data
      read(4,*,END=99998)(Y(I),I=1,NPTS)
      close (96,status='DELETE') 
C     close (96)
C     open (97,file='unit97.dat',status='unknown')
C        write(97,121)thmin,step,thmax,(y(i),i=1,npts)
c21      format(3f8.3,'RIGAKU DATA FILE',/,8(f7.0,1x))
C     close (97)
      return
99998 STOP 'END OF FILE unit=4'
99999 stop 'IS THE FILE NOT RIGAKU FORMATED?'
      end
c
C suboutine (QPAINIT) to compute mass fractions before starting the refinement
c
      subroutine qpainit
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      INTEGER PTR,fondo
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      COMMON/ALLP/FINAL(nfinal,2),ILOC
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/G3/COND,IORD1,IORD2,TH,NUM 
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,PLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94
      COMMON/G2/S1,S2,ss2,S3,ss4,D1,D2,D4,R1,R2,R3,r2nobk,r3nobk
      COMMON/F1/SMM(MSZ,MSZ),V(MSZ)
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),xmas(16)
      COMMON/MULTIP/TMASSA(99),MLTPHASE,xmltp(99),MURT(NATS),SAQF(99)
     *,wtis(99)
      common/dirCV/ DCSM(6,6),DCV(6)
      DIMENSION VOL(99), W(99), XMASS(99)
      DIMENSION FR(99), FRP(99)  
      COMMON/VOLUME/VOLi(99),GCOM(99)
      COMMON/BKGSCALE/SCABKG(99)
      write(6,6996)
6996  format(/'        >>> QPA before starting the refinement <<<',/,
     *        '        >>                                    <<<')      
        DO 10 IP=1,NPHASE
      CALL DIRECT(DCSM,DCV,IP)
          tmassa(ip)=0.   
          IOF=0
          IF(IP.GT.1) THEN
            DO 44340 IIPHAS=2,IP
44340         IOF = IOF + NATOM(IIPHAS-1)
          END IF
          N=NATOM(IP)
          do  428 i=1,n 
          icoco=ptr(i+iof)
          tmassa(IP) = tmassa(IP) + xl(i+iof,5)*xmas(icoco)*xmltp(ip)
C      print '(a,f12.5)','xl(i+iof)=', xl(i+iof,5)
428         continue
C      print '(a,f12.5, I2)', ' tmassa = ', tmassa(ip), ip  
      XFAC = 3.141592654 / 180.000000
      DCV(4) = XFAC * DCV(4)
      DCSM(4,4) = DCSM(4,4) * XFAC
      DCV(5) = XFAC * DCV(5)
      DCSM(5,5) = DCSM(5,5) * XFAC
      DCV(6) = DCV(6) * XFAC
      DCSM(6,6) = DCSM(6,6) * XFAC
c-----Calculations of VOLUME and SVZM (=W) for each phase
      ARGCOS= 1-(COS(DCV(4)))**2-(COS(DCV(5)))**2-(COS(DCV(6)))**2
     1 + 2 * (COS(DCV(4))) * (COS(DCV(5))) * (COS(DCV(6)))
      V0 = DCV(1) * DCV(2) * DCV(3)
      VOL(IP) = V0 * SQRT(ARGCOS)
C      print '(a,f12.5, I2)', ' tmassa = ', tmassa(ip), ip  
C      print '(a,f12.5)', ' argcos = ', argcos  
      VOSQ = 0.5*vol(ip)/argcos
      ARG1 = VOSQ*(2 * COS(DCV(4)) * SIN(DCV(4)) -
     1   2*SIN(DCV(4)) *COS(DCV(5)) *COS(DCV(6))) *DCSM(4,4)
      ARG2 = VOSQ*(2 * COS(DCV(5)) * SIN(DCV(5)) -
     1   2*SIN(DCV(5)) *COS(DCV(4)) *COS(DCV(6))) *DCSM(5,5)
      ARG3 = VOSQ*(2 * COS(DCV(6)) * SIN(DCV(6)) -
     1   2*SIN(DCV(6)) *COS(DCV(4)) *COS(DCV(5))) *DCSM(6,6)
c
C      print '(a,f12.5,i4)', 'SAQF(IP)= ',saqf(ip),ip
      W(IP) = PAR(IP,1) * tmassa(IP) * VOL(IP)/SAQF(IP)
C 
      WRITE(6,2705)IP, VOL(IP), tmassa(IP),1.66113*tmassa(ip)/vol(IP)
2705  FORMAT(' Volume('I2,')= ', F9.3, ' UCW= ',F7.2,
     *' U.C.Density = 'f7.3,' gr/cm^3')
10      CONTINUE
C ****** QUANTITATIVE ANALYSIS ***************
      write(6,6969)
6969  format(/)      
      WTOTAL = 0.000000
      DO 2707 I = 1, NPHASE
      WTOTAL = WTOTAL + W(I)
2707  CONTINUE
      DO 2708 I = 1, NPHASE
C      !cp nov 10 96
              XMASS(I) = 100. * W(I) / WTOTAL
2708  CONTINUE
      IINNMOL = 0
      DO 2713 I = 1, NPHASE
      IF (IINNMOL.EQ.1) GOTO 2713
      IF (NMOL(I).EQ.0) IINNMOL=1
2713    CONTINUE
      IF (IINNMOL.EQ.1) THEN
         DO 2714 I = 1, NPHASE
C ** printing results
         WRITE(6,2741) I, XMASS(I)
2741     FORMAT(' PHASE = ',I2,' => %MASS = ',F6.2,
     1     '%MOLAR = NOT COMPUTED')
2714     CONTINUE
      ELSE
C ****    CALCULATION OF MOLAR FRACTION  ****
       FT = 0.0000000
       DO 2709 I = 1, NPHASE
       FRP(I) = XMASS(I) * NMOL(I) / tmassa(I)
       FT = FT + FRP(I)
2709     CONTINUE
       DO 2710 I = 1, NPHASE
       FR(I) = 100. * FRP(I) / FT
C ** printing results
       WRITE(6,2711) I, XMASS(I), FR(I)  
2711     FORMAT(' PHASE = ',I2,' => %MASS = ',F6.2,2X,
     1     '%MOLAR = ',F6.2)
2710   CONTINUE
      ENDIF
      if(isphase.ne.0)then
       write(6,2717)
2717    format(/,' Considering Amorphous Content:')    
         sfIS=wtis(isphase)/xmass(isphase)
          if(sfIS.gt.1.)then
            write(6,2718)ISPHASE, wtis(isphase)
2718        format(' PROBLEM:Amount of Internal Standard (Phase #',i2,
     *             ') is less than the specified ' f6.2,'%.'/,
     *             ' Amorphous content not computed. Check ISWT in'
     *             ' line 11.2 for this phase')
            go to 2720
          end if
         do 2715 i=1,nphase
           write(6,2716)i,xmass(i)*sfIS  
2716       format(' PHASE = ',I2,' => %MASS = ',F6.2)
2715     continue   
       write(6,2719)100*(1.-sfIS)
2719       format(' AMORPHOUS  => %MASS = ',F6.2)       
      end if
2720  write(6,2712)
2712  format(1x,/)
      return
!      stop
      end
C      SUBROUTINE WRIT94 Gera uma saida para o 9411, pra poder Usar no ATOMS
      SUBROUTINE WRITE94(ISCALE,IDIF)
      include 'param.inc'
C      PARAMETER (idsz=2048 ,irs=512,NATS=64,msz=32 ,nov=128)
      REAL LAMDA
      INTEGER PTR,fondo
      COMMON/DC/SAVE(99,6),NSAVE
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
     * ,ULOR,VLOR,ZZZ,UC
      COMMON/CNTRLS/JOBTYP,NPROF,NPHASE,IOT,IPL,IPC,MAT,NXT,LST1,LST2,
     * LST3,IPL1,IPL2,IPLST,IPLOSS,IPLCAL,IPLPOL,
     * IPLCOM,IPLDIS,IPLAM,IPBIG,
     * MCYCLE,EPS,MAXS,INSTRM,MAXSX,MCYCLX,ICYRUN,
     * IPREF,IABSR,FONDO,IAS,iasym,sw,ibgd,idata,isphase,i2d94 
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
      CHARACTER*4 ATEXT(NATS),NTYP(NATS)
      COMMON/PARAC/ATEXT,NTYP 
      COMMON/JNK/ALOW(100),AHIGH(100),POS(100),BCK(100),
     * NATOM(99),NMOL(99),PREF(99,3),NBCKGD,NEXCRG,NSCAT
      CHARACTER SYMB(99,20)*1
      CHARACTER TITLE*70,PHSNM(99)*50
      COMMON/CHAR/SYMB,TITLE,PHSNM
      COMMON/COEFF/AC(10,16),POSI(30),SCAT(30),DFP(16),DFPP(16),XMAS(16)
      CHARACTER*4 NAM(16)
      COMMON/COEFC/NAM
      COMMON/MULTIP/TMASSA(99),MLTPHASE,xmltp(99),MURT(NATS),SAQF(99) 
     *,wtis(99)
      common/codebck/ibckcode
       REWIND 53
C line 1
      WRITE(53,1)TITLE
1     FORMAT(A70)
      JOBTYP=JOBTYP-1
C      NPROF=NPROF-1 
      INSTRM = INSTRM-1
C line 2.1
      WRITE(53,2)JOBTYP,NPROF,NPHASE,ibckcode,NEXCRG,NSCAT,INSTRM,IPREF,
     1          iabsr
2     FORMAT(9I4,21x,'LINE 2.1')
C line 3
      WRITE(53,14)IOT,IPL,IPC,MAT,NXT,LST1,LST2,LST3,IPL1,IPL2,IPLST
14    FORMAT(11i1,46X,'LINE 3')
C line 4
      WRITE(53,23)LAMDA,RATIO(2),BKPOS,WDT,CTHM,TMV,RLIM
23    FORMAT(3f8.5,5F8.4)
C line 5
      WRITE(53,30)MCYCLE,EPS,RELAX
30    FORMAT(1I4,5F4.2,33x,'CYCLS EPS RELAX P_CALC')
      IF(NBCKGD.LT.2)GOTO 120 
C line 6(*)
      WRITE(53,34)(POS(I),BCK(I),I=1,NBCKGD)
34    FORMAT(2F8.2)
120   IF(NEXCRG.LE.0)GOTO 122
C line 7(*)
      WRITE(53,341)(ALOW(I),AHIGH(I),I=1,NEXCRG)
341     format(2f8.2,41x,'EXCLUDED REGION')
122   IF(NSCAT.LE.0)GOTO 124
      DO 125 I=1,NSCAT
           IF (JOBTYP .EQ. 1.OR.JOBTYP.EQ.3) GOTO 1228
C line 8.1 XRD (*)
        WRITE(53,38)NAM(I),DFP(I),DFPP(I),XMAS(I),NSCAT
38        FORMAT(A4,3F8.4,29X,'SCATTERING SET ',I2)
        GOTO 126
C line 8.1 ND(*)
1228        WRITE(53,3838)NAM(I),DFP(I),XMAS(I),NSCAT
3838        FORMAT(A4,2F8.4,37X,'SCATTERING SET ',I2)
C line 8.2 XRD(*)
126       IF(JOBTYP.EQ.0.OR.JOBTYP.EQ.2)WRITE(53,39)(AC(J,I),J=1,9) 
39          FORMAT(9F8.5)
125     CONTINUE
124   CONTINUE
C line 9
      WRITE(53,43)MAXS
43    FORMAT(I8,49X,'PARAMS REFINED')
      N=0
      DO 4310 IIPHAS=1,NPHASE 
4310    N=N+NATOM(IIPHAS)
      DO 84 I=1,N
      DO 84 J=1,11
84        A(I,J)=SIGN(1.,A(I,J))*(FLOAT(10*LP(I,J))+ABS(A(I,J)))
      DO 87 I=1,NPHASE
      DO 88 J=1,6 
88        PAR(I,J+5)=SAVE(I,J)
      DO 87 J=1,27
87        APAR(I,J)=SIGN(1.,APAR(I,J))*(FLOAT(10*LPAR(I,J))+ABS(APAR(I,J
     *     )))
      DO 89 J=1,20
89      AGLB(J)=SIGN(1.,AGLB(J))*(FLOAT(10*LGLB(J))+ABS(AGLB(J)))
C line 10.1
      WRITE(53,45)GLB(1),AGLB(1),GLB(10),AGLB(10),GLB(11),Aglb(11)
45    FORMAT(6F8.4,9x,'ZER DISP TRANS + CODEWORS')
C new lines in the ICF !cp may 10 97 and (ibgd test) !cp jun 97
C      if (ibgd.eq.1) goto 4600 
C      if (iax.eq.0) then
C line 10.2
C         write(53,634)glb(20),glb(18),glb(19)
C634      format(3f8.4,33X,'AM MON1 MON2')
C line 10.21
C         write(53,6342)aglb(20),aglb(18),aglb(19)
C6342     format(3f8.4,33x,'CODEWORS')
C     LINES 10.3 AND 10.31
C4600   IF(NBCKGD.EQ.0)WRITE(53,46)(GLB(J),J=2,7),(AGLB(J),J=2,7)
       WRITE(53,46)0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
46    FORMAT(6f9.4,3x,'BACKGROUND',/,6F9.4,3X,'CODEWORDS') 
477   DO 81 K=1,NPHASE
      IOF=0
      IF(K.GT.1) THEN
        DO 4311 IIPHAS=2,K
4311        IOF = IOF + NATOM (IIPHAS-1)
      END IF
C line 11.1
      WRITE(53,465)PHSNM(K),K
465     FORMAT(A50,7X,'PHASE NUMBER ',I2)
C line 11.2
      WRITE(53,48)NATOM(K),NMOL(K), (PREF(K,I),I=1,3)
48      FORMAT(2I4,8X,3F4.1,29X,'#ATMS #FU PREFDIR')
      N=NATOM(K)
C line 11.3
      WRITE(53,460)(SYMB(K,I),I=1,20)
460     FORMAT(20A1,37X,'SPACE GROUP')
C Changing N to 'so' !cp oct 96 
C           do 3947 isof=1,n
C           xl(isof+iof,5)=xl(isof+iof,5)/murt(isof+iof)
C3947         continue
C
      WRITE(53,65)(ATEXT(I+IOF),NTYP(I+IOF),
     *     (XL(I+IOF,J),J=1, 5),(A(I+IOF,J),J=1, 5),
     *        (XL(I+IOF,J),J=6,11),(A(I+IOF,J),J=6,11),
     *            I=1,N)
65      FORMAT(A4,a4,8X,5F8.5,2x,'LBL NTYP x y z B So'/,
     *      16X,5F8.2,2x,'CODEWORDS'/,6F8.5,10X,'BETAS',/,
     *      6F8.2,10x,'CODEWORDS')
      WRITE(53,50) PAR(K,1),PAR(K,2),
     *                 APAR(K,1),APAR(K,2),
     *PAR(K,3),PAR(K,4),PAR(K,5),PAR(K,21),par(k,20),
     *         PAR(K,15),PAR(K,16),
     *APAR(K,3),APAR(K,4),APAR(K,5),APAR(K,21),apar(k,20),
     *         APAR(K,15),APAR(K,16),
     *              (PAR(K,I),I=6,11), (APAR(K,I),I=6,11), 
     *              PAR(K,12),PAR(K,13),PAR(K,14),
     *                 APAR(K,12),APAR(K,13),APAR(K,14),
     *              PAR(K,17),PAR(K,18),PAR(K,19),
     *                 APAR(K,17),APAR(K,18),APAR(K,19),
     1              par(k,24),par(k,25),par(k,26),
     1                 apar(k,24),apar(k,25),apar(k,26),
     1              par(k,27),
     1                 apar(k,27)
C !cp oct 96. End of FORMAT MODIFICATION
50      FORMAT(G8.3,F8.4,41X,'SCALE Bo(OVERALL)',/,2F8.2,/,
     *        7F8.5,1X,'U V W CT Z X Y',/,7F8.2,/,
     *        6F8.4,9X,'CELL PARAMETERS',/,6F8.2,/,
     *        3F8.5,33X,'PREF1 PREF2 R/RCF_ASYM',/,3F8.2,/,
     *        3F8.4,33X,'NA NB NC (MIX_PARAMS)',/,3F8.2,/,
     1        3F8.4,33X,'NA NB NC (HIGH SIDE)',/,3F8.2,/,
     1         f8.4,49x,'PEARSON ASYM.FACTOR',/,F8.2)
81      CONTINUE
      IF(IPL.NE.0)WRITE(5,191)ISCALE,IDIF
191   FORMAT(2I8,41X,'LINE PRINTER INFO')
C lines commented to avoid LINE 13 in the ICF    !cp jul 97
C      IF(IPL2.NE.0)WRITE(5,171)IFY,IFYC,IFM,IFD,IFB,EXPAND
C      IF(IPL2.EQ.0)WRITE(5,171)  2,   1,  1,  1,0,0.95
c171   FORMAT(5I1,3X,F8.4,41x,'CALCOMP INFO')
C                       JOBTYP  must be reproduced !!!
151   JOBTYP=JOBTYP+1
      RETURN
      END 
C subroutine to compute size&strain (NsizeStrain)
            subroutine size(k)
      include 'param.inc'
      real lamda             
      COMMON/PARAMS/GLB(20),AGLB(20),LGLB(20),PAR(99,30),APAR(99,30),
     * LPAR(99,30),XL(NATS,11),LP(NATS,11),A(NATS,11),
     *  PTR(NATS),RELAX(4),RATIO(2)
c      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,BKPOS,WDT
c     * ,ULOR,VLOR,ZZZ,UC
      COMMON/G1/THMIN,STEP,THMAX,U,V,W,LAMDA(2),TMV,CTHM,RLIM,SBX,WDT 
     * ,ULOR,VLOR,ZZZ,UC
      common/sizestrain/sizeG(15),strainG(15),sizeL(15),strainL(15)
     *,siz(15),strain(15),NsizeStrain
            rado=3.1415927/180.
            a7 = 2.69269
            b7 = 2.42843
            c7 = 4.47163
            d7 = 0.07842
      open (17,file='instr.dat',status='old')
      read(17,*)
      read(17,*)Ui,Zi,Xi,Yi
c na primeira linha de padrao.dat sera apenas para comentario
c a segunda linha ira ler U,Z,X,Y da FWHM do padrao
      close(17)
      Ua = par(k,3)
      Za = par(k,20)
      Xa = par(k,15)
      Ya = par(k,16)
      if (Zi .NE. 0.) then
        write(6,7108)
        write(*,7108)
7108   format( ' SIZE NOT COMPUTED * bad standard? - check Z'
     *' for standard.')
        goto 20
      end if
      if(Za .LT. 0. .OR. Ya .LT. Yi)then
            write(6,710)k
            write(*,710)k
        goto 20
      end if
      if(Ya .EQ. Yi .AND. Za .EQ.0.) then
            write(6,711)k
            write(*,711)k
        IsizeTest = 1
      else
        IsizeTest = 0  
      end if
710   format(' Phase ',i2,': Size-strain not Computed * check Z and'
     *,' Y for sample and standard')
711   format('Phase ',i2,': Infinite Size * check Z and Y values',
     *' for sample and standard')
c
C COMPUTE SIZE
c
C        COMPUTE SIZE BASED ON Gauss-COS Z parameter [sizeG(k)]
c
      if(Za.GT.0.) then
        hsg = sqrt(Za)*rado
        sizeG(k) =  lamda(1) / hsg 
      else  
        hsg = 0.
        sizeG(k) = 99999.        
      end if
c
C        COMPUTE SIZE BASED ON Loren-COS Y parameter [sizeL(k)]
c  par(k,16) = Y (Lorentz)
c
      if(Ya .gt. Yi) then
        hsl = (Ya - yi) * rado
        sizeL(k) =  lamda(1) / hsl  
      else 
        hsl = 0.
        sizeL(k) = 99999.
      end if  
c
C        COMPUTE WEIGHTED SIZE BASED on Gauss & Lorentz [siz(k)]
c
      if(IsizeTest.eq.0)then
              hs  = (hsg**5 + a7 * hsg**4 * hsl + b7 * hsg**3 * hsl**2
     !      + c7 * hsg**2 * hsl**3 + d7 * hsg * hsl**4 + hsl**5)**(0.2)
         siz(k) = lamda(1) / hs
      else
         siz(k) = 0.0   
      end if
c
C       STRAIN COMPUTATIONS
10    if(Ua .LT. Ui .OR. Xa .LT. Xi)then
          write(6,720)k
          write(*,720)k
          goto 20
      end if
      if(Xa .EQ. Xi .AND. Ua .EQ. Ui)  then
          write(6,721)k
          write(*,721)k
          IstrainTest = 1 
      else
          IstrainTest = 0
      end if
c      
720   format(' Phase ',i2, ': Size-Strain not computed * check U and X'
     *' for sample and standard')
721   format(' Phase ',i2, ': Strain not Present * check U and X values'
     *' for sample and standard')
c
C COMPUTE STRAIN
c
c compute strain based on Gauss-U * [strainG(k)]
c
      if(Ua .GT. Ui)  then
          hdG = sqrt( (Ua - ui) )*rado
          strainG(k) =  hdG / 2
      else
          hdG = 0.
          strainG(k) = 0.0
      end if
c
c compute strain based on Lorentz-X * [strainL(k)]
c
      if(Xa .gt. Xi)  then
          hdL = (Xa - xi) * rado
          strainL(k) =  hdL / 2
      else
         hdL = 0.0  
         strainL(k) = 0.0
      end if
c compute weighted strain: based on Gauss-U and Lorentz-X * [strain(k)]
C
      if (IstrainTest.eq.0)then
         hd  = (hdg**5 + a7 * hdg**4 * hdl + b7 * hdg**3 * hdl**2
     !     + c7 * hdg**2 * hdl**3 + d7 * hdg * hdl**4 + hdl**5)**(0.2)
         strain(k) = hd / 2
      else
         strain(k) = 0.
      end if
c      
      write(6,15)k,sizeG(k),strainG(k),k,sizeL(k), strainL(k),k,siz(k),
     *           strain(k)
15    format(4x,'Gauss(phase ',i2,'): Size =',g9.4,' Strain =',g8.3,/,
     *     2x,'Lorentz(phase ',i2,'): Size =',g9.4,' Strain =',g8.3,/,
     *    1x,'Weighted(phase ',i2,'): Size =',g9.4,' Strain =',g8.3)
20    return
      end
      
