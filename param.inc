C       PARAMETER (IDSZ= 8000,IRS=4096,NATS=128,MSZ= 64,NOV=512)
C       PARAMETER (IDSZ=20000,IRS=2048,NATS=128,MSZ=140,NOV=256)
C       parameter (idsz=24000,irs=2048,NATS=128,MSZ=154,nov=512)
C       PARAMETER (IDSZ=24000,IRS=8192,NATS=256,MSZ=256,NOV=8192)
c       PARAMETER (IDSZ=65530,IRS=12000,NATS=512,MSZ=256,NOV=8192)
c*****                                 Para a homepage size2k
c        PARAMETER (IDSZ=24000,IRS=2048,NATS=128,MSZ=154,NOV=512)
c*****                                 Para o homepage  size2ki
        PARAMETER (IDSZ=90000,IRS=20000,NATS=512,MSZ=256,NOV=15000)
        parameter (nfinal=6048)

c			 PARAMETER INCLUDE FILE
c	In order to get DBWS working with other parameters
c	 (e.g. IDSZ, IRS, NATS, MSZ, NOV) you can easily change 
c	desired PARAMETERs in this file and next recompile the
c	the source code using following commands:
c		fl *.for
c		link db1+db2+db3,dbws/e;
c
c     IDSZ   - maximum number of data points in tape 4
c     IRS    - maximum number of reflections in the pattern
c     NATS   - maximum number of atoms in the problem
c     MSZ    - maximum number of refinable parameters (matrix size)
c     NOV    - maximum number of Bragg reflections contributing to one step
c     nfinal - dimension of the table FINAL(nfinal) in COMMON/ALLP/
c		It must be equal to (#atoms*11 + #phases*27 + 11), at least,
c		where #atoms & #phases are actually numbers of atoms
c		and phases in the input control file for the current refinement.
c 	Previous declaration FINAL(8*MSZ,2) caused error in case of many atoms.

