#ifndef PARAM_INC_H
#define PARAM_INC_H



// IDSZ   - maximum number of data points in tape 4
const int IDSZ = 2*8000;  //90000;

// IRS    - maximum number of reflections in the pattern
const int IRS = 2*4096; // 20000;

// NATS   - maximum number of atoms in the problem
const int NATS = 2*128; //512;

// MSZ    - maximum number of refinable parameters (matrix size)
const int MSZ = 2*64;  //256;

// NOV    - maximum number of Bragg reflections contributing to one step
const int NOV = 4*512; //15000;



//     nfinal - dimension of the table FINAL(nfinal) in COMMON/ALLP/
//		It must be equal to (#atoms*11 + #phases*27 + 11), at least,
//		where #atoms & #phases are actually numbers of atoms
//		and phases in the input control file for the current refinement.
// 	Previous declaration FINAL(8*MSZ,2) caused error in case of many atoms.
const int NFINAL = 6048;


#endif // PARAM_INC_H
