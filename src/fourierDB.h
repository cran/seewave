/* ------- Short Term Fourier Transform on a decimal file return db matrix by: Alicia STOTZ (University of Evry) ------- */
# ifndef DEF_FOURIERDB
# define DEF_FOURIERDB
# include <stdio.h>
# include <stdlib.h>
# include <fftw3.h>
# include <math.h>
SEXP fourierDBC (SEXP matrixNorm);
# endif 

//Prototype of fourier function return norm Matrix
SEXP fourierNormC (int numberOfLine, float *buffer, int *ptrwl, double *ptrovlp);

//the argument: norm matrix.
SEXP fourierDBC (SEXP MATRIXNORM)
{
	SEXP matrixDB;//return matrix
	double *matrixdbptr;//pointer dB matrix
	double *matrixptr;//pointer norm matrix
	SEXP dimMatrix;//input matrix
	int nrow = 0, ncol = 0, i = 0, j = 0;


//Recovery of the size of the matrix
	dimMatrix = getAttrib(MATRIXNORM, R_DimSymbol);

//mension of the matrix are retracted into the variables: nrow and ncol
	nrow = INTEGER (dimMatrix)[0];
	ncol = INTEGER (dimMatrix)[1];

//Allocation memory
	PROTECT(MATRIXNORM = coerceVector( MATRIXNORM, REALSXP));
	matrixptr = REAL(MATRIXNORM);

//Memory allocation of out matrix 
	PROTECT(matrixDB = allocMatrix(REALSXP, nrow, ncol));
//matrix pointer real type
	matrixdbptr = REAL(matrixDB);

	for (j = 0; j < ncol; j++)
	{
		for (i = 0; i < nrow; i++)
		{
//Calculate dB matrix
			matrixdbptr[ i + nrow * j] = (20 * log10(matrixptr[i + nrow * j]));
		}
	}

UNPROTECT(2);
return matrixDB;

}
