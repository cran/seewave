/* ------- Short Term Fourier Transform on a decimal file with mean return + matrix by: Alicia STOTZ (University of Evry) ------- */
# ifndef DEF_FOURIERMEAN
# define DEF_FOURIERMEAN
# include <stdio.h>
# include <stdlib.h>
# include <fftw3.h>
# include <math.h>
SEXP fourierMeanC (SEXP MATRIX, int *ptrnorm, int *ptrdb);
# endif 

//the arguments: gorring matrix, pointer norm and pointer dB
SEXP fourierMeanC (SEXP MATRIX, int *ptrnorm, int *ptrdb)
{

	SEXP matrixMean;//return matrix
	double *matrixmeanptr;//pointer mean matrix
	double *matrixptr;//pointer  matrix
	SEXP dimMatrix;//input matrix
	int nrow = 0, ncol = 0, i = 0, j = 0;
	SEXP vecMean;
	double *vecmeanptr;
	SEXP List; //return list
	SEXP names;// arguments names
	double max = 0, tmp = 0, maxMatrice = 0;


//Recovery of the size of the matrix
	dimMatrix = getAttrib(MATRIX, R_DimSymbol);

//mension of the matrix are retracted into the variables: nrow and ncol
	nrow = INTEGER (dimMatrix)[0];
	ncol = INTEGER (dimMatrix)[1];

//Allocation memory
	PROTECT(MATRIX = coerceVector( MATRIX, REALSXP));
	matrixptr = REAL(MATRIX);

//Memory allocation of out matrix 
	PROTECT(matrixMean = allocMatrix(REALSXP, nrow, ncol));
//matrix pointer real type
	matrixmeanptr = REAL(matrixMean);
//Memory allocation return List
	PROTECT( List = allocVector (VECSXP, 2));
//Memory allocation argument names
	PROTECT( names = allocVector (STRSXP, 2));
//Memory allocation mean vector
	PROTECT ( vecMean = allocVector (REALSXP, nrow));
//mean vector pointer real type
	vecmeanptr = REAL(vecMean);

	SET_STRING_ELT(names, 0, mkChar ("mean"));
	SET_STRING_ELT (names, 1, mkChar ("amp"));
	setAttrib(List, R_NamesSymbol, names);

//Gorring matrix are take in return matrix
	for (j = 0; j < ncol; j++)
	{
		for (i = 0; i < nrow; i++)
		{
			matrixmeanptr[ i + nrow * j] = matrixptr[i + nrow * j];
			if (matrixptr[i + nrow * j] > maxMatrice) maxMatrice = matrixptr[i + nrow * j];
		}
	}


//Calculate the mean vector from matrix
	if ( *ptrnorm == 0 && *ptrdb == 0)
	{
//Calculate of line mean
		for (j = 0; j < nrow; j++)
		{
			for (i = 0; i < ncol; i++)
			{

				tmp = tmp + matrixptr[j+i *nrow];
			}

			vecmeanptr[j] = tmp /ncol;
			tmp = 0;
		}
	}

//Calculate the mean vector from matrix
	else if ( (*ptrnorm == 1 && *ptrdb == 0) || (*ptrnorm == 0 && *ptrdb == 1) )
	{
//Calculate of line mean
		for (j = 0; j < nrow; j++)
		{
			for (i = 0; i < ncol; i++)
			{
				
				tmp = tmp + matrixptr[j+i *nrow];
			}

			vecmeanptr[j] = tmp /ncol;
			if (vecmeanptr[j] > max)
				{
					max = vecmeanptr[j];
				}
			tmp = 0;
		}
//vector of line mean are normed
		for (j = 0; j < nrow; j++)
		{
			vecmeanptr[j] = vecmeanptr[j]/max;
		}
//Calculate of norm Matrix
			for (j = 0; j < ncol; j++)
			{
				for (i = 0; i < nrow; i++)
				{
					matrixmeanptr[ i + nrow * j] = matrixptr[i + nrow * j]/maxMatrice;
				}
			}
//vector of line mean are take in deciBel
		if ( *ptrdb == 1)
		{
			for (j = 0; j < nrow; j++)
			{
				vecmeanptr[j] =( 20 * log10 (vecmeanptr[j]));
			}
//Calculate dB matrix
			for (j = 0; j < ncol; j++)
			{
				for (i = 0; i < nrow; i++)
				{
					matrixmeanptr[ i + nrow * j] = (20 * log10(matrixptr[i + nrow * j]/maxMatrice));
				}
			}
		}
	}

//add SEXP in return List
	SET_VECTOR_ELT (List, 0, vecMean);
	SET_VECTOR_ELT (List, 1, matrixMean);


UNPROTECT(5);
return List;

}

