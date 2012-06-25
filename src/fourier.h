/* ------- Short Term Fourier Transform on a decimal file by: Alicia STOTZ (University of Evry) ------- */
# ifndef DEF_FOURIER
# define DEF_FOURIER
# include <stdio.h>
# include <stdlib.h>
# include <fftw3.h>
# include <math.h>
SEXP fourierC (int numberOfLine, float *buffer,int *ptrwl, double *ptrovlp);
# endif 


//the arguments: number of line, array containing the hexadecimal sound, the window size and the overlap
SEXP fourierC(int numberOfLine, float *buffer,int *ptrwl, double *ptrovlp)
{
	double *in; //Table of input
	fftw_complex *out;//Table of exit
	fftw_plan plan_forward; //FFT plan
	int i = 0, j = 0, l = 0;
	int windowsize = 0; //STFT window
	int nrow = 0, ncol = 0; //Number of line and row
	double mod = 0;
	SEXP matrix;//out matrix
	double *matrixptr;//matrix pointer
	double *hanning;//hanning window
	
	
windowsize = *ptrwl;

//Calculating the number of lines in the output matrix
		if ( *ptrovlp == 0) ncol = (numberOfLine / windowsize);
		else ncol = (numberOfLine / windowsize)/ (*ptrovlp / 100);
//Calculating the number of row in the output matrix
	nrow = (windowsize / 2);

//Memory allocation of the temporary array	
	in = fftw_malloc( sizeof (double) * numberOfLine );
	hanning = fftw_malloc( sizeof (double) * numberOfLine );

//Memory allocation out table of the FFT
	out = fftw_malloc (sizeof (fftw_complex) * numberOfLine);

	if (hanning == NULL || in == NULL || out == NULL) {
		Rprintf("Could not allocate memory for file\n");
		return (R_NilValue);
	}

//Memory allocation of out matrix 
	PROTECT(matrix = allocMatrix(REALSXP, nrow, ncol));
//matrix pointer real type
	matrixptr = REAL(matrix);


/**********	STFT	**********/
for ( j = 0; j < ncol; j++)
{
	for ( i = 0; i < windowsize; i++)
	{

//calculate the hanning window in order to norm the in array		
	hanning[i] = 0.5 - 0.5 * cos (2.0*M_PI*(double)i/(double)(windowsize-1));
//the value of i to i + 512 are placed in a temporary array
	in[i] = (double) buffer[i+l] * hanning[i];
	}

//Plan and real type 1D
	plan_forward = fftw_plan_dft_r2c_1d (windowsize, in, out, FFTW_ESTIMATE);
//Fft Execusion
	fftw_execute (plan_forward);
	
//release of plan
fftw_destroy_plan(plan_forward);

		for ( i = 0; i < nrow ; i++)
	{
//Calculate modulus
		mod = sqrt( (out[i][0]*out[i][0]) + (out[i][1]*out[i][1]) );
//Completion of the matrix with the calculated modulus
		matrixptr[i + nrow * j] = mod;
	}
//A step increase of 512
if ( *ptrovlp == 0) l = l + windowsize;
else l = l + (nrow * (*ptrovlp/100));
}

//release of memory space
fftw_free(in);
fftw_free(out);
fftw_free(hanning);

UNPROTECT(1);

return matrix;

}


