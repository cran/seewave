/* ------Reading sound file from the name returned in R and calculate STFT by: Alicia STOTZ (University of Evry) ------*/

# include <sndfile.h> //Library for write or read sound file
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <R.h>
# include <Rinternals.h>
# include <Rdefines.h>
# include <Rmath.h>
# include <Rversion.h>
# include <unistd.h>
//# include <malloc.h> // unnecessary and stops Mac CRAN checks
# include <fftw3.h> //Library for fftw calculate
# include <math.h>

//Other fourier function
# include "fourier.h"
# include "fourierNorm.h"
# include "fourierDB.h"
# include "fourierMean.h"


//Prototype of fourier function
SEXP fourierC (int numberOfLine,float *buffer, int *ptrwl, double *ptrovlp);
//Prototype of fourier function with mean vector return
SEXP fourierMeanC (SEXP MATRIX, int *ptrnorm, int *ptrdb);
//Prototype of fourier function return norm Matrix
SEXP fourierNormC ( int numberOfLine, float *buffer,int *ptrwl, double *ptrovlp);
//Prototype of fourier function return dB Matrix
SEXP fourierDBC (SEXP matrixNorm);

//Function that read the .WAV file.
SEXP stft ( SEXP FILENAME, SEXP OVLP, SEXP WL,SEXP MEAN, SEXP NORM ,SEXP DB, SEXP VERBOSE)
{	

	char *fileName[1];//Sound name imput in R
	int *ptrnorm;//norm pointer
	int *ptrdb;//DB pointer
	int *ptrwl;//window pointer
	int *ptrmean = NULL; //mean pointer
	int *ptrverbose; //verbose pointer
	double *ptrovlp;//ovlp pointer
	int i = 0;
	long numFrames;
	float *buffer;//FFT buffer
	SF_INFO sndInfo;//Header sound file
	SNDFILE *sndFile;//Sound file
	int numberOfLine = 0;//Number Of lin in the sound
	SEXP outMatrix;//out matrix
	SEXP retList;//return List



//Recovery of the file name returned in R.
	PROTECT(FILENAME = AS_CHARACTER(FILENAME));
	*fileName = R_alloc(strlen(CHAR(STRING_ELT(FILENAME, 0))), sizeof(char));

	strcpy(*fileName, CHAR(STRING_ELT(FILENAME, 0)));
 

//Recovery of the ovlp returned in R.
	PROTECT(OVLP = coerceVector (OVLP, REALSXP));
	ptrovlp = REAL(OVLP);

//Recovery of the window size returned in R.
	PROTECT(WL = coerceVector (WL, INTSXP));
	ptrwl = INTEGER(WL);

//Recovery if the user want the norm matrix
	PROTECT(NORM = coerceVector (NORM, LGLSXP));
	ptrnorm = LOGICAL(NORM);

//Recovery if the user want the DB matrix
	PROTECT(DB = coerceVector (DB, LGLSXP));
	ptrdb = LOGICAL(DB);

//Recovery if the user want the MEAN matrix
	PROTECT(MEAN = coerceVector (MEAN, LGLSXP));
	ptrmean = LOGICAL(MEAN);

//Recovery if the user want the sound file information
	PROTECT(VERBOSE = coerceVector (VERBOSE, LGLSXP));
	ptrverbose = LOGICAL(VERBOSE);

// Open sound file
	sndFile = sf_open(*fileName, SFM_READ, &sndInfo);
	if (sndFile == NULL) {
		Rprintf("Error reading the source file '%s': %s\n", *fileName, sf_strerror(sndFile));
		return (R_NilValue);
	}

/*
// Check format - 16bit PCM
	if (sndInfo.format != (SF_FORMAT_WAV | SF_FORMAT_PCM_16)) {
		Rprintf("Input should be 16bit Wav\n");
		sf_close(sndFile);
		return (R_NilValue);
	}*/

// Check channels - mono
	if (sndInfo.channels != 1) {
		Rprintf("Wrong number of channels\n");
		sf_close(sndFile);
		return (R_NilValue);
	}

// Allocate memory
	buffer = malloc(sndInfo.frames * sizeof(float));

	if (buffer == NULL) {
		Rprintf("Could not allocate memory for file\n");
		sf_close(sndFile);
		return (R_NilValue);
	}

// Load data
	numFrames = sf_readf_float(sndFile, buffer, sndInfo.frames);

	if ( *ptrverbose == 1)
	{
//Print some of the info, and figure out how much data to read.
		Rprintf("File name = \'%s\'\n",*fileName);
		Rprintf("Number of samples = %d\n", (int) sndInfo.frames);
		Rprintf("Sampling rate (Hertz) = %d\n", sndInfo.samplerate);
		/* Rprintf("Number of channels = %d\n",sndInfo.channels); */
		Rprintf("Duration (seconds) = %f\n", (float)numFrames/sndInfo.samplerate);
		Rprintf ("Window length (number of samples)= %d\n", *ptrwl);
		Rprintf ("Window overlap (percentage) = %f\n",*ptrovlp);
	}

	
numberOfLine = (int) sndInfo.frames;

 	if ( *ptrnorm == 0 && *ptrdb == 0 && *ptrmean == 0)
	{
//Call fourier.h function
		outMatrix = fourierC (numberOfLine, buffer, ptrwl, ptrovlp);
	}
	else if (  (*ptrmean == 1 && *ptrnorm == 0 && *ptrdb ==0) || ( *ptrmean == 1 && *ptrnorm == 1 && *ptrdb == 0) || ( *ptrmean == 1 && *ptrdb == 1 && *ptrnorm == 0) )
	{
//Call fourierMean.h
		retList = fourierMeanC (fourierC (numberOfLine, buffer, ptrwl, ptrovlp), ptrnorm, ptrdb);;
		sf_close(sndFile);
		free(buffer);
		UNPROTECT(7);
		return retList;
	}
	else if (  *ptrnorm == 1 && *ptrdb == 0)
	{
//Call fourierNorm.h function
		outMatrix = fourierNormC (numberOfLine, buffer, ptrwl, ptrovlp);	
	}
	else if ( *ptrnorm == 0 && *ptrdb == 1)
	{
//Call fourierDB.h function
		outMatrix = fourierDBC (fourierNormC (numberOfLine, buffer, ptrwl, ptrovlp));
	}
//Error
	else
	{
		Rprintf ("The arguments returned do not match");
		sf_close(sndFile);
		free(buffer);
		UNPROTECT(7);
		return (R_NilValue);
	}
		
		
// Check correct number of samples loaded
	if (numFrames != sndInfo.frames) {
		Rprintf("Could not read enough frames for source\n\n");
		sf_close(sndFile);
		free(buffer);
		UNPROTECT(7);
		return (R_NilValue);
	}

	sf_close(sndFile);

	free(buffer);
	
UNPROTECT(7);

return outMatrix;
}
