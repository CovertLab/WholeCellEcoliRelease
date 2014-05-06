#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "arrayobject.h"

#include "polymerize.h"
#include <math.h>
#include <gsl/gsl_randist.h>

#define TRUE 1
#define FALSE 0
#define MIN(a, b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a, b) ( ((a) > (b)) ? (a) : (b) )

/* #### Globals #################################### */

/* ==== Set up the methods table ====================== */
static PyMethodDef _polymerizeMethods[] = {
	{"polymerize", polymerize, METH_VARARGS},
	{NULL, NULL}     /* Sentinel - marks the end of this structure */
};

/* ==== Initialize the functions ====================== */
void init_polymerize()  {
	(void) Py_InitModule("_polymerize", _polymerizeMethods);
	import_array();  // Must be present for NumPy.  Called first after above line.
}

/* #### Vector Extensions ############################## */

static PyObject *polymerize(PyObject *self, PyObject *args)
{
	PyArrayObject *deficitAAs, *requiredAAs, *aaCounts;
	PyArrayObject *updatedAAs, *aasUsed;
	long int seed = 0;
	int nRows, nCols;
	double elngRate;

	if (!PyArg_ParseTuple(args, "dO!O!O!O!O!I", &elngRate,
		&PyArray_Type, &deficitAAs,
		&PyArray_Type, &requiredAAs,
		&PyArray_Type, &aaCounts,
		&PyArray_Type, &updatedAAs,
		&PyArray_Type, &aasUsed,
		&seed))
		return NULL;
	if (NULL == deficitAAs)
		return NULL;
	if (NULL == requiredAAs)
		return NULL;
	if (NULL == aaCounts)
		return NULL;
	if (NULL == updatedAAs)
		return NULL;
	if (NULL == aasUsed)
		return NULL;

	nRows = PyArray_DIM(deficitAAs, 0);
	nCols = PyArray_DIM(deficitAAs, 1);

	/* Numpy array data is getting passed in row-major order */
	_polymerize(
		elngRate, PyArray_DATA(deficitAAs), PyArray_DATA(requiredAAs),
		PyArray_DATA(aaCounts), PyArray_DATA(updatedAAs), PyArray_DATA(aasUsed),
		seed, nRows, nCols
		);

	Py_RETURN_NONE;
}

static void _polymerize(
	const double elngRate, long int *deficitAAs, long int const *requiredAAs,
	long int *aaCounts, long int *updatedAAs, long int *aasUsed, const long int seed,
	const int nRows, const int nCols
	)
{
	size_t K = nCols;
	const gsl_rng_type *T;
	gsl_rng * r;
	int results[K];
	double minResults[K];

	unsigned char stop = FALSE;
	unsigned long int i, j;

	/* Random number generator setup */
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	for(i = 0; i < (int) elngRate  && !stop; i++){
		for(j = 0; j < nRows && !stop; j++){

			if(allZero_long(aaCounts, nCols)){
				stop = TRUE;
				break;
			}
			if(allZero_long(deficitAAs + j * nCols, nCols))
				continue;

			long int selectedIdx;
			
			elementwiseMin((deficitAAs + j * nCols), aaCounts, minResults, K);

			if(allZero_double(minResults, K))
				continue;
			
			gsl_ran_multinomial(r, K, 1, minResults, (unsigned int *) results);

			selectedIdx = whereNonZero(results, K);

			if(selectedIdx < 0){
				fprintf(stderr, "Non-zero index not found!\n");
				return;
			}
			
			updatedAAs[j * nCols + selectedIdx] += 1;
			aaCounts[selectedIdx] -= 1;
			aasUsed[selectedIdx] += 1;
			deficitAAs[j * nCols + selectedIdx] -= 1;

			if(minResults[selectedIdx] == 0){
				fprintf(
					stderr,
					"Selected index [%ld] that had zero probability of being selected!\n",
					selectedIdx
					);
				return;
			}
		}

	}

	gsl_rng_free(r);

}

static long int whereNonZero(const int *array, long int len)
{
	long int i;
	for(i = 0; i < len; i++){
		if(array[i] != 0)
			return i;
	}
	return -1;
}

static unsigned char allZero_long(const long int *array, long int len)
{
	long int i;
	for(i = 0; i < len; i++){
		if(array[i] != 0)
			return FALSE;
	}
	return TRUE;
}

static unsigned char allZero_double(const double *array, long int len)
{
	long int i;
	for(i = 0; i < len; i++){
		if(array[i] != 0)
			return FALSE;
	}
	return TRUE;
}

static void elementwiseMin(const long int *a, const long int *b, double *result, long int len)
{
	long int i;
	for(i = 0; i < len; i++){
		result[i] = (double) MIN(a[i], b[i]);
		result[i] = MAX(result[i], 0.0);
	}
}
