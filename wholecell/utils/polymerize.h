/* Header to test of C modules for arrays for Python: C_test.c */

/* ==== Prototypes =================================== */

void init_polymerize(void);

static PyObject *polymerize(PyObject *self, PyObject *args);
static void _polymerize(
	const double elngRate, long int *deficitAAs, long int const *requiredAAs,
	long int *aaCounts, long int *updatedAAs, long int *aasUsed, const long int seed,
	const int nRows, const int nCols
	);

static long int whereNonZero(const int *array, long int len);
static unsigned char allZero_long(const long int *array, long int len);
static unsigned char allZero_double(const double *array, long int len);
static void elementwiseMin(
	const long int *a, const long int *b,
	double *result, long int len
	);

