#include "constants.h"

void testhello(void);
void calculateIntersections(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, int sizeP, int sizeQ);
void countIntersections(double *polyPX, double *polyPY, double *polyQX,  double *polyQY, int sizeP, int sizeQ, int *countNonDegenIntP, int *countNonDegenIntQ, double **intersectionsP, double **intersectionsQ, int **initLabelsP, int **initLabelsQ);
