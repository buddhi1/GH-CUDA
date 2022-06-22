#include "constants.h"

void calculateIntersections(
                double *polyPX, double *polyPY, double *polyQX,  double *polyQY, 
                int sizeP, int sizeQ);
void calculateIntersections(
                double *polyPX, double *polyPY, 
                double *polyQX,  double *polyQY, 
                int sizeP, int sizeQ, 
                int *countNonDegenIntP, int *countNonDegenIntQ, 
                double **intersectionsP, double **intersectionsQ, int **alphaValuesP, int **alphaValuesQ,
                int **initLabelsP, int **initLabelsQ,
                int **neighborMapP, int **neighborMapQ, int **neighborP, int **neighborQ);
void calculateIntersectionsMultipleComponents(
                double *polyPX, double *polyPY, 
                double *polyQX,  double *polyQY, 
                int sizeP, int sizeQ, int *sizesPP, int *sizesQQ, int sizePP, int sizeQQ,
                int *countNonDegenIntP, int *countNonDegenIntQ, 
                double **intersectionsP, double **intersectionsQ, int **alphaValuesP, int **alphaValuesQ,
                int **initLabelsP, int **initLabelsQ,
                int **neighborMapP, int **neighborMapQ, int **neighborP, int **neighborQ);
