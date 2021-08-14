// #include <Python.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

typedef struct{
    int size;
    double* mean;
    double* prevMean;

} cluster;

typedef struct Node Node;

struct Node{
    double* data;
    Node* next;
};

typedef struct{
    int rows;
    int columns;
    double ** data;
    
} matrix;

typedef struct{
    double value;
    int column;
} eigenVal;

enum goal{spk,wam,ddg,lnorm,jacobi};

#define true 1
#define false 0
#define MaxJacobiIter 100
#define epsilon 0.001


void testPrint(char arr[]);
double ** kmeansFunc(int k, int maxIter, int numOfPoints, int d, Node* pointsMatrix, cluster* centroids);
Node* getPoints(char* fileName, int* numOfPoints, int* finald);
// cluster* makeClusters(int k, int d, PyObject ** centroids);
// Node* getPointsMatrix(int d, int n, PyObject ** pointsMatrix);
void doKmeans(cluster* clusters, Node* points, int d, int k, int numOfPoints, int maxIteratiorns);
int continueLoop( int convergence, int count, int maxIter);
void iterateClusters(cluster* clusters, int length, int d);
int minCluster(double* point, cluster* clusters, int K, int d);
double calcDistance(double* p1, double* p2, int d);
void addPoint(cluster c, double* point, int d);
int checkCon(cluster* clusters, int k, int d);
Node* addNext(Node* node, double* point);
void printClusters(cluster* clusters, int k, int d);
void freeMemo(cluster* clusters, Node* points, int k, int numOfpoints);

cluster* makeClustersSp(int k, Node* points, int numOfPoints,int d);

matrix* formMatW(Node* points, int numOfPoints, int d);
matrix* formMatD(matrix* matW);
matrix* formMatLnorm(matrix* matD , matrix* matW);
matrix* minusRootMat(matrix* mat);
matrix* formMatI(int dimention);
matrix* pointsInMat(Node* points, int numOfPoints, int d);

matrix* mulMatrices(matrix* mat1, matrix* mat2, int free1, int free2);
matrix* addMatrices(matrix* mat1, matrix* mat2, int dec, int free1, int free2);
matrix* newMatrix(int rows, int columns);
void freeMatrix(matrix * m);

void jacobiAlg(matrix* mat, matrix* eigenValues, matrix* eigenVectors);
matrix* getRotationMatrixValues(matrix* mat, double *c, double *s, int *rowPivot, int *colPivot);
void findLargestCell(matrix* mat, int *row, int *col) ;
void calcNextJacobiMatrix(matrix* matA, matrix* matB, double c, double s, int i, int j);
int hasConvergence(matrix* matA, matrix* matB);
double calcOff(matrix* m);
matrix* calcInitialVectorsFromJacobi(matrix* eigenValues, matrix* eigenVectors);
int* eigenGapHeuristic(matrix* matA, int* k);
int findMaxGap(eigenVal* values, int length);
int isDiagonal(matrix* m);
int compareEigenVal(const void * a, const void * b);
