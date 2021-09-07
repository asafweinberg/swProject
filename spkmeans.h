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

//enum goal{spk,wam,ddg,lnorm,jacobi};

#define true 1
#define false 0
#define MaxJacobiIter 100
#define epsilon 0.000000000000001

/*-----------------------------------------------------------------------------*/
/* ******************************************** SPK Functions ******************************************** */
/*-----------------------------------------------------------------------------*/

/* -------------------------------------------- Run Flow Functions -------------------------------------------- */
void runWam(Node* points, int numOfPoints, int d);
void runDdG(Node* points, int numOfPoints, int d);
matrix * runLnorm(Node* points, int numOfPoints, int d, int printMat0);
matrix * runSpk(int k, Node* points, int numOfPoints, int d); //returns matrix T nxk
void runJacobi(Node* points, int numOfPoints, int d, int k);
double ** runMainFlow(int k, char* myGoal, char* fileName, int* finalK, int* numOfPoints); //return T or NULL


/* -------------------------------------------- goals Functions -------------------------------------------- */
matrix* formMatW(Node* points, int numOfPoints, int d);
matrix* formMatT(matrix * matU, int free1);
matrix* formMatD(matrix* matW, int freeW);
matrix* formMatLnorm(matrix* matD , matrix* matW , int freeD, int freeW);


/* -------------------------------------------- jacobi Functions ------------------------------------------- */
void jacobiAlg(matrix* lNorm, matrix** eigenValues, matrix** eigenVectors);
matrix* getRotationMatrixValues(matrix* mat, double *c, double *s, int *rowPivot, int *colPivot);
void findLargestCell(matrix* mat, int *row, int *col);
matrix* setPivotMatrix(int row, int col, double c, double s, int dim);
matrix* calcNextJacobiMatrix(matrix* matA, double c, double s, int i, int j);
int hasJacobiConvergence(matrix* matA, matrix* matB);
double calcOff(matrix* m);
matrix* calcInitialVectorsFromJacobi(matrix* eigenValues, matrix* eigenVectors, int initialK);
int* eigenGapHeuristic(matrix* matA, int* k);
int findMaxGap(eigenVal* values, int length);
int compareEigenVal(const void * a, const void * b);


/* -------------------------------------------- matrices Functions ------------------------------------------- */
matrix* minusRootMat(matrix* mat, int free1);
matrix* mulMatrices(matrix* mat1, matrix* mat2, int free1, int free2);
matrix* addMatrices(matrix* mat1, matrix* mat2, int dec, int free1, int free2);
matrix* formMatI(int dimention);
int isDiagonal(matrix* m);
matrix* newMatrix(int rows, int columns);
matrix* copyMatrix(matrix* m);
matrix* transposeMatrix(matrix* m, int free);

/*---------------------------------------------------------------------------------------------------*/
/* ******************************************** KMeans Functions ******************************************** */
/*---------------------------------------------------------------------------------------------------*/

/* -------------------------------------------- Data Convertions Functions ------------------------------------------- */
Node* getPointsFromT(double ** TDoubleArr, int d, int numOfpoints);
cluster * getClustersFromT(double ** TDoubleArr, int finalK);
matrix* pointsToMat(Node* points, int numOfPoints, int d);
double ** matToArr(matrix * m, int free1);


/* -------------------------------------------- Kmeans Algorithm Functions ------------------------------------------- */
Node* getPoints(char* fileName, int* numOfPoints, int* finald);
double ** kmeansFunc(int k, int maxIter, int numOfPoints, int d, Node* pointsMatrix, cluster* centroids);
void doKmeans(cluster* clusters, Node* points, int d, int k, int numOfPoints, int maxIteratiorns);


/* -------------------------------------------- Kmeans Helpers Functions ------------------------------------------- */
int getMinCluster(double* point, cluster* clusters, int K, int d);
double calcDistance(double* p1 , double* p2 ,int d);
void addPointToCluster(cluster c, double* point, int d);
int checkKmeansConvergence(cluster* clusters, int k, int d);
int continueLoop( int convergence, int count, int maxIter);
void iterateClusters(cluster* clusters, int length, int d);


/*---------------------------------------------------------------------------------------------------*/
/* ******************************************** General Functions ******************************************** */
/*---------------------------------------------------------------------------------------------------*/

/* -------------------------------------------- List Functions ------------------------------------------- */
Node* addNext(Node* node, double* point);
Node* addCurrentNext(Node* node, double* point);


/* -------------------------------------------- Free Memory Functions ------------------------------------------- */
void freeMatrix(matrix * m);
void free2DArray(double ** TDoubleArr, int d, int numOfpoints);
void freeKmeansMemory(cluster* clusters, Node* points, int k, int numOfpoints);
void freePointsList(Node * points, int numOfpoints);


/* -------------------------------------------- Print Functions ------------------------------------------- */
void printMatrix(matrix* A);
void printEigenfromMat(matrix* m);
void printMatOutput(matrix *m); //TODO check if it's according to rules and replace with printMatrix
void printEigenArr(eigenVal * arr, int length);
void printNodesList(Node * points, int d);
void printLst(Node* lst, int d); // TODO: duplicate method?
void printClusters(cluster* clusters, int k, int d);
void assertFunc(void* x);

