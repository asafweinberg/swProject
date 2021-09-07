#include "spkmeans.h"
#include <string.h>

/*-----------------------------------------------------------------------------*/
/* ******************************************** SPK Functions ******************************************** */
/*-----------------------------------------------------------------------------*/

/* -------------------------------------------- Run Flow Functions -------------------------------------------- */
void runWam(Node* points, int numOfPoints, int d)
{
    matrix * m = formMatW(points,numOfPoints,d);
    printMatrix(m);
    freeMatrix(m);
    freePointsList(points,numOfPoints);
}


void runDdG(Node* points, int numOfPoints, int d)
{
    matrix * m;
    m = formMatW(points,numOfPoints,d);
    m = formMatD(m, true);
    printMatrix(m);
    freeMatrix(m);
    freePointsList(points, numOfPoints);
}


matrix * runLnorm(Node* points, int numOfPoints, int d, int printMat0)
{
    matrix * m1, *m2, *m3;
    m1 = formMatW(points, numOfPoints, d);
    m2 = formMatD(m1, false);
    m3 = formMatLnorm(m2, m1, true, true);
    
    if (printMat0)
    {
        printMatrix(m3);
        freeMatrix(m3);
        freePointsList(points,numOfPoints);
        return NULL;
    }

    return m3;
}


matrix * runSpk(int k, Node* points, int numOfPoints, int d)
{   
    matrix * eigenValues, *eigenVectors;
    matrix * lNorm, * matT, *matU;

    if (numOfPoints <= k && k != 0)
    {
        printf("ERROR K>=N");
        assert(0);
    }

    lNorm = runLnorm(points, numOfPoints, d, false);
    jacobiAlg(lNorm, &eigenValues, &eigenVectors);
    matU = calcInitialVectorsFromJacobi(eigenValues, eigenVectors, k); 
    matT = formMatT(matU, true);
    freePointsList(points,numOfPoints);
    return matT; 
    
}


void runJacobi(Node* points, int numOfPoints, int d)
{
    matrix * eigenValues, *eigenVectors, *vectorsToPrint;
    matrix * mat;

    mat = pointsToMat(points, numOfPoints, d);
    jacobiAlg(mat, &eigenValues, &eigenVectors);

    printEigenfromMat(eigenValues);
    vectorsToPrint = transposeMatrix(eigenVectors, true);
    printMatrix(vectorsToPrint);

    freeMatrix(eigenValues);
    freeMatrix(eigenVectors); 
    freePointsList(points,numOfPoints);
}


double ** runMainFlow(int k, char* myGoal, char* fileName, int* finalK, int* numOfPoints)
{
    int d;
    Node* points;
    matrix* T;
    double** TDoubleArr;

    points = getPoints(fileName, numOfPoints, &d);

    if(!strcmp(myGoal, "spk"))
    {
        T = runSpk(k, points, *numOfPoints, d);
        *finalK = T->columns;
        *numOfPoints = T->rows;
        TDoubleArr = matToArr(T, true);
        return TDoubleArr;
    }

    if(!strcmp(myGoal, "wam"))
    {
        runWam(points, *numOfPoints, d);
        return NULL;
    }

    if(!strcmp(myGoal, "ddg"))
    {
        runDdG(points, *numOfPoints, d);
        return NULL;
    }

    if(!strcmp(myGoal, "lnorm"))
    {
        runLnorm(points, *numOfPoints, d, true);
        return NULL;
    }

    if(!strcmp(myGoal, "jacobi"))
    {
        runJacobi(points, *numOfPoints, d);
        return NULL;
    }

    return NULL;
}

/* -------------------------------------------- goals Functions -------------------------------------------- */

matrix* formMatW(Node* points, int numOfPoints, int d)
{
    matrix *Wmat, *pointsMat;
    int i, j;
    double dis;

    Wmat = newMatrix(numOfPoints, numOfPoints);
    pointsMat = pointsToMat(points, numOfPoints, d);

    for (i = 0; i < numOfPoints; i++)
    {
        for(j = 0 ; j < numOfPoints ; j++)
        {
            if(i != j)
            {
                dis = calcDistance(pointsMat->data[i], pointsMat->data[j], d);
                Wmat->data[i][j] = exp(-dis/2);
            }
        }
    }
    freeMatrix(pointsMat);
    return Wmat;
}


matrix* formMatT(matrix * matU, int free1)
{
    int i,j,l;
    double sumSquare;
    matrix * T = newMatrix(matU->rows , matU->columns);
    for (i = 0; i < T->rows; i++)
    {
        sumSquare = 0;
        for (l = 0; l < T->columns; l++)
        {
            sumSquare+=(matU->data[i][l]) * (matU->data[i][l]);
        }
        sumSquare = pow(sumSquare, 0.5);

        for (j = 0; j < T->columns; j++)
        {
            if (matU->data[i][j] == 0.0)
            {
                T->data[i][j] = 0.0;
            }
            else
            {
                T->data[i][j] = (matU->data[i][j]) / sumSquare;
            }
        }
    }

    if (free1)
    {
        freeMatrix(matU);
    }
    return T;
}


matrix* formMatD(matrix* matW, int freeW)
{
    int i, j, l;
    matrix* matD;
    double sumW;

    matD = newMatrix(matW->rows, matW->columns);
    
    for (i = 0; i < matD->rows; i++)
    {
        for (j = 0; j < matD->columns; j++)
        {
            if (i == j)
            {
                sumW = 0;
                for (l = 0; l < matW->columns; l++)
                {
                   sumW += matW->data[j][l];
                }
                matD->data[i][j] = sumW;
            }
            else
            {
                matD->data[i][j] = 0;
            }
        }
    }

    if(freeW)
    {
        freeMatrix(matW);
    }

    return matD;
}


matrix* formMatLnorm(matrix* matD , matrix* matW , int freeD, int freeW)
{
    matrix *lnormMat, *IMat, *rootD;

    IMat = formMatI(matD->rows);
    rootD = minusRootMat(matD, false);

    lnormMat = mulMatrices(rootD, matW, false, false);
    lnormMat = mulMatrices(lnormMat, rootD, true, true);
    lnormMat = addMatrices(IMat, lnormMat ,true, true, true);

    if(freeD)
    {
        freeMatrix(matD);
    }
    if(freeW)
    {
        freeMatrix(matW);
    }


    return lnormMat;
}

/* -------------------------------------------- jacobi Functions ------------------------------------------- */

void jacobiAlg(matrix* lNorm, matrix** eigenValues, matrix** eigenVectors)
{
    int iterations = 0, convergence = 0, rowPivot, colPivot;
    double c, s;
    matrix *matA, *matB, *pivot, *tempVectors;

    matA = lNorm;
    tempVectors = formMatI(matA -> rows);

    while (!isDiagonal(matA) && ++iterations <= MaxJacobiIter && !convergence) 
    {
        pivot = getRotationMatrixValues(matA, &c, &s, &rowPivot, &colPivot);
        matB = calcNextJacobiMatrix(matA, c, s, rowPivot, colPivot);

        tempVectors = mulMatrices(tempVectors, pivot, true, true);

        convergence = hasJacobiConvergence(matA, matB);
        freeMatrix(matA);
        matA = matB;
    }

    *eigenValues = matA;
    *eigenVectors = tempVectors;
}


matrix* getRotationMatrixValues(matrix* mat, double *c, double *s, int *rowPivot, int *colPivot)
{
    int signTheta;
    double **m, theta, t;

    m = mat -> data;

    findLargestCell(mat, rowPivot, colPivot);
    
    theta = (m[*colPivot][*colPivot] - m[*rowPivot][*rowPivot]) / (2 * m[*rowPivot][*colPivot]);
    signTheta = theta < 0? -1 : 1;
    t = signTheta / (fabs(theta) + pow((theta * theta + 1), 0.5));
    *c = 1 / pow((t * t + 1), 0.5);
    *s = (*c) * t;

    return setPivotMatrix(*rowPivot, *colPivot, *c, *s, mat -> rows);
}


void findLargestCell(matrix* mat, int *row, int *col) 
{
    int i, j;
    double big = 0, **m;

    m = mat -> data;

    for (i = 0; i < mat -> rows; i++)
    {
        for (j = 0; j < mat -> columns; j++)
        {
            if (i == j)
                continue;
            if (fabs(m[i][j]) > big) 
            {
                *row = i;
                *col = j;
                big = fabs(m[i][j]);
            }
        }
    }
}


matrix* setPivotMatrix(int row, int col, double c, double s, int dim)
{
    matrix* pivot;

    pivot = formMatI(dim);
    pivot->data[row][row] = c;
    pivot->data[row][col] = s;
    pivot->data[col][row] = -s;
    pivot->data[col][col] = c;

    return pivot;
}


matrix* calcNextJacobiMatrix(matrix* matA, double c, double s, int i, int j)
{
    int r;
    double **a, **b;
    matrix* matB;

    matB = copyMatrix(matA);

    a = matA -> data;
    b = matB -> data;

    for (r = 0; r < matA -> rows; r++)
    {
        if(r != i && r!= j)
        {
            b[r][i] = c*a[r][i] - s*a[r][j];
            b[i][r] = c*a[r][i] - s*a[r][j];

            b[r][j] = c*a[r][j] + s*a[r][i];
            b[j][r] = c*a[r][j] + s*a[r][i];
        }
    }

    b[i][i] = c*c*a[i][i] + s*s*a[j][j] - 2*s*c*a[i][j];
    b[j][j] = s*s*a[i][i] + c*c*a[j][j] + 2*s*c*a[i][j];
    b[i][j] = 0;
    b[j][i] = 0;

    return matB;
}


int hasJacobiConvergence(matrix* matA, matrix* matB)
{
    return fabs(calcOff(matA) - calcOff(matB)) <= epsilon;
}


double calcOff(matrix* m)
{
    int i, j;
    double off = 0.0;
    for (i = 0; i < m ->rows; i++)
    {
        for (j = 0; j < m->columns; j++)
        {
            if(i != j)
            {
                off += m->data[i][j];
            }
        }
        
    }
    return off;
}


matrix* calcInitialVectorsFromJacobi(matrix* eigenValues, matrix* eigenVectors, int initialK)
{
    int *vectorsIndices, k, i, j;
    matrix* finalVectors;

    k = initialK;

    vectorsIndices = eigenGapHeuristic(eigenValues, &k);
    finalVectors = newMatrix(eigenVectors->rows, k);

    for (j = 0; j < k; j++)
    {
        for (i = 0; i < eigenVectors->rows; i++)
        {
            finalVectors->data[i][j] = eigenVectors->data[i][vectorsIndices[j]];
        }
    }

    free(vectorsIndices);
    freeMatrix(eigenValues);
    freeMatrix(eigenVectors);

    return finalVectors;
}


int* eigenGapHeuristic(matrix* matA, int* k)
{
    eigenVal* values;
    int length, i, *vectorsIndices;

    length = matA->columns;
    values = (eigenVal*)calloc(length, sizeof(eigenVal));
    assertFunc(values);

    for (i = 0; i < length; i++)
    {
        values[i].column = i;
        values[i].value = matA->data[i][i];
    }

    qsort(values, length, sizeof(eigenVal), compareEigenVal);

    if (*k == 0)
    {
        *k = findMaxGap(values, length);
    }

    vectorsIndices = (int*)calloc(*k, sizeof(int));
    assertFunc(vectorsIndices);
    for (i = 0; i < *k; i++)
    {
        vectorsIndices[i] = values[i].column;
    }

    free(values);

    return vectorsIndices;
}


int findMaxGap(eigenVal* values, int length)
{
    double currDiff = 0.0, maxDiff = -1.0;
    int i, k = 0;

    for (i = 0; i < floor(length / 2); i++)
    {
        currDiff = fabs(values[i].value - values[i + 1].value);
        if(currDiff > maxDiff)
        {
            maxDiff = currDiff;
            k = i;
        }
    }
    return k + 1;
}


int compareEigenVal(const void * a, const void * b)
{
    eigenVal* aPtr , *bPtr;
    aPtr = (eigenVal*)a;
    bPtr = (eigenVal*)b;

    if (aPtr -> value == bPtr -> value )
    {
        return aPtr->column - bPtr->column;
    }

    return aPtr -> value < bPtr ->value ? -1 : 1 ;
}

/* -------------------------------------------- matrices Functions ------------------------------------------- */

matrix* minusRootMat(matrix* mat, int free1)
{
    int i;
    matrix* newMat;

    newMat = newMatrix(mat->rows, mat->columns);
    
    for(i = 0; i< mat->rows; i++)
    {   
        newMat->data[i][i] = pow(mat->data[i][i],-0.5);
    }

    if(free1)
    {
        freeMatrix(mat);
    }
    return newMat;
}


matrix* mulMatrices(matrix* mat1, matrix* mat2, int free1, int free2)
{
    int rows, columns, mat1Columns, i, j, l;
    matrix* mulMat;
    double sumElement;

    rows = mat1->rows;
    columns = mat2->columns;
    mat1Columns = mat1->columns;

    mulMat = newMatrix(rows, columns);
    
    for(i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            sumElement = 0;
            for(l = 0; l < mat1Columns; l++)
            {
                sumElement += (mat1->data[i][l])*(mat2->data[l][j]);
            }
            mulMat->data[i][j] = sumElement;
        }
    }

    if (free1)
    {
        freeMatrix(mat1);
    }

    if (free2)
    {
        freeMatrix(mat2);
    }
    return mulMat;
}


matrix* addMatrices(matrix* mat1, matrix* mat2, int dec, int free1, int free2)
{
    int rows, columns, i, j, sign;
    matrix* sumMat;

    rows = mat1->rows;
    columns = mat1->columns;
    sumMat = newMatrix(rows,columns);
    sign = dec? -1: 1;

    for(i = 0 ; i < rows ; i++)
    {
        for(j = 0; j < columns; j++)
        {
            sumMat->data[i][j] = (mat1->data)[i][j] + (sign * (mat2->data)[i][j]);
        }
    }

    if (free1)
    {
        freeMatrix(mat1);
    }

    if (free2)
    {
        freeMatrix(mat2);
    }
    return sumMat;

}


matrix* formMatI(int dimention)
{
    int i,j;
    matrix* matI;

    matI = newMatrix(dimention, dimention);

    for (i = 0 ; i < dimention ; i++)
    {   
        for (j = 0 ; j< dimention ; j++)
        {
            if (i == j)
            {
                matI->data[i][j] = 1;
            }
            else
            {
                matI->data[i][j] = 0;
            }
        }
    }
    return matI;
}


int isDiagonal(matrix* m)
{
    int i, j;
    for ( i = 0; i < m ->rows; i++)
    {
        for (j = 0; j < m->columns; j++)
        {
            if(i != j)
            {
                if (m->data[i][j] != 0)
                    return false;
            }
        }
        
    }
    return true;
}


matrix* newMatrix(int rows, int columns)
{
    matrix * m;
    int i;

    m = calloc(1, sizeof(matrix));
    assertFunc(m);
    m->rows = rows;
    m->columns = columns;
    m->data = calloc(rows,sizeof(double*));
    assertFunc(m->data);

    for(i = 0; i < rows; i++)
    {
        (m->data)[i] = calloc(columns,sizeof(double));
        assertFunc((m->data)[i]);
    }
    return m;
}


matrix* copyMatrix(matrix* m)
{
    matrix* newMat;
    int i,j;

    newMat = newMatrix(m->rows, m->columns);

    for (i = 0; i < m -> rows; i++)
    {
        for (j = 0; j < m->columns; j++)
        {
            newMat->data[i][j] = m->data[i][j];
        }
        
    }
    return newMat;
}

matrix* transposeMatrix(matrix* m, int free) 
{
    matrix* newMat;
    int i,j;

    newMat = newMatrix(m->rows, m->columns);

    for (i = 0; i < m -> rows; i++)
    {
        for (j = 0; j < m->columns; j++)
        {
            newMat->data[j][i] = m->data[i][j];
        }
        
    }

    if (free) 
        freeMatrix(m);

    return newMat;
}

/*---------------------------------------------------------------------------------------------------*/
/* ******************************************** KMeans Functions ******************************************** */
/*---------------------------------------------------------------------------------------------------*/

/* -------------------------------------------- Data Convertions Functions ------------------------------------------- */

Node* getPointsFromT(double ** TDoubleArr, int d, int numOfpoints)
{
    int i,j;
    double* point;
    Node* points, *current;
   

    points = (Node*)malloc(sizeof(Node));
    assertFunc(points);
    current = points;

    for (i=0 ; i<numOfpoints ; i++)
    {
       point = (double*)calloc(d,sizeof(double));
        assertFunc(point);

       for(j = 0; j < d; j++)
       { 
           point[j] = TDoubleArr[i][j];
       }
       current = addCurrentNext(current, point);
    }
    free2DArray(TDoubleArr, numOfpoints);
    return points;
}


cluster * getClustersFromT(double ** TDoubleArr, int finalK)
{
    cluster* clusters;
    int i, j, d;
    d = finalK;
    clusters = (cluster *)calloc(finalK, sizeof(cluster));
    assertFunc(clusters);

    for(i = 0; i < finalK; i++)
    {
        clusters[i].mean = (double*)calloc(d, sizeof(double));
        assertFunc(clusters[i].mean);
        clusters[i].prevMean = (double*)calloc(d, sizeof(double));
        assertFunc(clusters[i].prevMean);
        
        for(j = 0; j < d; j++)
        {   
            clusters[i].prevMean[j] = TDoubleArr[i][j];
        }
        clusters[i].size = 0;
    }
    return clusters;
}


matrix* pointsToMat(Node* points, int numOfPoints, int d)
{
    int i, j;
    matrix* pointsMat;
    Node* current;

    pointsMat = newMatrix(numOfPoints, d);
    current = points;

    for(i = 0 ; i < numOfPoints ; i++)
    {   
        for(j = 0 ; j < d ; j++)
        {
            pointsMat->data[i][j] = (current->data[j]);
        }
        current = current->next;
    }
    return pointsMat;
}


double ** matToArr(matrix * m, int free1)
{
    double ** arr;
    int i,j;

    arr = calloc(m->rows, sizeof(double*));
    assertFunc(arr);
    for(i = 0; i < m->rows; i++)
    {
        arr[i] = calloc(m->columns, sizeof(double));
        assertFunc(arr[i]);
        for(j = 0; j < m->columns; j++)
        {
            arr[i][j] = m->data[i][j];
        }
    }
    if(free1)
    {
        freeMatrix(m);
    }
    return arr;
}

/* -------------------------------------------- Kmeans Algorithm Functions ------------------------------------------- */

Node* getPoints(char* fileName, int* numOfPoints, int* finald) 
{
    FILE *myfile;
    double number;
    double* firstPoint = NULL;
    char c;
    int d=0, i;
    double* point;
    Node* points, *current;
    int numPoints = 0;
    
    myfile = fopen(fileName, "r");


    fscanf(myfile, "%lf%c", &number, &c);
    firstPoint = (double*)malloc(sizeof(double));
    assertFunc(firstPoint);
    firstPoint[0] = number;
    d++;
    if (c!='\n')
    {
        while (fscanf(myfile,"%lf%c", &number, &c) == 2)
        {
            d++;
            firstPoint = (double*)realloc(firstPoint, d*sizeof(double));
            assertFunc(firstPoint);
            firstPoint[d-1] = number;
            if(c == '\n' || c == EOF)
            {
                break;
            }
        }
    }

    points = (Node*)malloc(sizeof(Node));
    assertFunc(points);
    points->next = NULL;
    points->data = firstPoint;
    numPoints++;
    current = points;
    while (fscanf(myfile, "%lf%c", &number, &c) == 2)
    {
       point = (double*)calloc(d,sizeof(double));
       assertFunc(point);
       point[0] = number;
       for(i = 1; i < d; i++)
       {
           fscanf(myfile, "%lf%c", &number, &c); 
           point[i]=number;
       }
       
       current = addNext(current, point);
       numPoints++;
    }

    fclose(myfile);

    *numOfPoints = numPoints;
    *finald = d;
    return points;
}

double ** kmeansFunc(int k, int maxIter, int numOfPoints, int d, Node* pointsMatrix, cluster* centroids)
{
    Node* points;
    cluster* clusters;
    double ** clusterArr;
    int i,j;
    double* clus;
    
    points = pointsMatrix;

    clusters = centroids;
    doKmeans(clusters,points,d,k,numOfPoints, maxIter);

    clusterArr = (double**)calloc(k,sizeof(double*));
    assertFunc(clusterArr);
    for (i = 0; i < k; i++){
        clus = (double*)calloc(d,sizeof(double));
        assertFunc(clus);
        for (j = 0; j < d; j++)
        {
            clus[j] = clusters[i].prevMean[j];
        }
        clusterArr[i] = clus;
    }

    freeKmeansMemory(clusters,points,k,numOfPoints);
    return clusterArr;
    
}

void doKmeans(cluster* clusters, Node* points, int d, int k, int numOfPoints, int maxIteratiorns)
{
    Node* head=points;
    Node* current=head;
    int convergence=false, count = 0, min, i;
    
    while(continueLoop(convergence, count, maxIteratiorns))
    {   
        for (i = 0; i < numOfPoints; i++)
        {
            min = getMinCluster(current->data, clusters, k, d);
            addPointToCluster(clusters[min], current->data, d);
            clusters[min].size ++;
            current = current->next;
        }
        current = head;
        count++;
        convergence = checkKmeansConvergence(clusters, k, d);
        iterateClusters(clusters, k, d);
    }
}

/* -------------------------------------------- Kmeans Helpers Functions ------------------------------------------- */

int getMinCluster(double* point, cluster* clusters, int K, int d)
{
    int clusterIndex = 0, i;
    double minDist, distance;
    for (i = 0 ; i < K ; i++)
    {
        distance = calcDistance(clusters[i].prevMean, point, d);
        if(i == 0)
        {
            minDist = distance;
        }
        if (distance < minDist)
        {
            minDist = distance;
            clusterIndex = i;
        }
    }
    return clusterIndex;
}


double calcDistance(double* p1 , double* p2 ,int d)
{
    double dist=0;
    int i;
  
    for (i=0 ; i<d ; i++)
    {
        dist+=(p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    dist=pow(dist,0.5);
    return dist;
}


void addPointToCluster(cluster c, double* point, int d)
{
    int i;
    double new, currU;

    for(i=0 ; i<d ; i++)
    {
        currU = c.mean[i];
        new = ((currU*(c.size))+point[i])/((c.size)+1);
        c.mean[i] = new;
    }
}


int checkKmeansConvergence(cluster* clusters, int k, int d)
{
    int i, j;
    for(i = 0; i < k; i++)
    {
        if (clusters[i].size == 0)
            continue;
        for(j = 0; j < d; j++)
        {
            if (clusters[i].mean[j] != clusters[i].prevMean[j])
            {
                return false;
            }
        }
    }

    return true;
}


int continueLoop( int convergence, int count, int maxIter)
{
    return !convergence && count < maxIter; 
}


void iterateClusters(cluster* clusters, int length, int d)
{
    double* prev;
    int i;

    for (i = 0; i < length; i++)
    {
       if (clusters[i].size != 0) 
       {
            prev = clusters[i].prevMean;
            clusters[i].prevMean = clusters[i].mean;
            free(prev);
       }
       clusters[i].size = 0;
    }
    for (i = 0; i < length; i++)
    {
        clusters[i].mean = (double*)calloc(d, sizeof(double));
        assertFunc(clusters[i].mean);
    }
}

/*---------------------------------------------------------------------------------------------------*/
/* ******************************************** General Functions ******************************************** */
/*---------------------------------------------------------------------------------------------------*/

/* -------------------------------------------- List Functions ------------------------------------------- */

Node* addNext(Node* node, double* point) 
{
    Node* new;
    new = (Node*)calloc(1, sizeof(Node));
    assertFunc(new);


    new->next = NULL;
    new->data = point;

    node->next = new;
    return new;
}

Node* addCurrentNext(Node* node, double* point) 
{
    Node* new;
    new = (Node*)calloc(1,sizeof(Node));
    assertFunc(new);


    new->next = NULL;
    node->data = point; 

    node->next = new;
    return new;
}

/* -------------------------------------------- Free Memory Functions ------------------------------------------- */

void freeMatrix(matrix * m)
{
    int i;
    for(i = 0; i < m->rows; i++)
    {
        free(m->data[i]);
    }
    free(m->data);
    free(m);
}


void free2DArray(double ** TDoubleArr, int numOfpoints)
{
    int i;
    for (i = 0; i < numOfpoints; i++)
    {
        free(TDoubleArr[i]);
    }
    free(TDoubleArr);
}


void freeKmeansMemory(cluster* clusters, Node* points, int k, int numOfpoints)
{
    int i;

    for(i = 0 ; i < k ; i++)
    {
        free(clusters[i].mean);
        free(clusters[i].prevMean);
    }
    free(clusters);

    freePointsList(points, numOfpoints);
}


void freePointsList(Node * points, int numOfpoints)
{
    Node* iter,*prevNode;
    int i;

    iter = points;
    prevNode = NULL;

    for (i = 0; i < numOfpoints; i++)
    {
        free(iter->data);
        prevNode = iter;
        iter = prevNode->next;
        free(prevNode);
    }
    free(iter);
}

/* -------------------------------------------- Print Functions ------------------------------------------- */

void printMatrix(matrix* A)
{
    int i, j;
    for (i = 0; i < A->rows; i++)
    {
        for (j = 0; j < A->columns; j++)
        {
            if (j != A->columns-1)
                printf("%.4f,", A->data[i][j]);
            else
                printf("%.4f", A->data[i][j]);
        }
        if (i != A->rows-1) 
            printf("\n");
    }
}


void printEigenfromMat(matrix* m)
{
    int i;
    for(i=0 ; i<m->rows-1 ; i++)
    {
        printf("%.4f,",m->data[i][i]);  
    }
    printf("%.4f\n",m->data[i][i]);
}


void printEigenArr(eigenVal * arr, int length)
{
    int u;
    for (u=0 ; u<length ; u++)
    {
        printf("%.4f,",arr[u].value);
    }
    printf("\n");
}


void printClusters(cluster* clusters, int k, int d)
{
    int i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (j != d-1)
            {
                printf("%.4f,",clusters[i].prevMean[j]);
            }
            else
            {
              printf("%.4f",clusters[i].prevMean[j]);  
            }
        }
        if (i != k-1) 
            printf("\n");
    }
}

void assertFunc(void* x)
{
    if (!x)
    {
        printf("An Error Has Occured");
        exit(true);
    }
}


int main(int argc, char *argv[])
{   
    int k, d, numOfPoints, maxIter = 300, finalK;
    char *myGoal, *fileName;
    double ** TDoubleArr;
    Node* points;
    cluster* clusters;

    if (argc != 4)
    {
        assertFunc(NULL);
    }

    k = atoi(argv[1]);
    myGoal = argv[2];
    fileName = argv[3];

    TDoubleArr = runMainFlow(k, myGoal, fileName, &finalK, &numOfPoints);

    if(TDoubleArr == NULL)
    {
        return 0;
    }

    d = finalK;
    clusters = getClustersFromT(TDoubleArr,finalK);
    points = getPointsFromT(TDoubleArr,d,numOfPoints);

    doKmeans(clusters,points,d,finalK,numOfPoints, maxIter);
    printClusters(clusters, finalK, d);
    freeKmeansMemory(clusters,points,finalK,numOfPoints);
    return 0;
}





