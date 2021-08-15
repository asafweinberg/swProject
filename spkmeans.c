#include "spkmeans.h"
#include <string.h>

int main(int argc, char *argv[])
{
    //handle all goals
    int k,d, numOfPoints;
    char* myGoal;
    char* fileName;
    double ** TDoubleArr;
    int maxIter=300;
    // int maxIter = 300, numOfPoints, d;
    Node* points;
    cluster* clusters;
    

    k = atoi(argv[1]);
    myGoal = argv[2];
    fileName = argv[3];
  
  //  printf("%d",argc);
 //   printf("%d",maxIter);
    // k = 3;
    // fileName = "input_1.txt";
    // maxIter=600;
    
    TDoubleArr = runMainFlow(k,myGoal,fileName);

    if(TDoubleArr==NULL)
    {
        return 0;
    }
    
    // points = getPoints(fileName, &numOfPoints, &d);

    // printClusters(clusters, k, d);

    k=len(TDoubleArr[0]);
    d=k;
    numOfPoints=len(TDoubleArr);
    clusters=getClustersFromT(TDoubleArr);
    points=getPointsFromT(TDoubleArr);

    doKmeans(clusters,points,d,k,numOfPoints, maxIter);
    freeMemo(clusters,points,k,numOfPoints);
    return 0;
}

Node* getPointsFromT(double ** TDoubleArr)
{
    int d,i,j;
    double number;
    double* point;
    Node* points, *current;
    int numPoints=len(TDoubleArr);
    d=len(TDoubleArr[0]);

    points = (Node*)malloc(sizeof(Node));
    assert(points);
    current=points;

    for (i=0 ; i<numPoints ; i++)
    {
       point=(double*)calloc(d,sizeof(double));
       assert(point);

       for(j=0 ; j<d ; j++)
       { 
           point[j]=TDoubleArr[i][j];
       }
       current=addNext(current, point);
    }

    return points;
}

cluster * getClustersFromT(double ** TDoubleArr)
{
    cluster* clusters;
    int i,j,k,d;
    k=len(TDoubleArr[0]);
    d=k;
    clusters=(cluster *)calloc(k, sizeof(cluster));
    assert(clusters);

    for(i=0; i<k; i++)
    {
        clusters[i].mean = (double*)calloc(d, sizeof(double));
        assert(clusters[i].mean);
        clusters[i].prevMean = (double*)calloc(d,sizeof(double));
        assert(clusters[i].prevMean);
        
        for(j=0; j<d; j++)
        {   
            clusters[i].prevMean[j] = TDoubleArr[i][j];
        }
        clusters[i].size=0;
    }
    return clusters;
}

Node* getPoints(char* fileName, int* numOfPoints, int* finald) 
{
    FILE *myfile;
    double number;
    double* firstPoint = NULL;
    char c;
    int d=0;
    double* point;
    Node* points, *current;
    int numPoints=0;
    
    myfile = fopen(fileName, "r");

    fscanf(myfile, "%lf%c", &number, &c);

    //first number in first point
    fscanf(myfile, "%lf%c", &number, &c);
    firstPoint = (double*)malloc(sizeof(double));
    assert(firstPoint);
    firstPoint[0] = number;
    d++;
    if (c!='\n')
    {
        while (fscanf(myfile,"%lf%c", &number, &c) == 2)
        {
            d++;
            firstPoint = (double*)realloc(firstPoint, d*sizeof(double));
            assert(firstPoint);
            firstPoint[d-1] = number;
            if(c == '\n' || c == EOF)
            {
                break;
            }
        }
    }

    points = (Node*)malloc(sizeof(Node));
    assert(points);
    points->next=NULL;
    points->data=firstPoint;
    numPoints++;
    current=points;
    while (fscanf(myfile, "%lf%c", &number, &c) == 2)
    {
       point=(double*)calloc(d,sizeof(double));
       assert(point);
       point[0]=number;
       for(int i=1 ; i<d ; i++)
       {
           fscanf(myfile, "%lf%c", &number, &c); 
           point[i]=number;
       }
       
       current=addNext(current, point);
       numPoints++;
    }

    *numOfPoints=numPoints;
    *finald=d;
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
    assert(clusterArr);
    for (i=0; i<k; i++){
        clus = (double*)calloc(d,sizeof(double));
        assert(clus);
        for (j=0 ; j<d ; j++)
        {
            clus[j] = clusters[i].prevMean[j];
        }
        clusterArr[i] = clus;
    }

    freeMemo(clusters,points,k,numOfPoints);
    return clusterArr;
    
}
/*
Node* getPoints(int* numOfPoints, int* finald) 
{
    double number;
    double* firstPoint = NULL;
    char c;
    int d=0, i;
    double* point;
    Node* points, *current;
    int numPoints=0;
    

    scanf("%lf%c", &number, &c);
    firstPoint = (double*)calloc(1,sizeof(double));
    assert(firstPoint);
    firstPoint[0] = number;
    d++;
    if (c!='\n')
    {
        while (scanf("%lf%c", &number, &c) == 2)
        {
            d++;
            firstPoint = (double*)realloc(firstPoint, d*sizeof(double));
            assert(firstPoint);
            firstPoint[d-1] = number;
            if(c == '\n' || c == EOF)
            {
                break;
            }
        }
    }

    points = (Node*)calloc(1,sizeof(Node));
    assert(points);
    points->next=NULL;
    points->data=firstPoint;
    numPoints++;
    current=points;
    while (scanf("%lf%c", &number, &c) == 2)
    {
       point=(double*)calloc(d,sizeof(double));
       assert(point);
       point[0]=number;
       for(i=1 ; i<d ; i++)
       {
           scanf("%lf%c", &number, &c);
           point[i]=number;
       }
       
       current=addNext(current, point);
       numPoints++;
    }

    *numOfPoints=numPoints;
    *finald=d;
    return points;
}
*/
void doKmeans(cluster* clusters, Node* points, int d, int k, int numOfPoints, int maxIteratiorns)
{
    Node* head=points;
    Node* current=head;
    int convergence=false, count = 0, min, i;
    
    while(continueLoop(convergence, count, maxIteratiorns))
    {
        for (i = 0; i < numOfPoints; i++)
        {
            min = minCluster(current->data, clusters, k, d);
            addPoint(clusters[min], current->data, d);
            clusters[min].size ++;
            current = current->next;
        }
        current = head;
        count++;
        convergence = checkCon(clusters, k, d);
        iterateClusters(clusters, k, d);
    }
}


int minCluster(double* point, cluster* clusters, int K, int d){
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
            minDist=distance;
            clusterIndex=i;
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
    return dist;
}

void addPoint(cluster c, double* point, int d)
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

int checkCon(cluster* clusters, int k, int d)
{
    int i, j;
    for(i = 0; i < k; i++)
    {
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
       prev = clusters[i].prevMean;
       clusters[i].prevMean = clusters[i].mean;
        clusters[i].size = 0;
       free(prev);
    }
    for (i=0; i< length; i++)
    {
        clusters[i].mean = (double*)calloc(d, sizeof(double));
    }
}

Node* addNext(Node* node, double* point) 
{
    Node* new;
    new = (Node*)calloc(1,sizeof(Node));
    assert(new);


    new->next = NULL;
    node->data = point;

    node->next = new;
    return new;
}

void printClusters(cluster* clusters, int k, int d)
{
    int i, j;
    for(i=0 ; i<k ; i++)
    {
        for(j=0 ; j<d ; j++)
        {
            if(j!=d-1)
            {
                printf("%.4f,",clusters[i].prevMean[j]);
            }
            else
            {
              printf("%.4f\n",clusters[i].prevMean[j]);  
            }
        }
    }
}

void freeMemo(cluster* clusters, Node* points, int k, int numOfpoints)
{
    Node* iter=points ,*prevNode=NULL;
    int i;

    for(i=0 ; i<k ; i++)
    {
        free(clusters[i].mean);
        free(clusters[i].prevMean);
    }
    free(clusters);

    iter=points;
    for (i=0 ; i<numOfpoints ; i++)
    {
        free(iter->data);
        prevNode=iter;
        iter=prevNode->next;
        free(prevNode);
    }
}

void printLst(Node* lst, int d)
{
    int j;
    Node* curr = lst;

    while (curr->data != NULL)
    {
        for(int j=0 ; j<d ; j++)
        {
            if(j!=d-1)
            {
                printf("%.4f,", curr->data[j]);
            }
            else
            {
              printf("%.4f\n", curr->data[j]);  
            }
        }
        curr = curr->next;
    }
    
    
}

/*
cluster* makeClusters(int k, int d , PyObject ** centroids)
{
    cluster* clusters;
    PyListObject *item;
    PyFloatObject *pointDouble;
    int i,j;

    clusters=(cluster *)calloc(k, sizeof(cluster));
    assert(clusters);

    for(i=0; i<k; i++)
    {
        clusters[i].mean = (double*)calloc(d, sizeof(double));
        assert(clusters[i].mean);
        clusters[i].prevMean = (double*)calloc(d,sizeof(double));
        assert(clusters[i].prevMean);
        
        item = (PyListObject *) PyList_GetItem((PyObject *)centroids, i);
        for(j=0; j<d; j++)
        {   
            pointDouble = (PyFloatObject *) PyList_GET_ITEM(item,j);
            double value = PyFloat_AS_DOUBLE(pointDouble);
            clusters[i].prevMean[j] = value;
        }
        clusters[i].size=0;
    }
    return clusters;
}
*/
cluster* makeClustersSp(int k, Node* points, int numOfPoints,int d)
{
    // assumption: k is already calculated correctly if it was 0 before

}

matrix* formMatW(Node* points, int numOfPoints, int d)
{
    matrix *Wmat, *pointsMat;
    int i, j;
    Node* current;
    double dis;

    current = points;
    Wmat = newMatrix(numOfPoints, numOfPoints);
    pointsMat = pointsToMat(points, numOfPoints, d);

    for (i = 0; i < numOfPoints; i++)
    {
        for(j = 0 ; j < numOfPoints ; j++)
        {
            if(i != j)
            {
                dis = pow(calcDistance(pointsMat->data[i], pointsMat->data[j], d), 0.5);
                Wmat->data[i][j] = exp(-dis/2);
            }
        }
    }

    return Wmat;
}

matrix* formMatD(matrix* matW, int freeW)
{
    int i, j, l;
    matrix* matD;
    double sumW;

    matD=calloc(1,sizeof(matrix));
    matD->rows=matW->rows;
    matD->columns=matW->columns;
    matD->data=calloc(matW->rows, sizeof(double*));

    for (i = 0; i < matW->rows; i++)
    {
        matD->data[i] = calloc(matW->columns, sizeof(double));

        for (j = 0 ; j<matW->columns ; j++)
        {
            if (i==j)
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
        free(matW);
    }

    return matD;
}

matrix* formMatLnorm(matrix* matD , matrix* matW , int freeD, int freeW)
{
    matrix *lnormMat, *IMat, *rootD;

    IMat = formMatI(matD->rows);
    rootD = minusRootMat(matD, false);

    lnormMat = rootD;
    lnormMat = mulMatrices(lnormMat, matW, false, true); //dont free lnorm cause it points to rootD
    lnormMat = mulMatrices(lnormMat, rootD, true, true);
    lnormMat = addMatrices(IMat, lnormMat ,true, true, true);

    if(freeD)
    {
        free(matD);
    }
    if(freeW)
    {
        free(matW);
    }


    return lnormMat;
}

//only for diagonal matrix
matrix* minusRootMat(matrix* mat, int free1)
{
    int i,j;
    matrix* newMat;

    newMat = newMatrix(mat->rows, mat->columns);
    
    for(i = 0; i< mat->rows; i++)
    {   
        newMat->data[i][i] = pow(mat->data[i][i],-0.5);
    }

    if(free1)
    {
        free(mat);
    }
    return newMat;
}

matrix* formMatI(int dimention)
{
    int i,j;
    matrix* matI;

    matI = newMatrix(dimention, dimention);

    for(i = 0 ; i < dimention ; i++)
    {   
        for(j = 0 ; j< dimention ; j++)
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
            pointsMat->data[i][j] = current->data[j];
        }
        current = current->next;
    }
    return pointsMat;
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

matrix* addMatrices(matrix* mat1, matrix* mat2, int dec, int free1, int free2) // dec=1 -> dec 
{
    int rows, columns, i, j, sign;
    matrix* sumMat;

    assert((mat1->rows == mat2->rows) && (mat1->columns == mat2->columns));

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


void jacobiAlg(matrix* lNorm, matrix** eigenValues, matrix** eigenVectors)
{
    int iterations = 0, convergence = 0, rowPivot, colPivot;
    double c, s;
    matrix *matA, *matB, *pivot, *tempVectors;

    matA = lNorm;
    matB = newMatrix(matA -> rows, matA -> columns);
    tempVectors = formMatI(matA -> rows);
    // printMatrix(matA);

    while (!isDiagonal(matA) && ++iterations <= MaxJacobiIter && !convergence)
    {
        pivot = getRotationMatrixValues(matA, &c, &s, &rowPivot, &colPivot);
        // printMatrix(pivot);
        matB = calcNextJacobiMatrix(matA, c, s, rowPivot, colPivot); //TODO, make sure that matB changes
        tempVectors = mulMatrices(tempVectors, pivot, true, true);
        convergence = hasConvergence(matA, matB);
        freeMatrix(matA);
        matA = matB;
        // printMatrix(matA);
    }

    *eigenValues = matA;
    *eigenVectors = tempVectors;
}

matrix* getRotationMatrixValues(matrix* mat, double *c, double *s, int *rowPivot, int *colPivot)
{
    matrix *p;
    int signTheta;
    double **m, theta, t;

    m = mat -> data;

    p = newMatrix(mat -> rows, mat -> columns);
    findLargestCell(mat, rowPivot, colPivot);
    
    //TODO: Handle zeros of theta or cell and others

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

int hasConvergence(matrix* matA, matrix* matB)
{
    return (calcOff(matA) - calcOff(matB)) <= epsilon;
}

double calcOff(matrix* m)
{
    // improving running time possibility by multipling the sum of half by 2 because of simetry
    int i, j;
    double off = 0;
    for ( i = 0; i < m ->rows; i++)
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

    vectorsIndices = eigenGapHeuristic(eigenValues, &k);
    finalVectors = newMatrix(eigenVectors->rows, k);

    for (j = 0; j < k; j++)
    {
        for (i = 0; i < eigenVectors->rows; i++)
        {
            finalVectors->data[i][j] = eigenVectors->data[i][vectorsIndices[j]];
        }
    }        
    return finalVectors;
}

//TODO: handle initial k
// returns an array of the final starting vectors
int* eigenGapHeuristic(matrix* matA, int* k)
{
    eigenVal* values;
    int length, i, *vectorsIndices;

    length = matA->columns;
    values = (eigenVal*)calloc(length, sizeof(eigenVal));
    assert(values);

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
    for (i = 0; i < *k; i++)
    {
        vectorsIndices[i] = values[i].column;
    }

    //TODO: free values function
    free(values);

    return vectorsIndices;
}

int findMaxGap(eigenVal* values, int length)
{
    double currDiff, maxDiff = -1;
    int i, k;
    for (i = 0; i < (length-1) / 2; i++)
    {
        currDiff = fabs(values[i].value - values[i + 1].value); //TODO: check the indices
        if(currDiff > maxDiff)
        {
            maxDiff = currDiff;
            k = i;
        }
    }
    return k + 1;
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
    int i, j;

    m = calloc(1, sizeof(matrix));
    assert(m);
    m->rows = rows;
    m->columns = columns;
    m->data = calloc(rows,sizeof(double*));
    assert(m->data);

    for(i = 0; i < rows; i++)
    {
        (m->data)[i] = calloc(columns,sizeof(double));
        assert((m->data)[i]);
    }
    return m;
}

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

void printMatrix(matrix* A) {
    int i, j;
    printf("===================\n");
    for (i = 0; i < A->rows; i++)
    {
        for (j = 0; j < A->columns; j++)
        {
            printf("%.4f    ", A->data[i][j]);
        }
        printf("\n");
    }
}

int compareEigenVal(const void * a, const void * b)
{
    return ( ((eigenVal*)a)->value - ((eigenVal*)b)->value );
}

double ** runSpk(int k, Node* points, int numOfPoints, int d)
{   
    double** TDoubleArr;
    // cluster* clusters;

    if (numOfPoints <= k && k != 0)
    {
        printf("ERROR K>=N");
        assert(0);
    }

//TODO: call getLnorm and then call getTmatFromLnorm and return 

    return TDoubleArr; // ADD THE REST OF THE FLOW
    // clusters = makeClustersSp(k,points,numOfPoints,d);
    
}

void runWam(Node* points, int numOfPoints, int d)
{
    matrix * m = formMatW(points,numOfPoints,d);
    printMatrix(m);
}

void runDdG(Node* points, int numOfPoints, int d)
{
    matrix * m = formMatW(points,numOfPoints,d);
    m=formMatD(m,true); //may cause error
    printMatrix(m);
}

void runLnorm(Node* points, int numOfPoints, int d)
{
    matrix * m1,*m2,*m3;
    m1 = formMatW(points,numOfPoints,d);
    m2 = formMatD(m1,false);
    m3 = formMatLnorm(m2,m1,true,true);
    printMatrix(m3);
}

void runJacobi(Node* points, int numOfPoints, int d)
{
   //TODO
}

matrix* getTmatFromLnorm(matrix* lNorm, int k)
{
//TODO
}

double ** runMainFlow(int k, char* myGoal, char* fileName) //return T or NULL
{
    int maxIter = 300, numOfPoints, d;
    Node* points;
    cluster* clusters;

    points = getPoints(fileName, &numOfPoints, &d);

    if(!strcmp(myGoal,"spk")) //returns T
    {
        return runSpk(k,points,numOfPoints,d); //returns double **
    }


    if(!strcmp(myGoal,"wam"))
    {
        runWam(points,numOfPoints,d);
        return NULL;
    }

    if(!strcmp(myGoal,"ddg"))
    {
        runDdG(points,numOfPoints,d);
        return NULL;
    }

    if(!strcmp(myGoal,"lnorm"))
    {
        runLnorm(points,numOfPoints,d);
        return NULL;

    }

    if(!strcmp(myGoal,"jacobi"))
    {
        runJacobi(points,numOfPoints,d);
        return NULL;
    }
}


