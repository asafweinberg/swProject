
#include "spkmeans.h"


void testPrint(char arr[]);
double ** kmeansFunc(int k, int maxIter, int numOfPoints, int d, Node* pointsMatrix, cluster* centroids);
Node* getPoints(int* numOfPoints, int* finald); 
cluster* makeClusters(int k, int d, PyObject ** centroids);
Node* getPointsMatrix(int d, int n, PyObject ** pointsMatrix);
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


int main(int argc, char *argv[])
{
    //handle all goals
    int k , maxIter = 300, numOfPoints, d;
    char* fileName;
    Node* points;
    cluster* clusters;
    enum goal myGoal;

    k = atoi(argv[1]);
    myGoal = atoi(argv[2]);
    fileName = argv[3];
  
  //  printf("%d",argc);
 //   printf("%d",maxIter);
    // k = 3;
    // fileName = "input_1.txt";
    // maxIter=600;
    
    points = getPoints(fileName, &numOfPoints, &d);

    if (numOfPoints <= k && k != 0)
    {
        printf("ERROR K>=N");
        assert(0);
    }

   
    clusters = makeClustersSp(k,points,numOfPoints,d);
    
     
    // printClusters(clusters, k, d);
    doKmeans(clusters,points,d,k,numOfPoints, maxIter);
    freeMemo(clusters,points,k,numOfPoints);
    return 0;
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
    Wmat = calloc(1,sizeof(matrix));
    Wmat->rows = numOfPoints;
    Wmat->columns = numOfPoints;
    Wmat->data = calloc(numOfPoints,sizeof(double *));

    pointsMat = pointsInMat(points, numOfPoints, d);

    for (i = 0; i<numOfPoints; i++)
    {
        Wmat->data[i] = calloc(numOfPoints,sizeof(double));
        for(j = 0 ; j<numOfPoints ; j++)
        {
            dis = calcDistance(pointsMat[i],pointsMat[j]);  // needs to add root!
            Wmat->data[i][j]=exp(-dis/2);
        }
    }

    return Wmat;
}

matrix* formMatD(matrix* matW)
{
    int i,j,l;
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
                   sumW += matW[j][l];
                }
                matI[i][j] = sumW;
            }
            else
            {
                matI[i][j] = 0;
            }
        }
    }
    return matD;
}

matrix* formMatLnorm(matrix* matD , matrix* matW)
{
    matrix *lnormMat, *IMat;

    IMat = formMatI(matD->rows);
    lnormMat = minusRootMat(matD);
    lnormMat = mulMatrices(lnormMat,matW);
    lnormMat = mulMatrices(lnormMat,minusRootMat(matD));
    lnormMat = addMatrices(IMat, lnormMat ,1);
    return lnormMat;
}

matrix* minusRootMat(matrix* mat)
{
    int i,j;
    matrix* newMat;

    newMat = calloc(1,sizeof(matrix));
    newMat->rows = mat->rows;
    newMat->columns = mat->columns;
    newMat->data = calloc(mat->rows, sizeof(double*));
    
    for(i = 0; i< mat->rows; i++)
    {   
        newMat->data[i] = calloc(mat->columns,sizeof(double));
        for(j = 0; j< mat->columns; j++)
        {
            newMat[i][j] = pow(mat[i][j],-0.5);
        }
    }
    return newMat;
}

matrix* formMatI(int dimention)
{
    int i,j;
    matrix* matI;

    matI = calloc(1, sizeof(matrix));
    matI->rows = dimention;
    matI->rows = dimention;
    matI->data = calloc(dimention,sizeof(double*));


    for(i = 0 ; i< dimention ; i++)
    {   
        matI->data[i] = calloc(dimention,sizeof(double));
        for(j = 0 ; j< dimention ; j++)
        {
            if (i == j)
            {
                matI[i][j] = 1;
            }
            else
            {
                matI[i][j] = 0;
            }
        }
    }
    return matI;
}

matrix* pointsInMat(Node* points, int numOfPoints, int d)
{
    int i, j;
    matrix* pointsMat;
    Node* current;

    pointsMat = calloc(1, sizeof(matrix));
    pointsMat->rows = numOfPoints;
    pointsMat->rows = d;
    pointsMat->data = calloc(numOfPoints,sizeof(double*));
    current = points;

    for(i = 0 ; i< numOfPoints ; i++)
    {   
        pointsMat->data[i] = calloc(d,sizeof(double));
        for(j = 0 ; j< d ; i++)
        {
            pointsMat[i][j] = current->data[j];
        }
        current = current->next;
    }
    return pointsMat;
}

matrix* mulMatrices(matrix* mat1, matrix* mat2)
{
    int rows, columns, mat1Columns, i, j, l, sumElement;
    matrix* mulMat;

    rows = mat1->rows;
    columns = mat2->columns;
    mat1Columns = mat1->columns;

    mulMat = calloc(1, sizeof(matrix));
    assert(mulMat);
    mulMat->rows = rows;
    mulMat->columns = columns;
    mulMat->data = calloc(rows,sizeof(double *));
    assert(mulMat->data);
    for(i = 0 ; i<rows ; i++)
    {
        (mulMat->data)[i] = calloc(columns,sizeof(double));
        for (j = 0 ; j<columns ; j++)
        {
            sumElement = 0;
            for(l = 0 ; l<mat1Columns ; l++)
            {
                sumElement += (mat1->data[i][l])*(mat2->data[l][j]);
            }
            mulMat->data[i][j] = sumElement;
        }
    }
    return mulMat;
}

matrix* addMatrices(matrix* mat1, matrix* mat2, int dec) // dec=1 -> dec 
{
    // assumption: matrices equals in size
    int rows, columns, i, j;

    rows = mat1->rows;
    columns = mat1->columns;
    matrix* sumMat = calloc(1, sizeof(matrix));
    assert(sumMat);
    sumMat->rows = rows;
    sumMat->columns = columns;
    sumMat->data = calloc(rows,sizeof(double *));
    assert(sumMat->data);

    if (dec == 0)
    {
    for(i = 0 ; i<rows ; i++)
    {
        (sumMat->data)[i] = calloc(columns,sizeof(double));
        for(j = 0; j < columns; j++)
        {
            (sumMat->data)[i][j] = (mat1->data)[i][j]+(mat2->data)[i][j];
        }
    }
    }

    else
    {
        for(i = 0 ; i<rows ; i++)
    {
        (sumMat->data)[i] = calloc(columns,sizeof(double));
        for(j = 0; j<columns; j++)
        {
            (sumMat->data)[i][j] = (mat1->data)[i][j]-(mat2->data)[i][j];
        }
    }
    }

    return sumMat;

}


void jacobi(matrix* mat, matrix* eigenValues, matrix* eigenVectors)
{
    int iterations = 0, convergence = 0;
    double c, s;
    matrix *matA, *matB, *pivot, *tempVectors;

    matA = mat;
    matB = newMatrix(matA -> rows, matB -> columns); //TODO

    while (!isDiagonal(matA) && iterations < MaxJacobiIter && !convergence)
    {
        pivot = getRotationMatrixValues(matA, &c, &s); //TODO
        calcNextJacobiMatrix(matA, matB, c, s); //TODO, make sure that matB changes
        tempVectors = mulMatrices(tempVectors, pivot);
        convergence = hasConvergence(matA, matB);
        matA = matB;
        iterations++;
    }

    eigenValues = matA;
    eigenVectors = tempVectors;
}

matrix* getRotationMatrixValues(matrix* mat, double *c, double *s)
{
    matrix *p;
    int row, col, signTheta;
    double **m, theta, t;

    m = mat -> data;

    p = newMatrix(mat -> rows, mat -> columns);
    findlargestCell(mat, &row, &col);

    theta = (m[row][row] - m[col][col]) / (2 * m[row][col]);
    signTheta = abs(theta) / theta;
    t = signTheta / (abs(theta) + pow((theta * theta + 1), 0.5));
    *c = 1 / pow((t * t + 1), 0.5);
    *s = (*c) * t;
}

void findlargestCell(matrix* mat, int *row, int *col) 
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
            if (abs(m[i][j]) > big) 
            {
                *row = i;
                *col = j;
            }
        }
    }
    
}


void calcNextJacobiMatrix(matrix* matA, matrix* matB, double c, double s)
{

}

int hasConvergence(matrix* matA, matrix* matB)
{

}

// returns an array of the final starting vectors
int* eigenGapHeuristic(matrix* matA)
{

}

int isDiagonal(matrix* m)
{

}

