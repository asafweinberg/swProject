#include "../spkmeans.c"

/* test section */
void testMain(int isDebug);
void testMultiplyMatrixs(int isDebug);
double** pointsForTestE0();
void TestE0(int isDebug);
double** pointsForTestE1();
void TestE1(int isDebug);
void testJacobi(int isDebug);
void testEigen(int isDebug);
void testReadPoints_Input0();
void testReadPoints_Input1();

void setMatrixValue(matrix* A, int row, int col, double value);
int isMatrixEqual(matrix *A, matrix *B);
void freeMemPointsArr(double **pointsArr, int n);
double* createPointWithVals(int d, double *values);
void printPointsArr(double **pointArr, int n, int d);
void printPoint(double *point, int dim);

Node* pointsArrToList(double **pointsArr, int n);



# define isDoubleEqual(_a, _b) (fabs((_a) - (_b)) < 0.0001)

void testMain(int isDebug) {
    // testMultiplyMatrixs(isDebug);
    // testJacobi(isDebug);
    testEigen(isDebug);
    // TestE0(isDebug);
    // TestE1(isDebug);
}

/* ############# */
/* Tests section */
/* ############# */

void testMultiplyMatrixs(int isDebug) {
    matrix *A, *B, *C, *D;
    A = newMatrix(3, 3);
    setMatrixValue(A , 0 , 0 , 1.0);
    setMatrixValue(A , 1 , 0 , 2.0);
    setMatrixValue(A , 1 , 1 , 3.0);
    setMatrixValue(A , 2 , 0 , 4.0);
    setMatrixValue(A , 2 , 1 , 5.0);
    setMatrixValue(A , 2 , 2 , 6.0);
    
    B = newMatrix(3, 2);
    setMatrixValue(B , 0 , 0 , 1.0);
    setMatrixValue(B , 1 , 0 , 8.0);
    setMatrixValue(B , 2 , 0 , 7.0);
    setMatrixValue(B , 0 , 1 , 11.0);
    setMatrixValue(B , 1 , 1 , 6.0);
    setMatrixValue(B , 2 , 1 , 1.0);

    if (isDebug) {
        printf("Matrix A: \n");
        printMatrix(A);
        printf("Matrix B: \n");
        printMatrix(B);
    }
    
    C = mulMatrices(A, B, false, false);
    if (isDebug) {
        printf("Matrix C calculted: \n");
        printMatrix(C);
    }

    D = newMatrix(3, 2);
    setMatrixValue(D , 0 , 0 , 1.0);
    setMatrixValue(D , 1 , 0 , 26.0);
    setMatrixValue(D , 2 , 0 , 86.0);
    setMatrixValue(D , 0 , 1 , 11.0);
    setMatrixValue(D , 1 , 1 , 40.0);
    setMatrixValue(D , 2 , 1 , 80.0);

    if (isDebug) {
        printf("Matrix D result: \n");
        printMatrix(D);
    }

    (isMatrixEqual(C,D)) ?
        printf("'test Multiply Matrixs'\t\tresult: Great!\n") : 
        printf("'test Multiply Matrixs'\t\tresult: Problem!\n");

    freeMatrix(A);
    freeMatrix(B);
    freeMatrix(C);
    freeMatrix(D);
}

double** pointsForTestE0(int* d) {
    double **pointsArr;
    int numOfPoints = 5;
    *d = 4;
    int dim = *d;
    double pointVal1[4] = {0.1255,-0.4507,-0.232,-0.0987};
    double pointVal2[4] = {0.344,0.344,0.4419,-0.3662}; 
    double pointVal3[4] = {-0.1011,-0.2081,0.4794,-0.4699};
    double pointVal4[4] = {-0.3324,0.2877,0.3182,0.3166}; 
    double pointVal5[4] = {0.1958,-0.0248,0.0681,0.2088};

    pointsArr = (double**)calloc(numOfPoints, sizeof(double*));
    
    pointsArr[0] = createPointWithVals(dim, pointVal1);
    pointsArr[1] = createPointWithVals(dim, pointVal2);
    pointsArr[2] = createPointWithVals(dim, pointVal3);
    pointsArr[3] = createPointWithVals(dim, pointVal4);
    pointsArr[4] = createPointWithVals(dim, pointVal5);
    return pointsArr;
}

void TestE0(int isDebug) {
    double **pointsArr;
    matrix *W, *WA, *D, *DA, *L, *LA;
    int numOfPoints = 5, d; 
    Node * pointsLst;

    /* Genarate points arr as input */
    pointsArr = pointsForTestE0(&d);
    pointsLst = pointsArrToList(pointsArr, numOfPoints);
    if(isDebug)
    {
        printLst(pointsLst, d);
    }
    W = formMatW(pointsLst, numOfPoints, d);

    /* Calculate matrix W */
    WA = newMatrix(numOfPoints, numOfPoints);
    setMatrixValue(WA, 0 ,0, 0.0);
    setMatrixValue(WA, 0 ,1, 0.5776);
    setMatrixValue(WA, 1 ,0, 0.5776);
    setMatrixValue(WA, 0 ,2, 0.6478);
    setMatrixValue(WA, 2 ,0, 0.6478);
    setMatrixValue(WA, 0 ,3, 0.5743);
    setMatrixValue(WA, 3 ,0, 0.5743);
    setMatrixValue(WA, 0 ,4, 0.7375);
    setMatrixValue(WA, 4 ,0, 0.7375);
    setMatrixValue(WA, 1 ,1, 0.0);
    setMatrixValue(WA, 1 ,2, 0.6985);
    setMatrixValue(WA, 2 ,1, 0.6985);
    setMatrixValue(WA, 1 ,3, 0.6155);
    setMatrixValue(WA, 3 ,1, 0.6155);
    setMatrixValue(WA, 1 ,4, 0.6728);
    setMatrixValue(WA, 4 ,1, 0.6728);
    setMatrixValue(WA, 2 ,2, 0.0);
    setMatrixValue(WA, 2 ,3, 0.6152);
    setMatrixValue(WA, 3 ,2, 0.6152);
    setMatrixValue(WA, 2 ,4, 0.6483);
    setMatrixValue(WA, 4 ,2, 0.6483);
    setMatrixValue(WA, 3 ,3, 0.0);
    setMatrixValue(WA, 3 ,4, 0.7148);
    setMatrixValue(WA, 4 ,3, 0.7148);
    setMatrixValue(WA, 4 ,4, 0.0);

    if (isDebug == 1) {
        printf("\nTestE0 - points array: \n");
        printPointsArr(pointsArr, numOfPoints, d);
        printf("\nTestE0 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE0 - Matrix WA correct Matrix\n");
        printMatrix(WA);
    }
    
    (isMatrixEqual(W,WA)) ?
        printf("TestE0 - Matrix W\t\tresult: Great!\n") : 
        printf("TestE0 - Matrix W\t\tresult: Problem!\n");
    
    /* Calculate Matrix D */
    D = formMatD(W);    
    
    DA = newMatrix(numOfPoints, numOfPoints);
    setMatrixValue(DA, 0 ,0, 2.5372);
    setMatrixValue(DA, 1 ,1, 2.5644);
    setMatrixValue(DA, 2 ,2, 2.6098);
    setMatrixValue(DA, 3 ,3, 2.5199);
    setMatrixValue(DA, 4 ,4, 2.7733);

    if (isDebug == 1) {
        printf("\nTestE0 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE0 - Matrix D calc:\n");
        printMatrix(D);
        printf("\nTestE0 - Matrix A correct Matrix\n");
        printMatrix(DA);
    }
    
    (isMatrixEqual(D,DA)) ?
        printf("TestE0 - Matrix D\t\tresult: Great!\n") : 
        printf("TestE0 - Matrix D\t\tresult: Problem!\n");

    L = formMatLnorm(D,W);
    
    LA = newMatrix(numOfPoints, numOfPoints);
    setMatrixValue(LA, 0 ,0, 1.0);
    setMatrixValue(LA, 0 ,1, -0.2264);
    setMatrixValue(LA, 0 ,2, -0.2517);
    setMatrixValue(LA, 0 ,3, -0.2271);
    setMatrixValue(LA, 0 ,4, -0.2780);

    setMatrixValue(LA, 1 ,0, -0.2264);
    setMatrixValue(LA, 1 ,1, 1.0);
    setMatrixValue(LA, 1 ,2, -0.2700);
    setMatrixValue(LA, 1 ,3, -0.2421);
    setMatrixValue(LA, 1 ,4, -0.2523);

    setMatrixValue(LA, 2 ,0, -0.2517);
    setMatrixValue(LA, 2 ,1, -0.2700);
    setMatrixValue(LA, 2 ,2, 1.0);
    setMatrixValue(LA, 2 ,3, -0.2399);
    setMatrixValue(LA, 2 ,4, -0.2410);

    setMatrixValue(LA, 3 ,0, -0.2271);
    setMatrixValue(LA, 3 ,1, -0.2421);
    setMatrixValue(LA, 3 ,2, -0.2399);
    setMatrixValue(LA, 3 ,3, 1.0);
    setMatrixValue(LA, 3 ,4, -0.2704);

    setMatrixValue(LA, 4 ,0, -0.2780);
    setMatrixValue(LA, 4 ,1, -0.2523);
    setMatrixValue(LA, 4 ,2, -0.2410);
    setMatrixValue(LA, 4 ,3, -0.2704);
    setMatrixValue(LA, 4 ,4, 1.0);

    if (isDebug == 1) {
        printf("\nTestE0 - Matrix L calculated: \n");
        printMatrix(L);
        printf("\nTestE0 - Matrix LA correct Matrix\n");
        printMatrix(LA);
    }
    
    (isMatrixEqual(L,LA)) ?
        printf("TestE0 - Matrix L\t\tresult: Great!\n") : 
        printf("TestE0 - Matrix L\t\tresult: Problem!\n");

    freeMatrix(W);
    freeMatrix(WA);
    freeMatrix(D);
    freeMatrix(DA);
    freeMatrix(L);
    freeMatrix(LA);
    freeMemPointsArr(pointsArr, 3);
}

double** pointsForTestE1(int* dim) {
    double **pointsArr;
    int numOfPoints = 8;
    *dim = 5;
    double pointVal1[5] = {-0.1119,0.3605,0.2079,0.1336,0.0439};
    double pointVal2[5] = {-0.2852,0.3003,-0.0142,-0.0924,0.4535}; 
    double pointVal3[5] = {-0.1075,0.3295,0.4349,-0.4489,-0.4656};
    double pointVal4[5] = {-0.3084,0.1954,0.4023,-0.1842,0.0598}; 
    double pointVal5[5] = {0.378,0.0048,0.4656,-0.4093,0.2412};
    double pointVal6[5] = {-0.1625,0.1883,-0.0201,-0.0467,0.3151};
    double pointVal7[5] = {-0.4696,-0.2751,-0.4395,-0.3948,-0.0979};
    double pointVal8[5] = {-0.4219,0.4115,0.304,0.4548,0.1747};

    pointsArr = calloc(numOfPoints, sizeof(double*));
    assert(pointsArr != NULL); 
    
    pointsArr[0] = pointVal1;
    pointsArr[1] = pointVal2;
    pointsArr[2] = pointVal3;
    pointsArr[3] = pointVal4;
    pointsArr[4] = pointVal5;
    pointsArr[5] = pointVal6;
    pointsArr[6] = pointVal7;
    pointsArr[7] = pointVal8;
    return pointsArr;
}

void TestE1(int isDebug) {
    double **pointsArr;
    matrix *W, *WA, *D, *DA, *L, *LA;
    int numOfPoints = 8, d; 

    /* Genarate points arr as input */
    pointsArr = pointsForTestE1(&d);
    W = formMatW(pointsArrToList(pointsArr, numOfPoints), numOfPoints, d);

    /* Calculate Matrix W */
    WA = newMatrix(numOfPoints, numOfPoints);
    setMatrixValue(WA, 0 ,0, 0.0);
    setMatrixValue(WA, 0 ,1, 0.7598);
    setMatrixValue(WA, 0 ,2, 0.6679);
    setMatrixValue(WA, 0 ,3, 0.7975);
    setMatrixValue(WA, 0 ,4, 0.6455);
    setMatrixValue(WA, 0 ,5, 0.8041);
    setMatrixValue(WA, 0 ,6, 0.5717);
    setMatrixValue(WA, 0 ,7, 0.7875);

    setMatrixValue(WA, 1 ,1, 0.0);
    setMatrixValue(WA, 1 ,2, 0.5775);
    setMatrixValue(WA, 1 ,3, 0.7444);
    setMatrixValue(WA, 1 ,4, 0.6218);
    setMatrixValue(WA, 1 ,5, 0.8953);
    setMatrixValue(WA, 1 ,6, 0.6156);
    setMatrixValue(WA, 1 ,7, 0.6999);

    setMatrixValue(WA, 2 ,2, 0.0);
    setMatrixValue(WA, 2 ,3, 0.7273);
    setMatrixValue(WA, 2 ,4, 0.6318);
    setMatrixValue(WA, 2 ,5, 0.6063);
    setMatrixValue(WA, 2 ,6, 0.5535);
    setMatrixValue(WA, 2 ,7, 0.5594);

    setMatrixValue(WA, 3 ,3, 0.0);
    setMatrixValue(WA, 3 ,4, 0.6800);
    setMatrixValue(WA, 3 ,5, 0.7661);
    setMatrixValue(WA, 3 ,6, 0.6027);
    setMatrixValue(WA, 3 ,7, 0.7045);

    setMatrixValue(WA, 4 ,4, 0.0);
    setMatrixValue(WA, 4 ,5, 0.6584);
    setMatrixValue(WA, 4 ,6, 0.5180);
    setMatrixValue(WA, 4 ,7, 0.5331);

    setMatrixValue(WA, 5 ,5, 0.0);
    setMatrixValue(WA, 5 ,6, 0.6436);
    setMatrixValue(WA, 5 ,7, 0.7038);

    setMatrixValue(WA, 6 ,6, 0.0);
    setMatrixValue(WA, 6 ,7, 0.5091);

    setMatrixValue(WA, 7 ,7, 0.0);

    if (isDebug == 1) {
        printf("\nTestE1 - points array: \n");
        printPointsArr(pointsArr, numOfPoints, d);
        printf("\nTestE1 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE1 - Matrix A correct Matrix\n");
        printMatrix(WA);
    }
    
    (isMatrixEqual(W,WA)) ?
        printf("TestE1 - Matrix W\t\tresult: Great!\n") : 
        printf("TestE1 - Matrix W\t\tresult: Problem!\n");
    
    /* Calculate Matrix D */
    D = formMatD(W);    
    
    DA = newMatrix(numOfPoints, numOfPoints);
    setMatrixValue(DA, 0 ,0, 5.0340);
    setMatrixValue(DA, 1 ,1, 4.9143);
    setMatrixValue(DA, 2 ,2, 4.3239);
    setMatrixValue(DA, 3 ,3, 5.0225);
    setMatrixValue(DA, 4 ,4, 4.2886);
    setMatrixValue(DA, 5 ,5, 5.0778);
    setMatrixValue(DA, 6 ,6, 4.0143);
    setMatrixValue(DA, 7 ,7, 4.4974);

    if (isDebug == 1) {
        printf("\nTestE1 - Matrix W calculated: \n");
        printMatrix(W);
        printf("\nTestE1 - Matrix D calc:\n");
        printMatrix(D);
        printf("\nTestE1 - Matrix A correct Matrix\n");
        printMatrix(DA);
    }
    
    (isMatrixEqual(D,DA)) ?
        printf("TestE1 - Matrix D\t\tresult: Great!\n") : 
        printf("TestE1 - Matrix D\t\tresult: Problem!\n");

    L = formMatLnorm(W,D);
    
    LA = newMatrix(numOfPoints, numOfPoints);
    setMatrixValue(LA, 0 ,0, 1.0);
    setMatrixValue(LA, 0 ,1, -0.1528);
    setMatrixValue(LA, 0 ,2, -0.1432);
    setMatrixValue(LA, 0 ,3, -0.1586);
    setMatrixValue(LA, 0 ,4, -0.1389);
    setMatrixValue(LA, 0 ,5, -0.1590);
    setMatrixValue(LA, 0 ,6, -0.1272);
    setMatrixValue(LA, 0 ,7, -0.1655);

    setMatrixValue(LA, 1 ,1, 1.0);
    setMatrixValue(LA, 1 ,2, -0.1253);
    setMatrixValue(LA, 1 ,3, -0.1498);
    setMatrixValue(LA, 1 ,4, -0.1354);
    setMatrixValue(LA, 1 ,5, -0.1792);
    setMatrixValue(LA, 1 ,6, -0.1386);
    setMatrixValue(LA, 1 ,7, -0.1489);

    setMatrixValue(LA, 2 ,2, 1.0);
    setMatrixValue(LA, 2 ,3, -0.1561);
    setMatrixValue(LA, 2 ,4, -0.1467);
    setMatrixValue(LA, 2 ,5, -0.1294);
    setMatrixValue(LA, 2 ,6, -0.1329);
    setMatrixValue(LA, 2 ,7, -0.1269);

    setMatrixValue(LA, 3 ,3, 1.0);
    setMatrixValue(LA, 3 ,4, -0.1465);
    setMatrixValue(LA, 3 ,5, -0.1517);
    setMatrixValue(LA, 3 ,6, -0.1342);
    setMatrixValue(LA, 3 ,7, -0.1482);

    setMatrixValue(LA, 4 ,4, 1.0);
    setMatrixValue(LA, 4 ,5, -0.1411);
    setMatrixValue(LA, 4 ,6, -0.1248);
    setMatrixValue(LA, 4 ,7, -0.1214);

    setMatrixValue(LA, 5 ,5, 1.0);
    setMatrixValue(LA, 5 ,6, -0.1426);
    setMatrixValue(LA, 5 ,7, -0.1473);

    setMatrixValue(LA, 6 ,6, 1.0);
    setMatrixValue(LA, 6 ,7, -0.1198);

    setMatrixValue(LA, 7 ,7, 1.0);

    if (isDebug == 1) {
        printf("\nTestE1 - Matrix L calculated: \n");
        printMatrix(L);
        printf("\nTestE1 - Matrix LA correct Matrix\n");
        printMatrix(LA);
    }
    
    (isMatrixEqual(L,LA)) ?
        printf("TestE1 - Matrix L\t\tresult: Great!\n") : 
        printf("TestE1 - Matrix L\t\tresult: Problem!\n");

    freeMatrix(W);
    freeMatrix(WA);
    freeMatrix(D);
    freeMatrix(DA);
    freeMatrix(L);
    freeMatrix(LA);
    freeMemPointsArr(pointsArr, numOfPoints);
}

void testJacobi(int isDebug) {
    matrix *V, *expectedV, *A, *expectedA, *Vectors;
    int testResult;

    // V = newMatrix(3, 3);
    // Vectors = newMatrix(3, 3);

    A = newMatrix(3, 3);
    setMatrixValue(A, 0, 0, 3.0);
    setMatrixValue(A, 0, 1, 2.0);
    setMatrixValue(A, 0, 2, 4.0);

    setMatrixValue(A, 1, 0, 2.0);
    setMatrixValue(A, 1, 1, 0.0);
    setMatrixValue(A, 1, 2, 2.0);

    setMatrixValue(A, 2, 0, 4.0);
    setMatrixValue(A, 2, 1, 2.0);
    setMatrixValue(A, 2, 2, 3.0);

    expectedV = newMatrix(3, 3);
    setMatrixValue(expectedV, 0, 0, 1 / sqrt(2));
    setMatrixValue(expectedV, 0, 1, - 1 / (3 * sqrt(2)));
    setMatrixValue(expectedV, 0, 2, 2.0 / 3.0);
    setMatrixValue(expectedV, 1, 0, 0);
    setMatrixValue(expectedV, 1, 1, 4 / (3 * sqrt(2)));
    setMatrixValue(expectedV, 1, 2, 1.0 / 3.0);
    setMatrixValue(expectedV, 2, 0, - 1 / sqrt(2));
    setMatrixValue(expectedV, 2, 1, - 1 / (3 * sqrt(2)));
    setMatrixValue(expectedV, 2, 2, 2.0 / 3.0);

    expectedA = newMatrix(3, 3);
    setMatrixValue(expectedA, 0, 0, -1);
    setMatrixValue(expectedA, 1, 1, -1);
    setMatrixValue(expectedA, 2, 2, 8);

    jacobiAlg(A, &V, &Vectors);
    testResult = isMatrixEqual(V, expectedA) && isMatrixEqual(Vectors, expectedV);
    (testResult) ?
        printf("'test Jacobi'\t\t\tresult: Great!\n") :
        printf("'test Jacobi'\t\t\tresult: Problem!\n");
    if (isDebug == 1) {
        printf("Matrix A\n");
        printMatrix(V);
        printf("Expected Matrix A\n");
        printMatrix(expectedA);
        printf("Matrix V\n");
        printMatrix(Vectors);
        printf("Expected Matrix V\n");
        printMatrix(expectedV);
    }
}

void testEigen(int isDebug) {
    matrix *A , *U, *eigens, *expectedEigens1, *expectedEigens2, *V, *Vectors;
    int testResult;
    double eps = 0.0000001;
    int i, j, k;
    
    A = newMatrix(3, 3);
    setMatrixValue(A, 0, 0, 3.0);
    setMatrixValue(A, 0, 1, 2.0);
    setMatrixValue(A, 0, 2, 4.0);

    setMatrixValue(A, 1, 0, 2.0);
    setMatrixValue(A, 1, 1, 0.0);
    setMatrixValue(A, 1, 2, 2.0);

    setMatrixValue(A, 2, 0, 4.0);
    setMatrixValue(A, 2, 1, 2.0);
    setMatrixValue(A, 2, 2, 3.0);


    expectedEigens1 = newMatrix(3, 3);
    setMatrixValue(expectedEigens1, 0, 0, 1 / sqrt(2));
    setMatrixValue(expectedEigens1, 0, 1, 0);
    setMatrixValue(expectedEigens1, 0, 2, - 1 / sqrt(2));

    setMatrixValue(expectedEigens1, 1, 0, - 1 / (3 * sqrt(2)));
    setMatrixValue(expectedEigens1, 1, 1, 4 / (3 * sqrt(2)));
    setMatrixValue(expectedEigens1, 1, 2, - 1 / (3 * sqrt(2)));

    setMatrixValue(expectedEigens1, 2, 0, 2.0 / 3.0);
    setMatrixValue(expectedEigens1, 2, 1, 1.0 / 3.0);
    setMatrixValue(expectedEigens1, 2, 2, 2.0 / 3.0);


    expectedEigens2 = newMatrix(3, 3);
    setMatrixValue(expectedEigens2, 1, 0, 1 / sqrt(2));
    setMatrixValue(expectedEigens2, 1, 1, 0);
    setMatrixValue(expectedEigens2, 1, 2, - 1 / sqrt(2));

    setMatrixValue(expectedEigens2, 0, 0, - 1 / (3 * sqrt(2)));
    setMatrixValue(expectedEigens2, 0, 1, 4 / (3 * sqrt(2)));
    setMatrixValue(expectedEigens2, 0, 2, - 1 / (3 * sqrt(2)));

    setMatrixValue(expectedEigens2, 2, 0, 2.0 / 3.0);
    setMatrixValue(expectedEigens2, 2, 1, 1.0 / 3.0);
    setMatrixValue(expectedEigens2, 2, 2, 2.0 / 3.0);

    jacobiAlg(A, &V, &Vectors);
    printMatrix(V);
    printMatrix(Vectors);
    eigens = calcInitialVectorsFromJacobi(V, Vectors);
    testResult = true;
    printMatrix(eigens);
    for(i = 0; i < 3; i++) {
        for (j = 0; j < eigens->columns; j++)
        {
            if(abs(eigens->data[i][j] - expectedEigens1->data[i][j]) > epsilon &&
               abs(eigens->data[i][j] - expectedEigens2->data[i][j]) > epsilon) {
                testResult = false;
                break;
            }
        }
        if(!testResult)
            break;
    }

    // (testResult) ?
    //     printf("'test Eigens'\t\t\tresult: Great!\n") :
    //     printf("'test Eigens'\t\t\tresult: Problem!\n");
    // if (isDebug == 1) {
    //     printf("Eigens\n");
    // }

    // (eigengapGetK(eigensArr) == 1) ?
    //     printf("'test Eigengap Heuristic'\tresult: Great!\n") :
    //     printf("'test Eigengap Heuristic'\tresult: Problem!\n");

    // k = 1;
    // U = computeMatrixU(eigensArr, k);

    // assert(U->cols = k);
    // assert(U->rows = 3);
    // testResult = true;
    // MatrixIterRows(U, i) {
    //     MatrixIterCols(U, j) {
    //         if (!isDoubleEqual(getMatrixValue(U, i, j), getDataPointVal((eigensArr->arr)[j].vector, i))) {
    //             printf("im\n");
    //             testResult = false;
    //             break;
    //         }
    //     }
    // }
    (testResult) ?
        printf("'test Matrix U'\t\t\tresult: Great!\n") :
        printf("'test Matrix U'\t\t\tresult: Problem!\n");
}

// void testReadPoints_Input0() {
//     double **pointsArr, **pointsArrRes;
//     int i, testResult = true, argc = 3, numOfPoint = 5; 
//     char* argv[] = {"3", "spk", ".\\Test_files\\Test_files\\input_0.txt"}; /* set path! */
//     pointsArr = readPointsfromFile(argc, argv);
//     pointsArrRes = pointsForTestE0();

//     for (i = 0; i < numOfPoint; i++) {
//         if ( !isPointsEquel(pointsArr[i], pointsArrRes[i])) {
//             testResult = false;
//             break;
//         }
//     }

//     (testResult) ?
//         printf("TestE0 - Read Input\t\tresult: Great!\n") : 
//         printf("TestE0 - Read Input\t\tresult: Problem!\n");
// }

// void testReadPoints_Input1() {
//     double **pointsArr, **pointsArrRes;
//     int i, testResult = true, argc = 3, numOfPoint = 5; 
//     char* argv[] = {"3", "spk", ".\\Test_files\\Test_files\\input_1.txt"}; /* set path! */
//     pointsArr = readPointsfromFile(argc, argv);
//     pointsArrRes = pointsForTestE1();

//     for (i = 0; i < numOfPoint; i++) {
//         if ( !isPointsEquel(pointsArr[i], pointsArrRes[i])) {
//             testResult = false;
//             break;
//         }
//     }

//     (testResult) ?
//         printf("TestE1 - Read Input\t\tresult: Great!\n") : 
//         printf("TestE1 - Read Input\t\tresult: Problem!\n");
// }



void setMatrixValue(matrix* A, int row, int col, double value) {
    // assert(row < A -> rows && col < A -> columns);
    A->data[row][col] = value;
}

int isMatrixEqual(matrix *A, matrix *B) {
    int i, j;
    double eps = 0.0001; /* set to 4 digits after the dot */ 
    int isSymmetric;
    assert(A -> rows == B -> rows);
    assert(A -> columns == B -> columns);

    for (i = 0; i < A->rows; i++)
    {
        for (j = 0; j < B->columns; j++)
        {
            if (fabs(A->data[i][j] - B->data[i][j]) > epsilon)
            {
                return false;
            }
        }
    }

    return true;
}

void freeMemPointsArr(double **pointsArr, int n) {
    int i; 
    for (i = 0; i < n; i++) {
        free(pointsArr[i]);
    }
    free(pointsArr);
}

double* createPointWithVals(int d, double *values) {
    double *point;
    int i;    
    point = (double*)calloc(d, sizeof(double));
    for (i = 0; i < d; i++) {
        point[i] = values[i];
    }
    return point;
}

void printPoint(double *point, int dim) {
    int i;
    for (i = 0; i < dim; i++) {
        printf("%.4f", point[i]);
        if (i != dim - 1) {
            printf(",");
        }
    }
    printf("\n");
}

void printPointsArr(double **pointArr, int n, int d) {
    int i;
    for (i = 0; i < n; i++) {
        printPoint(pointArr[i], d);
    }
}

Node* pointsArrToList(double **pointsArr, int n)
{
    int i;
    Node* start, *curr, *temp;
    curr = (Node*)calloc(1, sizeof(Node));
    start = curr;
    printPointsArr(pointsArr, n, 4);
    printf("\n");
    for (i = 0; i < n; i++)
    {
        temp = addNext(curr, pointsArr[i]);
        curr = temp;
    }
    return start;
}

int main() {
    testMain(true);
    
    // testReadPoints(false, path);
    
    return 0;
}