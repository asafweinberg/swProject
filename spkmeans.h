#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
//#include <Python.h>

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

enum goal{spk,wam,ddg,lnorm,jacobi};

#define true 1
#define false 0