#define PY_SSIZE_T_CLEAN 
#include <Python.h>
#include "spkmeans.h"


static PyObject *runCFlow(PyObject *self, PyObject *args);
static PyObject *kmeansClustering(PyObject *self, PyObject *args);
cluster* makeClusters(int k, int d, PyObject ** centroids);
Node* getPointsMatrix(int d, int n, PyObject ** pointsMatrix);

static PyObject *runCFlow(PyObject *self, PyObject *args) 
{
    int k, i, j, finalK, numOfPoints;
    char* myGoal;
    char* fileName;
    double ** newPoints;
    PyObject  *wrapList, *miu = NULL;
    PyFloatObject *dNum = NULL;
    if (!PyArg_ParseTuple(args, "iss", &k, &myGoal, &fileName)) {
        return NULL;
    }

    newPoints = runMainFlow(k, myGoal, fileName, &finalK, &numOfPoints); // double ** 
    //can be merged to one func

    if (newPoints==NULL)
    {
        return Py_BuildValue("O", PyList_New(0)); //may cause error
    }

    wrapList = PyList_New(0);
    for (i=0 ; i<numOfPoints ; i++)
    {
        miu = PyList_New(0);    
        for (j=0 ; j<finalK ; j++)   //to check that point 0 exists
        {
            dNum= (PyFloatObject *)PyFloat_FromDouble(newPoints[i][j]);
            PyList_Append(miu, (PyObject*) dNum);
        }
        PyList_Append(wrapList, miu);
    }

    for (i=0; i<k; i++){
        free(newPoints[i]);
    }
    free(newPoints);

    return Py_BuildValue("O", wrapList); // returns to python the T mat
}

static PyObject *kmeansClustering(PyObject *self, PyObject *args) {
    int k, i, j, d, maxIterations, numberOfPoints;
    PyObject **points, **centroids, *wrapList, *miu = NULL;
    PyFloatObject *dNum = NULL;
    Node *nodesList;
    cluster* initialClusters;
    double** means;

    if (!PyArg_ParseTuple(args, "iiiiOO", &k, &maxIterations, &numberOfPoints, &d, &points, &centroids)) {
        return NULL;
    }

    nodesList = (Node*) getPointsMatrix(d, numberOfPoints, points);
    initialClusters = makeClusters(k, d, centroids);
    means = kmeansFunc(k, maxIterations, numberOfPoints, d, nodesList, initialClusters); //means are final centroids
    
    wrapList = PyList_New(0);
    for (i=0 ; i<k ; i++)
    {
        miu = PyList_New(0);
        for (j=0 ; j<d ; j++)
        {
            dNum= (PyFloatObject *)PyFloat_FromDouble(means[i][j]);
            PyList_Append(miu, (PyObject*) dNum);
        }
        PyList_Append(wrapList, miu);
    }

    for (i=0; i<k; i++){
        free(means[i]);
    }
    free(means);

    return Py_BuildValue("O", wrapList); // returns to python the final centroids
}

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

Node* getPointsMatrix(int d, int n, PyObject** pointsMatrix)
{
    int i,j;
    Node *points, *current;
    PyListObject *item;
    PyFloatObject *pointDouble;
    double* point;
    
    points = (Node*)calloc(1,sizeof(Node));
    assert(points);

    current = points;
    for (i=0; i<n; i++){
        item = (PyListObject *) PyList_GetItem((PyObject *) pointsMatrix , i);
        point=(double*)calloc(d,sizeof(double));

        for(j=0; j<d; j++){
            pointDouble= (PyFloatObject *)PyList_GET_ITEM(item,j);
            double value = PyFloat_AS_DOUBLE(pointDouble);
            point[j] = value;
        }
        current = addCurrentNext(current, point);
    }

    return points;
}



static PyMethodDef Module_Methods[] = {
    {
        "fit",      // name exposed to Python
        (PyCFunction) kmeansClustering, // C wrapper function
        METH_VARARGS,          // received variable args (but really just 1)
        "return clusters by kmeans method" // documentation
    },
    {
        "runCFlow",      // name exposed to Python
        (PyCFunction) runCFlow, // C wrapper function
        METH_VARARGS,          // received variable args (but really just 1)
        "perform the request action by goal" // documentation
    },
    {
        NULL, NULL, 0, NULL
    }
};

static struct PyModuleDef Moduledef = {
    PyModuleDef_HEAD_INIT,
    "mySpkmeansModule",     // name of module exposed to Python
    NULL, // module documentation
    -1,
    Module_Methods
};

PyMODINIT_FUNC PyInit_mySpkmeansModule(void) {
    PyObject *m;
    m = PyModule_Create(&Moduledef);
    if(!m) {
        return NULL;
    }
    return m;
}