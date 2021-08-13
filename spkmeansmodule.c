#include "spkmeans.h"


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
    means = kmeansFunc(k, maxIterations, numberOfPoints, d, nodesList, initialClusters);


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

    return Py_BuildValue("O", wrapList);
}

static PyMethodDef Module_Methods[] = {
    {
        "fit",      // name exposed to Python
        (PyCFunction) kmeansClustering, // C wrapper function
        METH_VARARGS,          // received variable args (but really just 1)
        "return clusters by kmeans method" // documentation
    }, {
        NULL, NULL, 0, NULL
    }
};

static struct PyModuleDef Moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",     // name of module exposed to Python
    NULL, // module documentation
    -1,
    Module_Methods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&Moduledef);
    if(!m) {
        return NULL;
    }
    return m;
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
        current = addNext(current, point);
    }

    return points;
}