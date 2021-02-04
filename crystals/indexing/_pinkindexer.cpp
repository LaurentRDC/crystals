/* Copyright (C) Laurent P. René de Cotret */
/* All rights reserved. */
/* This file is part of the crystals library. */
/* See LICENSE for copying conditions. */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// https://numpy.org/devdocs/reference/c-api/deprecations.html#deprecation-mechanism-npy-no-deprecated-api
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <Eigen/Dense>
#include <PinkIndexer.h>
#include <ExperimentSettings.h>

using namespace Eigen;
using namespace pinkIndexer;
using namespace std;

const char *index_pink_docstring = R""""(
Index reflections using pinkindexer.

Parameters
----------
peaks : ndarray, shape (N, 2)
    Peak locations on detector [m]
intensities : ndarray, shape (N,)
    Scatterign intensity for each peak in `peaks` [a.u.]
detector_distance : float
    Distance between the sample and detector [m]
beam_energy : float
    Radiation energy [eV]
divergence_angle : float
    Divergence angle [deg]
non_monochromaticity: float
    I don't know what that is yet.
detector_radius : float
    Detector radius [m]
reciprocal_lattice : ndarray, shape (3,3)
    Initial guess for the eciprocal lattice [1/A]

Raises
------
ValueError
    If the number of peaks does not match the intensities provided.
    If the dimensions of the reciprocal lattice are not adequate.
)"""";

static PyObject *PinkIndexerError = NULL;

static PyObject * index_pink(PyObject *self, PyObject *args, PyObject *kwargs) {
    float detectorDistance_m;
    float beamEenergy_eV;
    float divergenceAngle_deg;
    float nonMonochromaticity;
    float detectorRadius_m;
    PyArrayObject * py_peaks;
    PyArrayObject * py_intensities;
    PyArrayObject * py_sample_lattice;
    

    static char *kwlist[] = {"peaks", "intensities", "detector_distance", "beam_energy", "divergence_angle", "non_monochromaticity", "detector_radius", "reciprocal_lattice", NULL };

    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, "OOfffffO", kwlist,
        &py_peaks,
        &py_intensities, 
        &detectorDistance_m, 
        &beamEenergy_eV, 
        &divergenceAngle_deg, 
        &nonMonochromaticity, 
        &detectorRadius_m, 
        &py_sample_lattice))
        return NULL;

    // Cast data into Eigen types
    npy_intp num_peaks = PyArray_SHAPE(py_peaks)[0];
    npy_intp num_intensities = PyArray_SHAPE(py_intensities)[0];
    if (num_peaks != num_intensities) {
        PyErr_SetString(PyExc_ValueError, "Number of peaks does not match the intensities provided.");
    }
    double (*peaks_raw)[2] = (double (*)[2])PyArray_DATA(py_peaks);
    double *intensities_raw = (double *)PyArray_DATA(py_intensities);
    Matrix2Xf peaksOnDetector_m(2, num_peaks);
    RowVectorXf intensities(num_peaks);

    for (int i=0; i < num_peaks; i++) {
        intensities(i) = intensities_raw[i];
        peaksOnDetector_m(0, i) = peaks_raw[0][i];
        peaksOnDetector_m(1, i) = peaks_raw[1][i];
    }

    // Cast lattice into Eigen types
    npy_intp num_rows = PyArray_SHAPE(py_sample_lattice)[0];
    npy_intp num_cols = PyArray_SHAPE(py_sample_lattice)[1];
    if ((num_rows != 3) || (num_cols != 3)) {
        PyErr_SetString(PyExc_ValueError, "reciprocal_lattice is expected to be a 3x3 matrix");
        return NULL;
    };
    double (*lat)[3];
    lat = (double(*)[3])PyArray_DATA(py_sample_lattice);

    Vector3f aStar(lat[0][0], lat[0][1], lat[0][2]);
    Vector3f bStar(lat[1][0], lat[1][1], lat[1][2]);
    Vector3f cStar(lat[2][0], lat[2][1], lat[2][2]);
    Matrix3f basis;
    basis << aStar, bStar, cStar;
    Lattice sampleReciprocalLattice_1A(basis);

    static ExperimentSettings settings = ExperimentSettings(
        beamEenergy_eV, 
        detectorDistance_m, 
        detectorRadius_m, 
        divergenceAngle_deg, 
        nonMonochromaticity, 
        sampleReciprocalLattice_1A,
        0.01,                       // tolerance 
        2.528445006321113e-04       // reflectionRadius_1_per_A
    );

    PinkIndexer indexer(
        settings, 
        PinkIndexer::ConsideredPeaksCount::standard, 
        PinkIndexer::AngleResolution::standard, 
        PinkIndexer::RefinementType::firstFixedThenVariableLatticeParameters, 
        0.2 //maxResolutionForIndexing_1_per_A
    );

    Lattice indexedLattice;
    try {
        int threadCount = 6;
        Eigen::Array<bool, Eigen::Dynamic, 1> fittedPeaks;
        Vector2f centerShift;
        indexer.indexPattern(indexedLattice, centerShift, fittedPeaks, intensities, peaksOnDetector_m, threadCount);
    } 
    catch (exception &e)
    {
        PyErr_SetString(PinkIndexerError, e.what());
    }
    Matrix3f indexedBasis = indexedLattice.getBasis();

    // For some reason, creating a numpy array directly always resulted in
    // a segfault. Therefore, we return a list of rows instead.
    PyObject * py_indexed_lattice, * a1, * a2, * a3;
    py_indexed_lattice = PyList_New(3);
    a1 = PyList_New(3);
    a2 = PyList_New(3);
    a3 = PyList_New(3);

    for (int i=0; i <3; i++) {
        PyList_SetItem(a1, i, PyFloat_FromDouble((double)indexedBasis(0, i)));
        PyList_SetItem(a2, i, PyFloat_FromDouble((double)indexedBasis(1, i)));
        PyList_SetItem(a3, i, PyFloat_FromDouble((double)indexedBasis(2, i)));
    }
    PyList_SetItem(py_indexed_lattice, 0, a1);
    PyList_SetItem(py_indexed_lattice, 1, a2);
    PyList_SetItem(py_indexed_lattice, 2, a3);
    return py_indexed_lattice;
}

static PyMethodDef PinkIndexerMethods[] {
    {"index_pink", (PyCFunction)(void(*)(void))index_pink, METH_VARARGS | METH_KEYWORDS, index_pink_docstring},
    {NULL, NULL, 0, NULL}
};
    
static struct PyModuleDef pinkindexermodule = {
    PyModuleDef_HEAD_INIT,
    "pinkindexer",
    "Python interface to pinkindexer",
    0,
    PinkIndexerMethods,
};

PyMODINIT_FUNC PyInit__pinkindexer(void) {
    Eigen::initParallel();
    PyObject *module = PyModule_Create(&pinkindexermodule);

    PinkIndexerError = PyErr_NewException("_pinkindexer.PinkIndexerError", NULL, NULL);
    PyModule_AddObject(module, "PinkIndexerError", PinkIndexerError);

    return module;
}