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
    Scattering intensity for each peak in `peaks` [a.u.]
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
tolerance : float
    Relative tolerance which determines when a peak has been successfully fit.
reciprocal_lattice : ndarray, shape (3,3)
    Initial guess for the reciprocal lattice [1/A]
considered_peaks_count: int
    Controls the number of Bragg spots which are used to compute the indexing solution.
angle_resolution: int
    Set the resolution of the rotogram in terms of number of voxels $N$ 
    spanning $-\arctan \pi/4$ to $\arctan \pi/4$
refinement_type: int
    Determines which type of refinement to perform after initial indexing.
num_threads : int
    Number of threads to use for indexing.

Returns
-------
indexed : list of lists
    Reciprocal basis vectors for the indexed reflections.
num_indexed : int
    Number of peaks that were indexed successfully.

Raises
------
ValueError
    If the number of peaks does not match the intensities provided.
    If the dimensions of the reciprocal lattice are not adequate.
PinkIndexerError
    If indexing has failed.
)"""";

static PyObject *PinkIndexerError = NULL;

static PyObject * index_pink(PyObject *self, PyObject *args, PyObject *kwargs) {
    float detectorDistance_m;
    float beamEnergy_eV;
    float divergenceAngle_deg;
    float nonMonochromaticity;
    float detectorRadius_m;
    float tolerance;
    float reflectionRadius_1_per_A;
    int considered_peaks_count, angle_resolution, refinement_type;
    int num_threads;
    PyArrayObject * py_peaks;
    PyArrayObject * py_intensities;
    PyArrayObject * py_sample_lattice;
    

    static char *kwlist[] = { "peaks"
                            , "intensities"
                            , "detector_distance"
                            , "beam_energy"
                            , "divergence_angle"
                            , "non_monochromaticity"
                            , "detector_radius"
                            , "tolerance"
                            , "reflectionRadius_1_per_A"
                            , "reciprocal_lattice"
                            , "considered_peaks_count"
                            , "angle_resolution"
                            , "refinement_type"
                            , "num_threads"
                            , NULL 
                            };

    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, "OOfffffffOiiii", kwlist,
        &py_peaks,
        &py_intensities, 
        &detectorDistance_m, 
        &beamEnergy_eV, 
        &divergenceAngle_deg, 
        &nonMonochromaticity, 
        &detectorRadius_m,
        &tolerance, 
        &reflectionRadius_1_per_A,
        &py_sample_lattice,
        &considered_peaks_count,
        &angle_resolution,
        &refinement_type,
        &num_threads))
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
        // Note that pinkindexer expects peak positions
        // as columns, but peaks_raw is a (N, 2) array
        peaksOnDetector_m(0, i) = peaks_raw[i][0];
        peaksOnDetector_m(1, i) = peaks_raw[i][1];
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
    Lattice sampleReciprocalLattice_1A = Lattice(basis);

    static ExperimentSettings settings = ExperimentSettings(
        beamEnergy_eV, 
        detectorDistance_m, 
        detectorRadius_m, 
        divergenceAngle_deg, 
        nonMonochromaticity, 
        sampleReciprocalLattice_1A,
        tolerance,
        reflectionRadius_1_per_A
    );

    PinkIndexer::ConsideredPeaksCount consideredPeaksCount = static_cast<PinkIndexer::ConsideredPeaksCount>(considered_peaks_count);
    PinkIndexer::AngleResolution angleResolution = static_cast<PinkIndexer::AngleResolution>(angle_resolution);
    PinkIndexer::RefinementType refinementType = static_cast<PinkIndexer::RefinementType>(refinement_type);

    PinkIndexer indexer(
        settings, 
        consideredPeaksCount, 
        angleResolution, 
        refinementType, 
        0.2 //maxResolutionForIndexing_1_per_A
    );

    Lattice indexedLattice;
    int num_indexed;
    try {
        Eigen::Array<bool, Eigen::Dynamic, 1> fittedPeaks;
        Vector2f centerShift;
        // Returns the number of fitted peaks according to tolerance
        // Therefore, if result is 0, indexing has failed.
        num_indexed = indexer.indexPattern(indexedLattice, centerShift, fittedPeaks, intensities, peaksOnDetector_m, num_threads);
    } 
    catch (exception &e)
    {
        PyErr_SetString(PinkIndexerError, e.what());
        return NULL;
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
    return Py_BuildValue("(Oi)", py_indexed_lattice, num_indexed);
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