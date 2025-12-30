# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2021-06-07 21:31:28
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2021-06-30 14:42:48

import numpy
import h5py
import collections
import scipy.sparse as sp_sparse

FeatureBCMatrix = collections.namedtuple('FeatureBCMatrix', ['ids', 'names', 'barcodes', 'matrix'])


def read_10X_h5(filename):
    """Read 10X HDF5 files using h5py (replaces pytables for better compatibility)."""
    with h5py.File(filename, 'r') as f:
        if 'matrix' not in f:
            print("Matrix group does not exist in this file.")
            return None
        
        group = f['matrix']
        feature_group = group['features']
        
        # Read feature data
        ids = feature_group['id'][:]
        names = feature_group['name'][:]
        
        # Read barcode data
        barcodes = group['barcodes'][:]
        
        # Read sparse matrix data
        data = group['data'][:]
        indices = group['indices'][:]
        indptr = group['indptr'][:]
        shape = group['shape'][:]
        
        # Create sparse matrix
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        
        return FeatureBCMatrix(ids, names, barcodes, matrix)


def read_count(count_file, separator = "tab"):
    """Read count table as matrix."""

    if separator == "tab":
        sep = "\t"
    elif separator == "space":
        sep = " "
    elif separator == "comma":
        sep = ","
    else:
        raise Exception("Invalid separator!")

    infile = open(count_file, 'r').readlines()
    barcodes = infile[0].strip().split(sep)
    features = []
    matrix = []
    for line in infile[1:]:
        line = line.strip().split(sep)
        features.append(line[0])
        matrix.append([float(t) for t in line[1:]])
    if len(barcodes) == len(matrix[0]) + 1:
        barcodes = barcodes[1:]

    return {"matrix": matrix, "features": features, "barcodes": barcodes}


def write_10X_h5(filename, matrix, features, barcodes):
    """Write 10X HDF5 files, support both gene expression and peaks."""
    f = h5py.File(filename, 'w')
    datatype = "Gene"
    M = sp_sparse.csc_matrix(matrix, dtype=numpy.float32)
    B = numpy.array(barcodes, dtype='|S200')
    P = numpy.array(features, dtype='|S100')
    FT = numpy.array([datatype]*len(features), dtype='|S100')
    mat = f.create_group('matrix')
    mat.create_dataset('barcodes', data=B)
    mat.create_dataset('data', data=M.data)
    mat.create_dataset('indices', data=M.indices)
    mat.create_dataset('indptr', data=M.indptr)
    mat.create_dataset('shape', data=M.shape)
    fet = mat.create_group('features')
    fet.create_dataset('id', data=P)
    fet.create_dataset('name', data=P)
    f.close()
