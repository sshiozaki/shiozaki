/**
 * @file   hdf5_output.h
 * @brief  hdf5 output Header
 * @author T.Otani
 */
#include "hdf5_input.h"
#include <iostream>
#include <cstdlib>
using namespace std;
using namespace H5;

void HDF5input::read_1dINTARRAY(H5File &file,const string &dataName,INTARRAY1 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int n1 = (unsigned long)(dims_out[0]);

    data=Allocation::allocate1dINT(n1);

    dataset.read(&data[0],PredType::NATIVE_INT);
}

void HDF5input::read_2dINTARRAY(H5File &file,const string &dataName,INTARRAY2 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int n1 = (unsigned long)(dims_out[0]);
    int n2 = (unsigned long)(dims_out[1]);

    data=Allocation::allocate2dINT(n1,n2);

    dataset.read(&data[0][0],PredType::NATIVE_INT);
}

void HDF5input::read_3dINTARRAY(H5File &file,const string &dataName,INTARRAY3 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[3];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int n1 = (unsigned long)(dims_out[0]);
    int n2 = (unsigned long)(dims_out[1]);
    int n3 = (unsigned long)(dims_out[2]);

    data=Allocation::allocate3dINT(n1,n2,n3);

    dataset.read(&data[0][0][0],PredType::NATIVE_INT);
}

void HDF5input::read_4dINTARRAY(H5File &file,const string &dataName,INTARRAY4 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[4];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int n1 = (unsigned long)(dims_out[0]);
    int n2 = (unsigned long)(dims_out[1]);
    int n3 = (unsigned long)(dims_out[2]);
    int n4 = (unsigned long)(dims_out[3]);

    data=Allocation::allocate4dINT(n1,n2,n3,n4);

    dataset.read(&data[0][0][0][0],PredType::NATIVE_INT);
}

void HDF5input::read_1dDOUBLEARRAY(H5File &file,const string &dataName,DOUBLEARRAY1 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[1];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int NX = (unsigned long)(dims_out[0]);
    data=Allocation::allocate1dDOUBLE(NX);

    dataset.read(&data[0],PredType::NATIVE_DOUBLE);
}

void HDF5input::read_2dDOUBLEARRAY(H5File &file,const string &dataName,DOUBLEARRAY2 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int NX = (unsigned long)(dims_out[0]);
    int NY = (unsigned long)(dims_out[1]);
    data=Allocation::allocate2dDOUBLE(NX,NY);

    dataset.read(&data[0][0],PredType::NATIVE_DOUBLE);
}

void HDF5input::read_3dDOUBLEARRAY(H5File &file,const string &dataName,DOUBLEARRAY3 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[3];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int NX = (unsigned long)(dims_out[0]);
    int NY = (unsigned long)(dims_out[1]);
    int NZ = (unsigned long)(dims_out[2]);
    data=Allocation::allocate3dDOUBLE(NX,NY,NZ);

    dataset.read(&data[0][0][0],PredType::NATIVE_DOUBLE);
}

void HDF5input::read_4dDOUBLEARRAY(H5File &file,const string &dataName,DOUBLEARRAY4 &data)
{
    DataSet dataset = file.openDataSet(dataName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[4];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    int n1 = (unsigned long)(dims_out[0]);
    int n2 = (unsigned long)(dims_out[1]);
    int n3 = (unsigned long)(dims_out[2]);
    int n4 = (unsigned long)(dims_out[3]);

    data=Allocation::allocate4dDOUBLE(n1,n2,n3,n4);

    dataset.read(&data[0][0][0][0],PredType::NATIVE_DOUBLE);
}
