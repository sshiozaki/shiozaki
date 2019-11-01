/**
 * @file   hdf5_output.h
 * @brief  hdf5 output Header
 * @author T.Otani
 */
#include "hdf5_output.h"
#include <iostream>
#include <cstdlib>
using namespace std;
using namespace H5;

void HDF5output::write_1dINTARRAY(H5File &file,const string &dataName,const INTARRAY1 &data,const int size0)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[1];              // dataset dimensions
    dim[0] = size0;
    DataSpace dataspace( 1, dim );
    IntType datatype(PredType::STD_I32LE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0], PredType::STD_I32LE);
}

void HDF5output::write_2dINTARRAY(H5File &file,const string &dataName,const INTARRAY2 &data,const int size0,const int size1)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[2];              // dataset dimensions
    dim[0] = size0;
    dim[1] = size1;
    DataSpace dataspace( 2, dim );
    IntType datatype(PredType::STD_I32LE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0][0], PredType::STD_I32LE);
}

void HDF5output::write_3dINTARRAY(H5File &file,const string &dataName,const INTARRAY3 &data,const int size0,const int size1,const int size2)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[3];              // dataset dimensions
    dim[0] = size0;
    dim[1] = size1;
    dim[2] = size2;
    DataSpace dataspace( 3, dim );
    IntType datatype(PredType::STD_I32LE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0][0][0], PredType::STD_I32LE);
}

void HDF5output::write_4dINTARRAY(H5File &file,const string &dataName,const INTARRAY4 &data,const int size0,const int size1,const int size2,const int size3)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[4];              // dataset dimensions
    dim[0] = size0;
    dim[1] = size1;
    dim[2] = size2;
    dim[3] = size3;
    DataSpace dataspace( 4, dim );
    IntType datatype(PredType::STD_I32LE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0][0][0][0], PredType::STD_I32LE);
}

void HDF5output::write_1dDOUBLEARRAY(H5File &file,const string &dataName,const DOUBLEARRAY1 &data,const int size0)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[1];              // dataset dimensions
    dim[0] = size0;
    DataSpace dataspace( 1, dim );
    IntType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0], PredType::NATIVE_DOUBLE);
}

void HDF5output::write_2dDOUBLEARRAY(H5File &file,const string &dataName,const DOUBLEARRAY2 &data,const int size0,const int size1)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[2];              // dataset dimensions
    dim[0] = size0;
    dim[1] = size1;
    DataSpace dataspace( 2, dim );
    IntType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0][0], PredType::NATIVE_DOUBLE);
}

void HDF5output::write_3dDOUBLEARRAY(H5File &file,const string &dataName,const DOUBLEARRAY3 &data,const int size0,const int size1,const int size2)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[3];              // dataset dimensions
    dim[0] = size0;
    dim[1] = size1;
    dim[2] = size2;
    DataSpace dataspace( 3, dim );
    IntType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0][0][0], PredType::NATIVE_DOUBLE);
}

void HDF5output::write_4dDOUBLEARRAY(H5File &file,const string &dataName,const DOUBLEARRAY4 &data,const int size0,const int size1,const int size2,const int size3)
{
    H5std_string DATASET_NAME(dataName.c_str());

    hsize_t dim[4];              // dataset dimensions
    dim[0] = size0;
    dim[1] = size1;
    dim[2] = size2;
    dim[3] = size3;
    DataSpace dataspace( 4, dim );
    IntType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder( H5T_ORDER_LE );

    DataSet dataset = file.createDataSet(DATASET_NAME,datatype,dataspace);
    dataset.write( &data[0][0][0][0], PredType::NATIVE_DOUBLE);
}
