/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Author: Sanjib Sharma                                                     *
 * Copyright (c) 2012 Sanjib Sharma                                          *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of EBF.  The full EBF copyright notice, including       *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and COPYRIGHT, which can be found at the root           *
 * of the source code distribution tree.                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "ebf.h"




/**
 * \mainpage
 * \section Introduction
 * Ebf is the main class to do I/O on EBF files (filenames with .ebf extension).
 * An ebf file is a collection of data objects
 * with each object having a unique name (dataname) with in a given file.
 * The datanames follow the unix style pathname convention and must begin with "/" signifying
 * the root path with in the file. Data objects
 * are an array of primitive types like int, float, double, long etc. One can perform two types
 * of I/O operations. First is the bulk operation in which the entire data object is copied
 * from the file into memory and returned as an array (or scalar if it has only one element).
 * Second is the partial I/O operation in which a specified number of elements are read. This
 * allows more control over I/O operations. Bulk operations can be performed by a single function
 * call, but for partial I/O operations one has to perform
 * a three step process of open , read or write , and close. For reading prior to read
 * operation an array of suitable size also needs to be created.
 *
 * \section Example
 * A code showing the use of Ebf class is ebf_demo.c, click to view.
 *
 * \section Portability
 * The code has been tested with gcc compiler on 64 bit linux machine. For other platforms
 * it is strongly recommended to run the test program to check for compatibility. In general,
 * the code should run correctly wherever header inttypes.h and stdint.h is available.
 * If these are not available one should supply ones own. Basically one has to supply
 * the intX_t and uintX_t  types, X being 8,16,32 and 64, using typedef statement.
 * The types sizeof(float) and sizeof(double) must be 4 and 8 bytes.
 *
 *
 * <pre>
 * datatype                  bits    ebf_typecode   Function Suffix
 * -----------------------------------------------------------------
 * char                   is 8  bit, 1,             Char
 * int                    is 32 bit, 2              Int32
 * long long int          is 64 bit, 3              Int64
 * float                  is 32 bit, 4              Float32
 * double                 is 64 bit, 5              FLoat64
 * short int              is 16 bit, 6              Int32
 * signed char            is 8  bit, 9              Int8
 * unsigned char          is 8  bit, 10             UInt8
 * unsigned short int     is 16 bit, 11             UInt16
 * unsigned int           is 32 bit, 12             UInt32
 * unsigned long long int is 64 bit, 13             UInt64
 * </pre>
 *  On systems where the above assumption is not true, one might have to tweak the library
 *  to make it compatible.
 *
 * \section API
 * In the examples below int64_t is a 64 bit integer data type (long long int).
 * The unsigned version of this is euint64.
 *
 * \subsection Bulk
 * - void Ebf_WriteFloat32(const char *filename1,const char * dataname1, const efloat32* data, const char* mode1, const char* dataunit, int64_t dim1); \n
 *   Writes an array of count elements pointed to by pointer data.
 *   There are 11 versions of this function depending upon the data type of data e.g Ebf_WriteInt32(), Ebf_writeDouble() etc.
 *   The chosen function should match the data type of data. For appropriate suffixes see table above.
 *   @param filename1  The name of file in which to write data e.g. "test.ebf".
 *   @param dataname1  The name of the data object, this should begin with "/" , e.g. "/x".
 *   @param mode       The writing mode, "w" for fresh file and "a" to append to an existing file.
 *   @param data       An array or a scalar  of any one of the primitive data types.
 *   @param count      A 64 bit integer specifying the number of elements to be written.
 *   @param dataunit: Units of data e.g. 100 m/s.
 *
 *
 * - void Ebf_ReadFloat32(const char *filename1,const char *dataname1,float* data, int64_t count). \n
 *   Reads an array of count elements from an ebf file into memory pointed to by x.
 *   There are 11 versions of this function depending upon the data type of data.
 *   The chosen function should match the data type of data.
 *   @param filename1: (char *) The name of file from which to read data e.g. "test.ebf"
 *   @param dataname1: (char *) A valid name of a data object that exists in the file, e.g. "/x"
 *   @param data       : A pointer of any data type into which the data will be read.
 *   @param count   : (int64_t) a 64 bit integer specifying the number of elements to be read
 *               Use \n EbfDataInfo dinfo = Ebf_getDataInfo() \n to get number of elements
 *               as dinfo.elements.
 *
 * - void Ebf_WriteCString(const char *filename1,const char *dataname1,char* data, const char* mode1);
 *   Write a C-style null terminated string
 *
 *
 * - void Ebf_ReadCString(const char *filename1,const char *dataname1,char* data, size_t maxsize);
 *   Read c style null terminated string. If size of supplied arrays is small
 *   reads only part of it and issues a warning.
 *   @param filename1
 *   @param dataname1
 *   @param data
 *   @param maxsize maximum size of data array including null char
 *
 *
 * \subsection Partial
 *
 * - EbfFile efile = EbfFile_Create(). \n
 *   create an efile object
 *
 * - void EbfFile_OpenR(EbfFile* efile,const char* filename1,const char* dataname1). \n
 *   open an efile object for read operation and position the file pointer to the location of the first element.
 *   @param filename1: (String) The name of file on which to perform I/O
 *   @param dataname1: (String) The name of the data object, to be read or written.
 *
 * - void EbfFile_OpenW(EbfFile* efile,const char* filename1,const char* dataname1,const char* mode1,int datatype,const char* dataunit) \n
 *   open an efile object for write operations.
 *   @param mode1    : (String) "w" for writing a fresh file, "a"  for appending to an existing file.
 *               In general "r" is for reading.
 *   @param datatype: an integer type code specifying the intended destination datatype.
 *               See Ebf_TypeS for generating from a string.
 *   @param dataunit:(char *)  units of data e.g. 100 m/s.
 *
 * - void EbfFile_Open(EbfFile* efile,const char* filename1,const char* dataname1,const char* mode1,int datatype,const char* dataunit,int rank, int64_t *dim);
 *   General function to open an efile object. Useful when multidimensional data needs
 *   to be written. If mode is set to "r" the rest of the fields are not utlilized.
 *   @param rank The number of dimensions.
 *   @param dim  A pointer to an array specifying the dimensions. Note, the first
 *     dimension can be set to 0. When Ebf_close is invoked dim[0] adjusted
 *     to match the number of items written.
 *
 * - void EbfFile_wFloat32(EbfFile *efile,const float* data, int64_t count);
 *   writes count elements from memeory pointed to by x, to an open ebf file and increments the file pointer.
 *   There are 11 versions of this function one for each data type.
 *
 * - void EbfFile_rFloat32(EbfFile *efile,efloat32* data, int64_t data_size); \n
 *   Reads elements offset to offset+count, from file into memory pointed to by x  and advances the file pointer.
 *   There are also 11 versions of this function one for each data type.
 *   @param offset if set to -1 , it simply reads from the current position pointed to by file pointer.
 *
 *
 *
 *
 * - void EbfFile_Close(EbfFile *efile). \n
 *   Close the file and finish the write or read operation. After the object has been closed
 *   the same EbfFile object can be used for a fresh I/O. Note, for write operations
 *   the data is fully written to the file, only after the close() method has been called.
 *
 *
 * - void EbfFile_SaveTo(EbfFile *efile,const char* filename1,const char* mode1). \n
 *   This will first close the file. Subsequently, if the original file was opened in "w" mode
 *   it will  transfer the currently written data to the specified file with the specified mode and
 *   then delete the original file. Useful for copying data written un multiple files to one file.
 * .
 *
 * \subsection General
 * - int Ebf_GetDataInfo(const char * filename1,const char* dataname1,EbfDataInfo* dinfo). \n
 * Get the DataInfo structure, which can be queried for number of elements
 * data units and existence of data and so on.
 * @param filename1  name of file  e.g "check.ebf"
 * @param dataname1  data tag name within the file e.g "/x"
 * @param dinfo a structure in which the information about data will be put
 * @return an integer error code, if non zero then it means that either data
 *         item is not present or file is not present or problem loading the data item.
 * Among other things dinfo can be queried for
 *         - elements
 *         - data units
 *         .
 *
 * - int Ebf_ContainsKey(const char* filename1, const char* dataname1,int* ecode); \n
 *   A function to check existence of a data object in a file
 * @param filename1  name of file  e.g "check.ebf"
 * @param dataname1  data object name within the file e.g "/x"
 * @return true if item exists or else false
 *
 * - int Ebf_TypeS(const char* typestring). \n
 * Get the integer type code used by EBF to define data types.
 * @param typestring The name of data type. Valid values are \n
 * char, int8, int16, int32, int64, uint8, uint16, uint32, uint64, float32, float64
 * @return an integer giving the type code
 *
 * - void Ebf_Copy(const char* infile, const char* outfile, const char* mode1,const char* dataname,int *ecode); \n
 *  A function to copy a data object from one ebf file to another
 * @param infile source file from which to copy data
 * @param outfile destination file to which data is to be copied
 * @param mode1 I/O mode "w" or "a" , fresh write or append respectively
 * @param dataname the name of data to be copied
 *
 *
 *
 * @author Sanjib Sharma
 *
 */

int ebf_demo()
{
	/* generate some test data */
	int nsize;
	float x1[100];
	int32_t x2[100];
	double y1[100];
	int64_t y2[100];
	double y3[100];
	char mystr1[200];
	char mystr2[200];
	int64_t dims[8];
	int i;
	EbfDataInfo dinfo;
	/* create and initialize an efile structure/object */
	EbfFile efile = EbfFile_Create();

	nsize = 100;
	for (i = 0; i < nsize; ++i)
	{
		x1[i] = i;
		x2[i] = i;
	}
	/* write array x1 to file check.ebf */
	Ebf_WriteFloat32("check.ebf", "/x1", &x1[0], "w", "", nsize);
	Ebf_WriteInt32("check.ebf", "/x2", &x2[0], "a", "100 m/s", nsize);

	/* read data as long and double */
	/*  After openr(), one can allocate memory as y1=(double *)malloc(EbfFile_Elements(&efile)*8)*/
	EbfFile_OpenR(&efile, "check.ebf", "/x1");
	EbfFile_rFloat64(&efile, &y1[0], nsize);
	EbfFile_Close(&efile);
	EbfFile_OpenR(&efile, "check.ebf", "/x2");
	printf("%s %s \n", "Units are", efile.ebfh.unit);
	EbfFile_rInt64(&efile, &y2[0], nsize);
	EbfFile_Close(&efile);

	if (Ebf_ContainsKey_Info("check.ebf", "/x1", &dinfo))
	{
		/* one can allocate memory as y1=(double *)malloc(dinfo.elements*8)*/
		Ebf_ReadFloat64("check.ebf", "/x1", &y1[0], dinfo.elements);
	}

	/* Write a 2d array. dim[0] is adjusted to fit dataszie other dims are fixed.*/
	/* 80 elements starting from x1[20] written as 8x10 array*/
	dims[0] = 0;
	dims[1] = 10;
	EbfFile_Open(&efile, "check.ebf", "/x1", "w", Ebf_TypeS("float32"),
			"100 m/s", 2, dims);
	EbfFile_wFloat32(&efile, &x1[20], 80);
	EbfFile_Close(&efile);

	/* print the data */
	printf("%s \n", "Printing x1 x2 y1 y2");
	for (i = 90; i < nsize; ++i)
		printf("%f %f %f %f \n", (float) x1[i], (float) x2[i], (float) y1[i],
				(float) y2[i]);

	/* Write nsize elements of x1  and then of x2 as a double array*/
	EbfFile_OpenW(&efile, "check.ebf", "/test/x1x2_double", "a",
			Ebf_TypeS("double"), " 100 km/s");
	EbfFile_wFloat32(&efile, &x1[0], nsize);
	EbfFile_wInt32(&efile, &x2[0], nsize);
	EbfFile_Close(&efile);

	/*
	 Write elements 20 to 100 as 2x40 multi dimensional array
	 the first dimension can be set to zero and is automatically set,
	 depending upon the number of elements written, when close() method is invoked.
	 Also units are written. Can be set to blank string.
	 80 elements starting from x1[20] written as 2x40 array
	 */
	dims[0] = 0; /* can be set to anything is adjusted when EbfFile_Close() is called*/
	dims[1] = 40;
	EbfFile_Open(&efile, "check.ebf", "/test/x1_multi", "a",
			Ebf_TypeS("float32"), "100 m/s", 2, dims);
	EbfFile_wFloat32(&efile, &x1[20], 80);
	EbfFile_Close(&efile);

	/*Read all elements*/
	EbfFile_OpenR(&efile, "check.ebf", "/x1");
	/*  On can allocate memory  as double *y3=(int64_t *)malloc(EbfFile_Elements(&efile)*8)*/
	EbfFile_rFloat64(&efile, &y3[0], 80);
	EbfFile_Close(&efile);

	/*Read  elements 0 to 19 then next 40 and then next 40*/
	EbfFile_OpenR(&efile, "check.ebf", "/x1");
	EbfFile_rFloat64(&efile, &y3[0], 20);
	EbfFile_rFloat64(&efile, &y3[20], 40);
	EbfFile_rFloat64(&efile, &y3[60], 20);
	EbfFile_Close(&efile);

	/*Read  10 elements with an offset of 20 (x[20] to x[29])*/
	EbfFile_OpenR(&efile, "check.ebf", "/x1");
	EbfFile_Seek(&efile, 20);
	EbfFile_rFloat64(&efile, &y3[20], 10);
	EbfFile_Close(&efile);
	printf("%s %s \n", "Units are", efile.ebfh.unit);

	printf("%s \n", "printing values 20 to 30");
	for (i = 20; i < 30; ++i)
		printf("%d %e \n", i, y3[i]);

	/* write and read a string*/
	strcpy(mystr1, "Hello, World!");
	/* same as 		Ebf_WriteChar("test2.ebf", "/mystr",mystr1, "w","",strlen(mystr1));*/
	Ebf_WriteCString("check.ebf", "/mystr", mystr1, "w");
	Ebf_ReadCString("check.ebf", "/mystr", mystr2, 200);
	printf("%s \n", mystr2);

	return 0;
}





