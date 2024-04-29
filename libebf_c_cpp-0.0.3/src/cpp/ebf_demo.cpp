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

#include "ebfvector.hpp"

/**
 * \mainpage
 * \section Introduction
 * The libebf_cpp is main library to do I/O on EBF files (filenames with .ebf extension).
 * An ebf file is a collection of data objects
 * with each object having a unique name (dataname) with in a given file.
 * The datanames follow the unix style pathname convention and must begin with "/" signifying
 * the root path with in the file. Data objects
 * are an array of primitive types like int, float, double, long etc. One can perform two types
 * of I/O operations. First is the bulk operation in which the entire data object is copied
 * from the file into memory and returned as an array (or scalar if it has only one element).
 * Second is the partial I/O operation in which a specified number of elements are Read. This
 * allows more control over I/O operations. Bulk operations can be performed by a single function
 * call, but for partial I/O operations one has to perform
 * a three step process of Open , Read or Write , and Close. For reading prior to Read
 * operation an array of suitable size also needs to be created.
 *
 * \section  Naming-Conventions
 * We try to follow the google c++ style guide wherever possible.
 * Regular function names begin with Ebf_ with. First letter of each word in the function
 * is capitalized. Class names do not have underscore but otherwise follow same rule.
 * The first letter of first word in a class member is lowercase, for the rest of the words
 * the first letter is capitalized.
 *
 * \section Example
 * A code showing the use of the library is ebf_demo.cpp, click to view.
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
 * datatype                  bits    ebf_typecode
 * -----------------------------------------------------------------
 * char                   is 8  bit, 1
 * int32_t                is 32 bit, 2
 * int64_t                is 64 bit, 3
 * float                  is 32 bit, 4
 * double                 is 64 bit, 5
 * int16_t                is 16 bit, 6
 * int8_t                 is 8  bit, 9
 * uint8_t                is 8  bit, 10
 * uint16_t               is 16 bit, 11
 * uint32_t               is 32 bit, 12
 * unint64_t              is 64 bit, 13
 * </pre>
 *  On systems where the above assumption is not true, one might have to tweak the library
 *  to make it compatible.
 *
 * \section API
 * In the examples below int64_t is a 64 bit integer data type (long long int).
 * The unsigned version of this is uint64_t.
 *
 * \subsection Bulk-Read-Write
 * - template<class T2> void ebf::Write(const std::string& filename,
		const std::string& dataname, T2 * value1, const std::string& mode, std::string dataunit="", int64_t dim1=1,
		int64_t dim2 = 0, int64_t dim3 = 0, int64_t dim4 = 0)
 *   Writes an array of elements pointed to by pointer data. The last 5 arguments can have default values.
 *   @param filename1  The name of file in which to Write data e.g. "test.ebf".
 *   @param dataname1  The name of the data object, this should begin with "/" , e.g. "/x".
 *   @param data       An array or a scalar  of any one of the primitive data types.
 *   @param mode       The writing mode, "w" for fresh file and "a" to append to an existing file.
 *   @param dataunit:  Units of data e.g. 100 m/s.
 *   @param dim        A 64 bit integer specifying the number of elements to be written.
 *
 * - template<class T1, class T2> void ebf::WriteAs(const std::string& filename,
		const std::string& dataname, T2 * value1, const std::string& mode, std::string dataunit="",int64_t dim1=1,
		int64_t dim2 = 0, int64_t dim3 = 0, int64_t dim4 = 0) \n
 *   Writes an array of elements pointed to by pointer data. The last 5 arguments can have default values.
 *   T1 is the destination data type of the item to be written, while actual data type of item is T2.
 *   @param filename1  The name of file in which to Write data e.g. "test.ebf".
 *   @param dataname1  The name of the data object, this should begin with "/" , e.g. "/x".
 *   @param data       An array or a scalar  of any one of the primitive data types.
 *   @param mode       The writing mode, "w" for fresh file and "a" to append to an existing file.
 *   @param dataunit:  Units of data e.g. 100 m/s.
 *   @param dim        A 64 bit integer specifying the number of elements to be written.
 *
 * - template<class T1> void ebf::Read(const std::string& filename, const std::string& dataname,	T1 & value1) \n
 *   Reads an array of count elements from an ebf file into memory pointed to by x.
 *   @param filename1: The name of file from which to Read data e.g. "test.ebf"
 *   @param dataname1: A valid name of a data object that exists in the file, e.g. "/x"
 *   @param value1   : A container with resize function of any data type into which the data will be Read.
 *
 * - template<class T1> void ebf::Read(const std::string& filename, const std::string& dataname,T1* value1,int64_t ntot=1, int64_t offset1=-1)
 *   @param filename: The name of file from which to Read data e.g. "test.ebf"
 *   @param dataname: A valid name of a data object that exists in the file, e.g. "/x"
 *   @param value   : A pointer of any dtaa type into which the data will be Read.
 *   @param ntot           : The number of elements to Read
 *   @param offset1        : Offset in number of elements from where to start reading data from disc.
 *
 * - void ebf::WriteString(const std::string& filename, const std::string& dataname, const std::string &data, const std::string& mode);
 *   Write a c++ string
 *
 * - void ebf::ReadString(const std::string& filename, const std::string& dataname,std::string &mystr);
 *   Read a c++ string
 *
 * \subsection Partial-Read-Write
 *
 * - ebf::EbfFile efile. \n
 *   create an efile object
 * - void ebf::EbfFile::Open(const std::string& filename1, const std::string& dataname1, std::string mode1="r", int datatype=0, std::string dataunit="", int rank=1, int64_t* dims=NULL); \n
 *   Open a file and initialize file buffers to perform Read Write operations
 *   also positions the file to the location of the data.
 *   If a file is already Open from previous operation it is closed before opening a new buffer.
 *   @param filename1
 *   @param dataname1
 *   @param mode1     "w", "a" or "r"
 *   @param datatype the desired destination datatype
 *   @param dataunit units as a string
 *   @param rank     the rank or dimensionality of data
 *   @param dims     the shape or dimensionality of arrays being written.
 *   Depending upon the number of elements written, the value of dims[0] is adjusted such
 *   that the product dim[0]..dim[rank-1]= no of elements written. This is done the ebf file is closed.
 *
 * - template<class T1> void ebf::EbfFile::Write(T1* value1, int64_t data_size = 1) \n
 *   writes data_size elements from memeory pointed to by x, to an Open ebf file
 *   and increments the file pointer.
 *
 *
 * - template<class T1> void ebf::EbfFile::Read(T1 *value1, size_t data_size = 1) \n
 *   Reads data_size elements, from file into memory pointed to by x  and advances the file pointer.
 *
 *
 *
 *
 * - void ebf::EbfFile::Close(). \n
 *   Close the file and finish the Write or Read operation. After the object has been closed
 *   the same EbfFile object can be used for a fresh I/O. Note, for Write operations
 *   the data is fully written to the file, only after the Close() method has been called.
 *
 *
 * - void ebf::EbfFile::SaveTo(const std::string &filename1, const std::string& mode1); \n
 *   This will first Close the file. Subsequently, if the original file was opened in "w" mode
 *   it will  transfer the currently written data to the specified file with the specified mode and
 *   then delete the original file. Useful for copying data written un multiple files to one file.
 *
 * - int64_t ebf::EbfFile::headerpos()
 *   starting position of data-object in file (bytes)
 *
 * - int64_t ebf::EbfFile::datapos()
 *   starting position of data in file (bytes)
 *
 * - int64_t ebf::EbfFile::elements()
 *   total number of elements in data-object
 *
 * - int ebf::EbfFile::datatype()
 *   Data type of data in disc
 *
 * - int ebf::EbfFile::datasize()
 *   Size of data type in bytes
 *
 * - int ebf::EbfFile::capacity()
 *   total capacity in bytes that can be stored in current data-object
 *
 * - int ebf::EbfFile::rank()
 *   rank of data object
 *
 * - int64_t ebf::EbfFile::dim(int i)
 *   dimensions
 *
 * - std::string ebf::EbfFile::unit()
 *   unit string of data-object
 *
 * - std::string ebf::EbfFile::sdef()
 *   structure definition string
 *
 * - int64_t ebf::EbfFile::Index(int64_t i1,int64_t i2,int64_t i3)
 *   multi-dimensional Index to data
 *
 * - int64_t ebf::EbfFile::Index(int64_t i1,int64_t i2)
 *   multi-dimensional Index to data
 *
 * - void ebf::EbfFile::Seek(int64_t i)
 *   set to Read element no i
 *
 * - int64_t ebf::EbfFile::Tell()
 *   get the current element number
 *
 *\subsection General-Functions
 *
 * - int ebfc::Ebf_GetDataInfo(const char * filename1,const char* dataname1,EbfDataInfo* dinfo). \n
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
 *
 * - int ebf_impl1::TypeS(const std::string &typestring); \n
 * Get the integer type code used by EBF to define data types.
 * @param typestring The name of data type. Valid values are \n
 * char, int8, int16, int32, int64, uint8, uint16, uint32, uint64, float32, float64
 * @return an integer giving the type code
 *
 * - void ebf_impl1::ebfc::Copy(const char* infile, const char* outfile, const char* mode1,const char* dataname,int *ecode) \n
 *  A function to copy a data object from one ebf file to another
 * @param infile source file from which to copy data
 * @param outfile destination file to which data is to be copied
 * @param mode1 I/O mode "w" or "a" , fresh Write or append respectively
 * @param dataname the name of data to be copied
 * @param ecode an error code to test if operation had any error. Non zero on error.
 *
 *
 * \subsection Container
 * A container class implemented as a template. The operators [] and () can be used to access
 * data. The operator () allows multi-dimensional access upto 3 dimensions. Thre init method
 * can be used to load a new file. The methods rank() and size() can be used to query the
 * dimensions of the data.
 *
 *- template<class T> class ebf::EbfVector
 * - ebf::EbfVector(const std::string &fname, const std::string& dataname, std::string mode="r")
 *   Creates an instance of the class and load a data from an ebf file.
 * - void ebf::EbfVector::init(const std::string &fname,const std::string& dataname,std::string mode="r")
 *   Load a new data from a file.
 * - T& ebf::EbfVector::operator [](int64_t i)
 *   access the data
 * - T& ebf::EbfVector::operator ()(int64_t i1,int64_t i2)
 *   access 2-dimensional data
 * - T& ebf::EbfVector::operator ()(int64_t i1,int64_t i2,int64_t i3)
 *   access 3-dimensional data
 * - size_t ebf::EbfVector::size()
 *   total number of data elements
 * - int ebf::EbfVector::rank()
 *   the rank or number of dimensions
 * - const int64_t& ebf::EbfVector::dim(int i)
 *   the dimensions
 *
 * @author Sanjib Sharma
 *
 */


int ebf_demo()
{
		using namespace std;
	  	// generate some test data
	  	vector<float>   x1(100);
	  	vector<int32_t> x2(100);
	  	vector<double>  y1(100);
	  	vector<int64_t> y2(100);
	  	vector<double>  y3(100);
		int64_t dims[8];
		// create and initialize an efile structure/object
		ebf::EbfFile efile;
		ebf::EbfDataInfo dinfo;


	  	for(size_t i=0;i<x1.size();++i)
	  	{
	  		x1[i]=i;
	  		x2[i]=i;
	  	}
	  	// Write array x1 to file check.ebf
		ebf::Write("check.ebf", "/x1",&x1[0],"w","",x1.size());
	  	// append to file check.ebf array x2 with unit m/s
		ebf::Write("check.ebf", "/x2", &x2[0],"a","100 m*s^{-1}",x1.size());
		// Read data as long and double
		//  After openr(), one can allocate memory as y1=(double *)malloc(efile.elements()*8)
		efile.Open("check.ebf","/x1");
		efile.Read(&y1[0],efile.elements());
		efile.Close();
		efile.Open("check.ebf","/x2");
		cout<<"Units are "<<efile.unit()<<endl;
		efile.Read( &y2[0],efile.elements());
		efile.Close();

		// Read data
		if(ebf::ContainsKey("check.ebf","/x1",dinfo))
		{
			// one can allocate memory as y1=(double *)malloc(dinfo.elements*8)
			ebf::Read("check.ebf","/x1", &y1[0],dinfo.elements);
		}



		// Write a 2d array. dim[0] is adjusted to fit dataszie other dims are fixed.
		// 80 elements starting from x1[20] written as 8x10 array
		dims[0]=0;
		dims[1]=10;
		efile.Open("check.ebf","/x1","w",ebf::TypeS("float32"),"100 m*s^{-1}",2,dims);
		efile.Write(&x1[20],80);
		efile.Close();



		// print the data */
		cout<<"Printing x1 x2 y1 y2"<<endl;
		for(size_t i=70;i<80;++i)
				cout<<x1[i]<<" "<<x2[i]<<" "<<y1[i]<<" "<<y2[i]<<endl;

		// Write nsize elements of x1  and then of x2 as a double array
		efile.Open("check.ebf","/test/x1x2_double","a",ebf::TypeS("double")," 100 km*s^{-1}");
		efile.Write(&x1[0],x1.size());
		efile.Write(&x2[0],x2.size());
		efile.Close();



	  	// Write elements 20 to 100 as 2x40 multi dimensional array
	  	// the first dimension can be set to zero and is automatically set,
	  	// depending upon the number of elements written, when Close() method is invoked.
		// Also units are written. Can be set to blank string.
		// 80 elements starting from x1[20] written as 2x40 array

		dims[0]=0;	/* can be set to anything is adjusted when efile.Close() is called*/
		dims[1]=40;
		efile.Open("check.ebf","/test/x1_multi","a",ebf::TypeS("float32"),"100 m*s^{-1}",2,dims);
		efile.Write(&x1[20],80);
		efile.Close();

		//Read all elements
		efile.Open("check.ebf","/x1");
		//  On can allocate memory  as double *y3=(int64_t *)malloc(efile.elements()*8)
		efile.Read(&y3[0],80);
		efile.Close();

		//Read  elements 0 to 19 then next 40 and then next 20
		efile.Open("check.ebf","/x1");
		efile.Read( &y3[0],20);
		efile.Read( &y3[20],40);
		efile.Read( &y3[60],20);
		efile.Close();

		//Read  10 elements with an offset of 20 (x[20] to x[29])
		efile.Open("check.ebf","/x1");
		efile.Seek(20);
		efile.Read( &y3[20],10);
		efile.Close();
		cout<<"Units are "<<efile.unit()<<endl;


		cout<<"printing values 20 to 30"<<endl;
		for(int i=20;i<30;++i)
			cout<<i<<" "<<y3[i]<<endl;


		// using the EbfVector container
		ebf::EbfVector<float> y5("check.ebf","/x1");
		cout<<"Rank="<<y5.rank()<<"size="<<y5.size()<<" dim(0)="<<y5.dim(0)<<" dim(1)="<<y5.dim(1)<<endl;
		cout<<y5[9]<<" "<<y5[19]<<endl;
		cout<<y5(0,9)<<" "<<y5(1,19)<<endl;

		// Write and Read a string
		string mystr2;
		vector<float> x123;
		// same as 		ebf::WriteChar("test2.ebf", "/mystr",mystr1, "w","",strlen(mystr1));
		ebf::WriteString("check.ebf","/mystr","Hello World!","w");
		ebf::ReadString("check.ebf","/mystr",mystr2);
		cout<<mystr2<<endl;

		return 0;
}



