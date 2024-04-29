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

#include <sstream>
#include <cstring>
#include <cstdlib>

#include "ebf_demo.hpp"
#include "ebfvector.hpp"
#include "timer.hpp"
//#include "ebfextra.hpp"
//#include "ebfstruct.h"
//using namespace ebf;
using ebf::efloat32;
using ebf::efloat64;
using namespace std;

#define EBFLANGOPTCPP

#ifdef EBFLANGOPTCPP
#define Ebf_WriteChar ebf::Write
#define Ebf_WriteFloat32 ebf::Write
#define Ebf_WriteFloat64 ebf::Write
#define Ebf_WriteInt8 ebf::Write
#define Ebf_WriteInt16 ebf::Write
#define Ebf_WriteInt32 ebf::Write
#define Ebf_WriteInt64 ebf::Write
#define Ebf_WriteUInt8 ebf::Write
#define Ebf_WriteUInt16 ebf::Write
#define Ebf_WriteUInt32 ebf::Write
#define Ebf_WriteUInt64 ebf::Write
#define Ebf_ReadChar ebf::Read
#define Ebf_ReadFloat32 ebf::Read
#define Ebf_ReadFloat64 ebf::Read
#define Ebf_ReadInt8 ebf::Read
#define Ebf_ReadInt16 ebf::Read
#define Ebf_ReadInt32 ebf::Read
#define Ebf_ReadInt64 ebf::Read
#define Ebf_ReadUInt8 ebf::Read
#define Ebf_ReadUInt16 ebf::Read
#define Ebf_ReadUInt32 ebf::Read
#define Ebf_ReadUInt64 ebf::Read
#endif

template<class T>
void linspace(std::vector<T>&x,double xmin,double xmax,int size1)
{
	if(size1<2)
	{
		std::cout<<"linspace: size should be greater than 1"<<std::endl;
		exit(1);
	}
	x.resize(size1);
	double temp=(xmax-xmin)/(x.size()-1);
	for(unsigned int i=0;i<x.size();++i)
		x[i]=xmin+i*temp;

}

template<class T>
void printv(T a,T b,const string& s)
{
	std::cout<<"Printing Vector "<<s<<std::endl;
	while(a!=b)
		std::cout<<*a++<<" "<<std::endl;
//		std::cout<<scientific<<*a<<std::endl;
	std::cout<<std::endl;
}
template<class T>
void printv(T a,T b)
{
	printv(a,b,"");
}
template<class T>
void printv(std::vector<T> &a,const string& s)
{
	printv(&a[0],&a[a.size()],s);
}


int check_dinfo()
{
	int status=0;
	int i;
	efloat32 x1[10];
	efloat32 y1[10];
	ebf::EbfDataInfo dinfo;
	printf("%s \n","Testing check_dinfo()");
	for(i=0;i<10;++i)
		x1[i]=i;
	Ebf_WriteFloat32("check1.ebf","/x1",x1,"w","km/s",10);
	if(ebf::ContainsKey("check1.ebf","/x1",dinfo))
	{
		Ebf_ReadFloat32("check1.ebf","/x1",&y1[0],10);
		if((dinfo.ecode==0)&&(strcmp(dinfo.unit,"km/s")==0)&&(dinfo.rank==1)&&(dinfo.dim[0]==10))
			status=1;
	}
	ebf::ContainsKey("check1.ebf","/x2",dinfo);
	if(ebf::ebfc::Ebf_ContainsKey_Info("check1.ebf","/x2",&dinfo))
	{
		if((dinfo.ecode==0)&&(strcmp(dinfo.unit,"")==0)&&(dinfo.rank==0)&&(dinfo.headerpos <0))
			status=1;
	}
	return status;
}

int check_types()
{
	char x1[10];
	int32_t x2[10];
	int64_t x3[10];
	efloat32 x4[10];
	efloat64 x5[10];
	int16_t x6[10];
	int8_t x9[10];
	uint8_t x10[10];
	uint16_t x11[10];
	uint32_t x12[10];
	uint64_t x13[10];
	int status1=1,status;

	printf("%s \n","Testing getting typecode from given data type");

/* test for  .cpp files only */
#ifdef EBFLANGOPTCPP
	status=0;
	if((ebf::TypeV(x1[0])==1)&&(ebf::TypeV(x2[0])==2)&&(ebf::TypeV(x3[0])==3)&&(ebf::TypeV(x4[0])==4)&&(ebf::TypeV(x5[0])==5)&&(ebf::TypeV(x6[0])==6)&&(ebf::TypeV(x9[0])==9)&&(ebf::TypeV(x10[0])==10)&&(ebf::TypeV(x11[0])==11)&&(ebf::TypeV(x12[0])==12)&&(ebf::TypeV(x13[0])==13))
		status=1;
	if(status==0)
	{
		printf("%s \n","EXIT(): function getType not working properly");
	}
	status1*=status;

	status=0;
	if((ebf::TypeT<char>()==1)&&(ebf::TypeT<int>()==2)&&(ebf::TypeT<int64_t>()==3)&&(ebf::TypeT<float>()==4)&&(ebf::TypeT<double>()==5)&&(ebf::TypeT<short int>()==6)&&(ebf::TypeT<signed char>()==9)&&(ebf::TypeT<unsigned char>()==10)&&(ebf::TypeT<unsigned short int>()==11)&&(ebf::TypeT<unsigned int>()==12)&&(ebf::TypeT<uint64_t>()==13))
		status=1;

	if(status==0)
	{
		printf("%s \n","EXIT(): function getType not working properly");
	}
	status1*=status;
#endif

/* test for .c and .cpp files */
	status=1;
	if(sizeof(efloat32)!=4)
		{printf("%s \n","float not 32 bit"); status=0;}
	if(sizeof(int32_t)!=4)
		{printf("%s \n","int not 32 bit"); status=0;}
	if(sizeof(efloat64)!=8)
		{printf("%s \n","double not 64 bit"); status=0;}
	if(sizeof(int64_t)!=8)
		{printf("%s \n","int8 not 64 bit"); status=0;}
	if(sizeof(int16_t)!=2)
		{printf("%s \n","short int not 16 bit"); status=0;}
	if(sizeof(int8_t)!=1)
		{printf("%s \n","signed char not 8 bit"); status=0;}
	if(sizeof(uint8_t)!=1)
		{printf("%s \n","unsigned char not 8 bit"); status=0;}
	if(sizeof(uint16_t)!=2)
		{printf("%s \n","unsigned short int not 16 bit"); status=0;}
	if(sizeof(uint32_t)!=4)
		{printf("%s \n","unsigned int not 32 bit"); status=0;}
	if(sizeof(uint64_t)!=8)
		{printf("%s \n","uint8 not 64 bit"); status=0;}

	if(status==0)
	{
		printf("%s \n","EXIT(): 32 bit and 64 bit data type compatibility problem");
		printf("%s \n","hint: redefine datatypes (such as int64_t,uint64_t) in Ebf.h");
	}
	status1*=status;

	if(status1==0)
	{
		printf("%s \n","Test FAILED");
	}
	else
		printf("%s \n","SUCCESS");

	return status1;
}




int check_exceptions()
{
	printf("%s \n","Testing exception: overwrite protection");
	printf("%s \n"," and trying to read past end of record");
	Timer t1,t2;
	int nsize=100;
	vector<double> y1,y2;
	vector<int> y3(nsize);
    linspace(y1,0.0,nsize-1,nsize);
   	int status=1;
   	int status1=1;
   	y2.resize(nsize);

//	printf("%s \n","Testing exceptions: non existent files");
	try	{status1++; int x; ebf::Read("xxx.ebf", "/x1", &x);	} catch (exception &e)	{		status++;	}

   	t1.start();
   	ebf::Write("check1.ebf","/x0",&y1[0],"w","",nsize);
   	ebf::Write("check1.ebf","/x1",&y1[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x2",&y1[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x3",&y1[0],"a","",nsize);
//   	 printf("%s \n","check overwrite protection");
   	try{status1++; ebf::Write("check1.ebf","/x2",&y1[0],"a","",nsize);} catch (exception &e){status++;}

   	ebf::Read("check1.ebf","/x1",&y2[0],1);
   	ebf::Write("check1.ebf","/x0",&y1[0],"w");
   	ebf::Write("check1.ebf","/x1",&y1[0],"a");
   	ebf::Write("check1.ebf","/x2",&y1[0],"a");
   	ebf::Write("check1.ebf","/x3",&y1[0],"a");
//   	 printf("%s \n","check overwrite protection");
   	try{status1++; ebf::Write("check1.ebf","/x2",&y1[0],"a");} catch (exception &e){status++;}
   	ebf::Read("check1.ebf","/x1",&y2[0],1);

//   	printf("%s \n","check overwrite protection with multiple writes in other files in between");
   	ebf::Write("check2.ebf","/x0",&y1[0],"w");
   	ebf::Write("check2.ebf","/x1",&y1[0],"a");
   	ebf::Write("check2.ebf","/x2",&y1[0],"a");
   	ebf::Write("check2.ebf","/x3",&y1[0],"a");
   	try{status1++; ebf::Write("check2.ebf","/x1",&y1[0],"a");} catch (exception &e){status++;}
   	try{status1++; ebf::Write("check1.ebf","/x2",&y1[0],"a");} catch (exception &e){status++;}
   	ebf::Read("check1.ebf","/x1",&y2[0],1);
   	try{status1++; ebf::Write("check2.ebf","/x1",&y1[0],"a");} catch (exception &e){status++;}
  	try{status1++; ebf::Write("check1.ebf","/x2",&y1[0],"a");} catch (exception &e){status++;}

   	ebf::Write("check3.ebf","/x0",&y1[0],"w");
   	ebf::Read("check2.ebf","/x1",&y2[0]);
   	ebf::Write("check3.ebf","/x1",&y1[0],"a");
   	ebf::Read("check2.ebf","/x1",&y2[0]);
   	ebf::Read("check1.ebf","/x1",&y2[0]);
   	ebf::Write("check3.ebf","/x2",&y1[0],"a");
   	ebf::Read("check1.ebf","/x1",&y2[0]);
   	ebf::Read("check2.ebf","/x1",&y2[0]);
   	ebf::Write("check3.ebf","/x3",&y1[0],"a");
   	ebf::Read("check2.ebf","/x1",&y2[0]);
   	try{status1++; ebf::Write("check3.ebf","/x2",&y1[0],"a");} catch (exception &e){status++;}

//   	printf("%s \n","check trying to Read past end of record");
   	ebf::Write("check1.ebf","/x0",&y1[0],"w","",nsize);
   	ebf::Write("check1.ebf","/x1",&y1[0],"a","",nsize);
   	try{status1++; y2.resize(nsize+1); ebf::Read("check2.ebf","/x1",&y2[0],y2.size());}   	catch (exception &e){status++;}
   	try{status1++; y3.resize(nsize+1); ebf::Read("check2.ebf","/x1",&y3[0],y3.size());} 	catch (exception &e){status++;}

//   	printf("%s \n","check trying to Read past end of record with sequential Read and with different data types");
   	ebf::EbfFile efile;
   	y2.resize(nsize+1);

   	efile.Open("check1.ebf","/x1");
   	efile.Read(&y2[0],10);
   	efile.Read(&y2[0],1);
   	efile.Read(&y2[0],9);
   	efile.Read(&y2[0],80);
   	try{status1++;  efile.Read(&y2[0],1);   	}   	catch (exception &e){status++;}
   	efile.Close();

   	efile.Open("check1.ebf","/x1");
   	efile.Read( &y3[0],10);
   	efile.Read( &y3[0],1);
   	efile.Read( &y3[0],9);
   	efile.Read( &y3[0], 80);
   	try{status1++;  efile.Read( &y3[0],1);   	}   	catch (exception &e){status++;}
   	efile.Close();


	if (status1 == status)
	{
		printf("%s \n", "Test SUCCESS ");
		return 1;
	}
	else
	{
		printf("%s \n", "Test FAILED ");
		return 0;
	}

	return status;
}

int check_correctness()
{
	// initialize data with appropriate numbers for each type
	int score=0;
	int count=0;
	int status=1;
	int nsize=128;
	vector<char> x1(nsize),y1(nsize);
    linspace(x1,0,127,nsize);
	vector<int> x2(nsize),y2(nsize);
    linspace(x2,-65636,-65636-127,nsize);
	vector<int64_t> x3(nsize),y3(nsize);
    linspace(x3,-4294967296,-4294967296-127,nsize);
	vector<float> x4(nsize),y4(nsize);
    linspace(x4,1.23e20,128.23e20,nsize);
	vector<double> x5(nsize),y5(nsize);
    linspace(x5,1.23456789e200,128.23456789e200,nsize);
	vector<short int> x6(nsize),y6(nsize);
    linspace(x6,-256,-256-127,nsize);
	vector<signed char> x9(nsize),y9(nsize);
    linspace(x9,-128,-1,nsize);
	vector<unsigned char> x10(nsize),y10(nsize);
    linspace(x10,128,255,nsize);
	vector<unsigned short int> x11(nsize),y11(nsize);
    linspace(x11,256,256+127,nsize);
	vector<unsigned int> x12(nsize),y12(nsize);
    linspace(x12,65636,65636+127,nsize);
	vector<uint64_t> x13(nsize),y13(nsize);
    linspace(x13,4294967296,4294967296+127,nsize);

   	ebf::Write("check2.ebf","/x1",&x1[0],"w","",nsize);

    // Write in original types
   	ebf::Write("check1.ebf","/x1",&x1[0],"w","",nsize);
   	ebf::Write("check1.ebf","/x2",&x2[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x3",&x3[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x4",&x4[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x5",&x5[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x6",&x6[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x9",&x9[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x10",&x10[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x11",&x11[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x12",&x12[0],"a","",nsize);
   	ebf::Write("check1.ebf","/x13",&x13[0],"a","",nsize);

   	ebf::Write("check2.ebf","/x2",&x2[0],"a","",nsize);

   	ebf::Write("check1.ebf","/single/x1",&x1[0],"a");
   	ebf::Write("check1.ebf","/single/x2",&x2[0],"a");
   	ebf::Write("check1.ebf","/single/x3",&x3[0],"a");
   	ebf::Write("check1.ebf","/single/x4",&x4[0],"a");
   	ebf::Write("check1.ebf","/single/x5",&x5[0],"a");
   	ebf::Write("check1.ebf","/single/x6",&x6[0],"a");
   	ebf::Write("check1.ebf","/single/x9",&x9[0],"a");
   	ebf::Write("check1.ebf","/single/x10",&x10[0],"a");
   	ebf::Write("check1.ebf","/single/x11",&x11[0],"a");
   	ebf::Write("check1.ebf","/single/x12",&x12[0],"a");
   	ebf::Write("check1.ebf","/single/x13",&x13[0],"a");


    // Write with conversion to double
   	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x1", &x1[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x2", &x2[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x3", &x3[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x4", &x4[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x5", &x5[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x6", &x6[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x9", &x9[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x10", &x10[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x11", &x11[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x12", &x12[0], "a");
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/single/x13", &x13[0], "a");

   	ebf::Write("check2.ebf","/x3",&x3[0],"a","",nsize);


    // Write with conversion to double but for scalars
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x1", &x1[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x2", &x2[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x3", &x3[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x4", &x4[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x5", &x5[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x6", &x6[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x9", &x9[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x10", &x10[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x11", &x11[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x12", &x12[0], "a", "", nsize);
	ebf::WriteAs<double>("check1.ebf", "/AsDouble/x13", &x13[0], "a", "", nsize);


    // check scalar and array Read without type conversion
	{
		int i=0;
		y1[i]=y2[i]=y3[i]=y4[i]=y5[i]=y6[i]=y9[i]=y10[i]=y11[i]=y12[i]=y13[i]=0;
	}
   	ebf::Write("check2.ebf","/x4",&x4[0], "a", "",nsize);
   	ebf::Read("check1.ebf","/x1",&y1[0],1);
   	ebf::Read("check1.ebf","/x2",&y2[0],1);
   	ebf::Read("check1.ebf","/x3",&y3[0],1);
   	ebf::Read("check1.ebf","/x4",&y4[0],1);
   	ebf::Read("check1.ebf","/x5",&y5[0],1);
   	ebf::Read("check1.ebf","/x6",&y6[0],1);
   	ebf::Read("check1.ebf","/x9",&y9[0],1);
   	ebf::Read("check1.ebf","/x10",&y10[0],1);
   	ebf::Read("check1.ebf","/x11",&y11[0],1);
   	ebf::Read("check1.ebf","/x12",&y12[0],1);
   	ebf::Read("check1.ebf","/x13",&y13[0],1);
	{
   		int i=0;
		if(x1[i]!=y1[i])	{ status*=0; }
		if(x2[i]!=y2[i])	{ status*=0; }
		if(x3[i]!=y3[i])	{ status*=0; }
		if(x4[i]!=y4[i])	{ status*=0; }
		if(x5[i]!=y5[i])	{ status*=0; }
		if(x6[i]!=y6[i])	{ status*=0; }
		if(x9[i]!=y9[i])	{ status*=0; }
		if(x10[i]!=y10[i])	{ status*=0; }
		if(x11[i]!=y11[i])	{ status*=0; }
		if(x12[i]!=y12[i])	{ status*=0; }
		if(x13[i]!=y13[i])	{ status*=0; }
	}

	{
		int i=0;
		y1[i]=y2[i]=y3[i]=y4[i]=y5[i]=y6[i]=y9[i]=y10[i]=y11[i]=y12[i]=y13[i]=0;
	}
   	ebf::Write("check2.ebf","/x5",&x5[0],"a", "",nsize);
   	ebf::Read("check1.ebf","/single/x1",&y1[0],1);
   	ebf::Read("check1.ebf","/single/x2",&y2[0],1);
   	ebf::Read("check1.ebf","/single/x3",&y3[0],1);
   	ebf::Read("check1.ebf","/single/x4",&y4[0],1);
   	ebf::Read("check1.ebf","/single/x5",&y5[0],1);
   	ebf::Read("check1.ebf","/single/x6",&y6[0],1);
   	ebf::Read("check1.ebf","/single/x9",&y9[0],1);
   	ebf::Read("check1.ebf","/single/x10",&y10[0],1);
   	ebf::Read("check1.ebf","/single/x11",&y11[0],1);
   	ebf::Read("check1.ebf","/single/x12",&y12[0],1);
   	ebf::Read("check1.ebf","/single/x13",&y13[0],1);
	{
   		int i=0;
		if(x1[i]!=y1[i])	{ status*=0; }
		if(x2[i]!=y2[i])	{ status*=0; }
		if(x3[i]!=y3[i])	{ status*=0; }
		if(x4[i]!=y4[i])	{ status*=0; }
		if(x5[i]!=y5[i])	{ status*=0; }
		if(x6[i]!=y6[i])	{ status*=0; }
		if(x9[i]!=y9[i])	{ status*=0; }
		if(x10[i]!=y10[i])	{ status*=0; }
		if(x11[i]!=y11[i])	{ status*=0; }
		if(x12[i]!=y12[i])	{ status*=0; }
		if(x13[i]!=y13[i])	{ status*=0; }
	}


	for(size_t i=0;i<y1.size();++i)
	{
		y1[i]=y2[i]=y3[i]=y4[i]=y5[i]=y6[i]=y9[i]=y10[i]=y11[i]=y12[i]=y13[i]=0;
	}
   	ebf::Write("check2.ebf","/x6",&x6[0],"a", "",nsize);
   	ebf::Read("check1.ebf","/x1",&y1[0],nsize);
   	ebf::Read("check1.ebf","/x2",&y2[0],nsize);
   	ebf::Read("check1.ebf","/x3",&y3[0],nsize);
   	ebf::Read("check1.ebf","/x4",&y4[0],nsize);
   	ebf::Read("check1.ebf","/x5",&y5[0],nsize);
   	ebf::Read("check1.ebf","/x6",&y6[0],nsize);
   	ebf::Read("check1.ebf","/x9",&y9[0],nsize);
   	ebf::Read("check1.ebf","/x10",&y10[0],nsize);
   	ebf::Read("check1.ebf","/x11",&y11[0],nsize);
   	ebf::Read("check1.ebf","/x12",&y12[0],nsize);
   	ebf::Read("check1.ebf","/x13",&y13[0],nsize);
	for(size_t i=0;i<y1.size();++i)
	{
		if(x1[i]!=y1[i])	{ status*=0; }
		if(x2[i]!=y2[i])	{ status*=0; }
		if(x3[i]!=y3[i])	{ status*=0; }
		if(x4[i]!=y4[i])	{ status*=0; }
		if(x5[i]!=y5[i])	{ status*=0; }
		if(x6[i]!=y6[i])	{ status*=0; }
		if(x9[i]!=y9[i])	{ status*=0; }
		if(x10[i]!=y10[i])	{ status*=0; }
		if(x11[i]!=y11[i])	{ status*=0; }
		if(x12[i]!=y12[i])	{ status*=0; }
		if(x13[i]!=y13[i])	{ status*=0; }
	}

	printf("%s \n","Testing:");
	printf("%s \n","correctness of read and write, same type to same type");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;

    // check scalar and array Read of type converted written data
   	ebf::Write("check2.ebf","/x9",&x9[0],"a", "",nsize);
	vector<string> mytags;
	mytags.push_back("/x1");	mytags.push_back("/x2");	mytags.push_back("/x3");	mytags.push_back("/x4");
	mytags.push_back("/x5");	mytags.push_back("/x6");	mytags.push_back("/x9");	mytags.push_back("/x10");
	mytags.push_back("/x11");	mytags.push_back("/x12");	mytags.push_back("/x13");
	for(size_t j=0;j<mytags.size();++j)
	{
		vector<double> yd2(nsize,0);
		ebf::Read("check1.ebf",mytags[j],&yd2[0],nsize);
		vector<double> yd1(nsize,0);
		ebf::Read("check1.ebf","/AsDouble/single"+mytags[j],&yd1[0],1);
		if(yd1[0]!=yd2[0])	{ status*=0; cout<<"AsDouble/single"+mytags[j]<<" "<<yd1[0]<<" "<<yd2[0]<<endl;}
		for(int i=0;i<nsize;++i)
		{
			yd1[i]=0;
		}
		ebf::Read("check1.ebf","/AsDouble"+mytags[j],&yd1[0],nsize);
		for(size_t i=0;i<y1.size();++i)
		{
			if(yd1[i]!=yd2[i])	{ status*=0;  cout<<"AsDouble"+mytags[j]<<" "<<yd1[0]<<" "<<yd2[0]<<endl;}
		}
	}


	printf("%s \n","Testing:");
	printf("%s \n","correctness of write, other types to double ");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;

   	ebf::Write("check2.ebf","/x10",&x10[0],"a", "",nsize);


	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x1",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x1[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x2",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x2[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x3",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x3[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x4",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x4[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x5",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x5[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x6",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x6[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x9",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x9[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x10",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x10[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x11",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x11[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x12",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x12[i])!=yd[i]) { status*=0;}}
	}
	{
		vector<double> yd(nsize,0);	   	ebf::Read("check1.ebf","/x13",&yd[0],nsize);
		for(size_t i=0;i<y1.size();++i) {if(double(x13[i])!=yd[i]) { status*=0;}}
	}
	printf("%s \n","Testing:");
	printf("%s \n","correctness of read, other type to double");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	count++;
	score+=status;


   	ebf::Write("check2.ebf","/x11",&x11[0],"a", "",nsize);
	{		double yd;	   	ebf::Read("check1.ebf","/x1",&yd,1);		if(double(x1[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x2",&yd,1);		if(double(x2[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x3",&yd,1);		if(double(x3[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x4",&yd,1);	    if(double(x4[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x5",&yd,1);		if(double(x5[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x6",&yd,1);		if(double(x6[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x9",&yd,1);		if(double(x9[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x10",&yd,1);		if(double(x10[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x11",&yd,1);		if(double(x11[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x12",&yd,1);		if(double(x12[0])!=yd) {status*=0;}	}
	{		double yd;	   	ebf::Read("check1.ebf","/x13",&yd,1);		if(double(x13[0])!=yd) {status*=0;}	}
   	ebf::Write("check2.ebf","/x12",&x12[0],"a", "",nsize);
   	ebf::Write("check2.ebf","/x13",&x13[0],"a", "",nsize);

	printf("%s \n","Testing:");
	printf("%s \n","correctness of scalar read, other type to double");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	count++;
	score+=status;

	for(size_t j=0;j<mytags.size();++j)
	{
		vector<double> yd2(nsize,0);
		vector<double> yd1(nsize,0);
		ebf::Read("check1.ebf",mytags[j],&yd1[0],nsize);
		ebf::Read("check2.ebf",mytags[j],&yd2[0],nsize);
		for(size_t i=0;i<y1.size();++i)
		{
			if(yd1[i]!=yd2[i])	{ status*=0;  cout<<"AsDouble"+mytags[j]<<" "<<yd1[0]<<" "<<yd2[0]<<endl;}
		}
	}
	printf("%s \n","Testing:");
	printf("%s \n","correctness of data written in between to a different file");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	count++;
	score+=status;


//	for(size_t j=0;j<mytags.size();++j)
//	{
//		vector<double> yd1;
//		vector<double> yd2;
//		Read(EBFDATADIR"master_test1.ebf",mytags[j],yd1);
//		Read(EBFDATADIR"master_test1_swap.ebf",mytags[j],yd2);
//		for(size_t i=0;i<yd1.size();++i)
//		{
//			if(yd1[i]!=yd2[i])	{ status*=0;  cout<<mytags[j]<<" "<<yd1[0]<<" "<<yd2[0]<<endl;}
//		}
//	}

	printf("%s \n","Testing:");
	printf("%s \n","endian swap ");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	count++;
	score+=status;

	if(score==count)
		return 1;
	else
		return 0;


}


int check_masterfile(const string &dirname)
{
	// initialize data with appropriate numbers for each type
	int status=1;
	int nsize=128;
	vector<char> x1(nsize),y1(nsize);
    linspace(x1,0,127,nsize);
	vector<int> x2(nsize),y2(nsize);
    linspace(x2,-65636,-65636-127,nsize);
	vector<int64_t> x3(nsize),y3(nsize);
    linspace(x3,-4294967296,-4294967296-127,nsize);
	vector<float> x4(nsize),y4(nsize);
    linspace(x4,1.23e20,128.23e20,nsize);
	vector<double> x5(nsize),y5(nsize);
    linspace(x5,1.23456789e200,128.23456789e200,nsize);
	vector<short int> x6(nsize),y6(nsize);
    linspace(x6,-256,-256-127,nsize);
	vector<signed char> x9(nsize),y9(nsize);
    linspace(x9,-128,-1,nsize);
	vector<unsigned char> x10(nsize),y10(nsize);
    linspace(x10,128,255,nsize);
	vector<unsigned short int> x11(nsize),y11(nsize);
    linspace(x11,256,256+127,nsize);
	vector<unsigned int> x12(nsize),y12(nsize);
    linspace(x12,65636,65636+127,nsize);
	vector<uint64_t> x13(nsize),y13(nsize);
    linspace(x13,4294967296,4294967296+127,nsize);


string filename1=dirname+"master_test1.ebf";
string filename2=dirname+"master_test1_swap.ebf";
status=1;
{
	vector<char> yd1,yd2;	   		ebf::Read(filename1,"/x1",yd1);	ebf::Read(filename2,"/x1",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x1",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x1[i])!=yd1[i]) { status*=0;}
		if((x1[i])!=yd2[i]) { status*=0;}
		if(double(x1[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in char read");
}
{
	vector<int32_t> yd1,yd2;	   		ebf::Read(filename1,"/x2",yd1);	ebf::Read(filename2,"/x2",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x2",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x2[i])!=yd1[i]) { status*=0;}
		if((x2[i])!=yd2[i]) { status*=0;}
		if(double(x2[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in int32 read");
}
{
	vector<int64_t> yd1,yd2;	   		ebf::Read(filename1,"/x3",yd1);	ebf::Read(filename2,"/x3",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x3",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x3[i])!=yd1[i]) { status*=0;}
		if((x3[i])!=yd2[i]) { status*=0;}
		if(double(x3[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in int64 read");
}
{
	vector<float> yd1,yd2;	   		ebf::Read(filename1,"/x4",yd1);	ebf::Read(filename2,"/x4",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x4",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x4[i])!=yd1[i]) { status*=0;}
		if((x4[i])!=yd2[i]) { status*=0;}
		if(double(x4[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in float32 read");
}
{
	vector<double> yd1,yd2;	   		ebf::Read(filename1,"/x5",yd1);	ebf::Read(filename2,"/x5",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x5",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x5[i])!=yd1[i]) { status*=0;}
		if((x5[i])!=yd2[i]) { status*=0;}
		if(double(x5[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in float64 read");
}
{
	vector<int16_t> yd1,yd2;	   		ebf::Read(filename1,"/x6",yd1);	ebf::Read(filename2,"/x6",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x6",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x6[i])!=yd1[i]) { status*=0;}
		if((x6[i])!=yd2[i]) { status*=0;}
		if(double(x6[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in int16 read");
}
{
	vector<int8_t> yd1,yd2;	   		ebf::Read(filename1,"/x9",yd1);	ebf::Read(filename2,"/x9",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x9",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x9[i])!=yd1[i]) { status*=0;}
		if((x9[i])!=yd2[i]) { status*=0;}
		if(double(x9[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in int8 read");
}
{
	vector<uint8_t> yd1,yd2;	   		ebf::Read(filename1,"/x10",yd1);	ebf::Read(filename2,"/x10",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x10",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x10[i])!=yd1[i]) { status*=0;}
		if((x10[i])!=yd2[i]) { status*=0;}
		if(double(x10[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in uint8 read");
}
{
	vector<uint16_t> yd1,yd2;	   		ebf::Read(filename1,"/x11",yd1);	ebf::Read(filename2,"/x11",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x11",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x11[i])!=yd1[i]) { status*=0;}
		if((x11[i])!=yd2[i]) { status*=0;}
		if(double(x11[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in uint16 read");
}
{
	vector<uint32_t> yd1,yd2;	   		ebf::Read(filename1,"/x12",yd1);	ebf::Read(filename2,"/x12",yd2);
	vector<double> yd3; 	ebf::Read(filename1,"/x12",yd3);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x12[i])!=yd1[i]) { status*=0;}
		if((x12[i])!=yd2[i]) { status*=0;}
		if(double(x12[i])!=yd3[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in uint32 read");
}
{
	vector<uint64_t> yd1,yd2;	   		ebf::Read(filename1,"/x13",yd1);	ebf::Read(filename2,"/x13",yd2);
	for(size_t i=0;i<yd1.size();++i)
	{
		if((x13[i])!=yd1[i]) { status*=0;}
		if((x13[i])!=yd2[i]) { status*=0;}
	}
	if(status==0)
		printf("%s \n","Problem in uint64 read");
}


vector<string> mytags;
mytags.push_back("/x1");	mytags.push_back("/x2");	mytags.push_back("/x3");	mytags.push_back("/x4");
mytags.push_back("/x5");	mytags.push_back("/x6");	mytags.push_back("/x9");	mytags.push_back("/x10");
mytags.push_back("/x11");	mytags.push_back("/x12");	mytags.push_back("/x13");
for(size_t j=0;j<mytags.size();++j)
{
	vector<double> yd1;
	vector<double> yd2;
	ebf::Read(filename1,mytags[j],yd1);
	ebf::Read(filename2,mytags[j],yd2);
	cout<<mytags[j]<<" "<<yd1.size()<<endl;
	for(size_t i=0;i<yd1.size();++i)
	{
		if(yd1[i]!=yd2[i])	{ status*=0;  cout<<mytags[j]<<" "<<yd1[0]<<" "<<yd2[0]<<endl;}
	}
}



printf("%s \n","Testing: master_test1.ebf");
if (status == 1)
	printf("%s \n", "Test SUCCESS ");
else
	printf("%s \n", "Test FAILED ");

return status;
}




int check_ebfcopy()
{
	int count=0;
	int score=0;
	int status=1;
	int ecode=0;
	vector<double> y1,y2,y3,y4;
	size_t nsize=128;

	vector<char> x1(nsize);
    linspace(x1,0,127,nsize);
	vector<int> x2(nsize);
    linspace(x2,-65636,-65636-127,nsize);
	vector<int64_t> x3(nsize);
    linspace(x3,-4294967296,-4294967296-127,nsize);
	vector<float> x4(nsize);
    linspace(x4,1.234*1e-30,(1.234+127)*1e-30,nsize);
	vector<double> x5(nsize),y5(nsize);
    linspace(x5,1.23456789012*1e-300,(1.23456789012+127)*1e-300,nsize);
	vector<short int> x6(nsize);
    linspace(x6,-256,-256-127,nsize);
	vector<signed char> x9(nsize);
    linspace(x9,-128,-1,nsize);
	vector<unsigned char> x10(nsize);
    linspace(x10,128,255,nsize);
	vector<unsigned short int> x11(nsize);
    linspace(x11,256,256+127,nsize);
	vector<unsigned int> x12(nsize);
    linspace(x12,65636,65636+127,nsize);
	vector<uint64_t> x13(nsize);
    linspace(x13,4294967296,4294967296+127,nsize);


	printf("%s \n","Testing ebfcopy ");

   	ebf::Write("check1.ebf","/x1",&x1[0],"w", "",x1.size());
   	ebf::Write("check1.ebf","/x2",&x2[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x3",&x3[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x4",&x4[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x5",&x5[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x6",&x6[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x9",&x9[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x10",&x10[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x11",&x11[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x12",&x12[0],"a", "",x1.size());
   	ebf::Write("check1.ebf","/x13",&x13[0],"a", "",x1.size());

   	ebf::Write("check2.ebf","/test/x1",&x1[0],"w", "",x1.size());
   	ebf::Write("check2.ebf","/test/x2",&x2[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x3",&x3[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x4",&x4[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x5",&x5[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x6",&x6[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x9",&x9[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x10",&x10[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x11",&x11[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x12",&x12[0],"a", "",x1.size());
   	ebf::Write("check2.ebf","/test/x13",&x13[0],"a", "",x1.size());

   	ecode=ebf::Copy("check1.ebf","check3.ebf","w","");
   	ecode=ebf::Copy("check2.ebf","check3.ebf","a","");


	vector<string> mytags;
	mytags.push_back("/x1");
	mytags.push_back("/x2");
	mytags.push_back("/x3");
	mytags.push_back("/x4");
	mytags.push_back("/x5");
	mytags.push_back("/x6");
	mytags.push_back("/x9");
	mytags.push_back("/x10");
	mytags.push_back("/x11");
	mytags.push_back("/x12");
	mytags.push_back("/x13");
	string s=mytags[0]+" "+mytags[1]+" "+mytags[2];
   	ecode=ebf::Copy("check1.ebf","check4.ebf","w",s.c_str());
	for(size_t j=3;j<mytags.size();++j)
	{
	   	ecode=ebf::Copy("check1.ebf","check4.ebf","a",mytags[j].c_str());
	}

	printf("%s \n","testing exception in ebfcopy");
	int status1;
	try{status1=0; ecode=ebf::Copy("check1.ebf","check1.ebf","a","/x1");} 	catch (exception &e){status1=1;}
	status*=status1;
	try{status1=0; ecode=ebf::Copy("check1.ebf","check3.ebf","a","/x1");} catch (exception &e){status1=1;}
	status*=status1;
	try{status1=0; ecode=ebf::Copy("check1.ebf","check3.ebf","a","");} catch (exception &e){status1=1;}
	status*=status1;
	try{status1=0; ecode=ebf::Copy("check1.ebf","check3.ebf","a","/abc");} catch (exception &e){status1=1;}
	status*=status1;
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;

	for(size_t j=0;j<mytags.size();++j)
	{
		y1.resize(nsize);
		y2.resize(nsize);
		y3.resize(nsize);
		y4.resize(nsize);
		for(size_t i=0;i<nsize;++i)
		{
			y1[i]=0.0;			y2[i]=0.0;			y3[i]=0.0;			y4[i]=0.0;
		}
		ebf::Read("check1.ebf",mytags[j],&y1[0],nsize);
		ebf::Read("check3.ebf",mytags[j],&y2[0],nsize);
		ebf::Read("check3.ebf","/test"+mytags[j],&y3[0],nsize);
		ebf::Read("check4.ebf",mytags[j],&y4[0],nsize);



		if((y1.size() != nsize)||(y2.size()!= nsize)||(y3.size()!= nsize)||(y4.size()!= nsize))
		{
			status*=0;
			cout<<"Size mismatch tag="<<mytags[j]<<" "<<nsize<<" "<<y1.size()<<" "<<y2.size()<<" "<<y3.size()<<" "<<y4.size()<<endl;
			break;
		}
		for(size_t i=0;i<y1.size();++i)
		{
			if((y1[i] != y2[i])||(y1[i] != y3[i])||(y1[i] != y4[i]))
			{
				status*=0;
				cout<<"value mismatch tag="<<mytags[j]<<" "<<i<<" "<<y1[i]<<" "<<y2[i]<<" "<<y3[i]<<" "<<y4[i]<<endl;
				break;
			}
		}


		if(status==0)
			break;
	}

	if(status==0)
		cout<<"Checking: "<<setw(16)<<left<<ebf::TypeT<double>()<<" status="<<status<<endl;

	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	count++;
	score+=status;
	printf("%d %d \n", score,status);

	status=ecode;
	if(score==count)
		return 1;
	else
		return 0;


}


int test_ebf_locate_speed()
{

	Timer t1,t2;
	size_t nsize=10000;
	vector<double> y1;
	vector<int64_t> yi(nsize);
	vector<double> yd(nsize,0);
    linspace(y1,0.0,nsize-1,nsize);
	printf("%s \n","Testing speed to read and write a single item");

	vector<string> labels;
   	for(size_t j=0;j<y1.size();++j)
   	{
   		stringstream mystr;
   		mystr<<"/x"<<j;
   		labels.push_back(mystr.str());
   	}

   	t1.start();
   	ebf::Write("check1.ebf","/x0",&y1[0],"w");
   	for(size_t j=1;j<y1.size();++j)
   	{
   	   	ebf::Write("check1.ebf",labels[j],&y1[j],"a");
   	}
   	t1.printC("Item writing speed",y1.size());

    t1.start();
   	for(size_t j=0;j<y1.size();++j)
   	{
   	   	ebf::Read("check1.ebf",labels[j],&yd[j],1);
   		if(yd[j]!=double(j))
   		{
   			cout<<"Error "<<j<<" "<<yd[j]<<endl;
   			exit(1);
   		}
   	}
   	t1.printC("Item reading speed",y1.size());
//   	cout<<ebf::ebfutils::FileSize("check.ebf")<<endl;


//   	t1.start();
//   	{
//   		vector<string> key;
//   		vector<int64_t> value;
//   		ebf::EbfHT ebfht;
//   		ebfht.getKeyVal("check.ebf",key,value);
//   	}
//   	t1.printC("getKeyVal speed",y1.size());

return 0;
}



int test_ebf_speed()
{
	Timer t1,t2;
	int nsize=10000000;
	int lsize=10;
	vector<double> y1;
	vector<int64_t> yi(nsize,0);
	vector<double> yd(nsize,0);
    linspace(y1,-128,128,nsize);
   	ebf::EbfFile efile;
	printf("%s \n","test_ebf1() ");
	printf("%s \n","------------------------------------------------------");

   	ebf::Write("check2.ebf","/x1",&y1[0],"w", "",y1.size());

	vector<float> y2;
    linspace(y2,-128,128,nsize);
   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   	   	ebf::Write("check1.ebf","/x1",&y2[0],"w", "",y2.size());
   	}
   	t1.printC("efile write float ",y2.size()*sizeof(y2[0])*lsize);
   	t1.printC("efile write float operations",lsize*y2.size());

	vector<float> yf(nsize);
   	t1.start();
   	yf.resize(nsize);
   	for(int j=0;j<lsize;++j)
   	{
   	   	ebf::Read("check1.ebf","/x1",&yf[0],yf.size());
   	}
   	t1.printC("efile read float ",yf.size()*sizeof(yf[0])*lsize);
   	t1.printC("efile read float  operations",lsize*yf.size());

   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   	   	ebf::Write("check1.ebf","/x1",&y1[0],"w", "",y1.size());
   	}
   	t1.printC("efile write double ",y1.size()*sizeof(y1[0])*lsize);
   	t1.printC("efile write double operations",lsize*y1.size());

   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   	   	ebf::Read("check1.ebf","/x1",&yd[0],yd.size());
   	}
   	t1.printC("efile read double",yd.size()*sizeof(yd[0])*lsize);
   	t1.printC("efile read double operations",lsize*yd.size());


   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   		efile.Open("check2.ebf","/x1");
   		for(size_t i=0;i<y1.size();++i)
   			efile.Read(&yi[i],1);
   		efile.Close();
   	}
   	t1.printC("efile read single gen",y1.size()*8*lsize);


   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   		efile.Open("check2.ebf","/x1");
   		for(size_t i=0;i<y1.size();++i)
   			efile.Read(&yd[i],1);
   		efile.Close();
   	}
   	t1.printC("efile read single gen native",y1.size()*8*lsize);



   	ebf::EbfVector<int64_t> ev1;
   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   		ev1.init("check2.ebf","/x1");
   		for(size_t i=0;i<y1.size();++i)
   			yd[i]=ev1[i];
   	}
   	t1.printC("efile ebfvector seq read",y1.size()*8*lsize);


   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   		ebf::Read("check2.ebf","/x1",&yi[0],nsize);
   	}
   	t1.printC("efile read bulk",y1.size()*8*lsize);

// This should give compile error
   	//   	Read("check2.ebf","/x1",yd[0]);

   	t1.start();
   	for(int j=0;j<lsize;++j)
   	{
   		ebf::Read("check2.ebf","/x1",&yd[0],nsize);
   	}
   	t1.printC("efile read bulk native",y1.size()*8*lsize);

return 0;
}

void check_checksum()
{
	int nsize=10;
	vector<float> x1(nsize);
	vector<double> x2(nsize);
	int64_t checksum=0;

	ebf::Write("check1.ebf","/x1",&x1[0],"w","",x1.size());
	ebf::Write("check1.ebf","/x2",&x2[0],"a","",x2.size());
	ebf::Write("check1.ebf","/single/x1",&x1[0],"a");
	ebf::Write("check1.ebf","/single/x2",&x2[0],"a");
	ebf::Read("check1.ebf","/.ebf/info",&checksum);
	cout<<"checksum="<<checksum<<endl;
//	cout<<"(EBF, 0)    hash="<<Ebf_ebfhash("(EBF, 0) ",0)<<endl;
//	cout<<"(EBF, 1000) hash="<<Ebf_ebfhash("(EBF, 1000)",1000)<<endl;
}

void check_ascii_speed()
{
//	int nsize=10000000;
	int nsize=10000000;
	int nread=0;
	vector<float> x(nsize);
	vector<float> y(nsize);
	linspace(x,0,nsize-1,nsize);
	Timer t1;
	FILE* fp;

   	t1.start();
   	remove("check1.ascii");
	fp=fopen("check1.ascii","w");
	for(int i=0;i<nsize;++i)
		fprintf(fp,"%30f \n",x[i]);
	fclose(fp);
   	t1.printC("ascii write float",nsize);
   	t1.printC("ascii write float",nsize*4);

   	t1.start();
	fp=fopen("check1.ascii","r");
	for(int i=0;i<nsize;++i)
	{
		nread=fscanf(fp,"%f",&y[i]);
//		   	cout<<x[i]<<" "<<y[i]<<endl;
	}
	fclose(fp);
   	t1.printC("ascii read float",nsize);
   	t1.printC("ascii read float",nsize*4);
//   	cout<<y[nsize-1]<<endl;
   	if(nread!=1)
   		printf("%s \n"," ");


// buffering does not help time taken is to convert ascii to float
//   	vector<char> temp(24*nsize,32);
//   	t1.start();
//	for(int i=0;i<nsize;++i)
//		sprintf(&temp[i*24],"%f \n",x[i]);
//	fp=fopen("check.ascii","w");
//	fwrite(&temp[0],24*nsize,1,fp);
//	fclose(fp);
//   	t1.printC("ascii Write float",nsize);
//
//   	t1.start();
//	fp=fopen("check.ascii","r");
//	fread(&temp[0],24*nsize,1,fp);
//	fclose(fp);
//	for(int i=0;i<nsize;++i)
//	{
//		sscanf(&temp[i*24],"%f",&y[i]);
//	   	cout<<x[i]<<" "<<y[i]<<endl;
//	}
//   	t1.printC("ascii Read float",nsize);

}

void test_prog()
{
int lsize=100;
	int nsize=1000000;
vector<float> x1(nsize),x2(nsize),x3(nsize),x4(nsize);
linspace(x1,0,1,nsize);
linspace(x2,0,1,nsize);
linspace(x3,0,1,nsize);
linspace(x4,0,1,nsize);
ebf::Write("col4.ebf","/x1",&x1[0],"w", "",x1.size());
ebf::Write("col4.ebf","/x2",&x2[0],"a", "",x1.size());
ebf::Write("col4.ebf","/x3",&x3[0],"a", "",x1.size());
ebf::Write("col4.ebf","/x4",&x4[0],"a", "",x1.size());
Timer t1;
t1.start();
for (int j=0;j<lsize;++j)
{
	ebf::Read("col4.ebf","/x1",&x1[0],x1.size());
ebf::Read("col4.ebf","/x2",&x2[0],x1.size());
ebf::Read("col4.ebf","/x3",&x3[0],x1.size());
ebf::Read("col4.ebf","/x4",&x4[0],x1.size());
float temp=0.000001;
for(int i=0;i<nsize;++i)
{
	if((x1[i]<temp)&&(x2[i]<temp)&&(x3[i]<temp)&&(x3[i]<temp))
		cout<<x1[i]<<" "<<x2[i]<<" "<<x3[i]<<" "<<x4[i]<<endl;
}
}
t1.printC("sql query",4*nsize*100);

}

int check_ebfvector()
{
	ebf::EbfVector<int> y3,y4,y5;
	int nsize=10000,count,score;
	vector<double> x0(nsize);
    linspace(x0,0,nsize-1,nsize);
    ebf::EbfFile efile;
    int64_t dims[8];
    count=0;
    score=0;

    dims[0]=0;    dims[1]=2;    dims[2]=5;
    efile.Open("check1.ebf","/x1","w",3,"",3,&dims[0]);
    efile.Write(&x0[0],x0.size());
    efile.Close();
    y3.init("check1.ebf.tmp","/x1","w");
//    for(int i=0;i<nsize;++i)
//    	y3.push_back(x0[i]);
//    y3.flush();
//    y3.saveTo("check1.ebf","w");

//    Write("check1.ebf","/x1","w",&x0[0],x0.size());
    ebf::Write("check1.ebf","/x2",&x0[0],"a", "",x0.size());

    y4.init("check1.ebf","/x1");
    y5.init("check1.ebf","/x2");
    vector<int> ind;
    ind.push_back(1010);
    ind.push_back(1001);
    ind.push_back(10);
    ind.push_back(10000-1);
    ind.push_back(10000-10);
    ind.push_back(10);
    ind.push_back(0);
    ind.push_back(10000-1);
    ind.push_back(1001);
   int  status=1;
    for(size_t i=0;i<ind.size();++i)
    {
    	if((int(x0[ind[i]])!=y4[ind[i]])||(int(x0[ind[i]])!=y5[ind[i]]))
    		status=0;
    }
    printf("%s \n","Ebf vector read/write test");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;

    status=0;
    int x=0;
	try	{x=y4[nsize+1];} catch (exception &e){ status=1;}
	try	{x=y4[-1];} catch (exception &e){ status=1;}

    printf("%s %d \n","Ebf vector out of range check",x);
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;

    status=0;
	try	{x=y3[0];} catch (exception &e){ status=1;}
    printf("%s \n","Ebf vector read after write error check");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED");
	count++;
	score+=status;


	int i1=0;
	status=1;
    for(int i=0;i<10;++i)
        for(int j=0;j<2;++j)
            for(int k=0;k<5;++k)
            {
            	if((int(x0[i1])!=y4(i,j,k)))
            		status=0;
            	i1++;
            }

    printf("%s \n","Ebf vector multi dimension index test check");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;

	if(score==count)
		return 1;
	else
		return 0;
}

int check_partial_io_mismatch()
{
	ebf::EbfFile efile;
	int64_t dims[2];
	int status=0;
	double x=0;
	vector<double> y;
	dims[0]=0;
	dims[1]=16;
	efile.Open("check1.ebf","/x1","w",5,"",2,&dims[0]);
	efile.Write(&x,1);
	try{efile.Close();} catch(exception &e){status=1;}
	ebf::Read("check1.ebf","/x1",y);
	if ((y.size()==16)&&(status==1))
		return 1;
	else
	{
		printf("%s \n","Partial I/O check failed");
		return 0;
	}
}

int check_write_to_nonebf()
{
	int status=0;
	double x=0;
	FILE* fd=fopen("check1.ebf","wb");
	printf("%s \n","Checking write to non ebf");

	fwrite(&x,8,1,fd);
	fclose(fd);
	try{ebf::Write("check1.ebf","/x1",&x,"a");} catch(exception &e){status=1;}
	fd=fopen("check1.ebf","rb");
	fseek(fd,0,SEEK_END);
	long mypos=ftell(fd);
	fclose(fd);
	if ((mypos==8)&&(status==1))
		return 1;
	else
	{
		printf("%s %d \n","check write to non ebf file failed",(int)mypos);
		return 0;
	}
}

int check_ebftable()
{
	ebf::EbfDataInfo dinfo;
	int status=1;
	int nsize=10;
	int i,ecode;
	float x1[10];
	float y1[10];
	printf("%s \n","check_ebftable()");
	printf("%s \n","Testing swaped htable1");


	for(i=0;i<nsize;++i)
	{
		x1[i]=i;
		y1[i]=-1;
	}
	ecode=ebf::ebfc::EbfTable_InitSwap("check.ebf");
	ebf::Write("check.ebf","/x1",x1,"a","",nsize);
	ebf::Write("check.ebf","/x2",x1,"a","",nsize);
	ebf::Write("check.ebf","/x3",x1,"a","",nsize);
	ebf::Write("check.ebf","/x4",x1,"a","",nsize);
	ebf::Write("check.ebf","/x5",x1,"a","",nsize);
/*	EbfTable_Print("check.ebf"); */
	ebf::Write("check.ebf","/x6",x1,"a","",nsize);
	ebf::Write("check.ebf","/x7",x1,"a","",nsize);
	ebf::Write("check.ebf","/x8",x1,"a","",nsize);
	ebf::Write("check.ebf","/x9",x1,"a","",nsize);
/*		EbfTable.print("check.ebf"); */
	ecode=ebf::Read("check.ebf","/x1",y1,nsize);
	for (i=0;i<nsize;++i)
		if(x1[i]!=y1[i])
		{
			printf("%d %f %f \n",i,x1[i],y1[i]);
			status=0;
		}
	printf("%d %d\n",ecode,status);

/*		Ebf.info("check.ebf"); */
	printf("%s \n","Testing Ebf_Rename()");
	ecode=ebf::Rename("check.ebf","/x2","/x11");
	ebf::Info("check.ebf");
	if(ebf::ContainsKey("check.ebf","/x2",dinfo))
	{
		status=0;
	}
	if(ebf::ContainsKey("check.ebf","/x11",dinfo))
	{
		ebf::Read("check.ebf","/x11",y1,nsize);
		for (i=0;i<nsize;++i)
			if(x1[i]!=y1[i])
			{
				status=0;
			}
	}
	else
	{
		status=0;
	}


/*		Ebf.info("check.ebf"); */

	ecode=ebf::Rename("check.ebf","/x3","");
	if(ebf::ContainsKey("check.ebf","/x3",dinfo))
		status=0;
	if(ebf::ContainsKey("check.ebf","/.tr/x3.t0",dinfo))
	{
		ebf::Read("check.ebf","/.tr/x3.t0",&y1[0],nsize);
		for (i=0;i<nsize;++i)
			if(x1[i]!=y1[i])
				status=0;
	}
	else
	{
		status=0;
	}

	ebf::Write("check.ebf","/x3",x1,"a","",nsize);
	if(ebf::ContainsKey("check.ebf","/x3",dinfo))
	{
		ebf::Read("check.ebf","/x3",&y1[0],nsize);
		for (i=0;i<nsize;++i)
			if(x1[i]!=y1[i])
				status=0;
	}
	else
	{
		status=0;
	}


	ecode=ebf::Rename("check.ebf","/x3","");
	if(ebf::ContainsKey("check.ebf","/x3",dinfo))
		status=0;
	if(ebf::ContainsKey("check.ebf","/.tr/x3.t1",dinfo))
	{
		ebf::Read("check.ebf","/.tr/x3.t1",&y1[0],nsize);
		for (i=0;i<nsize;++i)
			if(x1[i]!=y1[i])
				status=0;
	}
	else
	{
		status=0;
	}

	if(status==0)
	{
		printf("%s \n","Test FAILED");
		return 0;
	}
	else
	{
		printf("%s \n","Test SUCCESS");
		return 1;
	}


}


int check_header256(const string& dirname)
{
    ebf::EbfFile efile;
    ebf::EbfDataInfo dinfo;
    int i,status;
    int32_t x[10],z[10];
    status=1;
    for(i=0;i<10;++i)
    	z[i]=i;
    x[0]=0;
	printf("%s \n","Testing header256.ebf ");
	string filename=dirname+"header256.ebf";
    if(ebf::ContainsKey(filename.c_str(),"/xvar",dinfo))
    {
    	if(dinfo.elements == 10)
    	{
        	ebf::Read(filename,"/xvar",x,10);
            for(i=0;i<10;++i)
            	if(x[i]!=z[i])
            		status=0;
    	}
    	else
    	{
    		status=0;
    	}
    }
    else
    {
    	status=0;
    }
    if(ebf::ContainsKey(filename.c_str(),"/yvar",dinfo))
    {
    	if(dinfo.elements == 10)
    	{
			ebf::Read(filename,"/yvar",x,10);
			for(i=0;i<10;++i)
				if((x[i]-10)!=z[i])
					status=0;
    	}
    	else
    	{
    		status=0;
    	}
    }
    else
    {
    	status=0;
    }

    if(status==1)
        printf("%s \n","Reading header256.ebf check PASSED");
    else
        printf("%s \n","Reading header256.ebf check FAILED");
    return status;
}

int check_scalar(const string& dirname)
{
    ebf::EbfDataInfo dinfo;
    int status;
    int64_t x1,y1; // is int64 to check non native or else could have been int32
    int64_t x2,y2;
    uint32_t x3,y3;
    status=1;
    x1=-65636;
    x2=-4294967296;
    x3=65636;
	printf("%s \n","Testing scalar.ebf ");
	string filename=dirname+"scalar.ebf";
    if(ebf::ContainsKey(filename,"/x2",dinfo))
    {
    	if(dinfo.elements == 1)
    	{
        	ebf::Read(filename,"/x2",&y1,1);
        	if(x1!=y1)
        		status=0;
    	}
    	else
    	{
    		status=0;
    	}
    }
    else
    {
    	status=0;
    }
    if(ebf::ContainsKey(filename,"/x3",dinfo))
    {
    	if(dinfo.elements == 1)
    	{
			ebf::Read(filename,"/x3",&y2,1);
        	if(x2!=y2)
        		status=0;
    	}
    	else
    	{
    		status=0;
    	}
    }
    else
    {
    	status=0;
    }
    if(ebf::ContainsKey(filename,"/x12",dinfo))
    {
    	if(dinfo.elements == 128)
    	{
			ebf::Read(filename,"/x12",&y3,1);
        	if(x3!=y3)
        		status=0;
    	}
    	else
    	{
    		status=0;
    	}
    }
    else
    {
    	status=0;
    }

    if(status==1)
        printf("%s \n","Reading scalar.ebf check PASSED");
    else
        printf("%s \n","Reading scalar.ebf check FAILED");
    return status;
}


int check_all(const string &dirname)
{
	int count=0;
	int score=0;

	ebf_demo();
	score+=check_types(); count++;
	score+=check_correctness(); count++;
	score+=check_ebftable(); count++;
	score+=check_ebfcopy(); count++;
	score+=check_exceptions(); count++;
	score+=check_ebfvector(); count++;
	score+=check_partial_io_mismatch(); count++;
	score+=check_write_to_nonebf(); count++;
	score+=check_dinfo(); count++;
	if(dirname.size() > 0)
	{
		cout<<dirname<<endl;
		score+=check_masterfile(dirname); count++;
		score+=check_header256(dirname); count++;
		score+=check_scalar(dirname); count++;
	}
//	check_checksum();
//	test_ebf_locate_speed();
//	test_ebf_speed();
//	check_ascii_speed();
//	test_prog();
	printf("%s %i %s %i \n","Score ",score,"/",count);
	if(score==count)
	{
		printf("%s \n","ebf_libc_test: All test PASSED");
		return 0;
	}
	else
	{
		printf("%s \n","ebf_libc_test: Some tests FAILED");
		return 1;
	}
}

struct mystruct1
{
	int x,y,z[4];
};
struct mystruct2
{
	int x,y,z1,z2,z3,z4,z5,z6,z7,z8;//,z[8];
//	mystruct1 point[2];
};

//int check_ebfndarray()
//{
//	int nsize=1000000;
//	vector<mystruct2> check1(nsize),check2(nsize);
//	vector<int> x(nsize);
//	vector<double> y1(nsize);
//	vector<double> y2(nsize);
//	check1[nsize-1].x=10;
//	check1[nsize-1].y=10;
////	check1[3].z[0]=10;
////	check1[3].z[1]=10;
////	check1[3].z[2]=10;
////	check1[3].z[3]=10;
////	check1[3].point[1].z[0]=10;
////	check1[3].point[1].z[1]=10;
////	check1[3].point[1].z[2]=10;
////	check1[3].point[1].z[3]=10;
//	ebfformat::EbfNDArray s1("check",&check1[0],nsize,sizeof(check1[0]),8);
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("x",&(check1[0].x)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("y",&(check1[0].y)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z1",&(check1[0].z1)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z2",&(check1[0].z2)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z3",&(check1[0].z3)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z4",&(check1[0].z4)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z5",&(check1[0].z5)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z6",&(check1[0].z6)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z7",&(check1[0].z7)));
//	s1._fieldarray.push_back(ebfformat::EbfNDArray("z8",&(check1[0].z8)));
////	s1._fieldarray.push_back(ebfformat::EbfNDArray("z",&(check1[0].z[0]),2,sizeof(check1[0].z[0]),8));
////	s1._fieldarray.push_back(ebfformat::EbfNDArray("point",&(check1[0].point[0]),8,sizeof(check1[0].point[0]),2));
////	s1._fieldarray[3]._fieldarray.push_back(ebfformat::EbfNDArray("x",&(check1[0].point[0].x),2,sizeof(check1[0].point[0].x),1));
////	s1._fieldarray[3]._fieldarray.push_back(ebfformat::EbfNDArray("y",&(check1[0].point[0].y),2,sizeof(check1[0].point[0].y),1));
////	s1._fieldarray[3]._fieldarray.push_back(ebfformat::EbfNDArray("z",&(check1[0].point[0].z[0]),2,sizeof(check1[0].point[0].z[0]),4));
//
//	ebfformat::EbfNDArray s2("check",&check2[0],nsize,sizeof(check2[0]),8);
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("x",&(check2[0].x)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("y",&(check2[0].y)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z1",&(check2[0].z1)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z2",&(check2[0].z2)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z3",&(check2[0].z3)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z4",&(check2[0].z4)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z5",&(check2[0].z5)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z6",&(check2[0].z6)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z7",&(check2[0].z7)));
//	s2._fieldarray.push_back(ebfformat::EbfNDArray("z8",&(check2[0].z8)));
////	s2._fieldarray.push_back(ebfformat::EbfNDArray("z",&(check2[0].z[0]),2,sizeof(check2[0].z[0]),8));
////	s2._fieldarray.push_back(ebfformat::EbfNDArray("point",&(check2[0].point[0]),8,sizeof(check2[0].point[0]),2));
////	s2._fieldarray[3]._fieldarray.push_back(ebfformat::EbfNDArray("x",&(check2[0].point[0].x),2,sizeof(check2[0].point[0].x),1));
////	s2._fieldarray[3]._fieldarray.push_back(ebfformat::EbfNDArray("y",&(check2[0].point[0].y),2,sizeof(check2[0].point[0].y),1));
////	s2._fieldarray[3]._fieldarray.push_back(ebfformat::EbfNDArray("z",&(check2[0].point[0].z[0]),2,sizeof(check2[0].point[0].z[0]),4));
//	Timer t1;
//	t1.start();
//	s1.CopyTo(s2);
//	t1.printC("struct speed",nsize*sizeof(check1[0]));
//	cout<<" "<<check1[nsize-1].x<<" "<<check2[nsize-1].y<<endl;
//	cout<<" "<<check1[nsize-1].y<<" "<<check2[nsize-1].y<<endl;
////	cout<<" "<<check1[3].z[0]<<" "<<check2[3].z[0]<<endl;
////	cout<<" "<<check1[3].z[1]<<" "<<check2[3].z[1]<<endl;
////	cout<<" "<<check1[3].z[2]<<" "<<check2[3].z[2]<<endl;
////	cout<<" "<<check1[3].z[3]<<" "<<check2[3].z[3]<<endl;
////	cout<<" "<<check1[3].point[1].z[0]<<" "<<check2[3].point[1].z[0]<<endl;
////	cout<<" "<<check1[3].point[1].z[1]<<" "<<check2[3].point[1].z[1]<<endl;
////	cout<<" "<<check1[3].point[1].z[2]<<" "<<check2[3].point[1].z[2]<<endl;
////	cout<<" "<<check1[3].point[1].z[3]<<" "<<check2[3].point[1].z[3]<<endl;
//
//	exit(1);
//}


//int GetTagNames1(const std::string& filename1, std::vector<std::string> &tagnames)
//{
//	using namespace std;
//	ebf::ebfc::EbfHeader ebfh;
//	tagnames.resize(0);
//	FILE* fd = fopen(filename1.c_str(), "rb");
//	if (fd != NULL)
//	{
//		fseek(fd, 0, SEEK_END);
//		int64_t loc = ftell(fd);
//		fseek(fd, 0, SEEK_SET);
//		while (ftell(fd) < loc)
//		{
//			ebf::ebfc::EbfHeader_read(&ebfh,fd);
//			fseek(fd,ebfh.capacity_,SEEK_CUR);
//			tagnames.push_back(ebfh.name);
//			if(ebfh.ecode!=0)
//			{
//				cout<<"Ebf error: in reading header"<<endl;
//				fclose(fd);
//				return 1;
//			}
//		}
//		fclose(fd);
//	}
//	else
//	{
//		cout<<"Ebf error: file not found"<<endl;
//	}
//	return 0;
//}

//int test_ebf_tagnames()
//{
//	vector<int> x1(10,1);
//	ebf::Write("check.ebf","/x1",&x1[0],"w","",x1.size());
//	ebf::Write("check.ebf","/mydata/x1",&x1[0],"a","",x1.size());
//	ebf::Write("check.ebf","/mydata/x2",&x1[0],"a","",x1.size());
//	ebf::Write("check.ebf","/mydata2/mydata/x1",&x1[0],"a","",x1.size());
//	ebf::Write("check.ebf","/mydata/x3",&x1[0],"a","",x1.size());
//	ebf::Write("check.ebf","/mydata/x4",&x1[0],"a","",x1.size());
//	vector<string> tagnames;
//	int ecode=GetTagNames2("check.ebf",tagnames);
//	cout<<ecode<<endl;
//	for(size_t i=0;i<tagnames.size();++i)
//		cout<<i<<" "<<tagnames[i]<<endl;
//
//	return 0;
//}

int main(int argc,char** argv)
{
	int err1,err2;
//	test_ebf_tagnames();

	std::string dirname="";
//	std::string dirname="/home/sharma/sw/share/ebf/";
	err1=ebf_demo();
	if(argc>1)
		dirname=argv[1];
	err2=check_all(dirname);
	if((err1==0)&&(err2==0))
		return 0;
	else
		return 1;
}
