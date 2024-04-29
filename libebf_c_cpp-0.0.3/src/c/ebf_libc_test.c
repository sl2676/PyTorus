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
#include "ebf_demo.h"
#include <string.h>

#define EBF_NSIZE 128

#define LINSPACE(x,xmin,xmax,size1,temp,i) \
{\
	temp=(xmax-xmin)/(size1-1);\
	if(size1<2)\
	{\
		printf("%s \n","linspace: size should be greater than 1");\
		exit(1);\
	}\
	for(i=0;i<size1;++i)\
		x[i]=xmin+i*temp;\
}



static void Ebf_Info(const char* filename)
{
	EbfHeader ebfh;
	FILE* fd = fopen(filename, "rb");
	int64_t loc;
	int i;
	printf("%-28s %9s %12s %-s \n","name","dtype","units","dims");
	printf("%s \n","--------------------------------------------------------------------------");
	if (fd != NULL)
	{
		fseek(fd, 0, SEEK_END);
		loc = ftell(fd);
		fseek(fd, 0, SEEK_SET);
		while (ftell(fd) < loc)
		{
			EbfHeader_read(&ebfh,fd);
			fseek(fd,ebfh.capacity_,SEEK_CUR);
			printf("%-28s  %9d %-s%10s%s ",ebfh.name,ebfh.datatype,"(",ebfh.unit,")");
			printf("%s "," (");
			for (i = 0; i < ebfh.rank - 1; ++i)
				printf("%d %s",(int)ebfh.dim[i],",");
			printf("%d %s %d \n",(int)ebfh.dim[ebfh.rank - 1],")",((int)ftell(fd)));
			if(ebfh.datatype==8)
				printf("%s \n %s \n","structure definition:",ebfh.sdef);
			if(ebfh.ecode!=0)
			{
				printf("%s \n","Ebf error: in reading header");
				exit(1);
			}
		}
		fclose(fd);
	}
	else
	{
		printf("%s \n","Ebf error: file not found");
	}


}


int check_dinfo()
{
	int status=0;
	int i;
	float x1[10];
	float y1[10];
	EbfDataInfo dinfo;
	printf("%s \n","Testing check_dinfo()");
	for(i=0;i<10;++i)
		x1[i]=i;
	Ebf_WriteFloat32("check.ebf","/x1",x1,"w","km/s",10);
	if(Ebf_ContainsKey_Info("check.ebf","/x1",&dinfo))
	{
		Ebf_ReadFloat32("check.ebf","/x1",&y1[0],10);
		if((dinfo.ecode==0)&&(strcmp(dinfo.unit,"km/s")==0)&&(dinfo.rank==1)&&(dinfo.dim[0]==10))
			status=1;
	}
	Ebf_ContainsKey_Info("check.ebf","/x2",&dinfo);
	if(Ebf_ContainsKey_Info("check.ebf","/x2",&dinfo))
	{
		if((dinfo.ecode==0)&&(strcmp(dinfo.unit,"")==0)&&(dinfo.rank==0)&&(dinfo.headerpos <0))
			status=1;
	}
	return status;
}

int check_types()
{
	int status1=1,status;
	printf("%s \n","Testing getting typecode from given data type");



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


int check_correctness()
{
	/* initialize data with appropriate numbers for each type */
	char mytags[11][100];
	char tempchar[512];
	int i,j;
	int score=0;
	int count=0;
	int status=1;
	int nsize=EBF_NSIZE;
	size_t ilsp;
	double xmin,xmax,templsp;
	char x1[EBF_NSIZE],y1[EBF_NSIZE];
	int32_t x2[EBF_NSIZE],y2[EBF_NSIZE];
	int64_t x3[EBF_NSIZE],y3[EBF_NSIZE];
	efloat32 x4[EBF_NSIZE],y4[EBF_NSIZE];
	efloat64 x5[EBF_NSIZE],y5[EBF_NSIZE];
	efloat64 yd[EBF_NSIZE],yd1[EBF_NSIZE],yd2[EBF_NSIZE];
	int16_t x6[EBF_NSIZE],y6[EBF_NSIZE];
	int8_t x9[EBF_NSIZE],y9[EBF_NSIZE];
	uint8_t x10[EBF_NSIZE],y10[EBF_NSIZE];
	uint16_t x11[EBF_NSIZE],y11[EBF_NSIZE];
	uint32_t x12[EBF_NSIZE],y12[EBF_NSIZE];
	uint64_t x13[EBF_NSIZE],y13[EBF_NSIZE];

	xmin=0; xmax=127;
    LINSPACE(x1,xmin,xmax,nsize,templsp,ilsp)
	xmin=-65636; xmax=-65636-127;
    LINSPACE(x2,xmin,xmax,nsize,templsp,ilsp);
	xmin=-4294967296; xmax=-4294967296-127;
    LINSPACE(x3,xmin,xmax,nsize,templsp,ilsp);
	xmin=1.23e20; xmax=128.23e20;
    LINSPACE(x4,xmin,xmax,nsize,templsp,ilsp);
	xmin=1.23456789e200; xmax=128.23456789e200;
    LINSPACE(x5,xmin,xmax,nsize,templsp,ilsp);
	xmin=-256; xmax=-256-127;
    LINSPACE(x6,xmin,xmax,nsize,templsp,ilsp);
	xmin=-128; xmax=-1;
    LINSPACE(x9,xmin,xmax,nsize,templsp,ilsp);
	xmin=128; xmax=255;
    LINSPACE(x10,xmin,xmax,nsize,templsp,ilsp);
	xmin=256; xmax=256+127;
    LINSPACE(x11,xmin,xmax,nsize,templsp,ilsp);
	xmin=65636; xmax=65636+127;
    LINSPACE(x12,xmin,xmax,nsize,templsp,ilsp);
	xmin=4294967296; xmax=4294967296+127;
    LINSPACE(x13,xmin,xmax,nsize,templsp,ilsp);

    /* write in original types */
   	Ebf_WriteChar("check2.ebf","/x1",&x1[0],"w","",nsize);
   	Ebf_WriteChar("check1.ebf","/x1",&x1[0],"w","",nsize);
   	Ebf_WriteInt32("check1.ebf","/x2",&x2[0],"a","",nsize);
   	Ebf_WriteInt64("check1.ebf","/x3",&x3[0],"a","",nsize);
   	Ebf_WriteFloat32("check1.ebf","/x4",&x4[0],"a","",nsize);
   	Ebf_WriteFloat64("check1.ebf","/x5",&x5[0],"a","",nsize);
   	Ebf_WriteInt16("check1.ebf","/x6",&x6[0],"a","",nsize);
   	Ebf_WriteInt8("check1.ebf","/x9",&x9[0],"a","",nsize);
   	Ebf_WriteUInt8("check1.ebf","/x10",&x10[0],"a","",nsize);
   	Ebf_WriteUInt16("check1.ebf","/x11",&x11[0],"a","",nsize);
   	Ebf_WriteUInt32("check1.ebf","/x12",&x12[0],"a","",nsize);
   	Ebf_WriteUInt64("check1.ebf","/x13",&x13[0],"a","",nsize);
   	Ebf_WriteInt32("check2.ebf","/x2",&x2[0],"a","",nsize);

   	Ebf_WriteChar("check1.ebf","/single/x1",&x1[0],"a","",1);
   	Ebf_WriteInt32("check1.ebf","/single/x2",&x2[0],"a","",1);
   	Ebf_WriteInt64("check1.ebf","/single/x3",&x3[0],"a","",1);
   	Ebf_WriteFloat32("check1.ebf","/single/x4",&x4[0],"a","",1);
   	Ebf_WriteFloat64("check1.ebf","/single/x5",&x5[0],"a","",1);
   	Ebf_WriteInt16("check1.ebf","/single/x6",&x6[0],"a","",1);
   	Ebf_WriteInt8("check1.ebf","/single/x9",&x9[0],"a","",1);
   	Ebf_WriteUInt8("check1.ebf","/single/x10",&x10[0],"a","",1);
   	Ebf_WriteUInt16("check1.ebf","/single/x11",&x11[0],"a","",1);
   	Ebf_WriteUInt32("check1.ebf","/single/x12",&x12[0],"a","",1);
   	Ebf_WriteUInt64("check1.ebf","/single/x13",&x13[0],"a","",1);
    /* write with conversion to double for scalars */
   	Ebf_WriteCharAs(5, "check1.ebf","/AsDouble/single/x1",&x1[0],"a","",1);
   	Ebf_WriteInt32As(5, "check1.ebf","/AsDouble/single/x2",&x2[0],"a","",1);
   	Ebf_WriteInt64As(5, "check1.ebf","/AsDouble/single/x3",&x3[0],"a","",1);
   	Ebf_WriteFloat32As(5, "check1.ebf","/AsDouble/single/x4",&x4[0],"a","",1);
   	Ebf_WriteFloat64As(5, "check1.ebf","/AsDouble/single/x5",&x5[0],"a","",1);
   	Ebf_WriteInt16As(5, "check1.ebf","/AsDouble/single/x6",&x6[0],"a","",1);
   	Ebf_WriteInt8As(5, "check1.ebf","/AsDouble/single/x9",&x9[0],"a","",1);
   	Ebf_WriteUInt8As(5, "check1.ebf","/AsDouble/single/x10",&x10[0],"a","",1);
   	Ebf_WriteUInt16As(5, "check1.ebf","/AsDouble/single/x11",&x11[0],"a","",1);
   	Ebf_WriteUInt32As(5, "check1.ebf","/AsDouble/single/x12",&x12[0],"a","",1);
   	Ebf_WriteUInt64As(5, "check1.ebf","/AsDouble/single/x13",&x13[0],"a","",1);
   	Ebf_WriteInt64("check2.ebf","/x3",&x3[0],"a","",nsize);

    /* write with conversion to double but for vector */
   	Ebf_WriteCharAs(5, "check1.ebf","/AsDouble/x1",&x1[0],"a","",nsize);
   	Ebf_WriteInt32As(5, "check1.ebf","/AsDouble/x2",&x2[0],"a","",nsize);
   	Ebf_WriteInt64As(5, "check1.ebf","/AsDouble/x3",&x3[0],"a","",nsize);
   	Ebf_WriteFloat32As(5, "check1.ebf","/AsDouble/x4",&x4[0],"a","",nsize);
   	Ebf_WriteFloat64As(5, "check1.ebf","/AsDouble/x5",&x5[0],"a","",nsize);
   	Ebf_WriteInt16As(5, "check1.ebf","/AsDouble/x6",&x6[0],"a","",nsize);
   	Ebf_WriteInt8As(5, "check1.ebf","/AsDouble/x9",&x9[0],"a","",nsize);
   	Ebf_WriteUInt8As(5, "check1.ebf","/AsDouble/x10",&x10[0],"a","",nsize);
   	Ebf_WriteUInt16As(5, "check1.ebf","/AsDouble/x11",&x11[0],"a","",nsize);
   	Ebf_WriteUInt32As(5, "check1.ebf","/AsDouble/x12",&x12[0],"a","",nsize);
   	Ebf_WriteUInt64As(5, "check1.ebf","/AsDouble/x13",&x13[0],"a","",nsize);

    /* check scalar and array read without type conversion */
	{
		i=0;
		y1[i]=y2[i]=y3[i]=y4[i]=y5[i]=y6[i]=y9[i]=y10[i]=y11[i]=y12[i]=y13[i]=0;
	}
   	Ebf_WriteFloat32("check2.ebf","/x4",&x4[0], "a", "",nsize);
   	Ebf_ReadChar("check1.ebf","/x1",&y1[0],1);
   	Ebf_ReadInt32("check1.ebf","/x2",&y2[0],1);
   	Ebf_ReadInt64("check1.ebf","/x3",&y3[0],1);
   	Ebf_ReadFloat32("check1.ebf","/x4",&y4[0],1);
   	Ebf_ReadFloat64("check1.ebf","/x5",&y5[0],1);
   	Ebf_ReadInt16("check1.ebf","/x6",&y6[0],1);
   	Ebf_ReadInt8("check1.ebf","/x9",&y9[0],1);
   	Ebf_ReadUInt8("check1.ebf","/x10",&y10[0],1);
   	Ebf_ReadUInt16("check1.ebf","/x11",&y11[0],1);
   	Ebf_ReadUInt32("check1.ebf","/x12",&y12[0],1);
   	Ebf_ReadUInt64("check1.ebf","/x13",&y13[0],1);
	{
   		i=0;
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
		i=0;
		y1[i]=y2[i]=y3[i]=y4[i]=y5[i]=y6[i]=y9[i]=y10[i]=y11[i]=y12[i]=y13[i]=0;
	}
   	Ebf_WriteFloat64("check2.ebf","/x5",&x5[0],"a", "",nsize);
   	Ebf_ReadChar("check1.ebf","/single/x1",&y1[0],1);
   	Ebf_ReadInt32("check1.ebf","/single/x2",&y2[0],1);
   	Ebf_ReadInt64("check1.ebf","/single/x3",&y3[0],1);
   	Ebf_ReadFloat32("check1.ebf","/single/x4",&y4[0],1);
   	Ebf_ReadFloat64("check1.ebf","/single/x5",&y5[0],1);
   	Ebf_ReadInt16("check1.ebf","/single/x6",&y6[0],1);
   	Ebf_ReadInt8("check1.ebf","/single/x9",&y9[0],1);
   	Ebf_ReadUInt8("check1.ebf","/single/x10",&y10[0],1);
   	Ebf_ReadUInt16("check1.ebf","/single/x11",&y11[0],1);
   	Ebf_ReadUInt32("check1.ebf","/single/x12",&y12[0],1);
   	Ebf_ReadUInt64("check1.ebf","/single/x13",&y13[0],1);
	{
   		i=0;
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


	for(i=0;i<nsize;++i)
	{
		y1[i]=y2[i]=y3[i]=y4[i]=y5[i]=y6[i]=y9[i]=y10[i]=y11[i]=y12[i]=y13[i]=0;
	}
   	Ebf_WriteInt16("check2.ebf","/x6",&x6[0],"a", "",nsize);
   	Ebf_ReadChar("check1.ebf","/x1",&y1[0],nsize);
   	Ebf_ReadInt32("check1.ebf","/x2",&y2[0],nsize);
   	Ebf_ReadInt64("check1.ebf","/x3",&y3[0],nsize);
   	Ebf_ReadFloat32("check1.ebf","/x4",&y4[0],nsize);
   	Ebf_ReadFloat64("check1.ebf","/x5",&y5[0],nsize);
   	Ebf_ReadInt16("check1.ebf","/x6",&y6[0],nsize);
   	Ebf_ReadInt8("check1.ebf","/x9",&y9[0],nsize);
   	Ebf_ReadUInt8("check1.ebf","/x10",&y10[0],nsize);
   	Ebf_ReadUInt16("check1.ebf","/x11",&y11[0],nsize);
   	Ebf_ReadUInt32("check1.ebf","/x12",&y12[0],nsize);
   	Ebf_ReadUInt64("check1.ebf","/x13",&y13[0],nsize);
	for(i=0;i<nsize;++i)
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

    /* check scalar and array read of type converted written data */
   	Ebf_WriteInt8("check2.ebf","/x9",&x9[0],"a", "",nsize);
   	strcpy(mytags[0],"/x1");   	strcpy(mytags[1],"/x2");   	strcpy(mytags[2],"/x3");   	strcpy(mytags[3],"/x4");
   	strcpy(mytags[4],"/x5");   	strcpy(mytags[5],"/x6");   	strcpy(mytags[6],"/x9");   	strcpy(mytags[7],"/x10");
   	strcpy(mytags[8],"/x11");   strcpy(mytags[9],"/x12");   strcpy(mytags[10],"/x13");

   	/*
	vector<string> mytags;
	mytags.push_back("/x1");	mytags.push_back("/x2");	mytags.push_back("/x3");	mytags.push_back("/x4");
	mytags.push_back("/x5");	mytags.push_back("/x6");	mytags.push_back("/x9");	mytags.push_back("/x10");
	mytags.push_back("/x11");	mytags.push_back("/x12");	mytags.push_back("/x13");
   	 */
	for(j=0;j<11;++j)
	{
		for(i=0;i<nsize;++i){			yd1[i]=0;		yd2[i]=0; }
		Ebf_ReadFloat64("check1.ebf",mytags[j],&yd2[0],nsize);
		strcpy(tempchar,"/AsDouble/single");		strcat(tempchar,mytags[j]);
		Ebf_ReadFloat64("check1.ebf",tempchar,&yd1[0],1);
		if(yd1[0]!=yd2[0])	{ status*=0; printf("%s %f %f \n",tempchar,yd1[0],yd2[0]);}
		for(i=0;i<nsize;++i)
			yd1[i]=0;
		strcpy(tempchar,"/asdouble");		strcat(tempchar,mytags[j]);
		Ebf_ReadFloat64("check1.ebf",tempchar,&yd1[0],nsize);
		for(i=0;i<nsize;++i)
		{
			if(yd1[i]!=yd2[i])	{ status*=0;  printf("%s %f %f \n",tempchar,yd1[0],yd2[0]);}
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

   	Ebf_WriteUInt8("check2.ebf","/x10",&x10[0],"a", "",nsize);


	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x1",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x1[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x2",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x2[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x3",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x3[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x4",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x4[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x5",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x5[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x6",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x6[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x9",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x9[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x10",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x10[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x11",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x11[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x12",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x12[i])!=yd[i]) { status*=0;}}
	}
	{
		for(i=0;i<nsize;++i){yd[i]=0;}   	Ebf_ReadFloat64("check1.ebf","/x13",&yd[0],nsize);
		for(i=0;i<nsize;++i) {if(((double)x13[i])!=yd[i]) { status*=0;}}
	}
	printf("%s \n","Testing:");
	printf("%s \n","correctness of read, other type to double");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	count++;
	score+=status;


   	Ebf_WriteUInt16("check2.ebf","/x11",&x11[0],"a", "",nsize);
	{	   	Ebf_ReadFloat64("check1.ebf","/x1",&yd[0],1);		if(((double)x1[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x2",&yd[0],1);		if(((double)x2[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x3",&yd[0],1);		if(((double)x3[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x4",&yd[0],1);	    if(((double)x4[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x5",&yd[0],1);		if(((double)x5[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x6",&yd[0],1);		if(((double)x6[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x9",&yd[0],1);		if(((double)x9[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x10",&yd[0],1);		if(((double)x10[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x11",&yd[0],1);		if(((double)x11[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x12",&yd[0],1);		if(((double)x12[0])!=yd[0]) {status*=0;}	}
	{	   	Ebf_ReadFloat64("check1.ebf","/x13",&yd[0],1);		if(((double)x13[0])!=yd[0]) {status*=0;}	}
   	Ebf_WriteUInt32("check2.ebf","/x12",&x12[0],"a", "",nsize);
   	Ebf_WriteUInt64("check2.ebf","/x13",&x13[0],"a", "",nsize);

	printf("%s \n","Testing:");
	printf("%s \n","correctness of scalar read, other type to double");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	count++;
	score+=status;

	for(j=0;j<11;++j)
	{
		for(i=0;i<nsize;++i){			yd1[i]=0;		yd2[i]=0; }
		Ebf_ReadFloat64("check1.ebf",mytags[j],&yd1[0],nsize);
		Ebf_ReadFloat64("check2.ebf",mytags[j],&yd2[0],nsize);
		for(i=0;i<nsize;++i)
		{
			if(yd1[i]!=yd2[i])	{ status*=0;  printf("%s %s %f %f \n","AsDouble",mytags[j],yd1[i],yd2[i]);}
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


int check_masterfile(const char* dirname)
{
	/* initialize data with appropriate numbers for each type */
	char mytags[11][100];
	char filename1[1204];
	char filename2[1204];
	int i,j;
	int status = 1;
	int nsize = EBF_NSIZE;
	size_t ilsp;
	double xmin, xmax, templsp;
	efloat64 yd1[EBF_NSIZE],yd2[EBF_NSIZE];
	char x1[EBF_NSIZE], y1[EBF_NSIZE], z1[EBF_NSIZE];
	int32_t x2[EBF_NSIZE], y2[EBF_NSIZE], z2[EBF_NSIZE];
	int64_t x3[EBF_NSIZE], y3[EBF_NSIZE], z3[EBF_NSIZE];
	efloat32 x4[EBF_NSIZE], y4[EBF_NSIZE], z4[EBF_NSIZE];
	efloat64 x5[EBF_NSIZE], y5[EBF_NSIZE], z5[EBF_NSIZE];
	int16_t x6[EBF_NSIZE], y6[EBF_NSIZE], z6[EBF_NSIZE];
	int8_t x9[EBF_NSIZE], y9[EBF_NSIZE], z9[EBF_NSIZE];
	uint8_t x10[EBF_NSIZE], y10[EBF_NSIZE], z10[EBF_NSIZE];
	uint16_t x11[EBF_NSIZE], y11[EBF_NSIZE], z11[EBF_NSIZE];
	uint32_t x12[EBF_NSIZE], y12[EBF_NSIZE], z12[EBF_NSIZE];
	uint64_t x13[EBF_NSIZE], y13[EBF_NSIZE], z13[EBF_NSIZE];


	if (1024 > (strlen(dirname)+strlen("/master_test1_swap.ebf")))
    {
    	strcpy(filename1,dirname);
    	strcat(filename1,"/master_test1.ebf");
    	strcpy(filename2,dirname);
    	strcat(filename2,"/master_test1_swap.ebf");
    }
    else
    {
    	printf("%s \n","Ebf Error string name too large for strcpy");
    	exit(1);
    }


	xmin = 0;
	xmax = 127;
	LINSPACE(x1, xmin, xmax, nsize, templsp, ilsp)
	xmin = -65636;
	xmax = -65636 - 127;
	LINSPACE(x2, xmin, xmax, nsize, templsp, ilsp);
	xmin = -4294967296;
	xmax = -4294967296 - 127;
	LINSPACE(x3, xmin, xmax, nsize, templsp, ilsp);
	xmin = 1.23e20;
	xmax = 128.23e20;
	LINSPACE(x4, xmin, xmax, nsize, templsp, ilsp);
	xmin = 1.23456789e200;
	xmax = 128.23456789e200;
	LINSPACE(x5, xmin, xmax, nsize, templsp, ilsp);
	xmin = -256;
	xmax = -256 - 127;
	LINSPACE(x6, xmin, xmax, nsize, templsp, ilsp);
	xmin = -128;
	xmax = -1;
	LINSPACE(x9, xmin, xmax, nsize, templsp, ilsp);
	xmin = 128;
	xmax = 255;
	LINSPACE(x10, xmin, xmax, nsize, templsp, ilsp);
	xmin = 256;
	xmax = 256 + 127;
	LINSPACE(x11, xmin, xmax, nsize, templsp, ilsp);
	xmin = 65636;
	xmax = 65636 + 127;
	LINSPACE(x12, xmin, xmax, nsize, templsp, ilsp);
	xmin = 4294967296;
	xmax = 4294967296 + 127;
	LINSPACE(x13, xmin, xmax, nsize, templsp, ilsp);
	Ebf_ReadChar(filename1, "/x1", &y1[0], nsize);
	Ebf_ReadInt32(filename1, "/x2", &y2[0], nsize);
	Ebf_ReadInt64(filename1, "/x3", &y3[0], nsize);
	Ebf_ReadFloat32(filename1, "/x4", &y4[0], nsize);
	Ebf_ReadFloat64(filename1, "/x5", &y5[0], nsize);
	Ebf_ReadInt16(filename1, "/x6", &y6[0], nsize);
	Ebf_ReadInt8(filename1, "/x9", &y9[0], nsize);
	Ebf_ReadUInt8(filename1, "/x10", &y10[0], nsize);
	Ebf_ReadUInt16(filename1, "/x11", &y11[0], nsize);
	Ebf_ReadUInt32(filename1, "/x12", &y12[0], nsize);
	Ebf_ReadUInt64(filename1, "/x13", &y13[0], nsize);

	Ebf_ReadChar(filename2, "/x1", &z1[0], nsize);
	Ebf_ReadInt32(filename2, "/x2", &z2[0], nsize);
	Ebf_ReadInt64(filename2, "/x3", &z3[0], nsize);
	Ebf_ReadFloat32(filename2, "/x4", &z4[0], nsize);
	Ebf_ReadFloat64(filename2, "/x5", &z5[0], nsize);
	Ebf_ReadInt16(filename2, "/x6", &z6[0], nsize);
	Ebf_ReadInt8(filename2, "/x9", &z9[0], nsize);
	Ebf_ReadUInt8(filename2, "/x10", &z10[0], nsize);
	Ebf_ReadUInt16(filename2, "/x11", &z11[0], nsize);
	Ebf_ReadUInt32(filename2, "/x12", &z12[0], nsize);
	Ebf_ReadUInt64(filename2, "/x13", &z13[0], nsize);

	status = 1;
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x1", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x1[i]) != y1[i])
			{
				status *= 0;
			}
			if ((x1[i]) != z1[i])
			{
				status *= 0;
			}
			if (((double) x1[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in char read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x2", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x2[i]) != y2[i])
			{
				status *= 0;
			}
			if ((x2[i]) != z2[i])
			{
				status *= 0;
			}
			if (((double) x2[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in int32 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x3", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x3[i]) != y3[i])
			{
				status *= 0;
			}
			if ((x3[i]) != z3[i])
			{
				status *= 0;
			}
			if (((double) x3[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in int64 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x4", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x4[i]) != y4[i])
			{
				status *= 0;
			}
			if ((x4[i]) != z4[i])
			{
				status *= 0;
			}
			if (((double) x4[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in float32 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x5", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x5[i]) != y5[i])
			{
				status *= 0;
			}
			if ((x5[i]) != z5[i])
			{
				status *= 0;
			}
			if (((double) x5[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in float64 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x6", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x6[i]) != y6[i])
			{
				status *= 0;
			}
			if ((x6[i]) != z6[i])
			{
				status *= 0;
			}
			if (((double) x6[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in int16 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x9", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x9[i]) != y9[i])
			{
				status *= 0;
			}
			if ((x9[i]) != z9[i])
			{
				status *= 0;
			}
			if (((double) x9[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in int8 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x10", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x10[i]) != y10[i])
			{
				status *= 0;
			}
			if ((x10[i]) != z10[i])
			{
				status *= 0;
			}
			if (((double) x10[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in uint8 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x11", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x11[i]) != y11[i])
			{
				status *= 0;
			}
			if ((x11[i]) != z11[i])
			{
				status *= 0;
			}
			if (((double) x11[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in uint16 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x12", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x12[i]) != y12[i])
			{
				status *= 0;
			}
			if ((x12[i]) != z12[i])
			{
				status *= 0;
			}
			if (((double) x12[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in uint32 read");
	}
	{
		for (i = 0; i < nsize; ++i)
		{
			yd1[i] = 0;
		}
		Ebf_ReadFloat64(filename1, "/x13", &yd1[0], nsize);
		for (i = 0; i < nsize; ++i)
		{
			if ((x13[i]) != y13[i])
			{
				status *= 0;
			}
			if ((x13[i]) != z13[i])
			{
				status *= 0;
			}
			if (((double) x13[i]) != yd1[i])
			{
				status *= 0;
			}
		}
		if (status == 0)
			printf("%s \n", "Problem in uint64 read");
	}


   	strcpy(mytags[0],"/x1");   	strcpy(mytags[1],"/x2");   	strcpy(mytags[2],"/x3");   	strcpy(mytags[3],"/x4");
   	strcpy(mytags[4],"/x5");   	strcpy(mytags[5],"/x6");   	strcpy(mytags[6],"/x9");   	strcpy(mytags[7],"/x10");
   	strcpy(mytags[8],"/x11");   strcpy(mytags[9],"/x12");   strcpy(mytags[10],"/x13");
	for(j=0;j<11;++j)
	{
		for(i=0;i<nsize;++i){			yd1[i]=0;		yd2[i]=0; }
		Ebf_ReadFloat64(filename1,mytags[j],&yd1[0],nsize);
		Ebf_ReadFloat64(filename2,mytags[j],&yd2[0],nsize);
		for(i=0;i<nsize;++i)
		{
			if(yd1[i]!=yd2[i])	{ status*=0;  printf(" %s %f %f \n",mytags[j],yd1[i],yd2[i]);}
		}
	}



	printf("%s \n", "Testing: master_test1.ebf");
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");

	return status;
}

int check_ebftable()
{
	EbfDataInfo dinfo;
	int status=1;
	int nsize=10;
	int i,ecode;
	float x1[10];
	float y1[10];
	printf("%s \n","check_ebftable()");
	printf("%s \n","Testing swaped htable");
	for(i=0;i<nsize;++i)
	{
		x1[i]=i;
	}
	ecode=EbfTable_InitSwap("check.ebf");
	Ebf_WriteFloat32("check.ebf","/x1",x1,"a","",nsize);
	Ebf_WriteFloat32("check.ebf","/x2",x1,"a","",nsize);
	Ebf_WriteFloat32("check.ebf","/x3",x1,"a","",nsize);
	Ebf_WriteFloat32("check.ebf","/x4",x1,"a","",nsize);
	Ebf_WriteFloat32("check.ebf","/x5",x1,"a","",nsize);
/*	EbfTable_Print("check.ebf"); */
	Ebf_WriteFloat32("check.ebf","/x6",x1,"a","",nsize);
	Ebf_WriteFloat32("check.ebf","/x7",x1,"a","",nsize);
	Ebf_WriteFloat32("check.ebf","/x8",x1,"a","",nsize);
	Ebf_WriteFloat32("check.ebf","/x9",x1,"a","",nsize);
/*		EbfTable.print("check.ebf"); */
	Ebf_ReadFloat32("check.ebf","/x1",y1,nsize);
	for (i=0;i<nsize;++i)
		if(x1[i]!=y1[i])
			status=0;

/*		Ebf.info("check.ebf"); */
	printf("%s \n","Testing Ebf_Rename()");
	ecode=Ebf_Rename("check.ebf","/x2","/x11");
	if(Ebf_ContainsKey_Info("check.ebf","/x2",&dinfo))
		status=0;
	if(Ebf_ContainsKey_Info("check.ebf","/x11",&dinfo))
	{
		Ebf_ReadFloat32("check.ebf","/x11",y1,nsize);
		for (i=0;i<nsize;++i)
			if(x1[i]!=y1[i])
				status=0;
	}
	else
	{
		status=0;
	}
/*		Ebf.info("check.ebf"); */

	ecode=Ebf_Rename("check.ebf","/x3","");
	if(Ebf_ContainsKey_Info("check.ebf","/x3",&dinfo))
		status=0;
	if(Ebf_ContainsKey_Info("check.ebf","/.tr/x3.t0",&dinfo))
	{
		Ebf_ReadFloat32("check.ebf","/.tr/x3.t0",&y1[0],nsize);
		for (i=0;i<nsize;++i)
			if(x1[i]!=y1[i])
				status=0;
	}
	else
	{
		status=0;
	}

	Ebf_WriteFloat32("check.ebf","/x3",x1,"a","",nsize);
	if(Ebf_ContainsKey_Info("check.ebf","/x3",&dinfo))
	{
		Ebf_ReadFloat32("check.ebf","/x3",&y1[0],nsize);
		for (i=0;i<nsize;++i)
			if(x1[i]!=y1[i])
				status=0;
	}
	else
	{
		status=0;
	}


	ecode=Ebf_Rename("check.ebf","/x3","");
	if(Ebf_ContainsKey_Info("check.ebf","/x3",&dinfo))
		status=0;
	if(Ebf_ContainsKey_Info("check.ebf","/.tr/x3.t1",&dinfo))
	{
		Ebf_ReadFloat32("check.ebf","/.tr/x3.t1",&y1[0],nsize);
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


int check_ebfcopy()
{
	char mytags[11][100];
	char tempchar[512];
	int count=0;
	int score=0;
	int ecode=0;
	EbfDataInfo dinfo1,dinfo2,dinfo3,dinfo4;
	int i,j;
	int status = 1,status1;
	int nsize = EBF_NSIZE;
	size_t ilsp;
	double xmin, xmax, templsp;
	char x1[EBF_NSIZE];
	int32_t x2[EBF_NSIZE];
	int64_t x3[EBF_NSIZE];
	efloat32 x4[EBF_NSIZE];
	efloat64 x5[EBF_NSIZE];
	int16_t x6[EBF_NSIZE];
	int8_t x9[EBF_NSIZE];
	uint8_t x10[EBF_NSIZE];
	uint16_t x11[EBF_NSIZE];
	uint32_t x12[EBF_NSIZE];
	uint64_t x13[EBF_NSIZE];
	efloat64 y1[EBF_NSIZE],y2[EBF_NSIZE],y3[EBF_NSIZE],y4[EBF_NSIZE];

	xmin = 0;
	xmax = 127;
	LINSPACE(x1, xmin, xmax, nsize, templsp, ilsp)
	xmin = -65636;
	xmax = -65636 - 127;
	LINSPACE(x2, xmin, xmax, nsize, templsp, ilsp);
	xmin = -4294967296;
	xmax = -4294967296 - 127;
	LINSPACE(x3, xmin, xmax, nsize, templsp, ilsp);
	xmin = 1.23e20;
	xmax = 128.23e20;
	LINSPACE(x4, xmin, xmax, nsize, templsp, ilsp);
	xmin = 1.23456789e200;
	xmax = 128.23456789e200;
	LINSPACE(x5, xmin, xmax, nsize, templsp, ilsp);
	xmin = -256;
	xmax = -256 - 127;
	LINSPACE(x6, xmin, xmax, nsize, templsp, ilsp);
	xmin = -128;
	xmax = -1;
	LINSPACE(x9, xmin, xmax, nsize, templsp, ilsp);
	xmin = 128;
	xmax = 255;
	LINSPACE(x10, xmin, xmax, nsize, templsp, ilsp);
	xmin = 256;
	xmax = 256 + 127;
	LINSPACE(x11, xmin, xmax, nsize, templsp, ilsp);
	xmin = 65636;
	xmax = 65636 + 127;
	LINSPACE(x12, xmin, xmax, nsize, templsp, ilsp);
	xmin = 4294967296;
	xmax = 4294967296 + 127;
	LINSPACE(x13, xmin, xmax, nsize, templsp, ilsp);

	printf("%s \n","Testing ebfcopy ");
   	Ebf_WriteChar("check1.ebf","/x1",&x1[0],"w","",nsize);
   	Ebf_Info("check1.ebf");
   	Ebf_WriteInt32("check1.ebf","/x2",&x2[0],"a","",nsize);
   	Ebf_WriteInt64("check1.ebf","/x3",&x3[0],"a","",nsize);
   	Ebf_WriteFloat32("check1.ebf","/x4",&x4[0],"a","",nsize);
   	Ebf_WriteFloat64("check1.ebf","/x5",&x5[0],"a","",nsize);
   	Ebf_WriteInt16("check1.ebf","/x6",&x6[0],"a","",nsize);
   	Ebf_WriteInt8("check1.ebf","/x9",&x9[0],"a","",nsize);
   	Ebf_WriteUInt8("check1.ebf","/x10",&x10[0],"a","",nsize);
   	Ebf_WriteUInt16("check1.ebf","/x11",&x11[0],"a","",nsize);
   	Ebf_WriteUInt32("check1.ebf","/x12",&x12[0],"a","",nsize);
   	Ebf_WriteUInt64("check1.ebf","/x13",&x13[0],"a","",nsize);


   	Ebf_WriteChar("check2.ebf","/test/x1",&x1[0],"w","",nsize);
   	Ebf_WriteInt32("check2.ebf","/test/x2",&x2[0],"a","",nsize);
   	Ebf_WriteInt64("check2.ebf","/test/x3",&x3[0],"a","",nsize);
   	Ebf_WriteFloat32("check2.ebf","/test/x4",&x4[0],"a","",nsize);
   	Ebf_WriteFloat64("check2.ebf","/test/x5",&x5[0],"a","",nsize);
   	Ebf_WriteInt16("check2.ebf","/test/x6",&x6[0],"a","",nsize);
   	Ebf_WriteInt8("check2.ebf","/test/x9",&x9[0],"a","",nsize);
   	Ebf_WriteUInt8("check2.ebf","/test/x10",&x10[0],"a","",nsize);
   	Ebf_WriteUInt16("check2.ebf","/test/x11",&x11[0],"a","",nsize);
   	Ebf_WriteUInt32("check2.ebf","/test/x12",&x12[0],"a","",nsize);
   	Ebf_WriteUInt64("check2.ebf","/test/x13",&x13[0],"a","",nsize);



   	ecode=Ebf_Copy("check1.ebf","check3.ebf","w","");
   	ecode=Ebf_Copy("check2.ebf","check3.ebf","a","");

   	strcpy(mytags[0],"/x1");   	strcpy(mytags[1],"/x2");   	strcpy(mytags[2],"/x3");   	strcpy(mytags[3],"/x4");
   	strcpy(mytags[4],"/x5");   	strcpy(mytags[5],"/x6");   	strcpy(mytags[6],"/x9");   	strcpy(mytags[7],"/x10");
   	strcpy(mytags[8],"/x11");   strcpy(mytags[9],"/x12");   strcpy(mytags[10],"/x13");
   	strcpy(tempchar,mytags[0]);
   	strcat(tempchar," ");
   	strcat(tempchar,mytags[1]);
   	strcat(tempchar," ");
   	strcat(tempchar,mytags[2]);
   	ecode=Ebf_Copy("check1.ebf","check4.ebf","w",tempchar);
	for(j=3;j<11;++j)
	{
	   	ecode=Ebf_Copy("check1.ebf","check4.ebf","a",mytags[j]);
	}


	printf("%s \n","testing exception in ebfcopy");
	{status1=0; ecode=Ebf_Copy("check1.ebf","check1.ebf","a","/x1");} 	if(ecode!=0){status1=1;}
	status*=status1;
	{status1=0; ecode=Ebf_Copy("check1.ebf","check3.ebf","a","/x1");} if(ecode!=0){status1=1;}
	status*=status1;
	{status1=0; ecode=Ebf_Copy("check1.ebf","check3.ebf","a","");} if(ecode!=0){status1=1;}
	status*=status1;
	{status1=0; ecode=Ebf_Copy("check1.ebf","check3.ebf","a","/abc");} if(ecode!=0){status1=1;}
	status*=status1;
	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;


	for(j=0;j<11;++j)
	{
		for(i=0;i<nsize;++i)
		{
			y1[i]=0.0;			y2[i]=0.0;			y3[i]=0.0;			y4[i]=0.0;
		}

		Ebf_ReadFloat64("check1.ebf",mytags[j],&y1[0],nsize);
		Ebf_ReadFloat64("check3.ebf",mytags[j],&y2[0],nsize);
		strcpy(tempchar,"/test");		strcpy(tempchar,mytags[j]);
		Ebf_ReadFloat64("check3.ebf",tempchar,&y3[0],nsize);
		Ebf_ReadFloat64("check4.ebf",mytags[j],&y4[0],nsize);

		Ebf_ContainsKey_Info("check1.ebf",mytags[j],&dinfo1);
		Ebf_ContainsKey_Info("check3.ebf",mytags[j],&dinfo2);
		Ebf_ContainsKey_Info("check3.ebf",tempchar,&dinfo3);
		Ebf_ContainsKey_Info("check4.ebf",mytags[j],&dinfo4);



		if((dinfo1.elements != nsize)||(dinfo1.elements!= nsize)||(dinfo1.elements!= nsize)||(dinfo1.elements!= nsize))
		{
			status*=0;
			printf("%s %s \n","Size mismatch tag=",mytags[j]);
			break;
		}
		for(i=0;i<nsize;++i)
		{
			if((y1[i] != y2[i])||(y1[i] != y3[i])||(y1[i] != y4[i]))
			{
				status*=0;
				printf("%s %s %d %f %f %f %f \n","value mismatch tag=",mytags[j],i,y1[i],y2[i],y3[i],y4[i]);
				break;
			}
		}


		if(status==0)
			break;
	}

	if(status==0)
		printf("%s %d \n","Checking: double status=",status);

	if (status == 1)
		printf("%s \n", "Test SUCCESS ");
	else
		printf("%s \n", "Test FAILED ");
	count++;
	score+=status;
	printf("%d %d \n", score,status);

	if(score==count)
		return 1;
	else
		return 0;


}


int check_exceptions()
{
   	EbfFile efile=EbfFile_Create();
	int nsize=100,ecode;
	double y1[100],y2[100];
   	int status=1;
   	int status1=1;
	int y3[100];
	size_t ilsp;
	double xmin, xmax, templsp;
	printf("%s \n","Testing exception: overwrite protection");
	printf("%s \n"," and trying to read past end of record");
	xmin = 0.0;	xmax = nsize-1;
	LINSPACE(y1, xmin, xmax, nsize, templsp, ilsp)

	printf("%s \n","Testing exceptions: non existent files");
	{status1++; ecode=Ebf_ReadFloat64("xxx.ebf", "/x1", &y1[0],1);	} if(ecode!=0)	{		status++;	}
  	 printf("%s \n","check overwrite protection");

   	Ebf_WriteFloat64("check1.ebf","/x0",&y1[0],"w","",nsize);
   	Ebf_WriteFloat64("check1.ebf","/x1",&y1[0],"a","",nsize);
   	Ebf_WriteFloat64("check1.ebf","/x2",&y1[0],"a","",nsize);
   	Ebf_WriteFloat64("check1.ebf","/x3",&y1[0],"a","",nsize);
/*   	 printf("%s \n","check overwrite protection"); */
   	{status1++; ecode=Ebf_WriteFloat64("check1.ebf","/x2",&y1[0],"a","",nsize);} if(ecode!=0){status++;}


   	Ebf_ReadFloat64("check1.ebf","/x1",&y2[0],1);
   	Ebf_WriteFloat64("check1.ebf","/x0",&y1[0],"w","",1);
   	Ebf_WriteFloat64("check1.ebf","/x1",&y1[0],"a","",1);
   	Ebf_WriteFloat64("check1.ebf","/x2",&y1[0],"a","",1);
   	Ebf_WriteFloat64("check1.ebf","/x3",&y1[0],"a","",1);
   	 printf("%s \n","check overwrite protection");
   	{status1++; ecode=Ebf_WriteFloat64("check1.ebf","/x2",&y1[0],"a","",1);} if(ecode!=0){status++;}

   	Ebf_ReadFloat64("check1.ebf","/x1",&y2[0],1);

/*   	printf("%s \n","check overwrite protection with multiple writes in other files in between"); */
   	Ebf_WriteFloat64("check2.ebf","/x0",&y1[0],"w","",1);
   	Ebf_WriteFloat64("check2.ebf","/x1",&y1[0],"a","",1);
   	Ebf_WriteFloat64("check2.ebf","/x2",&y1[0],"a","",1);
   	Ebf_WriteFloat64("check2.ebf","/x3",&y1[0],"a","",1);
   	{status1++; ecode=Ebf_WriteFloat64("check2.ebf","/x1",&y1[0],"a","",1);} if(ecode!=0){status++;}
   	{status1++; ecode=Ebf_WriteFloat64("check1.ebf","/x2",&y1[0],"a","",1);} if(ecode!=0){status++;}
   	Ebf_ReadFloat64("check1.ebf","/x1",&y2[0],1);
   	{status1++; ecode=Ebf_WriteFloat64("check2.ebf","/x1",&y1[0],"a","",1);} if(ecode!=0){status++;}
  	{status1++; ecode=Ebf_WriteFloat64("check1.ebf","/x2",&y1[0],"a","",1);} if(ecode!=0){status++;}

   	Ebf_WriteFloat64("check3.ebf","/x0",&y1[0],"w","",1);
   	Ebf_ReadFloat64("check2.ebf","/x1",&y2[0],1);
   	Ebf_WriteFloat64("check3.ebf","/x1",&y1[0],"a","",1);
   	Ebf_ReadFloat64("check2.ebf","/x1",&y2[0],1);
   	Ebf_ReadFloat64("check1.ebf","/x1",&y2[0],1);
   	Ebf_WriteFloat64("check3.ebf","/x2",&y1[0],"a","",1);
   	Ebf_ReadFloat64("check1.ebf","/x1",&y2[0],1);
   	Ebf_ReadFloat64("check2.ebf","/x1",&y2[0],1);
   	Ebf_WriteFloat64("check3.ebf","/x3",&y1[0],"a","",1);
   	Ebf_ReadFloat64("check2.ebf","/x1",&y2[0],1);
   	{status1++; ecode=Ebf_WriteFloat64("check3.ebf","/x2",&y1[0],"a","",1);} if(ecode!=0){status++;}

/*   	printf("%s \n","check trying to read past end of record"); */
   	Ebf_WriteFloat64("check1.ebf","/x0",&y1[0],"w","",nsize-1);
   	Ebf_WriteFloat64("check1.ebf","/x1",&y1[0],"a","",nsize-1);
   	{status1++; ecode=Ebf_ReadFloat64("check2.ebf","/x1",&y2[0],nsize);} if(ecode!=0){status++;}
   	{status1++; ecode=Ebf_ReadInt32("check2.ebf","/x1",&y3[0],nsize);} if(ecode!=0){status++;}

/*   	printf("%s \n","check trying to read past end of record with sequential read and with different data types"); */

   	EbfFile_OpenR(&efile,"check1.ebf","/x1");
   	EbfFile_rFloat64(&efile,&y2[0],10);
   	EbfFile_rFloat64(&efile,&y2[0],1);
   	EbfFile_rFloat64(&efile,&y2[0],9);
   	EbfFile_rFloat64(&efile,&y2[0],79);
   	{status1++;  EbfFile_rFloat64(&efile,&y2[0],1);   	} if(efile.ecode!=0){status++;}
    EbfFile_Close(&efile);

   	EbfFile_OpenR(&efile,"check1.ebf","/x1");
   	EbfFile_rInt32(&efile, &y3[0],10);
   	EbfFile_rInt32(&efile,&y3[0],1);
   	EbfFile_rInt32(&efile, &y3[0],9);
   	EbfFile_rInt32(&efile, &y3[0], 80);
   	{status1++;  EbfFile_rInt32(&efile,&y3[0],1); }  if(efile.ecode!=0) {status++;}
   	EbfFile_Close(&efile);


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


void check_checksum()
{
	int nsize=10;
	float x1[10];
	double x2[10];
	int64_t checksum=0;

	Ebf_WriteFloat32("check1.ebf","/x1",&x1[0],"w","",nsize);
	Ebf_WriteFloat64("check1.ebf","/x2",&x2[0],"a","",nsize);
	Ebf_WriteFloat32("check1.ebf","/single/x1",&x1[0],"a","",1);
	Ebf_WriteFloat64("check1.ebf","/single/x2",&x2[0],"a","",1);
	Ebf_ReadInt64("check1.ebf","/.ebf/info",&checksum,1);
	printf("%s %"PRId64"\n","checksum=",checksum);
/*
//	cout<<"(EBF, 0)    hash="<<Ebf_ebfhash("(EBF, 0) ",0)<<endl;
//	cout<<"(EBF, 1000) hash="<<Ebf_ebfhash("(EBF, 1000)",1000)<<endl;
 */
}

int check_partial_io_mismatch()
{
	EbfFile efile=EbfFile_Create();
	EbfDataInfo dinfo;
	int64_t dims[2];
	double x=0;
	dims[0]=0;
	dims[1]=16;
	EbfFile_Open(&efile,"check1.ebf","/x1","w",5,"",2,&dims[0]);
	EbfFile_wFloat64(&efile,&x,1);
	EbfFile_Close(&efile);
	Ebf_ContainsKey_Info("check1.ebf","/x1",&dinfo);
	if ((dinfo.elements==16)&&(efile.ecode!=0))
		return 1;
	else
	{
		printf("%s %d \n","Partial I/O check failed ",(int)dinfo.elements);
		return 0;
	}
}

int check_write_to_nonebf()
{
	double x=0;
	int ecode=0;
	int64_t mypos;
	FILE* fd=fopen("check1.ebf","wb");
	printf("%s \n","Checking write to non ebf");

	fwrite(&x,8,1,fd);
	fclose(fd);
	ecode=Ebf_WriteFloat64("check1.ebf","/x1",&x,"a","",1);
	fd=fopen("check1.ebf","rb");
	fseek(fd,0,SEEK_END);
	mypos=ftell(fd);
	fclose(fd);
	if ((mypos==8)&&(ecode!=0))
		return 1;
	else
	{
		printf("%s %d \n","check write to non ebf file failed",(int)mypos);
		return 0;
	}
}

int check_header256(const char* dirname)
{
    EbfDataInfo dinfo;
    int i,status;
    char filename[1024];
    int32_t x[10],z[10];
    status=1;
    for(i=0;i<10;++i)
    	z[i]=i;
    x[0]=0;
    if (1024 > (strlen(dirname)+strlen("header256.ebf")))
    {
    	strcpy(filename,dirname);
    	strcat(filename,"header256.ebf");
    }
    else
    {
    	printf("%s \n","Ebf Error string name too large for strcpy");
    	exit(1);
    }
	printf("%s \n","Testing header256.ebf ");
    if(Ebf_ContainsKey_Info(filename,"/xvar",&dinfo))
    {
    	if(dinfo.elements == 10)
    	{
        	Ebf_ReadInt32(filename,"/xvar",x,10);
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
    if(Ebf_ContainsKey_Info(filename,"/yvar",&dinfo))
    {
    	if(dinfo.elements == 10)
    	{
			Ebf_ReadInt32(filename,"/yvar",x,10);
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


#undef LINSPACE
#undef EBF_NSIZE


int check_scalar(const char* dirname)
{
    EbfDataInfo dinfo;
    char filename[1024];
    int status;
    int32_t x1,y1;
    int64_t x2,y2;
    uint32_t x3,y3;
    status=1;
    x1=-65636;
    x2=-4294967296;
    x3=65636;
	printf("%s \n","Testing scalar.ebf ");
    if (1024 > (strlen(dirname)+strlen("scalar.ebf")))
    {
    	strcpy(filename,dirname);
    	strcat(filename,"scalar.ebf");
    }
    else
    {
    	printf("%s \n","Ebf Error string name too large for strcpy");
    	exit(1);
    }
    if(Ebf_ContainsKey_Info(filename,"/x2",&dinfo))
    {
    	if(dinfo.elements == 1)
    	{
        	Ebf_ReadInt32(filename,"/x2",&y1,1);
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
    if(Ebf_ContainsKey_Info(filename,"/x3",&dinfo))
    {
    	if(dinfo.elements == 1)
    	{
			Ebf_ReadInt64(filename,"/x3",&y2,1);
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
    if(Ebf_ContainsKey_Info(filename,"/x12",&dinfo))
    {
    	if(dinfo.elements == 128)
    	{
			Ebf_ReadUInt32(filename,"/x12",&y3,1);
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


int check_all(const char* dirname)
{
	int count=0;
	int score=0;
	score+=check_dinfo(); count++;
	score+=check_types(); count++;
	score+=check_correctness(); count++;
	score+=check_ebftable(); count++;
	score+=check_ebfcopy(); count++;
	if(strlen(dirname)>0)
	{
		score+=check_masterfile(dirname); count++;
		score+=check_header256(dirname); count++;
		score+=check_scalar(dirname); count++;
	}
/*	score+=check_exceptions(); count++; */
/*	score+=check_partial_io_mismatch(); count++; */
/*	score+=check_write_to_nonebf(); count++; */
	check_checksum();
	/*
	score+=check_ebfvector(); count++;
//	test_ebf_locate_speed();
//	test_ebf_speed();
//	check_ascii_speed();
//	test_prog();
 */
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

int main(int argc,char** argv)
{
	int err1,err2;
	err1=ebf_demo();
	if(argc>1)
		err2=check_all(argv[1]);
	else
		err2=check_all("");
/*		err2=check_all("/home/sharma/sw/share/ebf/"); */

	if((err1==0)&&(err2==0))
		return 0;
	else
		return 1;
}


