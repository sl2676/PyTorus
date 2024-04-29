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

/*----------------------------------------------------------------------*/
#ifdef __cplusplus
#include<stdexcept>
#endif

/*-------------------------------------------------------------*/
static void ebf_tfunc_char_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}


static void ebf_tfunc_int32_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}


static void ebf_tfunc_int64_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}


static void ebf_tfunc_float32_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}


static void ebf_tfunc_float64_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}

static void ebf_tfunc_int16_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}


static void ebf_tfunc_int8_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}


static void ebf_tfunc_uint8_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}

static void ebf_tfunc_uint16_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}


static void ebf_tfunc_uint32_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}


static void ebf_tfunc_uint64_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}



/*-------------------------------------------------------------*/


static void ebfutils_error(const char *s)
{
	printf("%s \n", s);
#ifdef __cplusplus
	throw std::runtime_error("");
#else
	exit(1);
#endif
}

static void ebfutils_assert(int temp, const char* s)
{
	if (temp == 0)
	{
		printf("%s %s \n", "Assertion failed as ", s);
#ifdef __cplusplus
		throw std::runtime_error("");
#else
	exit(1);
#endif
	}
}


static int ebfutils_fileExists(const char* filename)
{
	/*	ifstream infile; */
	FILE* fp = fopen(filename, "rb");
	if (fp != NULL)
	{
		fclose(fp);
		return 1;
	}
	else
	{
		return 0;
	}
}

static void ebfutils_SwapEndian(void* addr, int datasize,int64_t nmemb)
{
	int64_t pattern[3];
	int64_t i, j = 0, k;
	char c;
	pattern[0]=datasize;
	pattern[1]=nmemb;
	pattern[2]=0;
	while (pattern[j] > 0)
	{
		for (k = 0; k < pattern[j + 1]; k++)
		{
			for (i = 0; i < pattern[j] / 2; i++)
			{
				c = *((char*) addr + i);
				*((char*) addr + i) = *((char*) addr + (pattern[j] - i - 1));
				*((char*) addr + (pattern[j] - i - 1)) = c;
			}
			addr = ((char *) addr) + pattern[j];
		}
		j += 2;
	}
}



static int ebfutils_TypeSize(int x)
{
	switch (x)
	{
	case 1:
		return 1;
		break;
	case 2:
		return 4;
		break;
	case 3:
		return 8;
		break;
	case 4:
		return 4;
		break;
	case 5:
		return 8;
		break;
	case 6:
		return 2;
		break;
	case 9:
		return 1;
		break;
	case 10:
		return 1;
		break;
	case 11:
		return 2;
		break;
	case 12:
		return 4;
		break;
	case 13:
		return 8;
		break;
	default:
		ebfutils_error(
				"Ebf Error from ebfutils_getDataSize()-unrecognized data type ");
		return 0;
	}
}

int Ebf_TypeS(const char* typestring)
{
	int etype = 0;
	if (strcmp(typestring, "char") == 0)
	{
		etype = 1;
	}
	else if (strcmp(typestring, "int32") == 0)
	{
		etype = 2;
	}
	else if (strcmp(typestring, "int64") == 0)
	{
		etype = 3;
	}
	else if (strcmp(typestring, "float32") == 0)
	{
		etype = 4;
	}
	else if (strcmp(typestring, "float") == 0)
	{
		etype = 4;
	}
	else if (strcmp(typestring, "float64") == 0)
	{
		etype = 5;
	}
	else if (strcmp(typestring, "double") == 0)
	{
		etype = 5;
	}
	else if (strcmp(typestring, "int16") == 0)
	{
		etype = 6;
	}
	else if (strcmp(typestring, "int8") == 0)
	{
		etype = 9;
	}
	else if (strcmp(typestring, "uint8") == 0)
	{
		etype = 10;
	}
	else if (strcmp(typestring, "uint16") == 0)
	{
		etype = 11;
	}
	else if (strcmp(typestring, "uint32") == 0)
	{
		etype = 12;
	}
	else if (strcmp(typestring, "uint64") == 0)
	{
		etype = 13;
	}
	else
	{
		ebfutils_error("Ebf Error from Ebf_etype(): un recognized data type");
	}
	return etype;

}


/*----------------------------------------------------------------------------------------------*/

void EbfFile_ReadGeneral(EbfFile *efile,char* value1,size_t data_size,int typecode)
{
	static void (*EBFTFUNCARRAY[14][14])(char*, char*,size_t,size_t,size_t)={
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_char_char,
			ebf_tfunc_char_int32,
			ebf_tfunc_char_int64,
			ebf_tfunc_char_float32,
			ebf_tfunc_char_float64,
			ebf_tfunc_char_int16,
			NULL,
			NULL,
			ebf_tfunc_char_int8,
			ebf_tfunc_char_uint8,
			ebf_tfunc_char_uint16,
			ebf_tfunc_char_uint32,
			ebf_tfunc_char_uint64},
			{NULL,
			ebf_tfunc_int32_char,
			ebf_tfunc_int32_int32,
			ebf_tfunc_int32_int64,
			ebf_tfunc_int32_float32,
			ebf_tfunc_int32_float64,
			ebf_tfunc_int32_int16,
			NULL,
			NULL,
			ebf_tfunc_int32_int8,
			ebf_tfunc_int32_uint8,
			ebf_tfunc_int32_uint16,
			ebf_tfunc_int32_uint32,
			ebf_tfunc_int32_uint64},
			{NULL,
			ebf_tfunc_int64_char,
			ebf_tfunc_int64_int32,
			ebf_tfunc_int64_int64,
			ebf_tfunc_int64_float32,
			ebf_tfunc_int64_float64,
			ebf_tfunc_int64_int16,
			NULL,
			NULL,
			ebf_tfunc_int64_int8,
			ebf_tfunc_int64_uint8,
			ebf_tfunc_int64_uint16,
			ebf_tfunc_int64_uint32,
			ebf_tfunc_int64_uint64},
			{NULL,
			ebf_tfunc_float32_char,
			ebf_tfunc_float32_int32,
			ebf_tfunc_float32_int64,
			ebf_tfunc_float32_float32,
			ebf_tfunc_float32_float64,
			ebf_tfunc_float32_int16,
			NULL,
			NULL,
			ebf_tfunc_float32_int8,
			ebf_tfunc_float32_uint8,
			ebf_tfunc_float32_uint16,
			ebf_tfunc_float32_uint32,
			ebf_tfunc_float32_uint64},
			{NULL,
			ebf_tfunc_float64_char,
			ebf_tfunc_float64_int32,
			ebf_tfunc_float64_int64,
			ebf_tfunc_float64_float32,
			ebf_tfunc_float64_float64,
			ebf_tfunc_float64_int16,
			NULL,
			NULL,
			ebf_tfunc_float64_int8,
			ebf_tfunc_float64_uint8,
			ebf_tfunc_float64_uint16,
			ebf_tfunc_float64_uint32,
			ebf_tfunc_float64_uint64},
			{NULL,
			ebf_tfunc_int16_char,
			ebf_tfunc_int16_int32,
			ebf_tfunc_int16_int64,
			ebf_tfunc_int16_float32,
			ebf_tfunc_int16_float64,
			ebf_tfunc_int16_int16,
			NULL,
			NULL,
			ebf_tfunc_int16_int8,
			ebf_tfunc_int16_uint8,
			ebf_tfunc_int16_uint16,
			ebf_tfunc_int16_uint32,
			ebf_tfunc_int16_uint64},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_int8_char,
			ebf_tfunc_int8_int32,
			ebf_tfunc_int8_int64,
			ebf_tfunc_int8_float32,
			ebf_tfunc_int8_float64,
			ebf_tfunc_int8_int16,
			NULL,
			NULL,
			ebf_tfunc_int8_int8,
			ebf_tfunc_int8_uint8,
			ebf_tfunc_int8_uint16,
			ebf_tfunc_int8_uint32,
			ebf_tfunc_int8_uint64},
			{NULL,
			ebf_tfunc_uint8_char,
			ebf_tfunc_uint8_int32,
			ebf_tfunc_uint8_int64,
			ebf_tfunc_uint8_float32,
			ebf_tfunc_uint8_float64,
			ebf_tfunc_uint8_int16,
			NULL,
			NULL,
			ebf_tfunc_uint8_int8,
			ebf_tfunc_uint8_uint8,
			ebf_tfunc_uint8_uint16,
			ebf_tfunc_uint8_uint32,
			ebf_tfunc_uint8_uint64},
			{NULL,
			ebf_tfunc_uint16_char,
			ebf_tfunc_uint16_int32,
			ebf_tfunc_uint16_int64,
			ebf_tfunc_uint16_float32,
			ebf_tfunc_uint16_float64,
			ebf_tfunc_uint16_int16,
			NULL,
			NULL,
			ebf_tfunc_uint16_int8,
			ebf_tfunc_uint16_uint8,
			ebf_tfunc_uint16_uint16,
			ebf_tfunc_uint16_uint32,
			ebf_tfunc_uint16_uint64},
			{NULL,
			ebf_tfunc_uint32_char,
			ebf_tfunc_uint32_int32,
			ebf_tfunc_uint32_int64,
			ebf_tfunc_uint32_float32,
			ebf_tfunc_uint32_float64,
			ebf_tfunc_uint32_int16,
			NULL,
			NULL,
			ebf_tfunc_uint32_int8,
			ebf_tfunc_uint32_uint8,
			ebf_tfunc_uint32_uint16,
			ebf_tfunc_uint32_uint32,
			ebf_tfunc_uint32_uint64},
			{NULL,
			ebf_tfunc_uint64_char,
			ebf_tfunc_uint64_int32,
			ebf_tfunc_uint64_int64,
			ebf_tfunc_uint64_float32,
			ebf_tfunc_uint64_float64,
			ebf_tfunc_uint64_int16,
			NULL,
			NULL,
			ebf_tfunc_uint64_int8,
			ebf_tfunc_uint64_uint8,
			ebf_tfunc_uint64_uint16,
			ebf_tfunc_uint64_uint32,
			ebf_tfunc_uint64_uint64}
	};
	size_t ic,size,size1;

	if((efile->mode != 1)||(efile->fp==NULL)||(efile->ecode!=0))
	{
		efile->ecode=1;
		ebfutils_error(
				"Ebf Error from Ebf::read()-EBF error:file not opened in read mode");
	}
	else
	{
		if (((data_size) * (efile->ebfh.datasize)
				+ ftell(efile->fp)) > efile->ebfh.datapos_end)
		{
			efile->ecode=1;
			ebfutils_error(
					"Ebf Error from Ebf::read()- trying to read past end of record");
		}
		else
		{
			if (efile->ebfh.datatype == typecode)
			{
					if (fread(&value1[0],(efile->ebfh.datasize) * data_size, 1,efile->fp) != 1)
					{
						efile->ecode=1;
						ebfutils_error("Ebf Error from Ebf::read()- fread");
					}
					if(efile->ebfh.flag_swap == 1)
						ebfutils_SwapEndian(&value1[0],efile->ebfh.datasize,data_size);
			}
			else
			{
				ic = 0;
				size1 = sizeof(efile->tbuf.bufchar)/8;
				size = data_size;

				while (ic < size)
				{
					if ((ic + size1) > size)
						size1 = size - ic;
					if (fread(&(efile->tbuf.bufchar[0]),(efile->ebfh.datasize)*size1,1, efile->fp) != 1)
					{
						efile->ecode=1;
						ebfutils_error("Ebf Error from Ebf::read()- fread");
					}
					if (efile->ebfh.flag_swap == 1)
						ebfutils_SwapEndian(efile->tbuf.bufchar,efile->ebfh.datasize,size1);
					EBFTFUNCARRAY[efile->ebfh.datatype][typecode](efile->tbuf.bufchar,value1+ic*ebfutils_TypeSize(typecode),size1,ebfutils_TypeSize(efile->ebfh.datatype),ebfutils_TypeSize(typecode));
					ic = ic + size1;
				}
			}
		}
	}


}



void EbfFile_WriteGeneral(EbfFile *efile,char* value1,size_t data_size,int typecode)
{
	static void (*EBFTFUNCARRAY[14][14])(char*, char*,size_t,size_t,size_t)={
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_char_char,
			ebf_tfunc_char_int32,
			ebf_tfunc_char_int64,
			ebf_tfunc_char_float32,
			ebf_tfunc_char_float64,
			ebf_tfunc_char_int16,
			NULL,
			NULL,
			ebf_tfunc_char_int8,
			ebf_tfunc_char_uint8,
			ebf_tfunc_char_uint16,
			ebf_tfunc_char_uint32,
			ebf_tfunc_char_uint64},
			{NULL,
			ebf_tfunc_int32_char,
			ebf_tfunc_int32_int32,
			ebf_tfunc_int32_int64,
			ebf_tfunc_int32_float32,
			ebf_tfunc_int32_float64,
			ebf_tfunc_int32_int16,
			NULL,
			NULL,
			ebf_tfunc_int32_int8,
			ebf_tfunc_int32_uint8,
			ebf_tfunc_int32_uint16,
			ebf_tfunc_int32_uint32,
			ebf_tfunc_int32_uint64},
			{NULL,
			ebf_tfunc_int64_char,
			ebf_tfunc_int64_int32,
			ebf_tfunc_int64_int64,
			ebf_tfunc_int64_float32,
			ebf_tfunc_int64_float64,
			ebf_tfunc_int64_int16,
			NULL,
			NULL,
			ebf_tfunc_int64_int8,
			ebf_tfunc_int64_uint8,
			ebf_tfunc_int64_uint16,
			ebf_tfunc_int64_uint32,
			ebf_tfunc_int64_uint64},
			{NULL,
			ebf_tfunc_float32_char,
			ebf_tfunc_float32_int32,
			ebf_tfunc_float32_int64,
			ebf_tfunc_float32_float32,
			ebf_tfunc_float32_float64,
			ebf_tfunc_float32_int16,
			NULL,
			NULL,
			ebf_tfunc_float32_int8,
			ebf_tfunc_float32_uint8,
			ebf_tfunc_float32_uint16,
			ebf_tfunc_float32_uint32,
			ebf_tfunc_float32_uint64},
			{NULL,
			ebf_tfunc_float64_char,
			ebf_tfunc_float64_int32,
			ebf_tfunc_float64_int64,
			ebf_tfunc_float64_float32,
			ebf_tfunc_float64_float64,
			ebf_tfunc_float64_int16,
			NULL,
			NULL,
			ebf_tfunc_float64_int8,
			ebf_tfunc_float64_uint8,
			ebf_tfunc_float64_uint16,
			ebf_tfunc_float64_uint32,
			ebf_tfunc_float64_uint64},
			{NULL,
			ebf_tfunc_int16_char,
			ebf_tfunc_int16_int32,
			ebf_tfunc_int16_int64,
			ebf_tfunc_int16_float32,
			ebf_tfunc_int16_float64,
			ebf_tfunc_int16_int16,
			NULL,
			NULL,
			ebf_tfunc_int16_int8,
			ebf_tfunc_int16_uint8,
			ebf_tfunc_int16_uint16,
			ebf_tfunc_int16_uint32,
			ebf_tfunc_int16_uint64},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_int8_char,
			ebf_tfunc_int8_int32,
			ebf_tfunc_int8_int64,
			ebf_tfunc_int8_float32,
			ebf_tfunc_int8_float64,
			ebf_tfunc_int8_int16,
			NULL,
			NULL,
			ebf_tfunc_int8_int8,
			ebf_tfunc_int8_uint8,
			ebf_tfunc_int8_uint16,
			ebf_tfunc_int8_uint32,
			ebf_tfunc_int8_uint64},
			{NULL,
			ebf_tfunc_uint8_char,
			ebf_tfunc_uint8_int32,
			ebf_tfunc_uint8_int64,
			ebf_tfunc_uint8_float32,
			ebf_tfunc_uint8_float64,
			ebf_tfunc_uint8_int16,
			NULL,
			NULL,
			ebf_tfunc_uint8_int8,
			ebf_tfunc_uint8_uint8,
			ebf_tfunc_uint8_uint16,
			ebf_tfunc_uint8_uint32,
			ebf_tfunc_uint8_uint64},
			{NULL,
			ebf_tfunc_uint16_char,
			ebf_tfunc_uint16_int32,
			ebf_tfunc_uint16_int64,
			ebf_tfunc_uint16_float32,
			ebf_tfunc_uint16_float64,
			ebf_tfunc_uint16_int16,
			NULL,
			NULL,
			ebf_tfunc_uint16_int8,
			ebf_tfunc_uint16_uint8,
			ebf_tfunc_uint16_uint16,
			ebf_tfunc_uint16_uint32,
			ebf_tfunc_uint16_uint64},
			{NULL,
			ebf_tfunc_uint32_char,
			ebf_tfunc_uint32_int32,
			ebf_tfunc_uint32_int64,
			ebf_tfunc_uint32_float32,
			ebf_tfunc_uint32_float64,
			ebf_tfunc_uint32_int16,
			NULL,
			NULL,
			ebf_tfunc_uint32_int8,
			ebf_tfunc_uint32_uint8,
			ebf_tfunc_uint32_uint16,
			ebf_tfunc_uint32_uint32,
			ebf_tfunc_uint32_uint64},
			{NULL,
			ebf_tfunc_uint64_char,
			ebf_tfunc_uint64_int32,
			ebf_tfunc_uint64_int64,
			ebf_tfunc_uint64_float32,
			ebf_tfunc_uint64_float64,
			ebf_tfunc_uint64_int16,
			NULL,
			NULL,
			ebf_tfunc_uint64_int8,
			ebf_tfunc_uint64_uint8,
			ebf_tfunc_uint64_uint16,
			ebf_tfunc_uint64_uint32,
			ebf_tfunc_uint64_uint64}
	};
	size_t ic,size,size1;


	if((efile->mode!=3)||(efile->fp==NULL)||(efile->ecode!=0))
	{
		efile->ecode=1;
		ebfutils_error("Ebf Error from Ebf::write()-EBF error: File not opened in write mode or ecode!=0");
	}
	else
	{
			if (efile->ebfh.datatype == typecode)
			{
				if(fwrite(&value1[0],efile->ebfh.datasize * data_size,1,efile->fp)!=1)
				{
					efile->ecode=1;
					ebfutils_error("Ebf Error from Ebf::write()- fwrite error");
				}
			}
			else
			{
				ic = 0;
				size1 = sizeof(efile->tbuf.bufchar)/8;
				size = data_size;
				while (ic < size)
				{
					if ((ic + size1) > size)
						size1 = size - ic;
					EBFTFUNCARRAY[typecode][efile->ebfh.datatype](value1+ic*ebfutils_TypeSize(typecode),efile->tbuf.bufchar,size1,ebfutils_TypeSize(typecode),ebfutils_TypeSize(efile->ebfh.datatype));
					if (fwrite(&(efile->tbuf.bufchar[0]), efile->ebfh.datasize * size1, 1, efile->fp)!= 1)
					{
						efile->ecode=1;
						ebfutils_error("Ebf Error from Ebf::write()-status not 1, file may be corrupted");
					}
					ic = ic + size1;
				}
			}
	}


}

void EbfFile_wChar(EbfFile *efile, const char* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,1);
}
void EbfFile_wInt32(EbfFile *efile, const int32_t* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,2);
}
void EbfFile_wInt64(EbfFile *efile, const int64_t* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,3);
}
void EbfFile_wFloat32(EbfFile *efile, const efloat32* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,4);
}
void EbfFile_wFloat64(EbfFile *efile, const efloat64* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,5);
}
void EbfFile_wInt16(EbfFile *efile, const short* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,6);
}
void EbfFile_wInt8(EbfFile *efile, const int8_t* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,9);
}
void EbfFile_wUInt8(EbfFile *efile, const uint8_t* value1,
		int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,10);
}
void EbfFile_wUInt16(EbfFile *efile, const uint16_t* value1,
		int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,11);
}
void EbfFile_wUInt32(EbfFile *efile, const uint32_t * value1,
		int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,12);
}
void EbfFile_wUInt64(EbfFile *efile, const uint64_t* value1, int64_t data_size)
{
	EbfFile_WriteGeneral(efile,(char *)value1,data_size,13);
}


void EbfFile_rChar(EbfFile *efile, char* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile,(char *)value1,data_size,1);
}
void EbfFile_rInt32(EbfFile *efile, int32_t* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,2);
}
void EbfFile_rInt64(EbfFile *efile, int64_t* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,3);
}
void EbfFile_rFloat32(EbfFile *efile, efloat32* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,4);
}
void EbfFile_rFloat64(EbfFile *efile, efloat64* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,5);
}
void EbfFile_rInt16(EbfFile *efile, int16_t * value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,6);
}
void EbfFile_rInt8(EbfFile *efile, int8_t* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,9);
}
void EbfFile_rUInt8(EbfFile *efile, uint8_t* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,10);
}

void EbfFile_rUInt16(EbfFile *efile, uint16_t * value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,11);
}
void EbfFile_rUInt32(EbfFile *efile, uint32_t* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,12);
}
void EbfFile_rUInt64(EbfFile *efile, uint64_t* value1, int64_t data_size)
{
	EbfFile_ReadGeneral(efile, (char *)value1,data_size,13);
}

void EbfFile_rVoid(EbfFile *efile, void* data, int64_t data_size)
{
	if (efile->mode != 1)
	{
		ebfutils_error(
				"Ebf Error from Ebf::readBytes-Reading not on or status is zero");
	}

	if ((efile->ebfh.datapos + data_size) > efile->ebfh.datapos_end)
		ebfutils_error("Ebf Error from Ebf::trying to read past end of record");

	if(fread(data,efile->ebfh.datasize*data_size, 1,efile->fp) != 1)
	{
		ebfutils_error("Ebf Error from Ebf::rVoid- fread error");
	}
	if(efile->ebfh.flag_swap == 1)
		ebfutils_SwapEndian(data,efile->ebfh.datasize,data_size);

	/*
	efile->patternData[0] = (efile->ebfh.datasize);
	efile->patternData[1] = data_size / (efile->ebfh.datasize);
	efile->patternData[2] = 0;
	ebfutils_my_fread1(data, &(efile->patternData[0]), 1, efile->ebfh.flag_swap,
			efile->fp);
			*/
}

int Ebf_ReadChar(const char *filename1, const char *dataname1, char* data,
		int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rChar(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadInt32(const char *filename1, const char *dataname1, int32_t* data,
		int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rInt32(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadInt64(const char *filename1, const char *dataname1, int64_t* data,
		int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rInt64(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadFloat32(const char *filename1, const char *dataname1, efloat32* data,
		int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rFloat32(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadFloat64(const char *filename1, const char *dataname1, efloat64* data,
		int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rFloat64(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadInt8(const char *filename1, const char *dataname1,
		int8_t* data, int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rInt8(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadUInt8(const char *filename1, const char *dataname1,
		uint8_t* data, int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rUInt8(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadInt16(const char *filename1, const char *dataname1,
		int16_t * data, int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rInt16(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadUInt16(const char *filename1, const char *dataname1,
		uint16_t * data, int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rUInt16(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadUInt32(const char *filename1, const char *dataname1,
		uint32_t* data, int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rUInt32(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_ReadUInt64(const char *filename1, const char *dataname1, uint64_t* data,
		int64_t ntot)
{
	EbfFile efile1 = EbfFile_Create();
	EbfFile_OpenR(&efile1, filename1, dataname1);
	if(efile1.ecode==0)
	{
		if (ntot > EbfFile_Elements(&efile1))
		{
			EbfFile_Close(&efile1);
			ebfutils_error(
					"Ebf Error from ebfread()::supplied size is larger than size of data on file");
			return 1;
		}
		EbfFile_rUInt64(&efile1, data, ntot);
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}


int Ebf_ReadCString(const char *filename1,const char *dataname1,char* data, size_t maxsize)
{
	size_t nsize;
	EbfFile efile1=EbfFile_Create();
	EbfFile_OpenR(&efile1,filename1,dataname1);
	data[0]=0;
	if(efile1.ecode==0)
	{
		nsize=EbfFile_Elements(&efile1);
		if ((nsize+1)>maxsize)
		{
			nsize=maxsize-1;
			printf("%s \n","Ebf Warning from Ebf_readCString()::supplied size is smaller than size of data on file");
		}
		EbfFile_rChar(&efile1,data,nsize);
		data[nsize]=0;
	}
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteCString(const char *filename1,const char *dataname1,char* data1, const char* mode1)
{
	return Ebf_WriteChar(filename1,dataname1,data1, mode1,"",strlen(data1));
}


int Ebf_WriteChar(const char *filename1, const char * dataname1,
		const char* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 1, dataunit,rank,dims);
	if(efile1.ecode==0)
		EbfFile_wChar(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt32(const char *filename1, const char * dataname1,
		const int32_t* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 2, dataunit,rank,dims);
	if(efile1.ecode==0)
		EbfFile_wInt32(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt64(const char *filename1, const char * dataname1,
		const int64_t* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 3, dataunit,rank,dims);
	EbfFile_wInt64(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteFloat32(const char *filename1, const char * dataname1,
		const efloat32* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 4, dataunit,rank,dims);
	EbfFile_wFloat32(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteFloat64(const char *filename1, const char * dataname1,
		const efloat64* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 5, dataunit,rank,dims);
	EbfFile_wFloat64(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt16(const char *filename1, const char * dataname1,
		const int16_t * data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 6, dataunit,rank,dims);
	EbfFile_wInt16(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt8(const char *filename1, const char * dataname1,
		const int8_t* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 9, dataunit,rank,dims);
	EbfFile_wInt8(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt8(const char *filename1, const char * dataname1,
		const uint8_t* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 10, dataunit,rank,dims);
	EbfFile_wUInt8(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt16(const char *filename1, const char * dataname1,
		const uint16_t * data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 11, dataunit,rank,dims);
	EbfFile_wUInt16(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt32(const char *filename1, const char * dataname1,
		const uint32_t* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 12, dataunit,rank,dims);
	EbfFile_wUInt32(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt64(const char *filename1, const char * dataname1,
		const uint64_t* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, 13, dataunit,rank,dims);
	EbfFile_wUInt64(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}


int Ebf_WriteCharAs(const int32_t typecode, const char *filename1, const char * dataname1,
		const char* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wChar(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt32As(const int32_t typecode,const char *filename1, const char * dataname1,
		const int32_t* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wInt32(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt64As(const int32_t typecode, const char *filename1, const char * dataname1,
		const int64_t* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wInt64(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteFloat32As(const int32_t typecode,const char *filename1, const char * dataname1,
		const efloat32* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wFloat32(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteFloat64As(const int32_t typecode,const char *filename1, const char * dataname1,
		const efloat64* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wFloat64(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt16As(const int32_t typecode,const char *filename1, const char * dataname1,
		const int16_t * data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wInt16(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteInt8As(const int32_t typecode,const char *filename1, const char * dataname1,
		const int8_t* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wInt8(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt8As(const int32_t typecode,const char *filename1, const char * dataname1,
		const uint8_t* data, const char* mode1,  const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wUInt8(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt16As(const int32_t typecode,const char *filename1, const char * dataname1,
		const uint16_t * data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wUInt16(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt32As(const int32_t typecode,const char *filename1, const char * dataname1,
		const uint32_t* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wUInt32(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}

int Ebf_WriteUInt64As(const int32_t typecode,const char *filename1, const char * dataname1,
		const uint64_t* data, const char* mode1, const char* dataunit, int64_t ntot)
{
	int64_t dims[8];
	int rank;
	EbfFile efile1 = EbfFile_Create();
	rank=1;
	dims[0]=ntot;
	EbfFile_Open(&efile1, filename1, dataname1, mode1, typecode, dataunit,rank,dims);
	EbfFile_wUInt64(&efile1, data, ntot);
	EbfFile_Close(&efile1);
	return efile1.ecode;
}




static void Efile_init(EbfFile* efile)
{
	efile->debug = 0;
	efile->fp = NULL;
	efile->patternData_size = 3;
	efile->magic = 314159;
	efile->mode = 0;
	efile->status = 0;
}

EbfFile EbfFile_Create()
{
	EbfFile efile;
	Efile_init(&efile);
	return efile;
}

/**
 * enquire header properties
 */
int64_t EbfFile_Elements(EbfFile *efile)
{
	return EbfHeader_elements(&(efile->ebfh));
}

int64_t EbfFile_Index3(EbfFile *efile, int64_t i1, int64_t i2, int64_t i3)
{
	if (efile->ebfh.rank == 3)
		return (i1 * (efile->ebfh.dim[1]) + i2) * (efile->ebfh.dim[2]) + i3;
	else
	{
		ebfutils_error("Ebf Error from index()::out of rank dim access");
	}
	return -1;
}
int64_t EbfFile_Index2(EbfFile *efile, int64_t i1, int64_t i2)
{
	if (efile->ebfh.rank == 2)
		return (i1 * (efile->ebfh.dim[1]) + i2);
	else
	{
		ebfutils_error("Ebf Error from index()::out of rank dim access");
	}
	return -1;
}

/**
 * set position
 */
void EbfFile_Seek(EbfFile *efile, int64_t offset)
{
	int64_t loc;
	if (efile->mode == 1)
	{
		if (offset >= 0)
		{
			loc=efile->ebfh.datapos+ offset * (efile->ebfh.datasize);
			if(loc>efile->ebfh.datapos_end)
			{
				ebfutils_error(
					"Ebf Error from Ebf::read()- trying to read past end of record");
			}

			fseek(efile->fp,loc,SEEK_SET);
		}
	}
	else
	{
		ebfutils_error(
				"Ebf Error from Ebf::position()- rewind allowed only in read mode");
	}
}




/**
 * get position
 */
int64_t EbfFile_Tell(EbfFile *efile)
{
	return (ftell(efile->fp) - efile->ebfh.datapos) / (efile->ebfh.datasize);
}

void EbfFile_OpenW(EbfFile *efile, const char* filename1, const char* dataname1,
		const char* mode1, int datatype, const char* dataunit)
{
	int64_t dims[1];
	dims[0] = 0;
	EbfFile_Open(efile, filename1, dataname1, mode1, datatype, dataunit, 1, &dims[0]);
}

void EbfFile_OpenR(EbfFile *efile, const char* filename1, const char* dataname1)
{
	EbfFile_Open(efile, filename1, dataname1, "r", 0, "", 0, NULL);
}

void EbfFile_Open(EbfFile *efile, const char* filename1, const char* dataname1,
		const char* mode1, int datatype, const char* dataunit,
		int rank, int64_t *dim1)
{
	int64_t location=-1;
	char mode2[10];
	char sdef[10];
	sdef[0]=0;
	efile->ecode=0;

	if (efile->magic != 314159)
		ebfutils_error("Efile not initliazed use EbfFile_create()");
	if (efile->fp != NULL)
	{
		EbfFile_Close(efile);
		ebfutils_error("File already open use EbfFile_close()");
	}

	if (efile->debug)
		printf("%s %s %s %s %s %s \n", "EbfFile::open() ", filename1, ":",
				dataname1, " mode=", mode1);


	ebfutils_assert(EBF_FILE_NAME_MAXSIZE > strlen(filename1),
			"filename_size>strlen(filename1)");
	strcpy(efile->filename, filename1);
	ebfutils_assert(strlen(mode1) < 10, "strlen(mode1)<10");
	strcpy(mode2, mode1);

	if (strcmp(mode1, "r") == 0)
	{
		efile->mode = 1;
		strcpy(mode2, "rb");
		efile->givenmode = 1;
	}
	else if (strcmp(mode1, "w") == 0)
	{
		efile->ecode=EbfTable_Init(efile->filename);
		efile->mode = 3;
		strcpy(mode2, "rb+");
		efile->givenmode = 2;
	}
	else if (strcmp(mode1, "a") == 0)
	{
		efile->mode = 3;
		strcpy(mode2, "rb+");
		efile->givenmode = 3;
	}
	else
	{
		printf("%s %s \n", "mode=", mode1);
		ebfutils_error("EBF Error from EbfFile::open(): Mode must be r w or a");
	}


	if ((efile->mode == 1) || (efile->mode == 3))
		location = EbfTable_Get(efile->filename, dataname1,&efile->ecode);

	if(efile->ecode==0)
		efile->fp = fopen(efile->filename, mode2);
	else
	{
		ebfutils_error("EBF Error from EBF::open(): probably not an ebf file");
	}


	if((efile->fp != NULL)&&(efile->ecode==0))
	{
		if (efile->mode == 3)
		{
			fseek(efile->fp, 0, SEEK_END);

			EbfHeader_set(&(efile->ebfh), dataname1, datatype,
					ebfutils_TypeSize(datatype), dataunit, sdef,rank, dim1);

			if(efile->ebfh.ecode!=0)
			{
				efile->ecode=efile->ebfh.ecode;
				printf("%s \n","EBF Error from ebfc::EbfHeader_set()");
			}
			if (location >= 0)
			{
				efile->ecode=1;
				printf("%s %s \n","EBF Error from EBF::open(): overwrite prevented as data object already exists OR not an ebf file ",dataname1);
			}

			if(efile->ecode==0)
				EbfHeader_write(&(efile->ebfh),efile->fp);

			if(efile->ebfh.ecode!=0)
			{
				efile->ecode=efile->ebfh.ecode;
				printf("%s \n","EBF Error from ebfc::EbfHeader_set()");
			}

		}
		else if (efile->mode == 1)
		{

			if (location >= 0)
			{
				fseek(efile->fp, location, SEEK_SET);
				EbfHeader_read(&(efile->ebfh), efile->fp);
			}
			else
			{
				efile->ecode=1;
				printf("%s %s %s %s\n","EBF Error from EBF::open(): Key not found:",efile->filename,":",dataname1);
			}

		}

	}
	else
	{
		efile->ecode=2;
		printf("%s %s \n","EBF Error from EBF::open(): can't open file:", efile->filename);
	}

	if(efile->ecode!=0)
	{
		if(efile->fp!=NULL)
		{
			fclose(efile->fp);
			efile->fp=NULL;
		}
		ebfutils_error("Ebf Error from Efile.open()");
	}

	if (efile->debug)
		printf("%s \n", "Efile::open() done");

}

void EbfFile_SaveTo(EbfFile *efile, const char* filename1, const char* mode1)
{
	int ecode=0;
	if (efile->givenmode == 2)
	{
		EbfFile_Close(efile);
		ecode=Ebf_Copy(efile->filename, filename1, mode1, efile->ebfh.name);
		if((ebfutils_fileExists(efile->filename) == 1)&&(ecode==0))
		{
			remove(efile->filename);
		}
		if(ecode!=0)
		{
			efile->ecode=ecode;
			ebfutils_error("EBF Error from EBF::saveTo(): while doing Ebf_Copy() ");
		}
	}
	else
	{
		EbfFile_Close(efile);
		efile->ecode=ecode;
		ebfutils_error(
				"EBF Error from EBF::saveTo():file must be opened in w mode for this operation");
	}
}

void EbfFile_Close(EbfFile *efile)
{
	int i = 0;
	int64_t temp;
	int64_t mypos;
	int ecode=0;
	int ecode_put=0;

	if (efile->fp != NULL)
	{
		if ((efile->mode == 2) || (efile->mode == 3))
		{
			if (ftell(efile->fp) > efile->ebfh.headerpos)
			{
				temp = 1;
				for (i = 1; i < efile->ebfh.rank; ++i)
					temp = temp * efile->ebfh.dim[i];

				   mypos=ftell(efile->fp);
			       if (temp*(efile->ebfh.datasize)*(efile->ebfh).dim[0] != (mypos-(efile->ebfh.datapos)))
			       {
			    	   if ((mypos-(efile->ebfh.datapos)) == 0)
			    	   {
			    		   efile->ebfh.rank=1;
			    		   efile->ebfh.dim[0]=0;
			    	   }
			    	   else if(temp==0)
			    	   {
			    		   efile->ebfh.rank=1;
						   efile->ebfh.dim[0]=(mypos-efile->ebfh.datapos)/efile->ebfh.datasize;
							ecode=1;
			    	   }
			    	   else
			    	   {
			    		   if (((mypos-(efile->ebfh.datapos))%(temp*efile->ebfh.datasize)) == 0)
			    		   {
			    			   efile->ebfh.dim[0]=(mypos-(efile->ebfh.datapos))/(temp*efile->ebfh.datasize);
			    		   }
			    		   else
			    		   {
			    			   int bufsize=(temp*efile->ebfh.datasize)-(mypos-(efile->ebfh.datapos))%(temp*efile->ebfh.datasize);
			    			   for(i=0;i<bufsize;++i)
			    				   fwrite(&temp,1,1,efile->fp);
								mypos=ftell(efile->fp);
								efile->ebfh.rank=1;
								efile->ebfh.dim[0]=(mypos-efile->ebfh.datapos)/efile->ebfh.datasize;
								ecode=1;
			    		   }
			    	   }
						efile->ebfh.capacity_=EbfHeader_elements(&(efile->ebfh)) * efile->ebfh.datasize;
						fseek(efile->fp, efile->ebfh.headerpos, SEEK_SET);
						EbfHeader_write(&(efile->ebfh),efile->fp);
			       }

					if ((mypos - efile->ebfh.datapos)!= efile->ebfh.capacity_)
						ebfutils_error("Ebf error:efile_close() capacity_< total data size");



			}
		}

		if (fclose(efile->fp) == 0)
			efile->fp = NULL;
		else
		{
			ebfutils_error(
					"Ebf Error from Efile.close()-error closing file fclose()");
		}

		if (efile->mode == 3)
		{
			EbfTable_Put(efile->filename, efile->ebfh.name,
					ecode_put=efile->ebfh.headerpos);
		}
		/* changed to after as correction done for corrupt file */
		if (ecode != 0)
		{
			printf("%s \n","Ebf error: Incorrect dimensions given, changing to 1 d and adjusting size");
			ebfutils_error("Ebf Error from Efile.close()-error closing file fclose()");
		}

	}
	else
	{
		ebfutils_error("trying to close an already closed efile");
	}
	efile->mode = 0;
	efile->givenmode = 0;
	if(ecode!=0)
		efile->ecode=ecode;
	if(ecode_put!=0)
		efile->ecode=ecode;
}






