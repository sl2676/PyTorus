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

#ifndef TIMER_H_
#define TIMER_H_
#include <ctime>
#include <sys/time.h>

//#include "sutils.h"

//#include<ctime>
//#include <iostream>

class Timer
{
private:
	clock_t startTime,stopTime;
//	time_t startTime1,stopTime1;
	double startTime1,stopTime1;
	timeval tim;
	double getTime();
public:
	Timer();
	~Timer();
	double elapsedTime();
	double currentTime();
	void start();
	void stop();
	void resume();
	void print();
	void print(const char *s,double n=0);
	void printC();
	void printC(const char *s,double n=0);
};

//ostream& operator<< (ostream& os, const Timer& timer1);

#endif /*TIMER_H_*/
