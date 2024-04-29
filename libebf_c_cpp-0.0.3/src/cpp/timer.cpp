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

#include <iostream>
#include <iomanip>
#include "timer.hpp"


Timer::Timer()
{
	start();
}

Timer::~Timer()
{
}

double Timer:: getTime()
{
	gettimeofday(&tim, NULL);
	return tim.tv_sec+(tim.tv_usec/1000000.0);
}

void Timer::start()
{
	startTime=clock();
	stopTime=clock();

	startTime1=getTime();
	stopTime1=getTime();

}

void Timer::stop()
{
	stopTime=clock();
	stopTime1=getTime();
}


void Timer::resume()
{
	startTime=clock()-(stopTime-startTime);
	startTime1=getTime()-(stopTime1-startTime1);
	stopTime=clock();
	stopTime1=getTime();
}


double Timer::elapsedTime()
{
	return double(stopTime-startTime)/CLOCKS_PER_SEC;
}

double Timer::currentTime()
{
	return double(clock()-startTime)/CLOCKS_PER_SEC;
}


void Timer::print()
{
	using namespace std;
	stop();
	cout<<left<<setw(12)<<" Time "<<setw(12)<<elapsedTime()<<endl;
}
void Timer::print(const char *s,double n)
{
	using namespace std;
	stop();
	if(n>0)
	{
		cout<<left<<setw(36)<<s<<setw(12)<<elapsedTime()<<"     (  "<<n*1e-6/elapsedTime()<<" Mops )"<<endl;
		cout<<left<<setw(36)<<s<<setw(12)<<difftime(stopTime1,startTime1)<<"     (  "<<n*1e-6/difftime(stopTime1,startTime1)<<" Mops )"<<endl;
	}
	else
	{
		cout<<left<<setw(36)<<s<<setw(12)<<elapsedTime()<<endl;
		cout<<left<<setw(36)<<s<<setw(12)<<difftime(stopTime1,startTime1)<<endl;
	}
}

void Timer::printC()
{
	using namespace std;
	cout<<left<<setw(12)<<" Time "<<setw(12)<<currentTime()<<endl;
	cout<<left<<setw(12)<<" Time "<<setw(12)<<(getTime()-startTime1)<<endl;
}

void Timer::printC(const char *s,double n)
{
	using namespace std;
	if(n>0)
		cout<<left<<setw(36)<<s<<setw(12)<<(getTime()-startTime1)<<"     (  "<<n*1e-6/(getTime()-startTime1)<<" Mops )"<<endl;
	else
	{
		cout<<left<<setw(36)<<s<<setw(12)<<currentTime()<<endl;
		cout<<left<<setw(12)<<" Time "<<setw(12)<<(getTime()-startTime1)<<endl;
	}
}

