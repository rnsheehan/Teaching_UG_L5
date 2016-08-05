#ifndef ATTACH_H
#define ATTACH_H

// This is the Attach library
// It contains a list of commonly used standard C++ libraries
// For the moment this is the minimum of libraries we require
// R. Sheehan 23 - 1 - 2013

// For more information go to http://www.cplusplus.com/

/*
cstdlib:
This header defines several general purpose functions, including dynamic memory management, random number generation, 
communication with the environment, integer arthmetic, searching, sorting and converting.

This library uses what are called streams to operate with physical devices such as keyboards, printers, terminals or with any other type of
files supported by the system. Streams are an abstraction to interact with these in an uniform way; 
All streams have similar properties independently of the individual characteristics of the physical media they are associated with.

iostream:
Header that defines the standard input/output stream objects

iomanip:
Header providing parametric manipulators, i.e. use it to change the parameters of the output stream

cmath:
declares a set of functions to compute common mathematical operations and transformations

using namespace std:
this tells MSVS that you want to use functions included in the standard library in your code
*/

#include <cstdlib> // this is equivalent to stdlib.h
#include <cstdio> // this is equivalent to stdio.h
#include <iostream> // this is C++ specific
#include <iomanip> // this is C++ specific

#include <cmath> // this is equivalent to math.h

#include <fstream> // this library enables us to read in data from files and to write data to files
#include <sstream> // this library is needed for file IO, specifically string stream objects
#include <string> // this library enabels us to manipulate strings

#include <vector> // array object from the standard template library

using namespace std; // enables you to output material to the console

// Constants used in the integral calculations
static const double EPS = (3.0e-12);

// definition of pi
static const double p = (atan(1.0)); // pi / 4
static const double PI = (4.0*p); // pi
static const double PI_2 = (2.0*p); // pi / 2

// Template function declarations

template <class T> T Signum(T a)
{
	//The sign operator
	T darg;
	
	return ( (darg=(a)) >= (T)(0) ? (T)(1) : -(T)(1) ); // Setting the Sign of zero to be 1
}

template <class T> T DSQR(T a)
{
	//Efficient squaring operator
	T darg;
	return ( (darg=(a)) == (T)(0) ? (T)(0) : darg*darg );
}

#include "Root_Finder.h"

#endif