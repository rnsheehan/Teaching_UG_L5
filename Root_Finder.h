#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

// Declare the root finder class
// R. Sheehan 8 - 3 - 2013

// normally I wouldn't reccommend the use of typedef, but in this situation it is the simplest solution
// we will use dfunc as a type for a pointer to a function that returns a double
// we can then pass functions to the root finding class

typedef double (*dfunc)(double x); // custom type for a pointer to a function of type double

// Data type for an interval [xlower, xupper]
// We will use this object to store the endpoints of subintervals known to contain a root of a function
class interval{
public:
	// Constructor
	interval(); 
	interval(double xl, double xu); 

	// Methods

	void set_xl_xu(double xl, double xu); 

	double get_x_lower();
	double get_x_upper();

private:
	double xlower; 
	double xupper; 
};

// root finding object
// with this class we will be able to locate bounds on the roots of a function over an interval [a, b]
// we can refine the location of those roots using either bisection mnethod, or newton-raphson method
// the implementation of the newton-raphson method does not require knowledge of the derivative of the function
// an accurate estimate of the derivative is provided by the derivative algorithm, it is computed using a 
// process known as Richardson extrapolation

class find_root{
public:
	// Constructor
	find_root(); 
	find_root(double (*func)(double), double a, double b, int ntest); 

	// Methods
	void set_interval(double a, double b, int ntest); // change the interval properties

	void bracket_roots(); // function for bracketing the roots
	void bisection_search(double toler); // bisection method root search
	void newton_raphson_search(double toler); // Newton-Raphson method root search

	double derivative(double x, double dx); // compute the value of the derivative of f at the point x by Richardson extrapolation

private:
	int nsub; // number of subintervals to search inside [a, b]
	int nroots; // the number of roots found on the interval [a, b]

	interval search_space; // this will define the interval [a, b]
	vector<interval> sub_intervals; // an array of sub-intervals known to contain a root

	dfunc the_func; // the function whose roots are being sought
	
	vector<double> the_roots; // vector to store the roots
};

#endif
