#ifndef ATTACH_H
#include "Attach.h"
#endif

// interval object
// Constructors
interval::interval()
{
	// Default constructor
	xlower = xupper = 0.0; 
}

interval::interval(double xl, double xu)
{
	// Constructor
	set_xl_xu(xl,xu); 
}

//Methods
void interval::set_xl_xu(double xl, double xu)
{
	// set the values of xlower and xupper
	// ensure that xl < xu
	xlower = min(xl, xu);

	xupper = max(xl, xu); 
}

double interval::get_x_lower()
{
	return xlower;
}

double interval::get_x_upper()
{
	return xupper;
}

// root search object
// Constructors
find_root::find_root()
{
	// Default constructor
	nsub = nroots = 0;
}

find_root::find_root(double (*func)(double),double a, double b, int ntest)
{
	// Constructor for single variable functions
	// search over ntest subintervals of [a, b] for the roots of func(x)

	nsub = ntest; 

	nroots = 0; 

	search_space.set_xl_xu(a, b); 

	the_func = func; 
}

// Methods

//void find_root::test_func()
//{
//	// method to ensure that the function under consideration has been passed correctly
//	// not required
//
//	cout<<"Function value  = "<<the_func(7.0)<<endl;
//}

void find_root::set_interval(double a, double b, int ntest)
{
	// change the interval properties

	nsub = ntest; 

	nroots = 0; 

	search_space.set_xl_xu(a, b); 

}

void find_root::bracket_roots()
{
	// Find the sub-intervals of [a, b] known to contain a root of the equation f(x) = 0

	double fl, fu; 
	double xl, xu; 
	double dx = ((search_space.get_x_upper() - search_space.get_x_lower()) / (nsub - 1)); 

	xl = search_space.get_x_lower(); 
	xu = xl + dx; 
	for(int i=1; i<=nsub; i++){
		
		// Evaluate the function on the sub-interval
		if(i==1){
			fl = Signum(the_func(xl)); 
		}
		else{
			fl = fu; // minimise the number of function evaluations
		}

		fu = Signum(the_func(xu));

		// Perform the bisection test
		if(fl*fu < 0.0){
			// the sub-interval contains a root so it is stored
			nroots++; 
			sub_intervals.push_back(interval(xl,xu)); 
		}

		// update the endpoints of the sub-interval
		xl = xu; 
		xu += dx; 
	}

	cout<<"The function contains "<<nroots<<" roots on the interval [ "<<search_space.get_x_lower()<<" , "<<search_space.get_x_upper()<<" ]\n"; 
	if(nroots > 0){
		cout<<"The roots are located in \n";
		for(int i=0; i<nroots; i++){
			cout<<"[ "<<sub_intervals[i].get_x_lower()<<" , "<<sub_intervals[i].get_x_upper()<<" ]\n";
		}
		cout<<"\n";
	}
}

void find_root::bisection_search(double toler)
{
	// Locate the roots of the function f over sub-intervals known to contain roots
	// to within a value of toler

	if(nroots > 0){
	
		int iter = 0, max_iter = 500; 
		bool cgt = false; 
		double r=0.0, xl, xu, dx, fl, fr; 

		for(int i = 0; i<nroots; i++){
			// Loop over each sub-interval and locate the root there in

			iter = 0; 
			cgt = false; 

			xl = sub_intervals[i].get_x_lower(); 
			xu = sub_intervals[i].get_x_upper(); 

			fl = Signum(the_func(xl)); 

			while(iter < max_iter){
			
				dx = 0.5*(xu-xl); 

				r = xl + dx; 
				
				if(fabs(dx) < toler){
					cout<<"Bisection algorithm has converged to within tolerance in "<<iter<<" iterations\n"; 
					the_roots.push_back(r); 
					cgt = true; 
					break; 
				}
			
				// Apply the bisection test
				// update the endpoints to the interval known to
				// contain the root
				fr = Signum(the_func(r)); 

				if(fl*fr > 0.0){
					xl = r; 
					fl = fr;
				}
				else{
					xu = r; 
				}

				iter++; 
			}

			if(cgt){
				cout<<"The root of f is located at "<<r<<"\n"; 
				cout<<"f ( "<<r<<" ) = "<<the_func(r)<<"\n\n";
			}
			else{
				cout<<"Bisection method failed to converge to a root\n\n"; 
			}

		}

	}
}

void find_root::newton_raphson_search(double toler)
{
	// Locate the roots of the function f over sub-intervals known to contain roots
	// to within a value of toler
	// This assumes that the bracketing algorithm has already been applied

	if(nroots > 0){
	
		int iter = 0, max_iter = 500; 
		bool cgt = false; 
		double rnew = 0.0, rold = 0.0, xl, xu, dx, fl, fr; 

		for(int i = 0; i<nroots; i++){
			// Loop over each sub-interval and locate the root there in

			iter = 0; 
			cgt = false; 

			xl = sub_intervals[i].get_x_lower(); 
			xu = sub_intervals[i].get_x_upper(); 

			rold = xl; 

			while(iter < max_iter){
			
				dx = the_func(rold) / derivative(rold,0.1); 

				rnew = rold - dx; 

				if(fabs(dx) < toler){
					cout<<"Newton-Raphson algorithm has converged to within tolerance in "<<iter<<" iterations\n"; 
					the_roots.push_back(rnew); 
					cgt = true; 
					break; 
				}
			
				// Ensure that N-R alg stays inside the original interval
				if(rnew > xu || rnew < xl){
					//cout<<"Newton-Raphson algorithm has jumped out of the search interval\n"; 
					rold = xl + 0.5*(xu-xl); 

					fl = Signum(the_func(xl)); 
					
					fr = Signum(the_func(rold)); 
					
					if(fl*fr > 0.0){
						xl = rold;
					}
					else{
						xu = rold; 
					}
				}
				else{
					rold = rnew; 
				}

				iter++; 
			}

			if(cgt){
				cout<<"The root of f is located at "<<rnew<<"\n"; 
				cout<<"f ( "<<rnew<<" ) = "<<the_func(rnew)<<"\n\n";
			}
			else{
				cout<<"Newton-Raphson method failed to converge to a root\n\n"; 
			}

		}

	}
}

double find_root::derivative(double x, double dx)
{
	// Use Richardson extrapolation to estimate the derivative of a
	// function at the point x
	// this implementation is based on the dfridr algorithm given in
	// NRinC by Press et al. 
	// this implementation uses less function evaluations than a
	// recursive implementation of Richardson extrapolation

	if(dx == 0.0){
		cout<<"Derivative cannot be computed with zero step size\n";
		return 1.0; // avoid division by zero
	}
	else{
		const int ntab = 10; 
		const double con = 1.4, con2 = DSQR(con); 
		const double big = std::numeric_limits<double>::max(); 
		const double safe = 2.0; 

		int i, j; 
		double err, errt, fac, hh, ans; 

		double a[ntab][ntab]; // keep the function values in a table

		hh = dx; 

		a[0][0] = (the_func(x+hh)-the_func(x-hh)) / (2.0*hh); // first approximation to f'(x)

		err = big; 

		for(i=1; i<ntab; i++){

			hh /= con; 
			
			a[0][i] = (the_func(x+hh)-the_func(x-hh)) / (2.0*hh); // approximation to f'(x) with smaller step size
			
			fac = con2; 
			
			// extrapolate the derivative to higher orders without extra function evaluations
			for(j=1; j<=i; j++){

				a[j][i] = (a[j-1][i]*fac-a[j-1][i-1]) / (fac-1.0); 
				
				fac = con2*fac; 
				
				errt = max(fabs(a[j][i]-a[j-1][i]), fabs(a[j][i]-a[j-1][i-1])); 
				
				// compute the new error with the error from the previous step
				if(errt <= err){
					err = errt; // err is an error estimate for the computed derivative
					ans = a[j][i]; // update the derivative estimate
				}

			}

			// if error has increased significantly stop
			if(fabs(a[i][i]-a[i-1][i-1]) >= safe*err){
				break; 
			}

		}

		//cout<<"The value of the derivative at x = "<<x<<" is "<<ans<<"\n"; 

		return ans; 
	}
}
