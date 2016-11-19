//This is the script that calculates the optimized gaussian elimination method for a specified tridiagonal matrix, also error calculations take place here
#include<iostream>
#include<cmath>
#include<time.h>
#include<fstream>
using namespace std;

void Error(double *v,double *v_exact,double *error,long n); //relative error

void V(double *v,double *x,long n); //exact solution values
void f(double *F,double *x,long n); //adds function f values to array F
void linespace(double *x,long n);  //adds x values to the x array for the function F

int main(){
    double start, finish;
    long n = 1E1;

    double *error = new double[n-2];   //array for errorvalues 


    double *x = new double[n+1];       //Discretized x values for the function
    double *F = new double[n+1];       //exact function values and f_tilde
    double *v = new double[n+1];       //numerical solution values
    double *v_exact = new double[n+1]; //exact solution values

    ////// The following functions fill up the x, F and v_exact arrays

    linespace(x,n);         //x values for x array, values from 0 to 1.

    f(F,x,n);               // assigning function values multiplied with the x_step value squared (h*h)

    V(v_exact, x, n);       //exact solution values for v added to array v_exact



    F[1] = -F[1]/2.;
    //we now have what we need to implement the forward substitution
    //Here begins the forward substitution

    start = clock();                      //clock value before gaussian elimination
    //forward substitution
    for(int i = 2; i < n-1; i++){                  //this loops gives F tilde values, a product of the forward substitution
        F[i] =  -(F[i] - F[i-1])/(double(i+1)/(i)); //4*n flops 
    }
    
    F[n-1] = (F[n-1] - F[n-2]);       //F[n]_tilde strays from general formula for F[i] tilde
    v[n-1] = F[n-1]/(double(n)/(n-1));         // the expression for v_(n) differs from the the general formula
                                    //- for other v_(i) in the back.subst. section 

    //Backward substitution
    for(int i = n-2; i > 0; i--){   // "i--" will reduce i by 1 for each loop iteration
        v[i] =  -(F[i]  -  v[i+1]/(double(i+1)/(i))); //2nflops
    }
    finish = clock();                    //clock value after gaussian elimination
    cout << (finish -start)/CLOCKS_PER_SEC<<endl;

    //error calculations
    Error(v, v_exact, error, n);

    ///Writing the result to file
    fstream outfile; //file we're writing to
    outfile.open("v.dat",ios::out);
    for (int i = 0; i<n+1; i++){
        outfile << v[i] << " " << v_exact[i] << endl;
    }
    outfile.close();

    fstream Outfile; //file we're writing the ERROR to
    Outfile.open("error.dat",ios::out);
    for (int i = 0; i<n-1; i++){
        Outfile << error[i] << endl;
    }
    Outfile.close();
    //freeing memory
    delete [] error;
    delete [] x;
    delete [] F;
    delete [] v;
    delete [] v_exact;
}   //end of main


void linespace(double *x,long n){ //discreet x values function
    for(long i = 0; i<n+1; i++){
        x[i] = double(i)/(n);
    }
}
void f(double *F,double *x,long n){ //function f
    double h = 1./(n);    
    for(long i = 0; i<n+1; i++){
        F[i] = 100.*exp(-10.*x[i])*h*h; //filling up the f array with F
    }

}

void V(double *v,double *x,long n){ //exact solution values
    double h = x[1];    
    for(long i = 0; i<n+1; i++){
        v[i] = 1. - (1.-exp(-10.))*x[i] -exp(-10.*x[i]); //filling up the v_exact array with F
    }
}
void Error(double *v,double *v_exact,double *error,long n){ //error values
    for(long i = 1; i<n; i++){
        error[i-1] = log10(sqrt(pow((v[i]/v_exact[i]-1),2))); //adding relative error to array
    }
}
