//This is the script that calculates the unoptimized gaussian elimination method for a general tridiagonal matrix
#include<iostream>
#include<cmath>
#include<time.h>
#include<fstream>
using namespace std;

    void V(double *v,double *x,long n); //exact solution values
    void f(double *F,double *x,long n); //adds function f values to array F
    void linespace(double *x,long n);  //adds x values to the x array for the function F

    int main(){
    double start, finish;

    long n = 10;


    double *a = new double[n+1];   //Constants in the upper tridiagonal matrix A
    double *b = new double[n+1];   //diagonal constants, later becomes b_tilde
    double *c = new double[n+1];   //lower diagonal

    double *x = new double[n+1];       //Discretized x values for the function
    double *F = new double[n+1];       //exact function values and f_tilde
    double *v = new double[n+1];       //numerical solution values
    double *v_exact = new double[n+1]; //exact solution values

    ////// The following functions fill up the x, F and v_exact arrays

    linespace(x,n);         //x values for x array, values from 0 to 1.

    f(F,x,n);               // assigning function values multiplied with the x_step value squared (h*h)

    V(v_exact, x, n);       //exact solution values for v added to array v_exact


    for(int i = 0; i < n+1; i++){       //assigning non diagonal and diagonal constants to arrays
     
        a[i] = -1.;
        b[i] = 2.;
        c[i] = -1.;
        }


    F[1] = a[1]/b[1] * F[1];
    b[1] = b[1];            // there are the initial values of F_1 tilde and b_1 tilde

    //we now have what we need to implement the forward substitution
    //Here begins the forward substitution
    double counter = 0;
    start = clock();       

    start = clock();       
    //forward substitution
    for(int i = 2; i < n-1; i++){              //this loops gives F and b tilde values
        b[i] = b[i] - a[i]/b[i-1] * c[i-1];
        F[i] =  a[i+1]/b[i] * (F[i] - F[i-1]); //6*n flops
    }


    b[n-1] = b[n-1] - a[n-1]/b[n-2] * c[n-2]; //b[n] tilde, doesn't stray from general formula but had -
                                    // to ble added seperately from the loop due to index difficulties.

    F[n-1] = (F[n-1] - F[n-2]);       //F[n]_tilde strays from general formula for F[i] tilde
    v[n-1] = F[n-1]/b[n-1];         // the expression for v_(n) differs from the the general formula
                                    //- for other v_(i) in the back.subst. section 


   /////Backward substitution
    for(int i = n-2; i > 0; i--){   // "i--" will reduce i by 1 for each loop iteration
        v[i] =  (F[i]  -  a[i+1]*c[i]*v[i+1]/b[i]) /a[i+1]; //5nflops
    }
    finish = clock();
    counter += (finish -start)/CLOCKS_PER_SEC;
    cout << (finish -start)/CLOCKS_PER_SEC <<" "<<counter <<endl;

    ///Writing the result to file
    fstream outfile; //file we're writing to
    outfile.open("v.dat",ios::out);
    for (int i = 0; i<n+1; i++){
        outfile << v[i] << " " << v_exact[i] << endl;
    }
    outfile.close();
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] x;
    delete [] F;
    delete [] v;
    delete [] v_exact;
} //end of main


void linespace(double *x,long n){ //discreet x values function
    for(long i = 0; i<n+1; i++){
        x[i] = double(i)/n;
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
