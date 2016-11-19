// This is the only file that uses aramdillo. It's library solves here the problem with LU decomposition.
#include<iostream>
#include<cmath>
#include<time.h>
#include<fstream>
#include<armadillo>
using namespace std;
using namespace arma;

void f(vec &F,vec &x,long n); //adds function f values to array F
void linespace(vec &x,long n);  //adds x values to the x array for the function F

int main(){
    double start, finish;

    long n = 10 -2;

    mat A = mat(n+1,n+1,fill::zeros);  //Matrix A for LU comparing
    mat P,L,U;                         //Permutation, lower and upper diagonal matrices for LU
    vec v_LU = vec(n+1);
    vec x = vec(n+1);       //Discretized x values for the function
    vec F = vec(n+1);
    ////// The following functions fill up the x, F and v_exact arrays

    linespace(x,n);         //x values for x array, values from 0 to 1.

    f(F,x,n);               // assigning function values multiplied with the x_step value squared (h*h)




    ////// Filling up the Matrix and solving problem with LU

    for(int i = 0; i<n; i++){         //The tridiagonal elements
        A(i,i) = 2.;                     //middle diagonal elements
        A(i,i+1) = -1.;                  //upper diagonal
        A(i+1,i) = -1.;                  //lower diagonal
    }
    A(n,n) = 2;
    start = clock();                      //clock value before gaussian elimination
    lu(L, U, P, A);                     //LU decomposition
    v_LU = solve(trimatu(U), solve(trimatl(L), P*F) );  //solution
    
    finish = clock();                    //clock value after gaussian elimination
    cout << (finish -start)/CLOCKS_PER_SEC<<endl;

    

    vec solution = vec(n+1+2);                //this guy we write to file since boundaries are ignored by armadillo.
    solution(0) = 0;
    solution(n+2) = 0;
    for(int i = 0; i<n+1; i++){
        solution(i+1) = v_LU(i);
    }
    ///Writing the result to file
    fstream outfile; //file we're writing to
    outfile.open("v_LU.dat",ios::out);
    for (int i = 0; i<n+3; i++){
        outfile << solution(i) << endl;
    }
    outfile.close();

} //end of main


void linespace(vec &x,long n){ //discreet x values function
    for(long i = 0; i<n+1; i++){
        x(i) = double(i+1)/(n);
    }
}
void f(vec &F,vec &x,long n){ //function f
    double h = 1./(n+1);    
    for(long i = 0; i<n+1; i++){
        F(i) = 100.*exp(-10.*x(i))*h*h; //filling up the f array with F
    }
}
