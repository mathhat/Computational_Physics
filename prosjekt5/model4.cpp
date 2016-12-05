using namespace std;
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <string>
#include <sstream>
using namespace arma;
fstream outfile;


void transactions(mat& C, vec& M,double exp,int Number_of_agents,double lambda,double alpha,double gamma);
void deal(vec& M,int N,double m0);
void write(vec& M,int N,string filename);

int main(int argc, char *argv[]){
	cout << argc << endl;
		if (argc < 5){
		cout << "You're gonna have to give me some command line arguments man" << endl;
		cout << "After the file name, write a number representing the 10 exponent for transactions."<< endl;
		cout << "Second CMD argument is gonna be lambda, the saving factor. Any number between 0 and 1 will do." << endl;
		cout << "Now, write another number between 0 and infinity, representing Alpha, This guy helps people trade with their equals." << endl;
		cout << "Lastly, write a number for gamma. Make it 0, 1, 2, 3 or 4" << endl;
		return 0;
		}
    stringstream s1;
    stringstream s2;
    stringstream s3;
	stringstream s4;
	int exponent;
	exponent = atoi(argv[1]); //The exponent gathered from the CMD decides
							  //- how many transactions we simulate (10^exponent)
	double lambda;
	lambda = atof(argv[2]); //The saving coefficient decides how rapid the cash disperses

	double alpha;
	alpha = atof(argv[3]); //alpha
	cout <<alpha <<endl;

	double gamma;
	gamma = atof(argv[4]); //alpha
	cout <<gamma <<endl;
	srand(2);
	int Number_of_agents = 500;  //Number of cash traders (agents)
	double m0 = 1; 				 //Start up cash for each agent

	s1 << exponent;
	s2 << lambda;
	s3 << alpha;
	s4 << gamma;
	string filename = "exp"+s1.str()+"lambda"+s2.str()+"alpha"+s3.str()+"gamma"+s4.str();
	outfile.open(filename,ios::out);

	for(int cycle = 0; cycle <1000; cycle ++){
		mat C = mat(Number_of_agents,Number_of_agents,fill::zeros);
		vec agents;  					  //Agents' money in a vector
		deal(agents,Number_of_agents,m0); //function gives each agent his starup moeny
		transactions(C,agents,exponent,Number_of_agents,lambda,alpha,gamma);	  //function where all transactions take place
		write(agents, Number_of_agents,filename);								  //write data to file
	}
} //end main






void deal(vec& M,int N,double m0){
	M = vec(N);
	for (int i = 0; i < N; i++)
		M(i) = m0;
}



void transactions(mat& C,vec& M,double exp,int Number_of_agents, double lambda,double alpha,double gamma){
	double eps;
	double md;
	int counter=0;
	int I;
	int J;
	int EXP =pow(10,exp);
	double randmax_inv = 1./RAND_MAX;
	while (counter < EXP){		//breaks after 10^exp transactions
		I = rand()%Number_of_agents;
		J = rand()%Number_of_agents;
		if (I-J){
			C(I,J) += 1;
			C(J,I) += 1;
			double r = rand()*randmax_inv*10;
			double p = pow(M(I)-M(J),-alpha) * pow(C(I,J)+1,gamma);
			if ( p*p > r*r){
				counter += 1;		//counts number of transactions
				eps = rand()*randmax_inv;
				md = (1.-lambda)*(eps*M(J)-(1.-eps)*M(I));
				M(I) = M(I)+md;
				M(J) = M(J)-md;
			}
		}
	}//while loops ends
}




void write(vec& M,int N,string filename){
	for (int i = 0; i < N; i++){
		outfile << setprecision(10) << M(i)<<endl;
	}
}


