using namespace std;
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <string>
using namespace arma;
fstream outfile;


void transactions(vec& M,double exp,int Number_of_agents,double lambda);
void deal(vec& M,int N,double m0);
void write(vec& M,int N,string filename);
void variance(int Number_of_agents,double exponent,double m0,double lambda);

int main(int argc, char *argv[]){
	cout << argc << endl;
		if (argc < 3){
		cout << "You're gonna have to give me some command line arguments man" << endl;
		cout << "After the file name, write a number representing the 10 exponent for transactions."<< endl;
		cout << "Second CMD argument is gonna be lambda, the saving factor. Any number between 0 and 1 will do." << endl;
		return 0;
	}
	int exponent;
	exponent = atoi(argv[1]); //The exponent gathered from the CMD decides
							  //- how many transactions we simulate (10^exponent)
	double lambda;
	lambda = atof(argv[2]); //The saving coefficient decides how rapid the cash disperses
	srand(2);
	int Number_of_agents = 500;  //Number of cash traders (agents)
	double m0 = 10; 				 //Start up cash for each agent	

    stringstream s1;
    stringstream s2;

	int exponent;
	exponent = atoi(argv[1]); //The exponent gathered from the CMD decides
							  //- how many transactions we simulate (10^exponent)
	double lambda;
	lambda = atof(argv[2]); //The saving coefficient decides how rapid the cash disperses
	srand(2);
	int Number_of_agents = 500;  //Number of cash traders (agents)
	double m0 = 10; 				 //Start up cash for each agent	
	s1 << exponent;
	s2 << lambda;

	string filename = "exp"+s1.str()+"lambda"+s2.str();
	outfile.open(filename,ios::out);
	for(int cycle = 0; cycle <1000; cycle ++){
		vec agents;  //Agents' money in a vector
		deal(agents,Number_of_agents,m0); //function gives each agent his starup moeny
	
		transactions(agents,exponent,Number_of_agents,lambda);	  //function where all transactions take place
		write(agents, Number_of_agents,filename);				//write data to file

	}	

	
} //end main






void deal(vec& M,int N,double m0){
	M = vec(N);
	for (int i = 0; i < N; i++)
		M(i) = m0;
}



void transactions(vec& M,double exp,int Number_of_agents, double lambda){

	double eps;
	double md;
	int counter=0;
	int I;
	int J;

	while (counter <pow(10,exp)){		//breaks after 10^exp transactions
		I = rand()%Number_of_agents;
		J = rand()%Number_of_agents;
		if (I-J){
			counter +=1;				//counts number of transactions
			eps = rand()/(1.*RAND_MAX);
			md = (1.-lambda)*(eps*M(J)-(1.-eps)*M(I));
			M(I) = M(I)+md;
			M(J) = M(J)-md;
		}
		
	}
}




void write(vec& M,int N,string filename){
	for (int i = 0; i < N; i++){
		outfile << setprecision(10) << M(i)<<endl;
	}
}




void variance(int Number_of_agents,double exponent,double m0,double lambda){
	fstream ofile;
	ofile.open("variance",ios::out);
	double variance;
	//variance
	for (double Exp = 0; Exp <=exponent;Exp+=0.1){
		variance=0;
		vec agents;  //Agents' money in a vector
		deal(agents,Number_of_agents,m0); //function gives each agent his starup moeny
		transactions(agents,Exp,Number_of_agents, lambda);	  //function where all transactions take place
		//write(agents, Number_of_agents,filename);
		for (int i = 0; i<Number_of_agents;i++){
			variance += (agents(i)-m0)*(agents(i)-m0)/Number_of_agents;
		}
		ofile << setprecision(10) << variance<<endl;
	}
}

