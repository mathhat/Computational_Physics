int numprocs, my_rank;
using namespace std;
#include <mpi.h>
#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>
#include <string>
using namespace arma;
ofstream ofile;

int magnet(mat& S, int L);
double ising_flip(mat& S, int L,double& E);
double energy(mat& S, double E, int L);
void Metropolis(mat& S,mat& S_old,int L,int macMCC,double Temperature,vec& LocalExpectationValues);
void WriteResultstoFile(double temperature, vec ExpectationValues);


int main(int nargs, char* args[]){

    mat S;
    mat S_old;
    int L;
    int maxMCC;
    double T;
    double T_end;
    double TempStep;
    //   MPI initializations
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    int seed=1*my_rank;
    srand(seed);
    if ((my_rank == 0)) {
        //filename= "data0.dat";
        L = 10;
        maxMCC = 70000;
        T = 1;
        T_end = 2.4;
        TempStep = 0.05;
    }
    // Declare new file name and add lattice size to file name, only master node opens file
    if (my_rank == 0) {
        string fileout = to_string(L+1);
        string argument = to_string(L+1);
        fileout.append(argument);
        ofile.open(fileout);
    }
    // broadcast to all nodes common variables since only master node reads from command line
    MPI_Bcast (&maxMCC, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&T_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&TempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Start Monte Carlo sampling by looping over the selected Temperatures
    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
    for (double Temperature = T; Temperature <= T_end; Temperature+=TempStep){
    //double Temperature = T;
    //for (maxMCC; maxMCC <= 100000; maxMCC*=1.2){
        vec LocalExpectationValues = zeros<mat>(5);
        // Start Monte Carlo computation and get local expectation values
        Metropolis(S,S_old,L,maxMCC,Temperature,LocalExpectationValues);
        // Find total average
        vec TotalExpectationValues = zeros<mat>(5);

        for( int i =0; i < 5; i++){
            MPI_Reduce(&LocalExpectationValues[i], &TotalExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        if ( my_rank == 0){
            ofstream ofile;
            WriteResultstoFile(Temperature, TotalExpectationValues);
        }
    }

    //if(my_rank == 0)  ofile.close();  // close output file

    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }
    // End MPI
    MPI_Finalize ();
    return 0;
}
void Metropolis(mat& S,mat& S_old,int L,int maxMCC,double Temperature,vec& LocalExpectationValues){
    int M_old = 0;				//Momentary magnetiztion
    int M_new = 0;				//Momentary magnetiztion
    int M_L=0;
    double M = 0; 				//Mean abs magnetization
    double M_sq = 0; 			//Mean squared magnetization

    double E = 0; 		//Total energy
    double E_mean = 0;      //Mean energy
    double E_mean_sq = 0;   //Mean squared energy

    double EnergySum = 0;		//Needed for calculating mean energy (expectation value)
    double EnergySum_L = 0;
    double EnergySum_sq = 0;	//For mean sqrd energy, needed for sigmaE and heat capacity

    double deltaE_new;			//Energy difference due to state transition
    double deltaE_old=0;			//Energy difference due to state transition

    double Cv;
    double chi;

    S = mat(L+1,L+1,fill::ones);
    for (int i =0; i <L;i++){
        for (int j = 0; j <L;j++){
            if (double(rand())/RAND_MAX < 0.5){
                S(i,j) =-1;}}}

    double E_new;
    double E_old;
    EnergySum = energy(S,E,L);
    EnergySum_sq = EnergySum*EnergySum/(L*L);

    M_old = ising_flip(S,L,E_old); // initial
    S_old=S;
    M += magnet(S,L);
    int accept = 0;
    for (int MCC = 0; MCC < maxMCC; MCC++) {
        for(int l = 0; l < L*L;l++){
            M_new = 2*ising_flip(S,L,E_new);
            deltaE_new =E_new - E_old;
            if (deltaE_new > 0){
                if (rand()*1./RAND_MAX < exp(-deltaE_new/Temperature)){
                    accept+=1;
                    EnergySum_L += deltaE_new;  // N*<E>
                    M_L += M_new;
                    S_old = S;
                }
                else{
                    M_L+=M_old;
                    EnergySum_L+=deltaE_old;
                }
            }
            else{
                EnergySum_L += deltaE_new;
                accept+=1;
                M_L += M_new;
                S_old = S;
            }
            E_old = E_new;
            M_old = M_new;
        }
        EnergySum += EnergySum_L;
        EnergySum_sq += EnergySum_L*EnergySum_L;
        EnergySum_L=0;
        if(M_L<0){M_L*=-1;}
        M+=M_L;
        M_sq += M_L*M_L;
        M_L=0;
    }
    E_mean = EnergySum / ((1.+maxMCC)*(L*L));
    E_mean_sq = EnergySum_sq / ((1.+maxMCC)*(L*L*L*L));

    M /=((1.+maxMCC)*L*L);
    M_sq /=((1.+maxMCC)*L*L*L*L);

    Cv = (E_mean_sq-E_mean*E_mean)/(Temperature*Temperature);
    chi =(M_sq -M*M)/Temperature;

    LocalExpectationValues(0) = accept;
    LocalExpectationValues(1) = maxMCC;
    LocalExpectationValues(2) = M_sq;
    LocalExpectationValues(3) = Cv;
    LocalExpectationValues(4) = chi;

    double E_mean_exact = -8. *(exp(8.)-exp(-8.))/(exp(8.)+exp(-8.)+6.);
    double M_mean_exact = 2.*(4*exp(8) + 8.)/(2.*exp(8.)+2.*exp(-8.)+12.);

    cout <<"exact E  = "	<< E_mean_exact<<endl;

    cout <<"exact M  = "	<<M_mean_exact<<endl;

    cout <<"Numeric E  = "	<<E_mean<<endl;

    cout<<"Numeric M = "<<setprecision(15) <<M << endl;
    cout<<"Numeric Expected Absolute Magnetization = " <<setprecision(15) <<M << endl;
    cout<<"Cycles = " <<setprecision(15) <<maxMCC<< endl;

    cout<<endl;

    //double Cv_exact=(128.*(exp(8./Temperature)+exp(-8./Temperature))/(2.*exp(8.)+2.*exp(-8.)+12.) - (8.*(-exp(8.)+exp(-8.))/(exp(8.)+exp(-8.)+6.))*(8.*(-exp(8.)+exp(-8.))/(exp(8.)+exp(-8.)+6.))) /(Temperature*Temperature);
    //double chi_exact =((2.*16*exp(8./Temperature) + 8.*2)/(2*exp(8./Temperature)+2*exp(-8./Temperature)+12.) - (2.*(4*exp(8) + 8.)/(2*exp(8)+2*exp(-8)+12)) * (2.*(4*exp(8) + 8.)/(2*exp(8)+2*exp(-8)+12))) / Temperature;
    //cout << "CV EXAC = "  << Cv_exact <<endl;
    //cout << "chi EXAC = "  << chi_exact <<endl;

}
//end of main



int magnet(mat& S,int L){
    int M_i = 0; //returns M_i to zero, or else it accumulates ALL spin values
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            M_i += S(i,j);
        }
    }
    if(M_i < 0){
        M_i = -M_i;
        return M_i;
    }
    else{
        return M_i;
    }
}


double ising_flip(mat& S, int L,double& E){
    int I;
    int J;
    int old_S;
    double max = RAND_MAX/double(L);
    I = (rand()-1)/max;
    J = (rand()-1)/max;
    old_S = S(I,J);
    S(I,J) *= -1;
    if ((I==0) or (J==0)){
        if ((I == 0) and (J!=0)){
            S(L,J) = S(I,J);
            E = 2*old_S*(S(I,J)*S(I+1,J)+S(I,J)*S(I,J+1)+S(I,J)*S(L,J)+S(I,J)*S(I,J-1));
        }
        else if ((J == 0) and (I!=0)){
            S(I,L) = S(I,J);
            E = 2*old_S*(S(I,J)*S(I+1,J)+S(I,J)*S(I,J+1)+S(I,J)*S(I-1,J)+S(I,J)*S(I,L));
        }
        else if ((J == 0) and (I==0)){
            S(I,L) = S(I,J);
            S(L,J) = S(I,J);
            E = 2*old_S*(S(I,J)*S(I+1,J)+S(I,J)*S(I,J+1)+S(I,J)*S(L,J)+S(I,J)*S(I,L));
        }
    }
    else{
        E = 2*old_S*(S(I,J)*S(I+1,J)+S(I,J)*S(I,J+1)+S(I,J)*S(I-1,J)+S(I,J)*S(I,J-1));
    }
    return S(I,J);
}

double energy(mat& S, double E, int L){
    E=0;
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            E += S(i,j)*S(i,j+1)+S(i,j)*S(i+1,j);
        }
    }
    return -E;
}
void WriteResultstoFile(double temperature, vec ExpectationValues)
{
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << ExpectationValues(0)/4.;
    ofile << setw(15) << setprecision(8) << ExpectationValues(3/4.);
    ofile << setw(15) << setprecision(8) << ExpectationValues(1)/4.;
    ofile << setw(15) << setprecision(8) << ExpectationValues(4)/4. << endl;
}
