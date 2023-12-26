#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Wavepacket.hpp"
using namespace std;


const std::vector<double> shakingfunctions {3.92699082, 3.92699082, 3.40339204, 0.,  3.66519143, 0., 3.14159265, 3.66519143, 3.92699082, 3.92699082, 3.66519143, 3.14159265, 3.14159265, 3.14159265, 3.14159265, 3.14159265, 2.35619449, 1.83259571, 
1.83259571,1.83259571,0.78539816, 3.40339204, 3.40339204, 2.87979327, 3.40339204, 3.40339204, 3.40339204, 3.40339204, 1.04719755, 1.04719755,0.78539816, 0.78539816};

const std::vector<double> shakingfunctions2p2e {3.40339204, 3.92699082, 1.83259571, 1.57079633, 2.0943951, 1.57079633, 2.61799388, 0., 3.92699082, 3.92699082,
2.0943951, 2.35619449, 2.35619449, 2.35619449, 2.35619449, 3.92699082, 3.92699082, 3.92699082, 2.61799388,
3.40339204, 3.92699082, 3.66519143, 3.66519143, 3.66519143, 3.66519143, 3.66519143, 3.66519143, 0.78539816,
0.78539816, 0.78539816, 0.78539816, 0.78539816};

const std::vector<double> shakingfunction2platt {1.83259571, 1.83259571, 1.83259571, 1.83259571, 1.83259571, 1.83259571, 1.04719755, 1.04719755, 1.04719755,0., 0., 0.52359878, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.26179939, 0., 0., 0., 0., 0., 0., 0., 0};

const int maxSteps = shakingfunctions.size();


int main()
{ 
    Eigen::VectorXd VictorResult(11);
    //VictorResult << 0.00794241, 0.02077387 ,0.02588339 ,0.05700539, 0.28713433, 0.17525071, 0.0383645 , 0.35900257, 0.00822938, 0.00634325, 0.01283137; //1pacc
    //VictorResult <<1.58957656e-03, 2.89464310e-02 ,1.29810901e-02, 1.78989933e-01, 2.83210388e-02, 4.13504639e-02, 5.59004786e-01,5.04655221e-02, 7.88692407e-02 ,1.94715680e-02 ,9.68698096e-06;//2p
    
    VictorResult << 1.94622104e-05, 9.39636337e-04 ,9.41312557e-03, 3.71801355e-02,
 3.47168680e-02, 3.00312460e-01 ,5.98330257e-01, 7.04380573e-03, 1.19118070e-02 ,1.12934387e-04, 5.02321998e-07; //multiparam lattice rough dt values
    double acceleration = 0.0;

    double omega = 11.5;

    std::ofstream out("Momentum.dat");
    std::ofstream out2("Position.dat");
    std::ofstream out3("CFI.dat");
    std::ofstream out4("Blochpop.dat");
    std::ofstream out5("ShakingFunction.dat");

    double latticedepth = 10.0;
    Wavepacket wp (acceleration, latticedepth);
    
    out << wp.momentum().transpose() << endl;
    out2 << wp.position().transpose() << endl;
    //out3 << environment.actualCFI() << '\t' << environment.latticeCFI() << '\t' <<environment.CFI_offdiagonal() << '\t' << environment.accQFI() << endl;
    out4 << wp.Blochpop().transpose() << endl;

    int sign = 1;
    for (int i=0; i< int(shakingfunction2platt.size()); i++){
        wp.step(sign *shakingfunction2platt[i], omega);
        out << wp.momentum().transpose() << endl;
        out2 << wp.position().transpose() << endl;
      // out3 << environment.actualCFI() << '\t' << environment.latticeCFI() << '\t' <<environment.CFI_offdiagonal() << '\t' << environment.accQFI() << endl;
        out4 << wp.Blochpop().transpose() << endl;
        out5 << shakingfunction2platt[i] << endl;
        sign = -1*sign;
    };

    std::cout<< wp.getLatticeDepth() <<std::endl;
    std::cout << VictorResult.transpose()<<std::endl;
    std::cout << wp.momentum().transpose()<<std::endl;
    std::cout<< (VictorResult-wp.momentum() ).transpose()<<std::endl;



    out.close();
    out2.close();
   // out3.close();
    out4.close();
    out5.close();
    exit(0);
};
