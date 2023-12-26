#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Wavepacket.hpp"
using namespace std;

#include <list>
#include <algorithm>

const std::vector<double> shakingfunctions {3.92699082, 3.92699082, 3.40339204, 0.,  3.66519143, 0., 3.14159265, 3.66519143, 3.92699082, 3.92699082, 3.66519143, 3.14159265, 3.14159265, 3.14159265, 3.14159265, 3.14159265, 2.35619449, 1.83259571, 
1.83259571,1.83259571,0.78539816, 3.40339204, 3.40339204, 2.87979327, 3.40339204, 3.40339204, 3.40339204, 3.40339204, 1.04719755, 1.04719755,0.78539816, 0.78539816};// acceleration 1p

const std::vector<double> shakingfunctions2p2e {3.40339204, 3.92699082, 1.83259571, 1.57079633, 2.0943951, 1.57079633, 2.61799388, 0., 3.92699082, 3.92699082,
2.0943951, 2.35619449, 2.35619449, 2.35619449, 2.35619449, 3.92699082, 3.92699082, 3.92699082, 2.61799388,
3.40339204, 3.92699082, 3.66519143, 3.66519143, 3.66519143, 3.66519143, 3.66519143, 3.66519143, 0.78539816,
0.78539816, 0.78539816, 0.78539816, 0.78539816}; // acceleration as 2 parameter


const std::vector<double> shakingfunctions2p2e_new {1.83259571, 0., 1.83259571, 2.87979327, 1.83259571, 1.83259571, 1.83259571, 3.40339204, 3.66519143,
3.40339204, 3.40339204, 3.14159265, 3.92699082, 3.92699082, 2.35619449, 2.35619449, 3.92699082, 3.92699082,
3.92699082, 3.66519143, 3.66519143, 3.66519143, 2.61799388, 3.66519143, 1.57079633, 1.57079633, 1.57079633,
1.04719755, 1.04719755, 1.04719755, 1.04719755, 1.57079633
};


const std::vector<double> shakingfunction2platt {1.83259571, 1.83259571, 1.83259571, 1.83259571, 1.83259571, 1.83259571, 1.04719755, 1.04719755, 1.04719755,0., 0., 0.52359878, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.26179939, 0., 0., 0., 0., 0., 0., 0., 0};

const int maxSteps = shakingfunctions.size();

const double omega = 11.5;

constexpr double ACC_NORM = 10.0/(1.9745e04*0.0054); // g divided by system units
constexpr double acclim = 0.1 * ACC_NORM;  // dacc is the epsilon in acceleration


int main()
{ 
    Eigen::VectorXd a_list = Eigen::VectorXd::LinSpaced(101,-0.0225,0.0225);
    Eigen::VectorXd V_list = Eigen::VectorXd::LinSpaced(51, 9.0, 11.0);

    std::vector<double> currentshaking = shakingfunctions2p2e_new;

  //  Eigen::VectorXd V_list(1) ;
  //  V_list <<10.0;

    std::ofstream aout("../Output/Acceleration.txt");
    aout << a_list<<endl;

    std::ofstream Vout("../Output/LatticeDepth.txt");
    Vout << V_list <<endl;
    
    std::ofstream out("../Output/AVIndex.txt");
    std::ofstream out2("../Output/MomentumProbability.txt");


    for (int aindex =0; aindex < a_list.size(); aindex++){
        for (int Vindex =0; Vindex < V_list.size() ; Vindex++){
            Wavepacket wp (a_list[aindex], V_list[Vindex]);
            int sign = 1;
            for (int i=0; i< int(shakingfunctions2p2e.size()); i++){
                wp.step(sign *shakingfunctions2p2e[i], omega);
                sign = -1*sign;
            };

            out << aindex <<'\t'<<Vindex << endl;
            out2 << wp.momentum().transpose()<<endl;
        };
    };

    exit(0);
};
