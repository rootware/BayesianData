#include <iostream>
#include <fstream>
#include "Environment_MhollandRL.hpp"

using namespace std;


const std::vector<double> ACC_LIST = {0.0, 0.005, -0.005, 0.01, -0.01, 0.02, -0.02, 0.03, -0.03, 0.04, -0.04};

const std::vector<double> shakingfunctions {3.92699082, 3.92699082, 3.40339204, 0.,  3.66519143, 0., 3.14159265, 3.66519143, 3.92699082, 3.92699082, 3.66519143, 3.14159265, 3.14159265, 3.14159265, 3.14159265, 3.14159265, 2.35619449, 1.83259571, 
1.83259571,1.83259571,0.78539816, 3.40339204, 3.40339204, 2.87979327, 3.40339204, 3.40339204, 3.40339204, 3.40339204, 1.04719755, 1.04719755,0.78539816, 0.78539816};

int main()
{ 


    double acceleration = 0.0;
    double latticedepth = 10.0;
    std::string filename = "action" + std::to_string(acceleration) + "maxsteps"+std::to_string(maxSteps)+".txt";

    std::fstream myfile(filename, std::ios_base::in);
    //std::ofstream out("Momentum.dat");
    std::ofstream out3("CFI.dat");
    //std::ofstream out4("Blochpop.dat");
    //std::ofstream out5("ShakingFunction.dat");



    Environment environment(std::vector<double>{acceleration, latticedepth}); 
   // out << environment.getWavepacket().momentum().transpose() << endl;
    //out2 << environment.getWavepacket().position().transpose() << endl;
    out3 << environment.actualCFI() << '\t' << environment.latticeCFI() << '\t' <<environment.CFI_offdiagonal() << '\t' << environment.accQFI() << endl;
    //out4 << environment.getWavepacket().Blochpop().transpose() << endl;

    int sign = 1;
    for (auto &p : shakingfunctions){
        environment.applyAction(p*sign);
       // out << environment.getWavepacket().momentum().transpose() << endl;
       // out2 << environment.getWavepacket().position().transpose() << endl;
        out3 << environment.actualCFI() << '\t' << environment.latticeCFI() << '\t' <<environment.CFI_offdiagonal() << '\t' << environment.accQFI() << endl;
        //out4 << environment.getWavepacket().Blochpop().transpose() << endl;
       // out5 << p << endl;
       sign*=-1;
    };

    Wavepacket wavepacket;
        


    out3.close();
    exit(0);
};
