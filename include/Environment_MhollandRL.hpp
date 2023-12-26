#pragma once

#include <cmath>
#include <iostream>
#include "Wavepacket.hpp"

// Toggle switches time symmetry? 
// Can keep constant for RL?


constexpr size_t maxSteps = 30 ;
constexpr double ACC_NORM = 10.0/(1.9745e04*0.0054); // g divided by system units
constexpr double dacc = 0.00005 * ACC_NORM;  // dacc is the epsilon in acceleration
constexpr  double dLatticeDepth = 0.1;
  

class Environment {

  Wavepacket wavepacket_retarded;
  Wavepacket wavepacket_central;
  Wavepacket wavepacket_accelerated;

  Wavepacket wavepacket_shallow;
  Wavepacket wavepacket_deep;



  double time_of_flight = M_PI/omega * maxSteps;


  int splitting = 4;

  // Normalize CFI using the QFI limit for Mach Zehnder interferometer
  // For derivation, see my Overleaf note
  double CFI_normalization = pow(2*splitting * time_of_flight*time_of_flight,2)/4 *2;  
  // The two at the end is because I want to compare result to twice the MZ full interferometer limit
  // Feel free to play with this

  double centralacc = 0.0* ACC_NORM; // central_acc is the central value for acceleration for CFI

  double centralLatticeDepth = 10; 

  std::vector<double> actionsDone;
  
  const double omega = 11.5;
public:

  // *************
  // The following are methods needed for doing RL on any Environment class 
  // *************

  void reset()
  {

    wavepacket_retarded.init();
    wavepacket_central.init();
    wavepacket_accelerated.init();

    wavepacket_deep.init();
    wavepacket_shallow.init();


    actionsDone.clear();
 
  }

  Environment(){
    initWavepackets();
    reset();
  }
  

  Environment(std::vector<double> PARAMS){
    centralLatticeDepth = PARAMS[1];
    centralacc = PARAMS[0] *ACC_NORM;

    initWavepackets();
    reset();
  }

  Feature state() const
  {
    Feature f = wavepacket_central.feature();//( (wavepacket_central.feature()*10).array().round() / 10.).matrix();

    // new : flooring
    f[f.size()-1] = stepCount();

    return f;
  }

  
  double fidelity() const
  {
    /**
    Fidelity calculates a value f that should approach 1 with successful training
    */
    double f=inverseCFIreward();

    return f;
  }

  size_t stepCount() const
  {
    return actionsDone.size();
  }

  bool done() const
  {
    return stepCount() >= maxSteps;
  }

  double reward() const
  {
    // Use fidelity to shape a reward function

   if (done()) {
        double fid = fidelity();

        return (fid/(1-fid));
    } else
      return 0;
  }

  void applyAction(double amplitude)
  { 

    //wavepacket.step(toggle.flip() == POS ? amplitude : -amplitude, omega);
    

    wavepacket_retarded.step(amplitude, omega);
    wavepacket_central.step(amplitude, omega);
    wavepacket_accelerated.step(amplitude, omega);

    wavepacket_deep.step(amplitude, omega);
    wavepacket_shallow.step(amplitude, omega);

    actionsDone.push_back(amplitude);
  
  }
  
  std::vector<double> actionRecord() const
  {
    return actionsDone;
  }


  // **************
  // The following are Environment methods written to calculate QFI and CFI
  // or methods that involve multiple copies of Wavepackets.
  // These methods do not directly interface with RL. 
  // Only reason they are public is for the Evolve and data recording reasons.
  // Includes getter and setter methods.
  // **************
  void initWavepackets(){

    wavepacket_accelerated.setParams(centralacc + dacc , centralLatticeDepth); // epsilon more accelerated
    wavepacket_retarded.setParams(centralacc - dacc, centralLatticeDepth); // epsilon less accelerated
    wavepacket_central.setParams(centralacc, centralLatticeDepth); //central value for which we want CFI

    wavepacket_shallow.setParams(centralacc, centralLatticeDepth - dLatticeDepth);
    wavepacket_deep.setParams(centralacc, centralLatticeDepth + dLatticeDepth);

  };

  Wavepacket getWavepacket(){
    return wavepacket_central;
  }

  void setAcceleration(double acceleration){
    this-> centralacc = acceleration;
    reset();
  }

  double getAcceleration() const {
    return centralacc;
  }


  // ***********************
  // The following are functions for calculating CFI, and QFI
  // ***********
  
  // This is the reward we are now using, and is called by fideltiy() function
  double inverseCFIreward() const {

    // I denotes CFI, and F denotes QFI
    double Iaa= this->actualCFI();
    double IaV= this->CFI_offdiagonal();
    double IVV= this->latticeCFI();
    double Faa = this->accQFI();


    //double f = (Iaa - pow(IaV,2)/IVV)/CFI_normalization;

    double f = (Iaa)/CFI_normalization;
    return f;
  }


  
  double actualCFI() const {
    // CFI in acceleration lol, i.e. I_aa. I called it "actual" CFI and don't want to refactor code
    Eigen::VectorXd derivative_acc = ( wavepacket_accelerated.momentum() - wavepacket_retarded.momentum() )/(2*dacc);
    double CFI_acc = (wavepacket_central.momentum().cwiseInverse()).dot( derivative_acc.cwiseProduct(derivative_acc) );
    return CFI_acc;
  }

  double latticeCFI() const {
    Eigen::VectorXd derivative_latt = ( wavepacket_deep.momentum() - wavepacket_shallow.momentum() )/(2*dLatticeDepth);
    double CFI_latt = (wavepacket_central.momentum().cwiseInverse()).dot( derivative_latt.cwiseProduct(derivative_latt) );
    return CFI_latt;
  }

  double CFI_offdiagonal() const {
    Eigen::VectorXd derivative_acc = ( wavepacket_accelerated.momentum() - wavepacket_retarded.momentum() )/(2*dacc);
    Eigen::VectorXd derivative_latt = ( wavepacket_deep.momentum() - wavepacket_shallow.momentum() )/(2*dLatticeDepth);
    double I_aV = (wavepacket_central.momentum().cwiseInverse()).dot( derivative_acc.cwiseProduct(derivative_latt) );
    return I_aV;
  };



    // to be filled
  Eigen::VectorXcd dPsi_acc() const{
    return ( wavepacket_accelerated.get_Psi() - wavepacket_retarded.get_Psi() )/(2*dacc);
   }

  Eigen::VectorXcd dPsi_V0() const {
    return ( wavepacket_deep.get_Psi() - wavepacket_shallow.get_Psi() )/(2*dLatticeDepth);
  }

  double accQFI() const{
    Eigen::VectorXcd dP_a = this->dPsi_acc();
    Eigen::VectorXcd Psi = wavepacket_central.get_Psi();
    std::complex<double> dada = dP_a.adjoint()*dP_a;
    std::complex<double> daP = dP_a.adjoint()*Psi;
    std::complex<double> Pda = Psi.adjoint()*dP_a;

    return 4.0* ( dada - daP*Pda).real();
  }

  // The following is code for calculating latticeQFI and offdiagonal QFI
  // I don't need it right now so bleh. Also need to optimize it.
  
  double latticeQFI() const {
    Eigen::VectorXcd dP_V = this->dPsi_V0();
    Eigen::VectorXcd Psi = wavepacket_central.get_Psi();

    std::complex<double> dVdV = dP_V.adjoint()*dP_V;
    std::complex<double> dVP = dP_V.adjoint()*Psi;
    std::complex<double> PdV = Psi.adjoint()*dP_V;

    return 4.0* ( dVdV - dVP*PdV).real();
  }

  double offdiagonalQFI() const {
    Eigen::VectorXcd dP_a= dPsi_acc();
    Eigen::VectorXcd dP_V = dPsi_V0();
    Eigen::VectorXcd Psi = wavepacket_central.get_Psi();

    std::complex<double> dadV = dP_a.adjoint()*dP_V;
    std::complex<double> daP = dP_a.adjoint()*Psi;
    std::complex<double> PdV = Psi.adjoint()*dP_V;


    return 4.0* ( dadV - daP*PdV).real();
  }


// This method is likely going to be deprecated soon.
  double CFI_ratio() const {
    /**
    // This method uses the definition of Classical Fisher Information Matrix 
    // see paper and OneNote notes for details


    // Derivative of P(p|a,V_0) w.r.t acceleration
    Eigen::VectorXd derivative_acc = ( wavepacket_accelerated.momentum() - wavepacket_retarded.momentum() )/(2*dacc);

    // Calculating I_(aa) = \sum_p [P(p|a,V_0)]^(-1) [\partial_a P(p|a,V_0) ]^2
    double CFI_acc = (wavepacket_central.momentum().cwiseInverse()).dot( derivative_acc.cwiseProduct(derivative_acc) );


    // Derivative of P(p|a,V_0) w.r.t lattice depth and I_(V_0,V_0)
    Eigen::VectorXd derivative_latt = ( wavepacket_deep.momentum() - wavepacket_shallow.momentum() )/(2*dLatticeDepth);
    double CFI_latt = (wavepacket_central.momentum().cwiseInverse()).dot( derivative_latt.cwiseProduct(derivative_latt) );

    // The off diagononal entry to CFIM: I_(a,V) is

    double I_aV = (wavepacket_central.momentum().cwiseInverse()).dot( derivative_acc.cwiseProduct(derivative_latt) );
    */

    /* Reward function is either
      1. I_(aa)
      2. I_(aa) - I_(V,V)
      2. I_(a,a) - |I_(a,V_0)| // This is preferred
    */
    double CFI_ratioed = this->actualCFI() - abs( this->CFI_offdiagonal() );
    return CFI_ratioed;
  }
  
};
