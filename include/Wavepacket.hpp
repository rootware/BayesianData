#pragma once

#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/FFT>

constexpr std::complex<double> I(0.0, 1.0);

constexpr int nMomentum = 11;
constexpr int nPosition = 120;

// constexpr int nStep = 20;

constexpr int nStep = 10000;

double mass = 1.0;
typedef Eigen::Matrix<double, nMomentum * 2 + 1 , 1> Feature;

class Wavepacket {

  double V0 = 10;
  double g = 0;
  Eigen::VectorXd p;
  double q;
  Eigen::MatrixXd H0;
  Eigen::MatrixXcd H1;
  Eigen::MatrixXcd H2;

  Eigen::VectorXcd w;
  Eigen::MatrixXcd v;

  Eigen::VectorXcd psi;

  
  Eigen::VectorXcd target;

public:

    Wavepacket()
  {
    this->g = 0.0;
    this->V0 = 10.0;
    std::cout<<"Should not be called" << std::endl;
    // init();
  }

  Wavepacket(double acceleration, double latticedepth){
    this->g = acceleration;
    this->V0 = latticedepth;

    double maxP = nMomentum - 1;
    p = Eigen::VectorXd::LinSpaced(nMomentum, -maxP, maxP);
    q = 0;
    H0 = p.cwiseAbs2().asDiagonal();
    H1 = Eigen::MatrixXcd::Zero(nMomentum, nMomentum); // sin
    Eigen::MatrixXcd ones = Eigen::MatrixXcd::Identity(nMomentum - 1, nMomentum - 1);
    // H1 sin after expanding kx and shaking function
    H1.block(0, 1, nMomentum - 1, nMomentum - 1) += V0 / 4 * I * ones;
    H1.block(1, 0, nMomentum - 1, nMomentum - 1) -= V0 / 4 * I * ones;//changed signs for V0 here
    H2 = Eigen::MatrixXcd::Zero(nMomentum, nMomentum); // cos
    H2.block(0, 1, nMomentum - 1, nMomentum - 1) += V0 / 4 * ones;
    H2.block(1, 0, nMomentum - 1, nMomentum - 1) += V0 / 4 * ones;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(H0 - H2);
    v = eigensolver.eigenvectors();


    groundState();
    target = v.col(3);
  }

  // This is not the most elegant way to re-initialize,
  // But I'll refactor the code later to be optimized
  // Eventually, I can just keep track of q and reset that,

  void init(){
    double maxP = nMomentum - 1;
    p = Eigen::VectorXd::LinSpaced(nMomentum, -maxP, maxP);
    q = 0;
    H0 = p.cwiseAbs2().asDiagonal();
    H1 = Eigen::MatrixXcd::Zero(nMomentum, nMomentum); // sin
    Eigen::MatrixXcd ones = Eigen::MatrixXcd::Identity(nMomentum - 1, nMomentum - 1);
    // H1 sin after expanding kx and shaking function
    H1.block(0, 1, nMomentum - 1, nMomentum - 1) += V0 / 4 * I * ones;
    H1.block(1, 0, nMomentum - 1, nMomentum - 1) -= V0 / 4 * I * ones;
    H2 = Eigen::MatrixXcd::Zero(nMomentum, nMomentum); // cos
    H2.block(0, 1, nMomentum - 1, nMomentum - 1) += V0 / 4 * ones;
    H2.block(1, 0, nMomentum - 1, nMomentum - 1) += V0 / 4 * ones;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(H0 - H2);
    v = eigensolver.eigenvectors();


    groundState();
    target = v.col(3);
    
  }

  void setParams(double acc, double latdepth){
    this->g = acc;
    this->V0 = latdepth;
    init();
  }


  Eigen::VectorXcd get_Psi() const{
    return psi;
  };

  void groundState()
  {
    
    psi = v.col(0);
    q = 0;

  }

  void reset(){

    init();
  }



  Eigen::VectorXd momentum() const
  {
    return psi.cwiseAbs2();
  }

  Eigen::VectorXd Blochpop() const{

    Eigen::VectorXd bloch(psi.cwiseAbs2());

    for (int i=0; i< psi.size();i++){
      bloch[i] = pow((psi.adjoint()*v.col(i)).norm(),2);
    }

    return bloch;
  }

  Eigen::VectorXd position() const
  {
    Eigen::FFT<double> fft;
    Eigen::VectorXcd result(psi.size());
    Eigen::VectorXcd pad = Eigen::VectorXcd::Zero(nPosition);
    pad.head(psi.size()) = psi;
    fft.fwd(result, pad);
    return result.cwiseAbs2();
  }
   
  double norm() const 
  {
    return psi.norm();
  }

  double fidelity() const
  {
    return (psi.adjoint() * target).norm();
  }

  Feature feature() const
  {
    Feature f;
    Eigen::VectorXcd result(psi);

    Eigen::VectorXcd dPsidt = -I * (H0-H2) * psi;
    Eigen::VectorXd dresult(psi.cwiseAbs2());

    for (long i = 0; i < psi.size() / 2; i++) {
      int idx = psi.size() - i - 1;
      result[i] = 0.5 * (psi[i] + psi[idx]);
      result[idx] = 0.5 * (psi[i] - psi[idx]);
    }

    for (long i = 0; i < psi.size() / 2; i++) {
      dresult[i] = 2.0 * ( conj(psi[i]) * dPsidt[i]).real();
    }
    

    
    f << result.cwiseAbs2(), dresult , 0.0;
    return f;
  }

  Eigen::VectorXcd rhs(const Eigen::VectorXcd& wavefn, double amplitude,
    double omega, double t)
  {
    double phi = amplitude * sin(omega * t);
    Eigen::MatrixXcd H = H0 + sin(phi) * H1 - cos(phi) * H2;
    return -I * H * wavefn;
  }

  void advancePsi(const double dt, const double amplitude, const double omega, const double t)
  {
    Eigen::VectorXcd k1 = rhs(psi, amplitude, omega, t);
    Eigen::VectorXcd k2 = rhs(psi + dt * k1 / 2, amplitude, omega, t + dt / 2);
    Eigen::VectorXcd k3 = rhs(psi + dt * k2 / 2, amplitude, omega, t + dt / 2);
    Eigen::VectorXcd k4 = rhs(psi + dt * k3, amplitude, omega, t + dt);
    psi += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    accelerate(-mass*g*dt);
  }

  void accelerate(double impulse)
  {
    q += impulse;
    Eigen::VectorXcd mom(p);
    mom.array() += q;
    H0 = mom.cwiseAbs2().asDiagonal();
  }

  void step(double amplitude, double omega)
  {
    double t = 0;
    double period = M_PI / omega; // half-cycle
    double dt = period / nStep;
    for (int n = 0; n < nStep; n++, t += dt)
      advancePsi(dt, amplitude, omega, t);
  }

  double getLatticeDepth(){
    return V0;
  }

  double getAcceleration(){
    return g;
  }

};
