#include <cmath>
#include <vector>
#include <iostream>
#include <cppad/cppad.hpp>

/*CALORICALLY PERFECT BOIVIN */

using namespace std;

using std::vector;
using CppAD::AD;
using CppAD::Value;
using CppAD::Var2Par;

extern "C" {
  void nspecies_(int& ns);
  void nreactions_(int& nr);
  void nelements_(int& ne);

  void getmolecularweights_(double* mw);
  void getmassfractions_(double* x, double* y);

  void getenthalpymass_(double* T, double* y, double* h);
  void getdensity_(double* T, double* p, double* y, double* rho);
  void gettemperature_(double* ei, double* y, double* T);
  void getpressure_(double *rho,double* T, double* y, double* pressure);

  void getsource_(double* p, double* h, double* T, double* y, double* fy);
  void getjacobian_(double* p, double* h, double* T, double* y, double* dfy);
  
  void getspecificheats_(double* T, double* s);
  void getenthalpiesofformation_(double* T, double* hf);
  void getinternalenergy_(double *rho,double *p, double* T, double* y, double* ei);

}

namespace mech {

  /***********************************************************************************/
  /* some constants                                                                  */
  /***********************************************************************************/
  double Tref =298.15; /*K*/
  double GasConstant = 8314.4621; /*J/(K kmol)*/
  double OneAtm      = 1.01325e5; /*kg/(m s^2)*/
  double OneThird    = 1.0/3.0;

  /***********************************************************************************/
  /* mech definitions                                                                */
  /***********************************************************************************/
  int kk = 9;
  int ii = 12;
  int mm = 3;

  vector<double> getMolecularWeights() {

    vector<double> mw(kk,0.0);
    mw[0] = 2.015880e+00; /*H2*/
    mw[1] = 1.007940e+00; /*H*/
    mw[2] = 3.199880e+01; /*O2*/
    mw[3] = 1.599940e+01; /*O*/
    mw[4] = 1.700734e+01; /*OH*/
    mw[5] = 3.300674e+01; /*HO2*/
    mw[6] = 3.401468e+01; /*H2O2*/
    mw[7] = 1.801528e+01; /*H2O*/
    mw[8] = 2.801348e+01; /*N2*/

    return(mw);
    
  };

  vector<double> getEnthalpyOfFormation() {

    vector<double> EF(kk,0.0);
    EF[0] = 0.0; /*H2*/
    EF[1] = 2.1799e+08; /*H*/
    EF[2] = 0.0; /*O2*/
    EF[3] = 1.557e+07; /*O*/
    EF[4] = 2.3135e+06; /*OH*/
    EF[5] = 3.8036e+05; /*HO2*/
    EF[6] = -3.8036e+05; /*H2O2*/
    EF[7] = -1.343e+07; /*H2O*/
    EF[8] = 0.0; /*N2*/

    return(EF);

  };

  /***********************************************************************************/
  /* thermodynamics                                                                  */
  /***********************************************************************************/
  template <class Type>
  void getSpecificHeats_R(Type& T, vector<Type>& cp0_R) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;

    if(tt0 < 1.0000e+03) {
      cp0_R[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 + 2.015720940000000e-08 * tt2 + -7.376117610000001e-12 * tt3;
    } else {
      cp0_R[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 + -1.795663940000000e-10 * tt2 + 2.002553760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 + 2.300816320000000e-18 * tt2 + -9.277323320000001e-22 * tt3;
    } else {
      cp0_R[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 + -4.735152350000000e-18 * tt2 + 4.981973570000000e-22 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[2] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 + -9.681295090000001e-09 * tt2 + 3.243728370000000e-12 * tt3;
    } else {
      cp0_R[2] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 + 2.094705550000000e-10 * tt2 + -2.167177940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[3] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 + -6.128066240000000e-09 * tt2 + 2.112659710000000e-12 * tt3;
    } else {
      cp0_R[3] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 + -1.001777990000000e-11 * tt2 + 1.228336910000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[4] = 4.125305610000000e+00 + -3.225449390000000e-03 * tt0 + 6.527646910000000e-06 * tt1 + -5.798536430000000e-09 * tt2 + 2.062373790000000e-12 * tt3;
    } else {
      cp0_R[4] = 2.864728860000000e+00 + 1.056504480000000e-03 * tt0 + -2.590827580000000e-07 * tt1 + 3.052186740000000e-11 * tt2 + -1.331958760000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[5] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 + -2.427638940000000e-08 * tt2 + 9.292251240000000e-12 * tt3;
    } else {
      cp0_R[5] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 + 1.142463700000000e-10 * tt2 + -1.079085350000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[6] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 + -2.157708130000000e-08 * tt2 + 8.624543630000000e-12 * tt3;
    } else {
      cp0_R[6] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 + 3.711859860000000e-10 * tt2 + -2.879083050000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[7] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 + -5.487970620000000e-09 * tt2 + 1.771978170000000e-12 * tt3;
    } else {
      cp0_R[7] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 + -9.704198700000000e-11 * tt2 + 1.682009920000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      cp0_R[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      cp0_R[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };


  };

  template <class Type>
  void getEnthalpies_RT(Type& T, vector<Type>& h0_RT) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;
    Type tt4 = 1.0 / T;

    if(tt0 < 1.0000e+03) {
      h0_RT[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 * 0.50 + -1.947815100000000e-05 * tt1 * OneThird + 2.015720940000000e-08 * tt2 * 0.25 + -7.376117610000001e-12 * tt3 * 0.20 + -9.179351730000000e+02 * tt4;
    } else {
      h0_RT[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 * 0.50 + 4.994567780000000e-07 * tt1 * OneThird + -1.795663940000000e-10 * tt2 * 0.25 + 2.002553760000000e-14 * tt3 * 0.20 + -9.501589220000000e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 * 0.50 + -1.995919640000000e-15 * tt1 * OneThird + 2.300816320000000e-18 * tt2 * 0.25 + -9.277323320000001e-22 * tt3 * 0.20 + 2.547365990000000e+04 * tt4;
    } else {
      h0_RT[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 * 0.50 + 1.615619480000000e-14 * tt1 * OneThird + -4.735152350000000e-18 * tt2 * 0.25 + 4.981973570000000e-22 * tt3 * 0.20 + 2.547365990000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[2] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 * 0.50 + 9.847302010000000e-06 * tt1 * OneThird + -9.681295090000001e-09 * tt2 * 0.25 + 3.243728370000000e-12 * tt3 * 0.20 + -1.063943560000000e+03 * tt4;
    } else {
      h0_RT[2] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 * 0.50 + -7.579666690000000e-07 * tt1 * OneThird + 2.094705550000000e-10 * tt2 * 0.25 + -2.167177940000000e-14 * tt3 * 0.20 + -1.088457720000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[3] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 * 0.50 + 6.643063960000000e-06 * tt1 * OneThird + -6.128066240000000e-09 * tt2 * 0.25 + 2.112659710000000e-12 * tt3 * 0.20 + 2.912225920000000e+04 * tt4;
    } else {
      h0_RT[3] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 * 0.50 + 4.194845890000000e-08 * tt1 * OneThird + -1.001777990000000e-11 * tt2 * 0.25 + 1.228336910000000e-15 * tt3 * 0.20 + 2.921757910000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[4] = 4.125305610000000e+00 + -3.225449390000000e-03 * tt0 * 0.50 + 6.527646910000000e-06 * tt1 * OneThird + -5.798536430000000e-09 * tt2 * 0.25 + 2.062373790000000e-12 * tt3 * 0.20 + 3.381538120000000e+03 * tt4;
    } else {
      h0_RT[4] = 2.864728860000000e+00 + 1.056504480000000e-03 * tt0 * 0.50 + -2.590827580000000e-07 * tt1 * OneThird + 3.052186740000000e-11 * tt2 * 0.25 + -1.331958760000000e-15 * tt3 * 0.20 + 3.718857740000000e+03 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[5] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 * 0.50 + 2.115828910000000e-05 * tt1 * OneThird + -2.427638940000000e-08 * tt2 * 0.25 + 9.292251240000000e-12 * tt3 * 0.20 + 2.948080400000000e+02 * tt4;
    } else {
      h0_RT[5] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 * 0.50 + -6.336581500000000e-07 * tt1 * OneThird + 1.142463700000000e-10 * tt2 * 0.25 + -1.079085350000000e-14 * tt3 * 0.20 + 1.118567130000000e+02 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[6] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 * 0.50 + 1.673357010000000e-05 * tt1 * OneThird + -2.157708130000000e-08 * tt2 * 0.25 + 8.624543630000000e-12 * tt3 * 0.20 + -1.770258210000000e+04 * tt4;
    } else {
      h0_RT[6] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 * 0.50 + -1.901392250000000e-06 * tt1 * OneThird + 3.711859860000000e-10 * tt2 * 0.25 + -2.879083050000000e-14 * tt3 * 0.20 + -1.786178770000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[7] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 * 0.50 + 6.520402110000000e-06 * tt1 * OneThird + -5.487970620000000e-09 * tt2 * 0.25 + 1.771978170000000e-12 * tt3 * 0.20 + -3.029372670000000e+04 * tt4;
    } else {
      h0_RT[7] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 * 0.50 + -1.640725180000000e-07 * tt1 * OneThird + -9.704198700000000e-11 * tt2 * 0.25 + 1.682009920000000e-14 * tt3 * 0.20 + -3.000429710000000e+04 * tt4;
    };

    if(tt0 < 1.0000e+03) {
      h0_RT[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 * 0.50 + -3.964481129109908e-06 * tt1 * OneThird + 5.642920880408571e-09 * tt2 * 0.25 + -2.445407041148433e-12 * tt3 * 0.20 + -1.020894198687962e+03 * tt4;
    } else {
      h0_RT[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 * 0.50 + -5.684761849244810e-07 * tt1 * OneThird + 1.009704225872734e-10 * tt2 * 0.25 + -6.753354387142974e-15 * tt3 * 0.20 + -9.227966980051905e+02 * tt4;
    };

  };

  template <class Type>
  void getEnthalpiesDerivatives(Type& T, vector<Type>& dh0dT) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;

    if(tt0 < 1.0000e+03) {
      dh0dT[0] = 2.344331120000000e+00 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 + 2.015720940000000e-08 * tt2 + -7.376117610000001e-12 * tt3;
    } else {
      dh0dT[0] = 3.337279200000000e+00 + -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 + -1.795663940000000e-10 * tt2 + 2.002553760000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[1] = 2.500000000000000e+00 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 + 2.300816320000000e-18 * tt2 + -9.277323320000001e-22 * tt3;
    } else {
      dh0dT[1] = 2.500000010000000e+00 + -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 + -4.735152350000000e-18 * tt2 + 4.981973570000000e-22 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[2] = 3.782456360000000e+00 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 + -9.681295090000001e-09 * tt2 + 3.243728370000000e-12 * tt3;
    } else {
      dh0dT[2] = 3.282537840000000e+00 + 1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 + 2.094705550000000e-10 * tt2 + -2.167177940000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[3] = 3.168267100000000e+00 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 + -6.128066240000000e-09 * tt2 + 2.112659710000000e-12 * tt3;
    } else {
      dh0dT[3] = 2.569420780000000e+00 + -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 + -1.001777990000000e-11 * tt2 + 1.228336910000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[4] = 4.125305610000000e+00 + -3.225449390000000e-03 * tt0 + 6.527646910000000e-06 * tt1 + -5.798536430000000e-09 * tt2 + 2.062373790000000e-12 * tt3;
    } else {
      dh0dT[4] = 2.864728860000000e+00 + 1.056504480000000e-03 * tt0 + -2.590827580000000e-07 * tt1 + 3.052186740000000e-11 * tt2 + -1.331958760000000e-15 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[5] = 4.301798010000000e+00 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 + -2.427638940000000e-08 * tt2 + 9.292251240000000e-12 * tt3;
    } else {
      dh0dT[5] = 4.017210900000000e+00 + 2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 + 1.142463700000000e-10 * tt2 + -1.079085350000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[6] = 4.276112690000000e+00 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 + -2.157708130000000e-08 * tt2 + 8.624543630000000e-12 * tt3;
    } else {
      dh0dT[6] = 4.165002850000000e+00 + 4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 + 3.711859860000000e-10 * tt2 + -2.879083050000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[7] = 4.198640560000000e+00 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 + -5.487970620000000e-09 * tt2 + 1.771978170000000e-12 * tt3;
    } else {
      dh0dT[7] = 3.033992490000000e+00 + 2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 + -9.704198700000000e-11 * tt2 + 1.682009920000000e-14 * tt3;
    };

    if(tt0 < 1.0000e+03) {
      dh0dT[8] = 3.298616281368291e+00 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 + 5.642920880408571e-09 * tt2 + -2.445407041148433e-12 * tt3;
    } else {
      dh0dT[8] = 2.926639911210682e+00 + 1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 + 1.009704225872734e-10 * tt2 + -6.753354387142974e-15 * tt3;
    };

  };

  template <class Type>
  void getEntropies_R(Type& T, vector<Type>& s0_R) {

    Type tt0 = T;
    Type tt1 = T * tt0;
    Type tt2 = T * tt1;
    Type tt3 = T * tt2;
    Type tt4 = 1.0 / T;
    Type tt5 = log(T);

    if(tt0 < 1.0000e+03) {
      s0_R[0] = 2.344331120000000e+00 * tt5 + 7.980520749999999e-03 * tt0 + -1.947815100000000e-05 * tt1 * 0.50 + 2.015720940000000e-08 * tt2 * OneThird + -7.376117610000001e-12 * tt3 * 0.25 + 6.830102380000000e-01;
    } else {
      s0_R[0] = 3.337279200000000e+00 * tt5 +  -4.940247310000000e-05 * tt0 + 4.994567780000000e-07 * tt1 * 0.50 + -1.795663940000000e-10 * tt2 * OneThird + 2.002553760000000e-14 * tt3 * 0.25 + -3.205023310000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[1] = 2.500000000000000e+00 * tt5 + 7.053328190000000e-13 * tt0 + -1.995919640000000e-15 * tt1 * 0.50 + 2.300816320000000e-18 * tt2 * OneThird + -9.277323320000001e-22 * tt3 * 0.25 + -4.466828530000000e-01;
    } else {
      s0_R[1] = 2.500000010000000e+00 * tt5 +  -2.308429730000000e-11 * tt0 + 1.615619480000000e-14 * tt1 * 0.50 + -4.735152350000000e-18 * tt2 * OneThird + 4.981973570000000e-22 * tt3 * 0.25 + -4.466829140000000e-01;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[2] = 3.782456360000000e+00 * tt5 + -2.996734160000000e-03 * tt0 + 9.847302010000000e-06 * tt1 * 0.50 + -9.681295090000001e-09 * tt2 * OneThird + 3.243728370000000e-12 * tt3 * 0.25 + 3.657675730000000e+00;
    } else {
      s0_R[2] = 3.282537840000000e+00 * tt5 +  1.483087540000000e-03 * tt0 + -7.579666690000000e-07 * tt1 * 0.50 + 2.094705550000000e-10 * tt2 * OneThird + -2.167177940000000e-14 * tt3 * 0.25 + 5.453231290000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[3] = 3.168267100000000e+00 * tt5 + -3.279318840000000e-03 * tt0 + 6.643063960000000e-06 * tt1 * 0.50 + -6.128066240000000e-09 * tt2 * OneThird + 2.112659710000000e-12 * tt3 * 0.25 + 2.051933460000000e+00;
    } else {
      s0_R[3] = 2.569420780000000e+00 * tt5 +  -8.597411370000000e-05 * tt0 + 4.194845890000000e-08 * tt1 * 0.50 + -1.001777990000000e-11 * tt2 * OneThird + 1.228336910000000e-15 * tt3 * 0.25 + 4.784338640000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[4] = 3.992015430000000e+00 * tt5 + -2.401317520000000e-03 * tt0 + 4.617938410000000e-06 * tt1 * 0.50 + -3.881133330000000e-09 * tt2 * OneThird + 1.364114700000000e-12 * tt3 * 0.25 + -1.039254580000000e-01;
    } else {
      s0_R[4] = 3.092887670000000e+00 * tt5 +  5.484297160000000e-04 * tt0 + 1.265052280000000e-07 * tt1 * 0.50 + -8.794615559999999e-11 * tt2 * OneThird + 1.174123760000000e-14 * tt3 * 0.25 + 4.476696100000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[5] = 4.301798010000000e+00 * tt5 + -4.749120510000000e-03 * tt0 + 2.115828910000000e-05 * tt1 * 0.50 + -2.427638940000000e-08 * tt2 * OneThird + 9.292251240000000e-12 * tt3 * 0.25 + 3.716662450000000e+00;
    } else {
      s0_R[5] = 4.017210900000000e+00 * tt5 +  2.239820130000000e-03 * tt0 + -6.336581500000000e-07 * tt1 * 0.50 + 1.142463700000000e-10 * tt2 * OneThird + -1.079085350000000e-14 * tt3 * 0.25 + 3.785102150000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[6] = 4.276112690000000e+00 * tt5 + -5.428224169999999e-04 * tt0 + 1.673357010000000e-05 * tt1 * 0.50 + -2.157708130000000e-08 * tt2 * OneThird + 8.624543630000000e-12 * tt3 * 0.25 + 3.435050740000000e+00;
    } else {
      s0_R[6] = 4.165002850000000e+00 * tt5 +  4.908316940000000e-03 * tt0 + -1.901392250000000e-06 * tt1 * 0.50 + 3.711859860000000e-10 * tt2 * OneThird + -2.879083050000000e-14 * tt3 * 0.25 + 2.916156620000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[7] = 4.198640560000000e+00 * tt5 + -2.036434100000000e-03 * tt0 + 6.520402110000000e-06 * tt1 * 0.50 + -5.487970620000000e-09 * tt2 * OneThird + 1.771978170000000e-12 * tt3 * 0.25 + -8.490322080000000e-01;
    } else {
      s0_R[7] = 3.033992490000000e+00 * tt5 +  2.176918040000000e-03 * tt0 + -1.640725180000000e-07 * tt1 * 0.50 + -9.704198700000000e-11 * tt2 * OneThird + 1.682009920000000e-14 * tt3 * 0.25 + 4.966770100000000e+00;
    };

    if(tt0 < 1.0000e+03) {
      s0_R[8] = 3.298616281368291e+00 * tt5 + 1.408708904146039e-03 * tt0 + -3.964481129109908e-06 * tt1 * 0.50 + 5.642920880408571e-09 * tt2 * OneThird + -2.445407041148433e-12 * tt3 * 0.25 + 3.950623591964732e+00;
    } else {
      s0_R[8] = 2.926639911210682e+00 * tt5 +  1.487977101178227e-03 * tt0 + -5.684761849244810e-07 * tt1 * 0.50 + 1.009704225872734e-10 * tt2 * OneThird + -6.753354387142974e-15 * tt3 * 0.25 + 5.980528055036107e+00;
    };

  };

  template <class Type>
  void getGibbsFunctions_RT(Type& p, Type& T,vector<Type>& g0_RT) {

    vector<Type> h0_RT(kk, 0.0);
    vector<Type> s0_R(kk,  0.0);

    getEnthalpies_RT(T, h0_RT);
    getEntropies_R(T, s0_R);
    for(int k = 0; k < kk; ++k) { g0_RT[k] = h0_RT[k] - s0_R[k]; }

  };

  template <class Type>
  void getEqConstants(Type& p, Type& T,vector<Type>& keqs) {

    double       p0 = OneAtm;
    Type         RT = GasConstant * T;
    Type         C0 = p0 / RT;
    //Type				 referenceTemp=Tref;
    vector<Type> g0_RT(kk, 0.0);

    getGibbsFunctions_RT(p, T, g0_RT);
    for(int k = 0; k < kk; ++k) { g0_RT[k] = exp(g0_RT[k]); }

    keqs[0] = (g0_RT[3] * g0_RT[4]);
    keqs[1] = (g0_RT[1] * g0_RT[4]);
    keqs[2] = (g0_RT[1] * g0_RT[7]);
    keqs[5] = (g0_RT[0] * g0_RT[2]);
    keqs[7] = (C0 * g0_RT[7]);
    keqs[8] = (C0 * g0_RT[0]);

    keqs[0] /= (g0_RT[1] * g0_RT[2]);
    keqs[1] /= (g0_RT[0] * g0_RT[3]);
    keqs[2] /= (g0_RT[0] * g0_RT[4]);
    keqs[5] /= (g0_RT[1] * g0_RT[5]);
    keqs[7] /= (g0_RT[1] * g0_RT[4]);
    keqs[8] /= (g0_RT[1] * g0_RT[1]);

    keqs[3] = 0.0;
    keqs[4] = 0.0;
    keqs[6] = 0.0;
    keqs[9] = 0.0;
    keqs[10] = 0.0;
    keqs[11] = 0.0;

  };

  /***********************************************************************************/
  /* kinetics                                                                        */
  /***********************************************************************************/
  template <class Type>
  void getFalloffRates(Type& T, vector<Type>& C, vector<Type>& kfwd) {

    int          nfall = 2; 
    Type         tlog  = log(T);
    Type         rt    = 1.0 / T;
    Type         lpr;
    Type         cc;
    Type         nn;
    Type         f1;
    Type         lgf;
    vector<Type> pr(nfall);
    vector<Type> khi(nfall);
    vector<Type> klo(nfall);
    vector<Type> work(nfall);

    khi[0] =  exp(3.859056134271e+01 + -1.4000e+00 * tlog);
    khi[1] =  exp(5.505747506612e+01 + -1.9000e+00 * tlog - 2.497094791014803e+04 * rt);

    klo[0] =  exp(1.535237777756e+01 + 4.4000e-01 * tlog);
    klo[1] =  exp(3.780453580568e+01 + -1.3900e+00 * tlog - 2.582728713141888e+04 * rt);

    pr[0] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.6e+01 * C[7] + 1.0e+00 * C[8]) * klo[0] / khi[0];
    pr[1] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 6.0e+00 * C[7] + 1.0e+00 * C[8]) * klo[1] / khi[1];

    work[0] =  0.5;
    work[1] =  (1.0 - 7.3500e-01) * exp(-T * 1.0638e-02) + 7.3500e-01 * exp(-T * 5.6948e-04) + exp(-5.1820e+03/ T);

    for(int r = 0; r < nfall; ++r) {
      lpr     =  log10(pr[r]);
      cc      = -0.40 - 0.67 * log10(work[r]);
      nn      =  0.75 - 1.27 * log10(work[r]);
      f1      =  (lpr + cc)/(nn - 0.14 * (lpr + cc));
      work[r] =  log10(work[r])/(1 + f1 * f1);
      work[r] =  pow(10.0, work[r]);
      work[r] =  (pr[r] * work[r])/(1 + pr[r]);
    };

    kfwd[3] = khi[0] * work[0];
    kfwd[11] = khi[1] * work[1];

  };


  template <class Type>
  void updateRateConstants(Type& p, Type& T, vector<Type>& C,
			   vector<Type>& kfwd, vector<Type>& krev) {

    Type         tlog = log(T);
    Type         rt   = 1.0 / T;
    vector<Type> keqs(ii, 0.0);

    getEqConstants(p, T, keqs);

    kfwd[0] =  exp(3.119206719853e+01  - 7.0000e-01 * tlog - 8.589852132466873e+03 * rt);
    kfwd[1] =  exp(3.923951576293e+00  + 2.6700e+00 * tlog - 3.165568582001234e+03 * rt);
    kfwd[2] =  exp(1.397251430677e+01  + 1.3000e+00 * tlog - 1.829342634203600e+03 * rt);
    kfwd[4] =  exp(2.498312483765e+01  - 1.479350059217902e+02 * rt);
    kfwd[5] =  exp(2.353266853231e+01  - 4.137369271308603e+02 * rt);
    kfwd[6] =  exp(2.408710743206e+01  + 2.501665140791248e+02 * rt);
    kfwd[7] =  exp(3.822765584902e+01  - 2.0000e+00 * tlog);
    kfwd[8] =  exp(2.789338538040e+01  - 1.0000e+00 * tlog);
    kfwd[9] =  exp(2.182852266833e+01  - 6.975797027206365e+02 * rt);
    kfwd[10] =  exp(1.890310689320e+01 + 6.1000e-01 * tlog - 1.204407438455940e+04 * rt);

    getFalloffRates(T, C, kfwd);

    kfwd[7] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[7];
    kfwd[8] = (2.5e+00 * C[0] + 1.0e+00 * C[1] + 1.0e+00 * C[2] + 1.0e+00 * C[3] + 1.0e+00 * C[4] + 1.0e+00 * C[5] + 1.0e+00 * C[6] + 1.2e+01 * C[7] + 1.0e+00 * C[8]) * kfwd[8];

    for(int i = 0; i < ii; ++i) { krev[i] = kfwd[i] * keqs[i]; }

  }

  template <class Type>
  void updateSpeciesSourceTerm(Type& p, Type& h, Type& T,vector<Type>& Y,vector<Type>& FY) {

    Type           rho;
    Type           W;
    vector<Type>   C(kk,    0.0);
    vector<Type>   kfwd(ii, 0.0);
    vector<Type>   krev(ii, 0.0);
    vector<Type>   Rfwd(ii, 0.0);
    vector<Type>   Rrev(ii, 0.0);
    vector<Type>   Rnet(ii, 0.0);
    vector<Type>   omega(kk,0.0);
    vector<double> mw = getMolecularWeights();

    //mech::getTemperature(h, Told, Y, T);
      
for(int k = 0; k < kk; ++k) { W += Y[k] / mw[k]; }
    W   = 1.0/W;
    rho = (p * W)/(GasConstant * T);
    for(int k = 0; k < kk; ++k) { C[k] = rho * Y[k] / mw[k]; }
    mech::updateRateConstants(p, T, C, kfwd, krev);

    Rfwd[0] = kfwd[0] * C[1] * C[2];
    Rfwd[1] = kfwd[1] * C[0] * C[3];
    Rfwd[2] = kfwd[2] * C[0] * C[4];
    Rfwd[3] = kfwd[3] * C[1] * C[2];
    Rfwd[4] = kfwd[4] * C[1] * C[5];
    Rfwd[5] = kfwd[5] * C[1] * C[5];
    Rfwd[6] = kfwd[6] * C[5] * C[4];
    Rfwd[7] = kfwd[7] * C[1] * C[4];
    Rfwd[8] = kfwd[8] * C[1] * C[1];
    Rfwd[9] = kfwd[9] * C[5] * C[5];
    Rfwd[10] = kfwd[10] * C[0] * C[5];
    Rfwd[11] = kfwd[11] * C[6];

    Rrev[0] = krev[0] * C[3] * C[4];
    Rrev[1] = krev[1] * C[1] * C[4];
    Rrev[2] = krev[2] * C[1] * C[7];
    Rrev[5] = krev[5] * C[0] * C[2];
    Rrev[7] = krev[7] * C[7];
    Rrev[8] = krev[8] * C[0];

    for(int i = 0; i < ii; ++i) { Rnet[i] = Rfwd[i] - Rrev[i]; }

    omega[0] =  + Rnet[8] + Rnet[5] - Rnet[1] - Rnet[2] - Rnet[10] ;
    omega[1] =  + Rnet[1] + Rnet[2] + Rnet[10] - Rnet[0] - Rnet[3] - Rnet[4] - Rnet[5]
      - Rnet[7] - Rnet[8] - Rnet[8] ;
    omega[2] =  + Rnet[5] + Rnet[6] + Rnet[9] - Rnet[0] - Rnet[3] ;
    omega[3] =  + Rnet[0] - Rnet[1] ;
    omega[4] =  + Rnet[0] + Rnet[1] + Rnet[4] + Rnet[4] + Rnet[11] + Rnet[11] - Rnet[2]
      - Rnet[6] - Rnet[7] ;
    omega[5] =  + Rnet[3] - Rnet[4] - Rnet[5] - Rnet[6] - Rnet[9] - Rnet[9] - Rnet[10] ;
    omega[6] =  + Rnet[9] + Rnet[10] - Rnet[11] ;
    omega[7] =  + Rnet[7] + Rnet[2] + Rnet[6] ;

    for(int k = 0; k < kk; ++k) { FY[k] = mw[k] * omega[k] / rho; }
    
  }

  template <class Type>
  vector<Type> transpose(vector<Type>& mat)
  {

    int nn = mat.size();
    int n  = sqrt(nn);
    
    for (int i = 0; i < n; i++)
      {
	for (int j = i+1; j < n; j++)
	  {
	    Type temp  = mat[i*n+j];
	    mat[i*n+j] = mat[j*n+i];
	    mat[j*n+i] = temp;
	  }
      }
    
    return(mat);

  }
  
} // namespace mech

void nspecies_(int& nspec) {

  nspec = mech::kk;
  
}

void nreactions_(int& nreac) {

  nreac = mech::ii;
  
}

void nelements_(int& nelem) {

  nelem = mech::mm;
  
}

void getmolecularweights_(double* mw) {

  vector<double> mwts = mech::getMolecularWeights();
  for(int k = 0; k < mech::kk; ++k) { mw[k] = mwts[k]; }
  
}

void getmassfractions_(double* x, double* y) {

  double         mmw = 0.0;
  vector<double> mw  = mech::getMolecularWeights();

  for(int k = 0; k < mech::kk; k++) {
    double xk = max(x[k], 0.0);
    mmw  += mw[k] * xk;
  }
  
  mmw = 1.0/mmw;
  for(int k = 0; k < mech::kk; ++k) { y[k] = max(mw[k] * x[k] * mmw, 0.0); }
  
}

void getenthalpymass_(double* T, double* y, double* h) {

  double         RT = mech::GasConstant * (*T);
  vector<double> mw = mech::getMolecularWeights();
  vector<double> hi(mech::kk, 0.0);
  vector<double> h0_RT(mech::kk);
  vector<double> cpi_R(mech::kk);

  mech::getEnthalpies_RT(mech::Tref, h0_RT);
  mech::getSpecificHeats_R(mech::Tref, cpi_R);
  *h = 0.0;
  for(int k = 0; k < mech::kk; ++k) {
    hi[k]  = RT*h0_RT[k]/mw[k]+cpi_R[k]*mech::GasConstant*(*T-mech::Tref)/mw[k];
    *h    += hi[k] * y[k];
  }
  
}

void getdensity_(double* T, double* p, double* y, double* rho) {

  double         RT  = mech::GasConstant * (*T);
  double         mmw = 0.0;
  vector<double> mw  = mech::getMolecularWeights();
  for(int k = 0; k < mech::kk; ++k) {
    mmw += y[k]/mw[k];
  }
  mmw = 1.0/mmw;

  *rho = (*p) * mmw / RT;
  
}

void getpressure_(double* rho, double* T, double *y, double* p){
	vector<double> Y(mech::kk, 0.0);
  double         mmw = 0.0;
  vector<double> mw  = mech::getMolecularWeights();
  for(int k = 0; k < mech::kk; ++k) {mmw += y[k]/mw[k];}
  mmw = 1.0/mmw;
  (*p)=(*rho)/mmw*mech::GasConstant*(*T);

}

void gettemperature_(double* ei, double* y, double* T) {

  vector<double> Y(mech::kk, 0.0);
  double         mmw = 0.0;
  double cp=0.0;
  vector<double> mw  = mech::getMolecularWeights();
  vector<double> cpi_R(mech::kk);
  mech::getSpecificHeats_R(mech::Tref, cpi_R);
  for(int k = 0; k < mech::kk; ++k) {Y[k] = y[k];}
  for(int k = 0; k < mech::kk; ++k) {cp+=Y[k]/mw[k]*cpi_R[k]*mech::GasConstant;} 
  for(int k = 0; k < mech::kk; ++k) {mmw += y[k]/mw[k];}
  mmw = 1.0/mmw;
  
//  getpressure_(rho,T,y,p)
  (*T)=(*ei+cp*mech::Tref)/(cp-mech::GasConstant/mmw);
//  (*T)=(*ei)/(cp-mech::GasConstant/mmw); 
 
}

void getinternalenergy_(double* rho,double *p, double* T, double* y, double* ei){

double internalEnergy=0.0;
vector<double> Y(mech::kk);
vector<double> cpi_R(mech::kk);
double cp=0.0;
vector<double> mw = mech::getMolecularWeights();
double         RT = mech::GasConstant * (*T);

mech::getSpecificHeats_R(mech::Tref, cpi_R);
internalEnergy=0.;
for(int k = 0; k < mech::kk; ++k) { Y[k] = y[k]; }
for(int k = 0; k < mech::kk; ++k) {cp+=Y[k]/mw[k]*cpi_R[k]*mech::GasConstant;}

//internalEnergy=cp*(*T);
internalEnergy=cp*(*T-mech::Tref);
internalEnergy-=(*p)/(*rho);
*ei=internalEnergy;

}

//void getspecificheats_(double* T, double* s){
//vector<double> cpi_R(mech::kk);
//vector<double> cpi(mech::kk);
//vector<double> mw = mech::getMolecularWeights();
//
//mech::getSpecificHeats_R(mech::Tref, cpi_R);
//for(int k = 0; k < mech::kk; ++k) {
//    cpi[k]  = mech::GasConstant * cpi_R[k] / mw[k];
//  }
//
// copy(cpi.begin(),cpi.end(),s);
//
//}

void getenthalpiesofformation_(double* T, double* hf){
vector<double> h0_RT(mech::kk);
vector<double> hif(mech::kk);
vector<double> mw = mech::getMolecularWeights();

mech::getEnthalpies_RT(mech::Tref, h0_RT);
for(int k = 0; k < mech::kk; ++k) {
 //   hif[k]  = mech::GasConstant/mw[k] *(h0_RT[k]*(*T)) ;
    hif[k]  = mech::GasConstant/mw[k] *(h0_RT[k]*(mech::Tref)) ;
  }
 copy(hif.begin(),hif.end(),hf);

}

void getsource_(double* p, double* h, double* T, double* y, double* fy) {

  vector<double> Y(mech::kk,  0.0);
  vector<double> FY(mech::kk, 0.0);

  for(int k = 0; k < mech::kk; ++k) { Y[k] = y[k]; }
  mech::updateSpeciesSourceTerm(*p, *h, *T, Y, FY);
  copy(FY.begin(), FY.end(), fy);

}

void getjacobian_(double* p, double* h, double* Told, double* y, double* dfy) {

  int                  kk = mech::kk;
  vector<double>       Y(kk, 0.0);
  vector<double>       FY(kk,0.0);
  vector<double>       DFY(kk*kk,0.0);
  vector< AD<double> > a_Y(kk, 0.0);
  vector< AD<double> > a_FY(kk,0.0);
  vector< AD<double> > a_p(1,0.0);
  vector< AD<double> > a_h(1,0.0);
  vector< AD<double> > a_T(1,0.0);
  AD<double> pressure;
  AD<double> enthalpy;
  AD<double> temperature;

  /* AD stuff */
  for(int k = 0; k < mech::kk; ++k) { Y[k] = y[k]; }
  copy(Y.begin(), Y.end(), a_Y.begin());

  /* AD Stuff */
  CppAD::Independent(a_Y); 
  pressure=*p;
  temperature=*Told;
  enthalpy=*h;
  mech::updateSpeciesSourceTerm(pressure,enthalpy,temperature, a_Y, a_FY);
  CppAD::ADFun<double> tapef(a_Y, a_FY);
  DFY  = tapef.Jacobian(Y);
  DFY  = mech::transpose(DFY);

  copy(DFY.begin(), DFY.end(), dfy);
   
}

