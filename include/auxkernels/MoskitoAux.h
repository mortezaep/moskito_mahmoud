#pragma once

#include "AuxKernel.h"

class MoskitoAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MoskitoAux(const InputParameters & parameters);

  virtual ~MoskitoAux() {}

protected:
  virtual Real computeValue();
  void ktf(const Real & t, Real tt[], Real kt[]);
  void shear_stress(const Real & rom, const Real & vm, const Real & roughness, const Real & diameter, const Real & mum, Real & s2);
  void fco2f(const Real & t, const Real & p, const Real & mfco2, Real & fco2);
  void fphf(const Real & t, const Real & ph, Real fph[]);
  void crtf(const Real & t, const Real tt[], const Real kt[], const Real fph[], const Real & fco2, const Real & s, Real & crt);
  void gep(const Real & t, const Real & p, const Real & u, const Real & ph, Real & cr7);

//  Real _value;         ///< The value being set for this kernel

  const VariableValue & _p; ///< Coupled variable
  const VariableValue & _m; ///< Coupled variable
  const MaterialProperty<Real> & _T;
  const MaterialProperty<Real> & _dia;
  const MaterialProperty<Real> & _rom;
  const MaterialProperty<Real> & _mu;

  const Real & _PH;
  const Real & _roughness;

  Real tt[2], kt[2], fph[2];

  Real s2,fco2,crt,cr7;

};
