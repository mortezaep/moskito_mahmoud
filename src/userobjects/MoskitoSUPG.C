/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*                                                                        */
/*  This program is free software: you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation, either version 3 of the License, or     */
/*  (at your option) any later version.                                   */
/*                                                                        */
/*  This program is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          */
/*  GNU General Public License for more details.                          */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>  */
/**************************************************************************/

#include "MoskitoSUPG.h"

registerMooseObject("MoskitoApp", MoskitoSUPG);

InputParameters
MoskitoSUPG::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  MooseEnum Method("optimal doubly_asymptotic critical transient_brooks transient_tezduyar", "optimal");
  params.addParam<MooseEnum>("supg_coeficient", Method,
        "The method for calculating SU/PG coefficent (tau)");

  return params;
}

MoskitoSUPG::MoskitoSUPG(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _method(getParam<MooseEnum>("supg_coeficient"))
{
}

void
MoskitoSUPG::execute(){}

void
MoskitoSUPG::initialize(){}

void
MoskitoSUPG::finalize(){}

void
MoskitoSUPG::SUPGCalculator(const Real & diff, const Real & dt, const Elem * ele, const RealVectorValue & v, RealVectorValue & SUPG_coeff, Real & alpha, Real & CrNr) const
{
  Real v_n = v.norm();
  Real h_n = ele->volume();

  if (v_n != 0.0)
  {
    if (diff == 0.0)
      alpha = std::numeric_limits<Real>::max();
    else
      alpha = 0.5 * v_n * h_n / diff;

    CrNr = v_n * dt / h_n;
    // SUPG_coeff = tau(alpha, diff, dt, v_n, h_n) * v;
    SUPG_coeff = tau(alpha, diff, dt, v_n, h_n) * v/v_n;
  }
  else
  {
    alpha = 0.0;
    SUPG_coeff.zero();
  }
}

Real
MoskitoSUPG::tau(const Real & alpha, const Real & diff, const Real & dt, const Real & v, const Real & h) const
{
  Real tau = 0.0;

  switch (_method)
  {
    case M::optimal:
      tau += Optimal(alpha) * h / (2.0 * v);
      break;
    case M::doubly_asymptotic:
      tau += DoublyAsymptotic(alpha) * h / (2.0 * v);
      break;
    case M::critical:
      tau += Critical(alpha) * h / (2.0 * v);
      break;
    case M::transient_brooks:
      // Brooks & Hughes 1982
      tau += Optimal(alpha) * h / (sqrt(15.0) * v);
      break;
    case M::transient_tezduyar:
      tau += Temporal(v, h, diff, dt);
      break;
  }

  return tau;
}

Real
MoskitoSUPG::Optimal(const Real & alpha) const
{
  Real s = 0.0;
  if (alpha < 0.01)
    s = alpha * (1.0 / 3.0 + alpha * alpha * (-1.0 / 45.0 + 18.0 / 8505.0 * alpha * alpha)); //taylor expansion
  else
    s = 1.0 / std::tanh(alpha) - 1.0 / alpha;
  return s;
}

Real
MoskitoSUPG::Temporal(const Real & v, const Real & h, const Real & diff, const Real & dt) const
{
  Real s1, s2, s3;
  // Tezduyar & Osawa 2000
  s1 = 2.0 * v / h;
  s2 = (dt!=0.0 ? 2.0 / dt : 0.0);
  s3 = 4.0 * diff / ( h * h );

  return (1.0 / sqrt(s1 * s1 + s2 * s2 + s3 * s3));
}

Real
MoskitoSUPG::DoublyAsymptotic(const Real & alpha) const
{
  Real s = 0.0;
  if (alpha <= 3.0)
    s = alpha / 3.0;
  else
    s = 1.0;
  return s;
}

Real
MoskitoSUPG::Critical(const Real & alpha) const
{
  Real s = 0.0;
  if (alpha <= 1.0)
    s = 0.0;
  else
    s = 1.0 - 1.0 / alpha;
  return s;
}
