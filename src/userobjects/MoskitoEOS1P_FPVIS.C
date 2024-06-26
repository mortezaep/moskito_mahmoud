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

#include "MoskitoEOS1P_FPVIS.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_FPVIS);

InputParameters
MoskitoEOS1P_FPVIS::validParams()
{
  InputParameters params = MoskitoViscosity1P::validParams();

  params.addRequiredParam<UserObjectName>("fp", "The name of the user object for fluid properties");
  // params.addRequiredParam<UserObjectName>("SinglePhase_fp",
  //         "The name of the FluidProperties UserObject");

  return params;
}

MoskitoEOS1P_FPVIS::MoskitoEOS1P_FPVIS(const InputParameters & parameters)
  : MoskitoViscosity1P(parameters),
    _fp_eos(getUserObject<SinglePhaseFluidProperties>("fp"))
    
{
}

Real
MoskitoEOS1P_FPVIS::mu(Real pressure, Real temperature) const
{
  return _fp_eos.mu_from_p_T(pressure, temperature);
}

void
MoskitoEOS1P_FPVIS::dmu_dpT(
    Real pressure, Real temperature, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  _fp_eos.rho_from_p_T(pressure, temperature, mu, dmu_dp, dmu_dT);
}
