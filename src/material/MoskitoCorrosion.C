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

#include "MoskitoCorrosion.h"

registerMooseObject("MoskitoApp", MoskitoCorrosion);

InputParameters
MoskitoCorrosion::validParams()
{
  InputParameters params = Material::validParams();

    params.addClassDescription("Material for claculating corrosion rate");
    params.addRequiredCoupledVar("pressure", "Fluid temperature (K)");
    params.addRequiredParam<Real>("PH", "PH of environment");

    return params;
}

MoskitoCorrosion::MoskitoCorrosion(const InputParameters & parameters)
  : Material(parameters),
    _p(coupledValue("pressure")),
    _crr(declareProperty<Real>("corrosion_rate")),
    _PH(getParam<Real>("PH"))
    //_velocity(getMaterialProperty<Real>("velocity")),

{
}

void
MoskitoCorrosion::computeQpProperties()
{

}

Real
MoskitoCorrosion::TemperatureWFinterface(const Real & Uto)
{
  //Real Twf = 0.0;
  //Twf += _Rto * Uto * _ft * _Tf[_qp] + _lambda_form * _Tform[_qp];
  //Twf /= _Rto * Uto * _ft + _lambda_form;

  return Uto;
}

void
MoskitoCorrosion::ktf(const Real & t, Real tt[], Real kt[])
{
  double ktlow, ktup, tlow, tup;

	if(t>=20 && t<=40)
	{
		tlow = 20;
		tup = 40;
		ktlow = 4.762;
		ktup = 8.927;
	}
	if(t>=40 && t<=60)
	{
		tlow = 40;
		tup = 60;
		ktlow = 8.927;
		ktup = 10.695;
	}
	if(t>=60 && t<=80)
	{
		tlow = 60;
		tup = 80;
		ktlow = 10.695;
		ktup = 9.949;
	}
	if(t>=80 && t<=90)
	{
		tlow = 80;
		tup = 90;
		ktlow = 9.949;
		ktup = 6.250;
	}
	if(t>=90 && t<=120)
	{
		tlow = 90;
		tup = 120;
		ktlow = 6.250;
		ktup = 7.770;
	}
	if(t>=120 && t<=150)
	{
		tlow = 120;
		tup = 150;
		ktlow = 7.770;
		ktup = 5.203;
	}

	tt[0] = tlow;
	tt[1] = tup;
	kt[0] = ktlow;
	kt[1] = ktup;

	return;
}
