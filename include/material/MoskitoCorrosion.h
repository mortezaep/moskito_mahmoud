/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*  Co-developed by Sebastian Held                                        */
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

#pragma once

#include "Material.h"
#include "Function.h"


class MoskitoCorrosion : public Material
{
public:
  static InputParameters validParams();

  MoskitoCorrosion(const InputParameters & parameters);
  virtual void computeQpProperties() override;

  Real TemperatureWFinterface(const Real & Uto);
  void ktf(const Real & t, Real tt[], Real kt[]);

protected:
  // Fluid temperature at the center of the pipe
  const VariableValue & _p;

  // Well thermal resistivity
  MaterialProperty<Real> & _crr;

  // formation properties
  const Real & _PH;

  // imported from other materials
  //const MaterialProperty<Real> & _Dti;

};
