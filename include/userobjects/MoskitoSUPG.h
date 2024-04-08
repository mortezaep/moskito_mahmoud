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

#pragma once

#include "GeneralUserObject.h"
#include "RankTwoTensor.h"


class MoskitoSUPG : public GeneralUserObject
{
public:
  static InputParameters validParams();
  
  MoskitoSUPG(const InputParameters & parameters);

  virtual void execute();
  virtual void initialize();
  virtual void finalize();

  void SUPGCalculator(const Real & diff, const Real & dt, const Elem * ele, const RealVectorValue & v, RealVectorValue & SUPG_coeff, Real & alpha, Real & CrNr) const;
  Real tau(const Real & alpha, const Real & diff, const Real & dt, const Real & v, const Real & h) const;

protected:
  Real Optimal(const Real & alpha) const;
  Real Temporal(const Real & v, const Real & h, const Real & diff, const Real & dt) const;
  Real DoublyAsymptotic(const Real & alpha) const;
  Real Critical(const Real & alpha) const;

  MooseEnum _method;
  enum M {optimal, doubly_asymptotic, critical, transient_brooks, transient_tezduyar};
};
