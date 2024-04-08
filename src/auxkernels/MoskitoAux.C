#include "MoskitoAux.h"

registerMooseObject("MoskitoApp", MoskitoAux);

InputParameters
MoskitoAux::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable");

  params.addRequiredParam<Real>("PH", "PH of liquid phase");
  params.addRequiredParam<Real>("roughness", "roughness again");

  return params;
}

MoskitoAux::MoskitoAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  _p(coupledValue("pressure")),
  _m(coupledValue("massrate")),
  _T(getMaterialProperty<Real>("temperature")),
  _dia(getMaterialProperty<Real>("well_diameter")),
  _rom(getMaterialProperty<Real>("density")),
  _mu(getMaterialProperty<Real>("viscosity")),
  _PH(getParam<Real>("PH")),
  _roughness(getParam<Real>("roughness"))

{
}

Real
MoskitoAux::computeValue()
{

    ktf(_T[_qp]-273.15, tt, kt);
	  fco2f(_T[_qp]-273.15, _p[_qp]/100000.0, 100.0, fco2);
	  fphf(_T[_qp]-273.15, _PH, fph);
    shear_stress(_rom[_qp], _m[_qp]/(3.14/4.0*_dia[_qp]*_dia[_qp]), _roughness, _dia[_qp], _mu[_qp], s2);
	  crtf( _T[_qp]-273.15, tt, kt, fph, fco2, s2, crt);


    //ktf(60.0, tt, kt);
	  //fco2f(60.0, 100.0, 100.0, fco2);
	  //fphf(60.0, 6.0, fph);
    //shear_stress(_rom[_qp], _m[_qp]/(3.14/4.0*_dia[_qp]*_dia[_qp]), _roughness, _dia[_qp], _mu[_qp], s2);
	  //crtf(60.0, tt, kt, fph, fco2, 50.0, crt);

    //gep(_T[_qp]-273.15, _p[_qp]/1000000.0, _m[_qp]/(3.14/4.0*_dia[_qp]*_dia[_qp]), _PH, cr7);
    gep(_T[_qp]-273.15, _p[_qp]/1000000.0, _m[_qp]/_rom[_qp]/(3.14/4.0*_dia[_qp]*_dia[_qp]), _PH, cr7);

    return cr7;

}

void
MoskitoAux::ktf(const Real & t, Real tt[], Real kt[])
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

void
MoskitoAux::shear_stress(const Real & rom, const Real & vm, const Real & roughness, const Real & diameter, const Real & mum, Real & s2)
{
  double re, f;

	re = rom*vm*diameter/mum;

	f = 0.001375*(1.0 + pow((20000.0*roughness/diameter + 1000000.0/re), 0.33));

	s2 = 0.5*rom*f*vm*vm;

	return;
}

void
MoskitoAux::fco2f(const Real & t, const Real & p, const Real & mfco2, Real & fco2)
{
  double a, pbar, pco2, tk;

	pbar = p;
	tk = t+273.15;

	if(pbar <= 250)
	{
		a = pow(10.0, pbar*(0.0031 - 1.4/tk));
	}
	if(pbar > 250)
	{
		a = pow(10.0, 250*(0.0031 - 1.4/tk));
	}

	pco2 = pbar*mfco2;

	fco2 = pco2*a/100.0;

	return;
}

void
MoskitoAux::fphf(const Real & t, const Real & ph, Real fph[])
{
  double fphlow, fphup;

	if(t>=20 && t<=40 && ph>=3.5 && ph<=4.6)
	{
		fphlow = 2.0676 - 0.2309*ph;
		fphup = 2.0676 - 0.2309*ph;
	}
	if(t>=20 && t<=40 && ph>=4.6 && ph<=6.5)
	{
		fphlow = 5.1885 - 1.2353*ph + 0.0708*ph*ph;
		fphup = 5.1885 - 1.2353*ph + 0.0708*ph*ph;
	}
	if(t>=40 && t<=60 && ph>=3.5 && ph<=4.6)
	{
		fphlow = 2.0676 - 0.2309*ph;
		fphup = 1.836 - 0.1818*ph;
	}
	if(t>=40 && t<=60 && ph>=4.6 && ph<=6.5)
	{
		fphlow = 5.1885 - 1.2353*ph + 0.0708*ph*ph;
		fphup = 15.444 - 6.1291*ph + 0.8204*ph*ph - 0.0371*ph*ph*ph;
	}
	if(t>=60 && t<=80 && ph>=3.5 && ph<=4.6)
	{
		fphlow = 1.836 - 0.1818*ph;
		fphup = 2.6727 - 0.3636*ph;
	}
	if(t>=60 && t<=80 && ph>=4.6 && ph<=6.5)
	{
		fphlow = 15.444 - 6.1291*ph + 0.8204*ph*ph - 0.0371*ph*ph*ph;
		fphup = 331.68*exp(-1.2618*ph);
	}
	if(t>=80 && t<=90 && ph>=3.5 && ph<=4.57)
	{
		fphlow = 2.6727 - 0.3636*ph;
		fphup = 3.1355 - 0.4673*ph;
	}
	if(t>=80 && t<=90 && ph>=4.57 && ph<=4.6)
	{
		fphlow = 2.6727 - 0.3636*ph;
		fphup = 21254*exp(-2.1811*ph);
	}
	if(t>=80 && t<=90 && ph>=4.6 && ph<=5.62)
	{
		fphlow = 331.68*exp(-1.2618*ph);
		fphup = 21254*exp(-2.1811*ph);
	}
	if(t>=80 && t<=90 && ph>=5.62 && ph<=6.5)
	{
		fphlow = 331.68*exp(-1.2618*ph);
		fphup = 0.4014 - 0.0538*ph;
	}
	if(t>=90 && t<=120 && ph>=3.5 && ph<=4.3)
	{
		fphlow = 3.1355 - 0.4673*ph;
		fphup = 1.5375 - 0.125*ph;
	}
	if(t>=90 && t<=120 && ph>=4.3 && ph<=4.57)
	{
		fphlow = 3.1355 - 0.4673*ph;
		fphup = 5.9757 - 1.157*ph;
	}
	if(t>=90 && t<=120 && ph>=4.57 && ph<=5.0)
	{
		fphlow = 21254*exp(-2.1811*ph);
		fphup = 5.9757 - 1.157*ph;
	}
	if(t>=90 && t<=120 && ph>=5.0 && ph<=5.62)
	{
		fphlow = 21254*exp(-2.1811*ph);
		fphup = 0.546125 - 0.071225*ph;
	}
	if(t>=90 && t<=120 && ph>=5.62 && ph<=6.5)
	{
		fphlow = 0.4014 - 0.0538*ph;
		fphup = 0.546125 - 0.071225*ph;
	}

	fph[0] = fphlow;
	fph[1] = fphup;

	return;
}

void
MoskitoAux::crtf(const Real & t, const Real tt[], const Real kt[], const Real fph[], const Real & fco2, const Real & s, Real & crt)
{
  double craux[2], kt_final, fph_final;
	int i;

	kt_final = kt[0] + (kt[1] - kt[0]) / (tt[1] - tt[0]) * (t - tt[0]);
	fph_final = fph[0] + (fph[1] - fph[0]) / (tt[1] - tt[0]) * (t - tt[0]);
	crt = kt_final * pow(fco2, 0.62) * pow(s/19.0, 0.146 + 0.0324 * log10(fco2)) * fph_final;

	return;
}

void
MoskitoAux::gep(const Real & t, const Real & p, const Real & u, const Real & ph, Real & cr7)
{

    double d[4];

    d[0] = t;
    d[1] = p;
    d[2] = u;
    d[3] = ph;

    const double G1C7 = -6.37490076288078;
    const double G1C9 = 9.71954790027811;
    const double G2C6 = -22.4713541949907;
    const double G2C8 = 112.371930589555;
    const double G2C1 = 16.9811572846352;
    const double G2C3 = 7.47999511703848;
    const double G3C0 = 12.3071233835765;
    const double G3C6 = -58.1855518566808;
    const double G3C7 = -4.79938154426978;
    const double G3C5 = 2.20796920024472;
    const double G4C2 = 7.59265536537209;
    const double G4C1 = 10.6737329569532;
    const double G5C6 = -5.9717883032398;
    const double G5C8 = 4.35994554792827;
    const double G5C9 = 8.28125857955026;
    const double G6C1 = 3.41839318569503;
    const double G6C6 = -3.72394969086718;
    const double G6C3 = 7.62807399214605;

    double y = 0.0;

    y = ((d[0]+((G1C9+G1C9)*G1C7))-d[1]);
    y *= (((d[1]-G2C8)-(G2C8-d[0]))-((G2C1*G2C3)-G2C6));
    y *= (((G3C0+G3C0)+(d[1]+G3C6))*((G3C7-G3C5)+d[3]));
    y *= (((d[1]-d[3])-(d[2]+d[3]))-((G4C2-G4C1)-d[1]));
    y *= (((G5C6+d[0])-(G5C8*d[1]))-((G5C9-d[2])-d[1]));
    y *= (((d[2]+d[1])+(G6C6*G6C3))-(d[3]*G6C1));

    cr7 = - 0.00000013*y + 508.116;
    cr7 = cr7 / 100;

    if(cr7<=0 || cr7>=35)
    {
      cr7 = 0;
    }

}
