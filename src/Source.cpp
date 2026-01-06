#include "Source.h"
#include <cmath>

double KSphere::getIPlus(uint it, uint ionType, double z0, double z1, double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand)
{
    double K_grid = K / ion.Z;
    if (ion.ionType == ionType && z >= z0 && z < z1 && r >= r0 && r < r1 && K_grid >= K0 && K_grid < K1)
    {
        return (I[it]+I[it-1])/(2.*VR3*Vv3InBand*coeffVolume*2.);
    }
    else
        return 0.0;
}

Source::Source(const Ion &ion, double z, double r, const darray &I) : ion(ion), I(I), z(z), r(r), coeffVolume(pow(ion.Z/ion.M, 1.5))
{
}

double KTheta::getIPlus(uint it, uint ionType, double z0, double z1, double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand)
{
    double K_grid = K / ion.Z;
    double mu_grid = K_grid * sin2Theta;
    if (ion.ionType == ionType && z >= z0 && z < z1 && r >= r0 && r < r1 && K_grid >= K0 && K_grid < K1 && mu_grid >= mub0 && mu_grid < mub1)
    {
        return (I[it]+I[it-1])/(2.*VR3*Vv3*coeffVolume*2.);
    }
    else
        return 0.0;
}

double KMaxwell::getIPlus(uint it, uint ionType, double z0, double z1, double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand)
{
    if (ion.ionType == ionType && z >= z0 && z < z1 && r >= r0 && r < r1)
    {
        
        double n = (I[it]+I[it-1]) / (2. * VR3);
        double V = 0.;

        if (K0 > mub0)
            V -= pow(K0 - mub0, 1.5);
        if (K0 > mub1)
            V += pow(K0 - mub1, 1.5);
        if (K1 > mub0)
            V += pow(K1 - mub0, 1.5);
        if (K1 > mub1)
            V -= pow(K1 - mub1, 1.5);

        const double dE = (K1-K0)/(NE-1);
        double integral = 0.;
        for (uint k = 0; k < NE-1; k++)
        {
            double E = K0 + dE*k;

            double expE1 = exp(-E*ion.Z/T);
            double expE2 = exp(-(E+dE)*ion.Z/T);
            if (E > mub0)
            {
                integral += (sqrt(E-mub0)*expE1 + sqrt(E+dE-mub0)*expE2)/2.;
            }
            if (E > mub1)
            {
                integral -= (sqrt(E-mub1)*expE1 + sqrt(E+dE-mub1)*expE2)/2.;
            }
        }
        integral *= dE;
        V *= 2./3.;
        return n*pow(ion.M/2./M_PI/T, 1.5)*integral/V;
    }
    else
        return 0.0;
}
