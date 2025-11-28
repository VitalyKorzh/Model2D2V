#ifndef __SOURCE_H__
#define __SOURCE_H__

#include <vector>
#include <cmath>
#include "Ion.h"

typedef unsigned uint;
typedef std::vector <double> darray;

class Source
{
protected:
    int type;
    Ion ion;
    darray I;
    double z;
    double r;

    double coeffVolume;

public:
    Source(const Ion &ion, double z, double r, const darray &I);

    virtual double getIPlus(uint it, uint ionType, double z0, double z1, double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand)=0;
    virtual double getIMinus(uint it, uint ionType, double z0, double z1, double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) 
    {
        return getIPlus(it, ionType, z0, z1, r0, r1, VR3, K0, K1, mub0, mub1, Vv3, Vv3InBand);
    }

    uint getType() const { return type; }
    const Ion & getIon() const { return ion; }
    double getR() const { return r; }
    double getZ() const { return z; }

    virtual double getEnergy() const=0;

    virtual ~Source() {}
};

class KSpehere : public Source {
private:
    double K;
public:

    KSpehere(const Ion &ion, double z, double r, const darray &I, double K) : Source(ion, z, r, I), K(K) {
        type = 0;
    }

    double getEnergy() const override { return K; }

    double getIPlus(uint it, uint ionType, double z0, double z1, double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) override;

};

class KMaxwell : public Source {
private:
    double T;
    uint NE;
public:

    KMaxwell(const Ion &ion, double z, double r, const darray &I, double T, uint NE) : Source(ion, z, r, I), T(T), NE(NE) {
        type = 1;
    }

    double getIPlus(uint it, uint ionType, double z0, double z1, double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) override;
    double getEnergy() const override { return T; }
};


class KTheta : public Source
{
protected:
    double K;
    double sin2Theta;
    double theta;
public:

    KTheta(const Ion &ion, double z, double r, const darray &I, double K, double theta) : Source(ion, z, r, I), K(K), theta(theta)
    {
        type=2;
        double sinTheta = sin(theta);
        sin2Theta = sinTheta*sinTheta;
    }

    double getIPlus(uint it, uint ionType, double z0, double z1,double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) override;
    double getEnergy() const override { return K; }
    double getTheta() const { return theta; }

};

class KThetaPositive : public KTheta
{
public:
    KThetaPositive(const Ion &ion, double z, double r, const darray &I, double K, double theta) : KTheta(ion, z, r, I, K, theta) {
        type = 3;
    }

    double getIPlus(uint it, uint ionType, double z0, double z1,double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) override
    {
        return KTheta::getIPlus(it, ionType, z0, z1, r0, r1, VR3, K0, K1, mub0, mub1, Vv3, Vv3InBand) * 2.;
    }

    double getIMinus(uint it, uint ionType, double z0, double z1,double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) override
    {
        return 0.;
    }
};


class KThetaNegative : public KTheta
{
public:
    KThetaNegative(const Ion &ion, double z, double r, const darray &I, double K, double theta) : KTheta(ion, z, r, I, K, theta) {
        type=4;
    }

    double getIMinus(uint it, uint ionType, double z0, double z1,double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) override
    {
        return KTheta::getIPlus(it, ionType, z0, z1, r0, r1, VR3, K0, K1, mub0, mub1, Vv3, Vv3InBand) * 2.;
    }

    double getIPlus(uint it, uint ionType, double z0, double z1,double r0, double r1, double VR3, double K0, double K1, double mub0, double mub1, double Vv3, double Vv3InBand) override
    {
        return 0.;
    }
};

#endif