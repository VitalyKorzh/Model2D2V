#ifndef __ION_H__
#define __ION_H__


typedef unsigned uint;

struct Ion
{
public:
    double Z;
    double M;
    uint ionType;

    Ion(double Z=1, double M=1, uint ionType=0) : Z(Z), M(M), ionType(ionType) {
    }
};

#endif