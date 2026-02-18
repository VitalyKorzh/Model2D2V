#ifndef __POTENTIAL_H__
#define __POTENTIAL_H__


#include "linear_math/Vector.h"
#include "Ion.h"
#include "Mesh.h"

class Potential
{
private:
    linear_math::Vector h;
    linear_math::Vector g;


    void countG(const linear_math::Vector &f, const Ion &ion, const Mesh &mesh);
    void countH(const Ion &ion, const Mesh &mesh);

public:
    Potential(const linear_math::Vector &f, const Ion &ion, const Mesh &mesh);

    const linear_math::Vector & getH() const { return h; };
    const linear_math::Vector & getG() const { return g; };
};



#endif