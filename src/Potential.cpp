#include "Potential.h"

void Potential::countG(const linear_math::Vector &f, const Ion &ion, const Mesh &mesh)
{
}

void Potential::countH(const Ion &ion, const Mesh &mesh)
{
}

Potential::Potential(const linear_math::Vector &f, const Ion &ion, const Mesh &mesh)
{
    countG(f, ion, mesh);
    countH(ion, mesh);
}