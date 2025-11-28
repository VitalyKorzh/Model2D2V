#ifndef __SOLVER_H__
#define __SOLVER_H__


#include <vector>
#include <ostream>
#include <fstream>
#include <iostream>
#include <cmath>
#include "Mesh.h"
#include "InputReader.h"

typedef std::vector<double> darray;
typedef std::vector <short> sharray;
typedef std::vector <uint> uiarray;
typedef unsigned uint;

class MatrixCSRtimer;

class Solver
{
private:

    std::ostream &os;

    const InputReader reader;
    const std::list <Ion> &ions;
    const Mesh &mesh;

    darray time;
    const darray &z;
    uint nz;
    const darray &r;
    uint nr;
    const darray &Bvac;
    darray Bz;
    darray phi;
    double Bmin;

    double normaN;
    double normaB;
    double normaE;

    sharray cellStatus;

    uiarray numberCellsInBand;
    uiarray numberCells;

    darray volume3R;
    darray volume3V;
    darray volume3VinBand;

    std::list<linear_math::Vector> f_plus;
    std::list<linear_math::Vector> f_minus;

    std::list <linear_math::Vector> na;

    double A1;

    void printTime() const;

    double delta_z_i(uint iz) { return z[iz+1] - z[iz]; }
    double z_i(uint iz) { return (z[iz]+z[iz+1])/2.; }
    double delta_r_i(uint ir) { return r[ir+1]-r[ir]; }
    double r_i(uint ir) { return (r[ir+1] + r[ir])/2; }

    void countVolume3R();
    void countVolume3V();
    void countNa();

    double getVolume3R(uint iz, uint ir) const { return volume3R[iz*nr+ir]; }
    double getVolume3V(uint iz, uint ir, uint cell, const Ion &ion) { return volume3V[(iz*nr+ir)*mesh.getNumberCells()+cell]*pow(ion.Z/ion.M, 1.5); }

    void printMesh() const;
    void printBPhi() const;
    void printF() const;

    void printParameterFromTime(const std::list <std::pair<std::string, double>> &list, std::string name) const;

    void printSource() const;

    void print(uint it, double t1, double tau, uint iterations) const;

    MatrixCSRtimer createMatrix(double tau) const;

    uint countF(linear_math::Vector &fPlus, linear_math::Vector &fMinus, const linear_math::Vector &ni, const Ion &ion, linear_math::SolverLinear<MatrixCSRtimer> &solver, uint it, double tau) const;

public:
    Solver(std::istream &in=std::cin, std::ostream &os=std::cout);

    const InputReader & getReader()  const { return reader; }
    bool isReadSuccess() const { return reader.isWork(); }

    void printStartInfo() const;

    void count();


    ~Solver();

};

#endif