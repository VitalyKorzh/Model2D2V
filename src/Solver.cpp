#include "Solver.h"
#include "PhysicValues.h"
#include "TimeProfiler.h"
#include "Potential.h"

#include <cmath>
#include <algorithm>

using namespace linear_math;

class MatrixCSRtimer : public MatrixCSR
{
private:
public:
    MatrixCSRtimer() : MatrixCSR() {}
    MatrixCSRtimer(const darray &values, const std::vector <uint> &col_indicies, const std::vector <uint> &row_ptr, uint h, uint w) : MatrixCSR(values, col_indicies, row_ptr, h, w) {};
    Vector operator*(const Vector &v) const override 
    {
        TimeProfiler t_mult("multiple matrix to vector");
        return MatrixCSR::operator*(v);
    }
    Vector transposeMult(const Vector &v) const override {
        TimeProfiler t_mult("multiple matrix to vector transpose");
        return MatrixCSR::transposeMult(v);
    }
};

void Solver::printMesh() const
{
    os << "# Mesh for draw:\n";

    for (uint iz = 0; iz < nz; iz++)
    {
        for (uint ir = 0; ir < nr; ir++)
        {
            uint cell = 0;
            os << "number_band=" << mesh.getNumberBand() << " number_cells=" << numberCells[iz*nr+ir] << "\n";
            for (uint iE = 0; iE < mesh.getNumberBand(); iE++)
            {
                os << numberCellsInBand[(iz*nr+ir)*mesh.getNumberBand()+iE] << " " << mesh.getBandMesh()[iE] - phi[iz*nr+ir] << " " << mesh.getBandMesh()[iE+1] - phi[iz*nr+ir] << " ";
                for (uint j = mesh.getCellsIndex(iE); j+1 < mesh.getCellsIndex(iE+1); j++)
                {
                    uint index = (iz*nr+ir)*mesh.getNumberCells() + cell;
                    if (cellStatus[index] != 0)
                    {
                        if (j == mesh.getCellsIndex(iE))
                            os << mesh.getCells()[j]/mesh.getRefB()*Bz[iz*nr+ir] << " ";
                        os << mesh.getCells()[j+1]/mesh.getRefB()*Bz[iz*nr+ir] << " " << volume3V[index] << " ";
                    }
                    cell++;
                }
                os << "\n";
            }
        }
    }


}

void Solver::printBPhi() const
{
    os << "# axis fields:\n";
    for (uint iz = 0; iz < nz; iz++)
    {
        for (uint ir = 0; ir < nr; ir++)
        {
            os << Bvac[iz] << " " << Bz[iz*nr+ir] << " " << phi[iz*nr+ir] << " ";
        }
        os << "\n";
    }

}

void Solver::printF() const
{
    os << "# f-count:\n";
    auto iplus = f_plus.begin();
    auto iminus = f_minus.begin();
    uint iion = 0;

    while (iplus != f_plus.end())
    {
        os << "ion " << iion << "\n";
        for (uint iz = 0; iz < nz; iz++)
        {
            for (uint ir = 0; ir < nr; ir++)
            {
                uint cell = 0;
                for (uint iE = 0; iE < mesh.getNumberBand(); iE++)
                {
                    for (uint j = mesh.getCellsIndex(iE); j+1 < mesh.getCellsIndex(iE+1); j++)
                    {
                        uint index = (iz*nr+ir)*mesh.getNumberCells()+cell;
                        if (cellStatus[index] != 0)
                        os << (*iplus)(index) << " " << (*iminus)(index) << " ";
                        cell++;
                    }
                    os << "\n";
                }
            }
        }
        iplus++;
        iminus++;
        iion++;
    }
}

void Solver::printParameterFromTime(const std::list<std::pair<std::string, double>> &list, std::string name) const
{
    os << "# \t" << name;
    
    auto it = list.begin();
    std::string approx_type = it->first;
    it++;
    if (it->first == "")
        os << "=" << it->second << "\n";
    else
    {
        os << " " << approx_type << list.size()-1 << "\n";
        while (it != list.end())
        {
            os << "# \t\t" << it->first << " " << it->second << "\n";
            it++;
        }
        
    }
}

void Solver::printSource() const
{
    auto is = reader.I_sources_list.begin();

    for (auto source : reader.sources)
    {

        os << "# source\n";
        if (source->getType() == 0)
            os << "# \tK-sphere\n";
        else if (source->getType() == 1)
            os << "# \tK-maxwell\n";
        else if (source->getType() == 2)
            os << "# \tK-theta\n";
        else if (source->getType() == 3)
            os << "# \tK-theta-positive\n";
        else if (source->getType() == 4)
            os << "# \tK-theta-negative\n";

        os << "# \t\tion " << source->getIon().ionType << "\n";
        os << "# \t\tr\n";
        os << "# \t\t\tvalue " << source->getR() << "\n";
        os << "# \t\tz\n";
        os << "# \t\t\tvalue " << source->getZ() << "\n";
        if (source->getType() != 1)
            os << "# \t\tK ";
        else
            os << "# \t\tT ";
        os << source->getEnergy() << "\n"; 
        if (source->getType() != 0 && source->getType() != 1) {
            os << "# \t\ttheta " << ((KTheta*) source)->getTheta() << "\n";
        }

        printParameterFromTime(*is, "I");

        is++;
    }
    os << "#\n";
}

void Solver::print(uint it, double t1, double tau, uint iterations) const
{
    os << "# it=" << it << " t=" << t1 << " dt=" << tau << "\n"; 

    printMesh();
    printBPhi();
    printF();

    os << "# iter=" << iterations << "\n";
    os << "#\n";
}

MatrixCSRtimer Solver::createMatrix(double tau) const
{
    const uint n = nr*nz*mesh.getNumberCells();
    uiarray col_indicies;
    uiarray row_ptr(n+1);
    darray values;

    col_indicies.reserve(n);
    values.reserve(n);

    row_ptr[0] = 0;

    for (uint i = 0; i < n; i++)
    {
        values.push_back(1.);
        col_indicies.push_back(i);
        row_ptr[i+1] = row_ptr[i]+1;
    }

    return MatrixCSRtimer(values, col_indicies, row_ptr, n, n);
}

uint Solver::countF(Vector &fPlus, Vector &fMinus, const Vector &ni, const Ion &ion, linear_math::SolverLinear<MatrixCSRtimer> &solver, uint it, double tau) const
{
    const uint fSize = nr*nz*mesh.getNumberCells(); 
    Vector Fplus(fSize, 0.);
    Vector Fminus(fSize, 0.);

    for (uint iz=0; iz < nz; iz++)
    {
        for (uint ir=0; ir < nr; ir++)
        {
            uint cell = 0;
            uint index = (iz*nr+ir)*mesh.getNumberCells();
            for (uint iE = 0; iE < mesh.getNumberBand(); iE++) 
            {
                double K1 = mesh.getBandMesh()[iE] - phi[iz*nr+ir];
                double K2 = mesh.getBandMesh()[iE+1] - phi[iz*nr+ir];
                for (uint j = mesh.getCellsIndex(iE); j+1<mesh.getCellsIndex(iE+1); j++)
                {
                    double mub1 = mesh.getCells()[j] * Bz[iz*nr+ir]/mesh.getRefB();
                    double mub2 = mesh.getCells()[j+1] * Bz[iz*nr+ir]/mesh.getRefB();
                    if (cellStatus[index+cell] != 0)
                    {
                        Fplus(index+cell) = fPlus(index+cell);
                        Fminus(index+cell) = fMinus(index+cell);

                        for (auto source : reader.sources)
                        {
                            Fplus(index+cell) += tau*source->getIPlus(it, ion.ionType, z[iz], z[iz+1], r[ir], r[ir+1], volume3R[nr*iz+ir], K1, K2, mub1, mub2, 
                                volume3V[(ir+iz*nr)*mesh.getNumberCells()+cell], volume3VinBand[(ir+iz*nr)*mesh.getNumberBand()+iE]);
                            Fminus(index+cell) += tau*source->getIMinus(it, ion.ionType, z[iz], z[iz+1], r[ir], r[ir+1], volume3R[nr*iz+ir], K1, K2, mub1, mub2,
                                volume3V[(ir+iz*nr)*mesh.getNumberCells()+cell], volume3VinBand[(ir+iz*nr)*mesh.getNumberBand()+iE]);
                        }

                    }
                    cell++;
                }
            }
        }
    }

    uint iterations = 0;

    fPlus = solver.solveBiCGStable(Fplus, reader.f_step_epsilon[it], fPlus);
    iterations += solver.getIterations();
    fMinus = solver.solveBiCGStable(Fminus, reader.f_step_epsilon[it], fMinus);
    iterations += solver.getIterations();

    return iterations;
}

Solver::Solver(std::istream &in, std::ostream &os) : os(os), reader(in), ions(reader.ions), mesh(reader.mesh), time(reader.getTime()),
                                                     z(reader.zArray), nz(reader.nz), r(reader.rArray), nr(reader.nr), Bvac(reader.Bvac), Bz(nz * nr, 0), phi(nz * nr, 0), normaN(reader.normaDensity),
                                                     normaB(reader.normaMagneticField), normaE(reader.normaEnergy), volume3R(nr * nz, 0), volume3V(mesh.getNumberCells() * nr * nz, 0)
{

    A1 = sqrt(2./PhysicValues::MP*normaE*PhysicValues::EV_TO_ERG);

    //Ereal = Za*E*normaEnergy*EV_TO_ERG
    //Breal = normaB*B
    //mureal = Za*mu*normaEnergy*EV_TO_ERG/normaB
    //phireal = phi/|e|*normaeEnergy*EV_TO_ERG
    //freal = normaN*f*(MP/normaE/EV_TO_ERG)^(3/2)

    Bmin = 1.;
    if (!Bvac.empty())
        Bmin = *std::min_element(Bvac.begin(), Bvac.end());

    for (uint iz = 0; iz < nz; iz++)
    {
        for (uint ir = 0; ir < nr; ir++)
        {
            Bz[iz*nr+ir] = Bvac[iz];
        }
    }

    f_plus = reader.f_plus0;
    f_minus = reader.f_minus0;

    for (uint i = 0; i < ions.size(); i++)
    {
        na.push_back(darray(nr*nz, 0.));
    }

    os.precision(reader.precision);
    os << std::scientific;
}

void Solver::printTime() const {
    os << "#\n# time-points\n";
    os << "#    n=" << time.size() << "\n";
    for (const auto &it : reader.time) 
        os << "#    " << it.second << " " << it.first << "\n";
}

void Solver::countVolume3R()
{
    for (uint iz=0; iz < nz; iz++)
    {
        double dz = delta_z_i(iz);
        for (uint ir=0; ir < nr; ir++)
        {
            double dr = delta_r_i(ir);
            double r  = r_i(ir);
            volume3R[iz*nr+ir] = 2.*M_PI*dz*dr*r;
        }
    }
}

void Solver::countVolume3V()
{
    cellStatus.resize(mesh.getNumberCells()*nz*nr, 0);
    numberCells.resize(nz*nr, 0);
    numberCellsInBand.resize(nz*nr*mesh.getNumberBand(), 0);
    volume3VinBand.resize(nz*nr*mesh.getNumberBand());

    const double volumeCoeff = 4. * M_PI / 3. * sqrt(2.); 

    for (uint iz = 0; iz < nz; iz++)
    {
        for (uint ir = 0; ir < nr; ir++)
        {
            double B_kl = Bz[iz*nr+ir];
            double phi_kl = phi[iz*nr+ir];

            uint cell = 0;
            
            numberCells[iz*nr+ir] = 0;
            
            for (uint iE = 0; iE < mesh.getNumberBand(); iE++)
            {
                numberCellsInBand[(iz*nr+ir)*mesh.getNumberBand()+iE] = 0;
                volume3VinBand[(iz*nr+ir)*mesh.getNumberBand()+iE] = 0;

                double K0_kl = mesh.getBandMesh()[iE] - phi_kl;
                double K1_kl = mesh.getBandMesh()[iE+1] - phi_kl;
                
                for (uint j = mesh.getCellsIndex(iE); j+1 < mesh.getCellsIndex(iE+1); j++)
                {
                    double volume = 0;
                    double muB0_kl = mesh.getCells()[j]/mesh.getRefB()*B_kl;
                    double muB1_kl = mesh.getCells()[j+1]/mesh.getRefB()*B_kl;
                    
                    short type = 0;
                    if (K0_kl > muB0_kl)
                    {
                        volume -= pow(K0_kl-muB0_kl, 1.5);
                        type++;
                    }
                    if (K1_kl > muB0_kl)
                    {
                        volume += pow(K1_kl-muB0_kl, 1.5);
                        type++;
                    }
                    if (K0_kl > muB1_kl)
                    {
                        volume += pow(K0_kl-muB1_kl, 1.5);
                        type++;
                    }
                    if (K1_kl > muB1_kl)
                    {
                        volume -= pow(K1_kl-muB1_kl, 1.5);
                        type++;
                    }
                    
                    volume3V[(ir+iz*nr)*mesh.getNumberCells()+cell] = volumeCoeff*volume;
                    if (type != 0)
                    {
                        numberCellsInBand[(iz*nr+ir)*mesh.getNumberBand()+iE]++;
                        numberCells[iz*nr+ir]++;
                        volume3VinBand[(ir+iz*nr)*mesh.getNumberBand()+iE] += volume3V[(ir+iz*nr)*mesh.getNumberCells()+cell];
                    }

                    cellStatus[(ir+iz*nr)*mesh.getNumberCells()+cell] = type;

                    cell++;
                }
            }

        }
    }

}

void Solver::countNa()
{
    auto if_plus = f_plus.begin();
    auto if_minus = f_minus.begin();
    auto ina = na.begin();

    for (const Ion &ion : ions)
    {
        for (uint iz = 0; iz < nz; iz++)
        {
            for (uint ir = 0; ir < nr; ir++)
            {
                uint index = iz*nr+ir; 
                (*ina)(index) = 0.;
                for (uint i = 0; i < mesh.getNumberCells(); i++)
                {
                    uint index_f = index*mesh.getNumberCells()+i;
                    (*ina)(index) += getVolume3V(iz, ir, i, ion)*((*if_plus)(index_f) + (*if_minus)(index_f));
                }
            }
        }

        if_plus++;
        if_minus++;
        ina++;
    }
}

void Solver::printStartInfo() const
{
    os << "# precision=" << reader.precision << "\n";
    os << "# maxwell-integrate-point=" << reader.maxwellIntegratePoints << "\n";
    os << "# normaN=" << normaN << "\n";
    os << "# normaE=" << normaE << "\n";
    os << "# normaB=" << normaB << "\n#\n";
    
    for (const Ion &ion : ions) {
        os << "# ion " << ion.ionType << "\n";
        os << "# \tZ=" << ion.Z << "\n";
        os << "# \tM=" << ion.M << "\n";
    }
    os << "#\n";

    os << "# number t step=" << time.size() << "\n";
    printTime();
    os << "# z-axis\n# \tn " << nz << "\n";
    for (const double & z0 : z)
        os << "# \t\t" << z0 << "\n";
    os << "# r-axis\n# \tn " << nr << "\n";
    for (const double & r0 : r)
        os << "# \t\t" << r0 << "\n";
    os << "# Bvac\n";
    for (const double & B0 : Bvac)
        os << "# \t" << B0 << "\n";
    os << "#\n";


    printSource();
    
    os << "# f-step\n";

    printParameterFromTime(reader.f_step_epsilon_list, "epsilon");

    os << "# \tzero-epsilon=" <<
    reader.zero_epsilon << "\n#\tthreads=" << reader.threads << "\n#\tlin-iter=" <<
    reader.limit_lin_iterations << "\n";


    if (reader.parameters.type == PreconditionerType::ILUT)
    {
        os << "# \tILUT\n";
        os << "# \t\tp " << reader.parameters.p_max << "\n";
        os << "# \t\ttol " << reader.parameters.threshold << "\n";
        os << "# \t\tcoeff " << reader.parameters.coeff_diagonal << "\n";
        os << "# \t\tshift " << reader.parameters.shift_diagonal << "\n";
        if (reader.parameters.fastType)
            os << "# \t\t[fast]\n";
        if (reader.parameters.isAbsoluteThreshold)
            os << "# \t\t[abs]\n";
        if (reader.parameters.fastTypeQueue)
            os << "# \t\t[queue]\n";
        os << "# \tILUT end\n";
    }
    else if (reader.parameters.type == PreconditionerType::ILU_K)
    {
        os << "# \tILUK\n";
        os << "# \t\tp " << reader.parameters.p_max << "\n";
        os << "# \t\tcoeff " << reader.parameters.coeff_diagonal << "\n";
        os << "# \t\tshift " << reader.parameters.shift_diagonal << "\n";
        os << "# \tILUK end\n";
    }

    os << "# f-step end\n";

    os << "#\n";
}

void Solver::count()
{   
    double t0 = time.front();

    {
        countVolume3R();
        countVolume3V();
        {
            TimeProfiler t_print("time print");
            print(0, 0, 0, 0);   
        }
    }

    SolverLinear <MatrixCSRtimer> solver(reader.limit_lin_iterations, reader.zero_epsilon);

    {
        TimeProfiler t_matrix("create matrix");
        solver.setMatrix(createMatrix(0));
        solver.setNumberThreads(reader.threads);
        solver.setPreconditionerType(reader.parameters);
    }

    uint iterations = 0;
    
    for (uint it = 1; it < time.size(); it++)
    {
        double t1 = time[it];
        double tau = t1 - t0;
        t0 = t1;

        //подсчет объем и концентрации
        {
            countVolume3R();
            countVolume3V();
            countNa();
        }

        //подсчет функции распределение на новом временном слое
        {
            auto ion = ions.begin();
            auto ifp = f_plus.begin();
            auto ifm = f_minus.begin();
            auto ina = na.begin();
            while (ion != ions.end())
            {
                TimeProfiler t_countF("count f");
                iterations = countF(*ifp, *ifm, *ina, *ion, solver, it, tau);

                ion++;
                ifp++;
                ifm++;
                ina++;
            }
        }

        //вывод информации
        {
            TimeProfiler t_print("time print");
            print(it, t1, tau, iterations);
        }
    }
}

Solver::~Solver()
{
    TimeProfiler::print(os);
}
