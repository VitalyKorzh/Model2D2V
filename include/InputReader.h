#ifndef __INPUT_READER_H__
#define __INPUT_READER_H__

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <istream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <map>
#include <list>
#include <unordered_map>
#include "Ion.h"
#include "Mesh.h"
#include "StringReader.h"
#include "Source.h"
#include "linear_math/SolverLinear.h"

class Solver;

class InputReader
{
private:  
    friend class Solver;

    Mesh mesh;

    std::map <double, std::string> time;
    darray timeArray;
    std::list <Ion> ions;

    double time_merge;

    double normaEnergy; // normaEnergy*E = eV
    double normaDensity; // normaDensity*n = cm^-3
    double normaMagneticField; //normaMagneticField * B = Gs

    uint maxwellIntegratePoints;
    
    uint precision;

    bool work;
    uint numberLine;
    
    std::string error_message;

    darray zArray;
    uint nz;
    darray rArray;
    uint nr;
    darray Bvac;

    uint threads;
    uint limit_lin_iterations;
    double zero_epsilon;
    std::list <std::pair<std::string, double>> f_step_epsilon_list;
    darray f_step_epsilon;
    linear_math::ILUParameters parameters;

    std::list <std::list<std::pair<std::string, double>>> I_sources_list;
    std::list<Source*> sources;
    std::list<linear_math::Vector> f_plus0;
    std::list<linear_math::Vector> f_minus0;

    double lin(double t1, double v1, double t2, double v2, double t) const;
    double log(double t1, double v1, double t2, double v2, double) const;

    uint countSpace(std::string line) const;

    void errorConfigConstNumberPar(std::string part1, const std::vector <std::string> PAR_NAMES, const bool *array, const uint N_STEP);

    auto & getline(std::istream &in, std::string &line, bool formatLine=false, bool ignoreEqual=false) 
    {
        numberLine++;
        std::getline(in, line);
        if (formatLine)
            line = StringReader::formatLine(line);
        if (ignoreEqual)
            std::replace(line.begin(), line.end(), '=', ' ');
        return in;
    }

    bool isComment(const std::string &line) const 
    {
        if (line.empty())
            return true;

        uint i = 0;
        while (line[i] == ' ' || line[i] == '\t') 
        {
            i++;
            if (i == line.size())
                return false;
        }
        return line[i] == '#'; 
    }

    std::string readWord(std::string line) 
    {
        std::istringstream iss(line);
        std::string word;
        iss >> word;
        return word;
    }

    void createMaxwell(uint ionType, uint ir, double ni, double Ti);

    bool readComponent(std::istream &in);

    bool readParameterFromTime(std::istream &in, std::string line, std::string name, darray &array, std::list<std::pair<std::string, double>> &array_list);

    bool readIon(std::istream &in);
    bool readInitial(std::istream &in);
    bool readPosition(uint &index, double &val, const std::string &line, std::string name, const darray &array, const uint n);

    bool readBandIndex(std::istream &in, uint &index);
    void skip(std::istream &in, std::string &line, bool ignoreEqual=false);
    bool addSequence(std::istream &in);
    bool addBandArray(std::istream &in);
    bool addCellsArray(std::istream &in, uint index);
    bool readBase(std::istream &in, base &base, const std::string nameBase="bottom");
    bool addCellsQuad(std::istream &in);
    bool readMesh(std::istream &in);
    bool readTime(std::string &line);
    
    bool readILUT(std::istream &in);
    bool readILUK(std::istream &in);
    bool readFStep(std::istream &in);

    bool readKSphere(std::istream &in);
    bool readKTheta(std::istream &in, std::string name="");
    bool readKMaxwell(std::istream &in);
    bool readSource(std::istream &in);

    std::map <double, std::string>::const_iterator findTime_ptr(std::string t_ptr) const;
    bool shouldEmplace(double t) const;
    void errorMessage(std::string error);

    inline void readKey(const std::string &line, std::string keyName, bool &key);

    bool verify();


    bool checkArray(bool *array, const uint N_PAR)
    {

        for (uint i = 0; i < N_PAR; i++)
        {
            if (!array[i])
                return false;
        }

        return true;
    }

public:
    InputReader(std::istream &in=std::cin);

    const Mesh & getMesh() const { return mesh; }
    darray getTime() const {
        darray t;
        t.reserve(time.size()); 
        for (const auto &it : time)
            t.push_back(it.first);
        return t; 
    }

    uint getPrecision() const { return precision; }

    bool isWork() const { return work; }
    std::string getError() const { return error_message; }


    ~InputReader();
};


#endif