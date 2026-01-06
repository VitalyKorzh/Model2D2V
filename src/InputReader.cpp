#include "InputReader.h"
#include "PhysicValues.h"
#include "Source.h"
#include <cmath>


using namespace linear_math;

std::map <double, std::string>::const_iterator InputReader::findTime_ptr(std::string t_ptr) const
{
    auto it = std::find_if(time.begin(), time.end(), 
        [&t_ptr](const std::pair<double, std::string>& pair) { return pair.second == t_ptr; });
    return it;
}

bool InputReader::shouldEmplace(double t) const
{
    if (time.empty())
        return true;

    auto it = time.lower_bound(t);

    if (it != time.end() && std::abs(it->first - t) <= time_merge)
        return false;
        
    if (it != time.begin() && std::abs((--it)->first - t) <= time_merge) {
        return false;
    }

    return true;
}

void InputReader::errorMessage(std::string error)
{
    error_message = "# " + std::to_string(numberLine) + ": " + error + "\n";
}

uint InputReader::countSpace(std::string line) const
{
    uint space = 0;
    while (line[space] == ' ')
        space++;
    return space;
}

void InputReader::errorConfigConstNumberPar(std::string part1, const std::vector<std::string> PAR_NAMES, const bool *array, const uint N_STEP)
{
    std::string part2 = "";
    for (uint ik = 0; ik < N_STEP; ik++) 
    {
        if (!array[ik]) 
        {
            if (part2 != "")
                part2 += ", ";
            
            part2 = part2 + PAR_NAMES[ik];   
        }
    }
    part2 = part2 + "]";
    errorMessage(part1+part2);
}

void InputReader::createMaxwell(uint ionType, uint ir, double ni, double Ti)
{
    auto ifp = f_plus0.begin();
    auto ifm = f_minus0.begin();
    auto ion = ions.begin();

    for (uint i = 0; i < ionType; i++)
    {
        ifp++;
        ifm++;
        ion++;
    }

    for (uint iz = 0; iz < nz; iz++)
    {
        uint index = (iz*nr+ir)*mesh.getNumberCells();

        uint cell = 0;
        for (uint iE = 0; iE < mesh.getNumberBand(); iE++)
        {
            double E1 = mesh.getBandMesh()[iE];
            double E2 = mesh.getBandMesh()[iE+1];

            for (uint j = mesh.getCellsIndex(iE); j+1 < mesh.getCellsIndex(iE+1); j++)
            {

                double muB1 = mesh.getCells()[j]*Bvac[iz]/mesh.getRefB();
                double muB2 = mesh.getCells()[j+1]*Bvac[iz]/mesh.getRefB();

                double V = 0.;

                if (E1 > muB1)
                    V -= pow(E1 - muB1, 1.5);
                if (E1 > muB2)
                    V += pow(E1 - muB2, 1.5);
                if (E2 > muB1)
                    V += pow(E2 - muB1, 1.5);
                if (E2 > muB2)
                    V -= pow(E2 - muB2, 1.5);

                const uint NE = maxwellIntegratePoints;
                const double dE = (E2-E1)/(NE-1);
                double integral = 0.;
                for (uint k = 0; k < NE-1; k++)
                {
                    double E = E1 + dE*k;

                    double expE1 = exp(-E*ion->Z/Ti);
                    double expE2 = exp(-(E+dE)*ion->Z/Ti);
                    if (E > muB1)
                    {
                        integral += (sqrt(E-muB1)*expE1 + sqrt(E+dE-muB1)*expE2)/2.;
                    }
                    if (E > muB2)
                    {
                        integral -= (sqrt(E-muB2)*expE1 + sqrt(E+dE-muB2)*expE2)/2.;
                    }
                }
                integral *= dE;
                V *= 2./3.;
                if (V != 0)
                {
                    (*ifp)(index+cell) = ni*pow(ion->M/2./M_PI/Ti, 1.5)*integral/V;
                    (*ifm)(index+cell) = ni*pow(ion->M/2./M_PI/Ti, 1.5)*integral/V;
                }
                cell++;
            }
        }
    }
}

bool InputReader::readComponent(std::istream &in)
{
    std::string line;

    uint ionType = 0;
    uint ir = 0;
    double ni = 0;
    double Ti = 0;
    const uint N_PAR = 4;


    bool array[] = {false, false, false, false};

    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line, true);
        array[0] = StringReader::getUnsignedParameter(line, "ion ", ionType) | array[0];
        if (line.find("r") != std::string::npos)
        {
            skip(in, line);
            double r;
            if (!(array[1] = readPosition(ir, r, line, "r", rArray, nr)))
                return false;
        }
        array[2] = StringReader::getDoubleParameter(line, "T ", Ti) | array[2];
        array[3] = StringReader::getDoubleParameter(line, "n ", ni) | array[3];
    }

    if (checkArray(array, N_PAR))
    {
        if (ir >= nr)
        {
            errorMessage("не правильная точка по радиусу [ir < nr]");
            return false;
        }

        if (ionType >= ions.size())
        {
            errorMessage("не существующий тип иона [ion < ionN]");
            return false;
        }

        createMaxwell(ionType, ir, ni, Ti);

    }
    else {
        errorConfigConstNumberPar(
            "не указаны все параметры component [",
            {"ion", "ir", "n", "T"},
            array,
            N_PAR
        );
        return false;
    }

    if (in.fail())
    {
        errorMessage("не удалось прочитать component");
        return false;
    }

    return true;
}

double InputReader::lin(double t1, double v1, double t2, double v2, double t) const
{
    double a = (v2 - v1) / (t2 - t1);
    double b = v1 - a * t1;
    return a * t + b;
}

double InputReader::log(double t1, double v1, double t2, double v2, double t) const
{
    if (t1 <= 0. || t2 <= 0. || t <= 0.)
        return lin(t1, v1, t2, v2, t);

    double a = (v2 - v1) / (std::log(t2) - std::log(t1));
    double b = v1 - a * std::log(t1);
    return a * std::log(t) + b;
}

bool InputReader::readParameterFromTime(std::istream &in, std::string line, std::string name, darray &array, std::list<std::pair<std::string, double>> &array_list)
{
    array.resize(time.size(), 0);
    if (line.find(name + "=") != std::string::npos)
    {
        double val = 0;
        if (StringReader::getDoubleParameter(line, name+"=", val))
        {
            for (uint it = 0; it < time.size(); it++)
                array[it] = val;
            array_list.emplace_back("", 0); // без аппроксимации
            array_list.emplace_back("", val);
        }
        else
        {
            errorMessage("не удалось прочитать " + name);
            return false;
        }
    }
    else if (line.find(name +" ") != std::string::npos)
    {
        uint N = 0;

        uint approxType = 0;
        if (line.find("lin") != std::string::npos)
            approxType = 1;
        else if (line.find("log") != std::string::npos)
            approxType = 2;
        
        if (approxType == 1)
            name += " lin";
        if (approxType == 2)
            name += " log";

        if (approxType == 1)
            array_list.emplace_back("lin ", 0);
        else if (approxType == 2)
            array_list.emplace_back("log ", 0);
        else
            array_list.emplace_back("", 0);
        
        
        if (StringReader::getUnsignedParameter(line, name+" ", N))
        {
            if (N == 0 || N > time.size())
            {
                errorMessage("не указано правильное число точек по времени");
                return false;
            }

            if (timeArray.empty())
                timeArray = getTime();

            uint index_prev = 0;
            for (uint it = 0; it < N; it++)
            {
                skip(in, line);
                std::istringstream iss(line);
                std::string time_ref;
                double val;
                iss >> time_ref >> val;
                if (iss.fail())
                {
                    errorMessage("не удалось прочитать " + name);
                    return false;
                }
                auto i = findTime_ptr(time_ref); 
                if (i == time.end())
                {
                    errorMessage("не удалось найти time pointer " + time_ref);
                    return false;
                }

                double t_it = i->first;
                uint index = timeArray.size();
                for (uint k = 0; k < timeArray.size(); k++)
                {
                    if (t_it == timeArray[k])
                    {
                        index = k;
                        break;
                    }
                }

                if (index >= timeArray.size())
                {
                    errorMessage("не найден точка по времени\n");
                    return false;
                }

                if (index <= index_prev && index_prev != 0)
                {
                    errorMessage("точки по времени указаны не по возрастанию");
                    return false;
                }

                if (approxType == 1)
                {
                    if (index == 0)
                        array[0] = val;
                    else 
                    {
                        for (uint l = index_prev; l <= index; l++)
                            array[l] = lin(timeArray[index_prev], array[index_prev], timeArray[index], val, timeArray[l]);
                    }
                }
                else if (approxType == 2)
                {
                    if (index == 0)
                        array[0] = val;
                        else 
                    {
                        for (uint l = index_prev; l <= index; l++)
                        array[l] = log(timeArray[index_prev], array[index_prev], timeArray[index], val, timeArray[l]);
                    }
                }

                for (uint l = index; l < timeArray.size(); l++)
                    array[l] = val;
                
                array_list.emplace_back(time_ref, val);

                index_prev = index;
            }
            
        }
        else
        {
            errorMessage("не удалось прочитать " + name);
            return false;   
        }
    }


    //for (uint i = 0; i < timeArray.size(); i++)
    //    std::cout << i << "\t" << time.find(timeArray[i])->second << "\t" << timeArray[i] << "\t" << array[i] << "\n";

    return true;
}

bool InputReader::readIon(std::istream &in)
{
    std::string line;

    double Z = 1.;
    double M = 1.;
    const uint N_PAR = 2;

    bool array[] = {false, false};

    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line, true);
        array[0] = StringReader::getDoubleParameter(line, "Z ", Z) | array[0];
        array[1] = StringReader::getDoubleParameter(line, "M ", M) | array[1];
    }

    if (!checkArray(array, N_PAR))
    {
        errorConfigConstNumberPar("не указаны все параметры иона [",
            {"Z=", "M="},
            array,
            N_PAR
        );
        return false;
    }
    
    if (M <= 0.)
    {
        errorMessage("не правильная масса иона [M > 0]");
        return false;
    }

    if (in.fail())
    {
        errorMessage("не удалось прочитать ion");
        return false;
    }
    
    ions.push_back(Ion(Z, M, ions.size()));
    return true;
}

bool InputReader::readInitial(std::istream &in)
{
    std::string line;

    if (!getline(in, line, true))
        return false;


    for (uint i = 0; i < ions.size(); i++)
    {
        f_plus0.push_back(darray(nr*nz*mesh.getNumberCells(), 0.));
        f_minus0.push_back(darray(nr*nz*mesh.getNumberCells(), 0.));
    }


    while (line.find("initial-plasma end") == std::string::npos)
    {
        if (!line.empty() && !isComment(line))
        {
            if (line.find("component") != std::string::npos)
            {
                if (!readComponent(in))
                    return false;
            }
        }

        skip(in, line);
        if (in.fail())
        {
            errorMessage("не найдено initial-plasma end");
            return false;
        }
    }
    

    return true;
}

bool InputReader::readPosition(uint &index, double &val, const std::string &line, std::string name, const darray &array, const uint n)
{
    if (StringReader::getDoubleParameter(line, "value ", val))
    {
        for (uint k = 0; k < n; k++)
        {
            if (array[k] <= val  && array[k+1] > val)
            {
                index = k;
                return true;
            }   
        }

        errorMessage(name + " вышел за границу");
        return false;
    }
    else if (StringReader::getUnsignedParameter(line, "index ", index))
    {
        if (index >= n)
        {
            errorMessage(name + " вышел за границу");
            return false;
        }

        val = (array[index]+array[index+1])/2.;

    }
    else {
        errorMessage("не удалось прочитать " + name);
        return false;
    }
    
    return true;
}

bool InputReader::readBandIndex(std::istream &in, uint &index)
{
    std::string line;
    skip(in, line);

    if (line.find("band") != std::string::npos) 
    {
        skip(in, line, true);
        if (line.find("E ") != std::string::npos) 
        {
            double value;
            StringReader::getDoubleParameter(line, "E ", value);
            index = mesh.findIndex(value);
        }
        else if (!StringReader::getUnsignedParameter(line, "index ", index))
        {
            errorMessage("не указан index или E");
            return false;
        }
    }
    else 
    {
        errorMessage("не задан band");
        return false;
    }

    if (index >= mesh.getNumberBand()) 
    {
        errorMessage("индекс вышел за границу band");
        return false;
    }
    return true;
}

void InputReader::skip(std::istream &in, std::string &line, bool ignoreEqual)
{
    getline(in, line, true, ignoreEqual);
    while ((line == " " || line == "" || isComment(line)) && !in.fail())
    {
        getline(in, line, true, ignoreEqual);
    }

    if (in.fail())
        line = "";
}

bool InputReader::addSequence(std::istream &in)
{   
    std::string line;             
    double min = 0;
    double max = 0;
    uint split = 0;
    double ratio = 1.;
    bool exp = false;
    const uint N_PAR = 3;

    bool array[] = {false, false, false};

    skip(in, line, true);
    
    while (line.find("E-axis end") == std::string::npos)
    {
        array[0] = StringReader::getDoubleParameter(line, "Emin ", min) | array[0];
        array[1] = StringReader::getDoubleParameter(line, "Emax ", max) | array[1];
        array[2] = StringReader::getUnsignedParameter(line, "split ", split) | array[2];
        StringReader::getDoubleParameter(line, "ratio ", ratio);
        readKey(line, "exp", exp);
        skip(in, line, true);
        if (in.fail())
        {
            errorMessage("не найдено закрытие E-axis end");
            return false;
        }
    }

    if (checkArray(array, N_PAR))
        mesh.addBandMesh(min, max, split, ratio, exp);
    else {
        errorConfigConstNumberPar(
            "не указаны все параметры разбиения sequence [",
            {"Emin", "Emax", "split"},
            array,
            N_PAR
        );
        return false;
    }

    if (min >= max)
    {
        errorMessage("не правильные граница [Emin < Emax]");
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать sequence");
        return false;
    }

    return true;
}

bool InputReader::addBandArray(std::istream &in)
{
    std::string line;
    skip(in, line);
    uint number = 0;
    if (!StringReader::getUnsignedParameter(line, "n=", number)) 
    {
        errorMessage("не указано количество границ [n=]");
        return false;
    }
    std::vector <double> band;
    band.reserve(number);
    
    for (uint i = 0; i < number; i++)
    {
        double value;
        in >> value;
        band.push_back(value);
    }
    
    
    if (!in.fail())
    mesh.addBandMesh(band);
    else {
        errorMessage("ошибки при чтения разбиения band через array");
        return false;
    }
    
    while (line.find("E-axis end") == std::string::npos || isComment(line))
    {
        if (!getline(in, line, true)) 
        {
            errorMessage("не найдено закрытие E-axis end");
            return false;
        }
    }
    
    return true;
}

bool InputReader::addCellsArray(std::istream &in, uint index)
{
    std::string line;
    skip(in, line);

    uint number = 0;
    if (!StringReader::getUnsignedParameter(line, "n=", number)) 
    {
        errorMessage("не указано число интервалов [n=]");
        return false;
    }

    darray cells;
    cells.reserve(number);

    for (uint i = 0; i < number; i++)
    {
        double value;
        in >> value;
        cells.push_back(value);
    }

    if (!in.fail())
    mesh.addCells(index, cells);
    else{
        errorMessage("ошибки при чтения разбиения cells через array");
        return false;
    }

    while (line.find("Mu-axis end") == std::string::npos || isComment(line))
    {
        if (!getline(in, line, true)) 
        {
            errorMessage("не найдено закрытие Mu-axis end");
            return false;
        }
    }    


    return true;
}

bool InputReader::readBase(std::istream &in, base &base, const std::string nameBase)
{
    uint index;
    if (!readBandIndex(in, index))
        return false;

    std::string line;

    double min = 0;
    double max = 0;
    const uint N_PAR = 2;

    bool array[] = {false, false};

    for (uint i = 0; i < N_PAR; i++) {
        skip(in, line, true);
        array[0] = StringReader::getDoubleParameter(line, "min ", min) | array[0];
        array[1] = StringReader::getDoubleParameter(line, "max ", max) | array[1];
    }

    if (checkArray(array, N_PAR))
    {
        base.index = index;
        base.min = min;
        base.max = max;
    }
    else 
    {
        errorConfigConstNumberPar(
            "не указаны все параметры разбиения "+nameBase + " [",
            {"min", "max"},
            array,
            N_PAR
        );
        return false;
    }

    if (min < 0.)
    {
        errorMessage("не правильные границы [min >= 0]");
        return false;
    }

    if (min >= max)
    {
        errorMessage("не правильные граница [min < max]");
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать base");
        return false;
    }

    return true;
}

bool InputReader::addCellsQuad(std::istream &in)
{
    std::string line;
    skip(in, line, true);

    const uint N_PAR = 3;
    bool array[] = {false, false, false};

    uint split = 0;
    double ratio = 1.;
    bool exp = false;

    base bottom;
    base top;

    while (line.find("Mu-axis end") == std::string::npos)
    {

        StringReader::getDoubleParameter(line, "ratio ", ratio);
        array[0] = StringReader::getUnsignedParameter(line, "split ", split) | array[0];
        readKey(line, "exp", exp);
        if (line.find("bottom") != std::string::npos) 
        {
            array[1] = true;
            if (!readBase(in, bottom, "bottom"))
                return false;
        }
        else if(line.find("top") != std::string::npos) 
        {
            array[2] = true;
            if (!readBase(in, top, "top"))
                return false;
        }


        skip(in, line, true);
        if (in.fail())
        {
            errorMessage("не найдено закрытие Mu-axis end");
            return false;
        }
    }

    if (checkArray(array, N_PAR))
        mesh.addCells(quad(bottom, top, split, ratio, exp));
    else 
    {
        errorConfigConstNumberPar(
            "не указаны все параметры разбиения quad [",
            {"split", "bottom", "top"},
            array,
            N_PAR
        );
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать quad");
        return false;
    }

    return true;
}

bool InputReader::readMesh(std::istream &in)
{
    std::string line;
    if (!getline(in, line, true))
        return false;

    double epsilon = 0.;
    double refB = 1.;

    bool findBand = false;
    bool findCells = false;

    while(line.find("mesh end") == std::string::npos)
    {
        if (!line.empty() && !isComment(line))
        {
            if(!findBand && StringReader::getDoubleParameter(line, "epsilon=", epsilon) && epsilon >= 0.) 
            {
                mesh.setEpsilon(epsilon);
            }

            if (!findBand && StringReader::getDoubleParameter(line, "refB=", refB) && refB > 0.) 
            {
                mesh.setRefB(refB);
            }

            if (line.find("E-axis") != std::string::npos && line.find("E-axis end") == std::string::npos) 
            {
                findBand = true;
                skip(in, line);

                if (line.find("sequence") != std::string::npos) 
                {
                    if (!addSequence(in))
                        return false;
                }
                else if (line.find("array") != std::string::npos) 
                {
                    if (!addBandArray(in))
                        return false;
                }
                else {
                    errorMessage("не известный тип разбиения y-axis [sequence, array]");
                    return false;
                }
            }
            else if (findBand && line.find("Mu-axis") != std::string::npos && line.find("Mu-axis end") == std::string::npos)
            {
                findCells = true;
                skip(in, line);
                if (line.find("array") != std::string::npos) 
                {
                    uint index;
                    if (!readBandIndex(in, index) || !addCellsArray(in, index))
                        return false;
                }
                else if (line.find("quad") != std::string::npos) 
                {
                    if (!addCellsQuad(in))
                        return false;
                }
                else 
                {
                    errorMessage("не известный тип разбиения Mu-axis [quad, array]");
                    return false;   
                }
            }
            else if (line.find("z-axis") != std::string::npos) 
            {
                skip(in, line, true);
                if (!StringReader::getUnsignedParameter(line, "n ", nz)) 
                {
                    errorMessage("не указано число разбиений по z");
                    return false;
                }
                if (nz == 0) 
                {
                    errorMessage("число разбиений должно n[>=1]");
                    return false;
                }
                zArray.clear();
                zArray.reserve(nz+1);
                for (uint i = 0; i < nz+1; i++) 
                {
                    double val;
                    in >> val;
                    if (i != 0 && val <= zArray.back()) 
                    {
                        errorMessage("сетка по z задается по возрастанию");
                        return false;
                    }
                    zArray.push_back(val);
                }
                if (in.fail())
                {
                    errorMessage("не удалось прочитать разбиение по z");
                    return false;
                }
            }
            else if (line.find("Bvac") != std::string::npos && nz > 0)
            {
                Bvac.clear();
                Bvac.resize(nz);
                for (uint i = 0; i < nz; i++)
                {
                    double val;
                    in >> val;
                    if (val <= 0)
                    {
                        errorMessage("не правильное значение магнитного поля [Bvac > 0]");
                        return false;
                    }
                    Bvac[i] = val;
                }
                
                if (in.fail())
                {
                    errorMessage("не удалось прочитать Bvac");
                    return false;
                }

            }
            else if (line.find("r-axis") != std::string::npos)
            {
                skip(in, line, true);
                if (!StringReader::getUnsignedParameter(line, "n ", nr)) 
                {
                    errorMessage("не указано число разбиений по r");
                    return false;
                }
                if (nr == 0) 
                {
                    errorMessage("число разбиений должно n[>=1]");
                    return false;
                }
                rArray.clear();
                rArray.reserve(nr+1);
                for (uint i = 0; i < nr+1; i++) 
                {
                    double val;
                    in >> val;
                    if (val < 0)
                    {
                        errorMessage("не правильное значение границы [r>=0]");
                        return false;
                    }

                    if (i != 0 && val <= rArray.back()) 
                    {
                        errorMessage("сетка по r задается по возрастанию");
                        return false;
                    }
                    rArray.push_back(val);
                }
                if (in.fail())
                {
                    errorMessage("не удалось прочитать разбиение по r");
                    return false;
                } 
            }
                
            if(in.fail()) {
                errorMessage("ошибки при чтения данных");
                return false;
            }
        }
        skip(in, line);
        if (in.fail()) 
        {
            errorMessage("не найдено закрытие mesh end");
            return false;
        }
    }

    if (!findBand) 
    {
        errorMessage("mesh не указаны E-axis");
        return false;   
    }

    if (!findCells) 
    {
        errorMessage("mesh не указаны Mu-axis");
        return false; 
    }

    if (zArray.size() == 0)
    {
        errorMessage("z-axis не указан");
        return false;
    }
    if (rArray.size() == 0)
    {
        errorMessage("r-axis не указан");
        return false;
    }
    if (Bvac.size() == 0)
    {
        errorMessage("Bvac не указан");
        return false;
    }
    
    return true;
}

bool InputReader::readTime(std::string &line)
{
    std::string t_ptr = readWord(line.substr(line.find("time ")+5));

    if (line.find("set ") != std::string::npos)
    {
        if (findTime_ptr(t_ptr) != time.end())
        {
            errorMessage(t_ptr + " уже существует");
            return false;
        }

        double t = 0;
        StringReader::getDoubleParameter(line, "set ", t);


        if (shouldEmplace(t))
            time.emplace(t, t_ptr);
    }
    else if (line.find("shift ") != std::string::npos) 
    {
        if (findTime_ptr(t_ptr) != time.end())
        {
            errorMessage(t_ptr + " уже существует");
            return false;
        }

        std::string t_p = readWord(line.substr(line.find("shift ")+6));

        auto it  = findTime_ptr(t_p);
        double t = 0;

        if (it == time.end())
        {
            errorMessage(t_p + " не найден");
            return false;
        }
        else
            t = it->first;

        double dt = 0.;
        StringReader::getDoubleParameter(line, t_p+" ", dt);
        t += dt;
        if (shouldEmplace(t))
            time.emplace(t, t_ptr);
    }
    else if (line.find("even ") != std::string::npos) 
    {
        std::istringstream iss(line.substr(line.find("even ")+5));
        std::string word;

        std::string t_ptr_start;
        std::string t_ptr_end;

        uint n = 0;
        iss >> t_ptr_start >> t_ptr_end >> n;

        if (iss.fail())
        {
            errorMessage("не удалось прочитать min max num");
            return false;
        }

        auto it_start = findTime_ptr(t_ptr_start);
        auto it_end = findTime_ptr(t_ptr_end);

        if (it_start == time.end())
        {
            errorMessage(t_ptr_start + " не найден");
            return false;
        }

        if (it_end == time.end())
        {
            errorMessage(t_ptr_end + " не найден");
            return false;
        }

        double min = it_start->first;
        double max = it_end->first;

        if (min >= max) {
            errorMessage("min >= max при разбиение time even");
            return false;
        }

        double dt = (max - min) / (n+1);
        for (uint i = 0; i < n; i++) 
        {
            double t = min + dt * (i+1);
            std::string t_ptr_full = t_ptr + std::to_string(i+1);
            if (findTime_ptr(t_ptr_full) != time.end()) 
            {
                errorMessage(t_ptr_full + " уже существует");
                return false;
            }

            if (shouldEmplace(t))
                time.emplace(t, t_ptr_full);
        }

    }
    else if (line.find("tau ") != std::string::npos) 
    {
        std::istringstream iss(line.substr(line.find("tau ")+4));
        std::string word;
        double min;
        double dt;
        uint n;
        iss >> min >> dt >> n;
        if (iss.fail())
        {
            errorMessage("не удалось прочитать min dt num");
            return false;
        }
        if (dt <= 0.) {
            errorMessage("min >= max при разбиение time tau");
            return false;
        }

        for (uint i = 0; i < n; i++) 
        {
            double t = min + dt * i;
            std::string t_ptr_full = t_ptr + std::to_string(i+1);
            if (findTime_ptr(t_ptr_full) != time.end()) 
            {
                errorMessage(t_ptr_full + " уже существует");
                return false;
            }

            if (shouldEmplace(t))
                time.emplace(t, t_ptr_full);
        }

    }
    else if (line.find("ratio ") != std::string::npos) 
    {
        if (findTime_ptr(t_ptr) != time.end())
        {
            errorMessage(t_ptr + " уже существует");
            return false;
        }

        std::string pt1 = readWord(line.substr(line.find("ratio ")+6));
        auto it1 = findTime_ptr(pt1);

        if (it1 == time.end()) 
        {
            errorMessage(pt1 + " не найден");
            return false;
        }

        std::string find_name = "ratio " + pt1;
        std::string pt2 = readWord(line.substr(line.find(find_name)+find_name.size()));
        auto it2 = findTime_ptr(pt2);

        if (it2 == time.end()) 
        {
            errorMessage(pt2 + " не найден");
            return false;
        }

        find_name += " " + pt2 + "";
        std::string iss_line = line.substr(line.find(find_name+" ")+find_name.size());
        std::istringstream iss(iss_line);

        double coeff = 1.;

        iss >> coeff;

        if (iss.fail())
        {
            errorMessage("не удалось прочитать coeff");
            return false;
        }

        if (coeff >= 1. || coeff <= 0.)
        {
            errorMessage("coeff от 0 до 1");
            return false;
        }

        double t = coeff * it1->first + (1. - coeff)*it2->first;

        if (shouldEmplace(t))
            time.emplace(t, t_ptr);

    }
    else {
        errorMessage("указан не правильный способ задание time");
        return false;
    }

    return true;
}

bool InputReader::readKSphere(std::istream &in)
{
    std::string line;

    uint ionType;
    double z;
    double r;
    double K;
    darray I;
    std::list <std::pair<std::string, double>> I_list;

    const uint N_PAR = 5;

    bool array[] = {false, false, false, false, false};


    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line);

        array[0] = StringReader::getUnsignedParameter(line, "ion ", ionType) | array[0];
        array[0] = StringReader::getUnsignedParameter(line, "ion=", ionType) | array[0];

        if (line.find("z") != std::string::npos && line.find("iz") == std::string::npos)
        {
            skip(in, line);
            uint iz;
            if (!(array[1] = readPosition(iz, z, line, "z", zArray, nz)))
                return false;
        }
        
        if (line.find("r") != std::string::npos)
        {
            skip(in, line);
            uint ir;
            if (!(array[2] = readPosition(ir, r, line, "r", rArray, nr)))
                return false;
        }

        array[3] = StringReader::getDoubleParameter(line, "K ", K) | array[3];
        array[3] = StringReader::getDoubleParameter(line, "K=", K) | array[3];


        if (line.find("I") != std::string::npos)
        {
            if (!readParameterFromTime(in, line, "I", I, I_list))
                return false;
            array[4] = true;
        }

    }

    if (checkArray(array, N_PAR))
    {
        if (z < zArray.front() || z >= zArray.back())
        {
            errorMessage("z вышел за границу");
            return false;
        }
        if (r < 0 || r >= rArray.back())
        {
            errorMessage("r вышел за границу");
            return false;
        }
        if (K < 0)
        {
            errorMessage("энергия K [>=0]");
            return false;
        }
        if (ionType >= ions.size())
        {
            errorMessage("не найден тип иона");
            return false;
        }
        for (uint i = 0; i < time.size(); i++)
        {
            if (I[i] < 0)
            {
                errorMessage("интенсивность I [>=0]");
                return false;
            }
        }

        auto ion = ions.begin();
        for (uint i = 0; i < ionType; i++)
            ion++;

        for (uint i = 0; i < time.size(); i++)
        {
            I[i] *= PhysicValues::E_A_TO_P_TO_S/normaDensity;
        }

        sources.push_back(new KSphere(*ion, z, r, I, K));
        I_sources_list.push_back(I_list);

    }
    else
    {
        errorConfigConstNumberPar(
            "не указаны все параметры source K-sphere [",
            {"ion", "z", "r", "K", "I"},
            array,
            N_PAR
        );
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать source K-sphere");
        return false;
    }

    return true;
}

bool InputReader::readKTheta(std::istream &in, std::string name)
{
    std::string line;

    uint ionType;
    double z;
    double r;
    double K;
    double theta;
    darray I;
    std::list <std::pair<std::string, double>> I_list;

    const uint N_PAR = 6;

    bool array[] = {false, false, false, false, false, false};


    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line);

        array[0] = StringReader::getUnsignedParameter(line, "ion ", ionType) | array[0];
        array[0] = StringReader::getUnsignedParameter(line, "ion=", ionType) | array[0];

        if (line.find("z") != std::string::npos && line.find("iz") == std::string::npos)
        {
            skip(in, line);
            uint iz;
            if (!(array[1] = readPosition(iz, z, line, "z", zArray, nz)))
                return false;
        }
        
        if (line.find("r") != std::string::npos)
        {
            skip(in, line);
            uint ir;
            if (!(array[2] = readPosition(ir, r, line, "r", rArray, nr)))
                return false;
        }

        array[3] = StringReader::getDoubleParameter(line, "K ", K) | array[3];
        array[3] = StringReader::getDoubleParameter(line, "K=", K) | array[3];


        if (line.find("I") != std::string::npos)
        {
            if (!readParameterFromTime(in, line, "I", I, I_list))
                return false;
            array[4] = true;
        }

        array[5] = StringReader::getDoubleParameter(line, "theta ", theta) | array[5];
        array[5] = StringReader::getDoubleParameter(line, "theta=", theta) | array[5];

    }

    if (checkArray(array, N_PAR))
    {
        if (z < zArray.front() || z >= zArray.back())
        {
            errorMessage("z вышел за границу");
            return false;
        }
        if (r < 0 || r >= rArray.back())
        {
            errorMessage("r вышел за границу");
            return false;
        }
        if (K < 0)
        {
            errorMessage("энергия K [>=0]");
            return false;
        }
        if (ionType >= ions.size())
        {
            errorMessage("не найден тип иона");
            return false;
        }
        if (theta < 0 || theta > 90)
        {
            errorMessage("угол theta вышел за пределы [0<=theta<=90]");
            return false;
        }
        for (uint i = 0; i < time.size(); i++)
        {
            if (I[i] < 0)
            {
                errorMessage("интенсивность I [>=0]");
                return false;
            }
        }

        auto ion = ions.begin();
        for (uint i = 0; i < ionType; i++)
            ion++;

        for (uint i = 0; i < time.size(); i++)
        {
            I[i] *= PhysicValues::E_A_TO_P_TO_S/normaDensity;
        }

        if (name == "")
            sources.push_back(new KTheta(*ion, z, r, I, K, theta*M_PI/180.));
        else if (name == "-positive")
            sources.push_back(new KThetaPositive(*ion, z, r, I, K, theta*M_PI/180.));
        else if (name == "-negative")
            sources.push_back(new KThetaNegative(*ion, z, r, I, K, theta*M_PI/180.));
        else
        {
            errorMessage("не верная наименование источника");
            return false;
        }
        I_sources_list.push_back(I_list);

    }
    else
    {
        errorConfigConstNumberPar(
            "не указаны все параметры source K-theta" + name + " [",
            {"ion", "z", "r", "K", "I", "theta"},
            array,
            N_PAR
        );
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать source K-theta" + name);
        return false;
    }

    return true;
}

bool InputReader::readKMaxwell(std::istream &in)
{
    std::string line;

    uint ionType;
    double z;
    double r;
    double T;
    darray I;
    std::list <std::pair<std::string, double>> I_list;

    const uint N_PAR = 5;

    bool array[] = {false, false, false, false, false};


    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line);

        array[0] = StringReader::getUnsignedParameter(line, "ion ", ionType) | array[0];
        array[0] = StringReader::getUnsignedParameter(line, "ion=", ionType) | array[0];

        if (line.find("z") != std::string::npos && line.find("iz") == std::string::npos)
        {
            skip(in, line);
            uint iz;
            if (!(array[1] = readPosition(iz, z, line, "z", zArray, nz)))
                return false;
        }
        
        if (line.find("r") != std::string::npos)
        {
            skip(in, line);
            uint ir;
            if (!(array[2] = readPosition(ir, r, line, "r", rArray, nr)))
                return false;
        }

        array[3] = StringReader::getDoubleParameter(line, "T ", T) | array[3];
        array[3] = StringReader::getDoubleParameter(line, "T=", T) | array[3];


        if (line.find("I") != std::string::npos)
        {
            if (!readParameterFromTime(in, line, "I", I, I_list))
                return false;
            array[4] = true;
        }

    }

    if (checkArray(array, N_PAR))
    {
        if (z < zArray.front() || z >= zArray.back())
        {
            errorMessage("z вышел за границу");
            return false;
        }
        if (r < 0 || r >= rArray.back())
        {
            errorMessage("r вышел за границу");
            return false;
        }
        if (T <= 0)
        {
            errorMessage("энергия T [>0]");
            return false;
        }
        if (ionType >= ions.size())
        {
            errorMessage("не найден тип иона");
            return false;
        }
        for (uint i = 0; i < time.size(); i++)
        {
            if (I[i] < 0)
            {
                errorMessage("интенсивность I [>=0]");
                return false;
            }
        }

        auto ion = ions.begin();
        for (uint i = 0; i < ionType; i++)
            ion++;

        for (uint i = 0; i < time.size(); i++)
        {
            I[i] *= PhysicValues::E_A_TO_P_TO_S/normaDensity;
        }

        sources.push_back(new KMaxwell(*ion, z, r, I, T, maxwellIntegratePoints));
        I_sources_list.push_back(I_list);

    }
    else
    {
        errorConfigConstNumberPar(
            "не указаны все параметры source K-maxwell [",
            {"ion", "z", "r", "T", "I"},
            array,
            N_PAR
        );
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать source K-maxwell");
        return false;
    }

    return true;
}

bool InputReader::readSource(std::istream &in)
{

    std::string line;

    skip(in, line);

    if (line.find("K-sphere") != std::string::npos)
    {
        if (!readKSphere(in))
            return false;
    }
    else if (line.find("K-theta-positive") != std::string::npos)
    {
        if (!readKTheta(in, "-positive"))
        return false;
    }
    else if (line.find("K-theta-negative") != std::string::npos)
    {
        if (!readKTheta(in, "-negative"))
        return false;
    }
    else if (line.find("K-theta") != std::string::npos)
    {
        if (!readKTheta(in))
            return false;
    }
    else if (line.find("K-maxwell") != std::string::npos)
    {
        if (!readKMaxwell(in))
            return false;
    }
    else 
    {
        errorMessage("тип источника не известен\n");
        return false;
    }


    return true;
}

inline void InputReader::readKey(const std::string &line, std::string keyName, bool &key)
{
    key = (line.find(keyName) != std::string::npos) | key;
}

bool InputReader::verify()
{
    bool isWork = true;
    if (time_merge < 0.)
    {
        std::cerr << "# time-merge >= 0.\n";
        isWork = false;
    }

    if (time.size() == 0)
    {
        std::cerr << "# не указаны точки по времени\n";
        isWork = false;
    }

    return isWork;
}

bool InputReader::readFStep(std::istream &in)
{
    f_step_epsilon.resize(time.size(), 1e-10);
    std::string line;

    if (!getline(in, line, true))
        return false;

    while (line.find("f-step end") == std::string::npos)
    {
        if (!line.empty() && !isComment(line))
        {
            StringReader::getDoubleParameter(line, "zero-epsilon=", zero_epsilon);
            StringReader::getDoubleParameter(line, "zero-epsilon ", zero_epsilon);

            if (line.find("epsilon") !=std::string::npos && line.find("zero-epsilon") == std::string::npos)
            {
                f_step_epsilon_list.clear();
                if (!readParameterFromTime(in, line, "epsilon", f_step_epsilon, f_step_epsilon_list))
                    return false;
            }

            StringReader::getUnsignedParameter(line, "threads=", threads);
            StringReader::getUnsignedParameter(line, "threads ", threads);
            StringReader::getUnsignedParameter(line, "lin-iter=", limit_lin_iterations);
            StringReader::getUnsignedParameter(line, "lin-iter ", limit_lin_iterations);

            if (line.find("ILUT") != std::string::npos && line.find("ILUT end") == std::string::npos)
            {
                if (!readILUT(in))
                    return false;
            }
            else if (line.find("ILUK") != std::string::npos && line.find("ILUK end") == std::string::npos) 
            {
                if (!readILUK(in))
                    return false;
            }
        }


        skip(in, line);
        if (in.fail()) 
        {
            errorMessage("не найдено закрытие f-step end");
            return false;
        }
    }


    if (threads < 1)
    {
        errorMessage("не правильное число потоков threads[>=1]");
        return false;
    }

    for (uint it = 0; it < time.size(); it++)
    {
        if (f_step_epsilon[it] <= 0.)
        {
            errorMessage("не правильное epsilon[>0]");
            return false;
        }
    }

    return true;
    
}

bool InputReader::readILUK(std::istream &in)
{
    std::string line;
    skip(in, line, true);

    uint k_max = 0;
    double coeff = 1.;
    double shift = 0.;
    while (line.find("ILUK end") == std::string::npos)
    {
        StringReader::getUnsignedParameter(line, "k ", k_max);
        StringReader::getDoubleParameter(line, "coeff ", coeff);
        StringReader::getDoubleParameter(line, "shift ", shift);

        skip(in, line, true);
        if (in.fail())
        {
            errorMessage("не найдено закрытие ILUT end");
            return false;
        }
    }

    parameters.type = linear_math::PreconditionerType::ILU_K;
    parameters.k_max = k_max;
    parameters.shift_diagonal = shift;
    parameters.coeff_diagonal = coeff;

    return true;
}

bool InputReader::readILUT(std::istream &in)
{
    std::string line;
    skip(in, line, true);

    uint p_max = 10;
    double tol = 1e-12;
    double coeff = 1.;
    double shift = 0.;
    bool fast = false;
    bool absolute = false;
    bool pre_sorted = false;

    while (line.find("ILUT end") == std::string::npos)
    {
        StringReader::getUnsignedParameter(line, "p ", p_max);
        StringReader::getDoubleParameter(line, "tol ", tol);
        StringReader::getDoubleParameter(line, "coeff ", coeff);
        StringReader::getDoubleParameter(line, "shift ", shift);
        readKey(line, "fast", fast);
        readKey(line, "abs", absolute);
        readKey(line, "queue", pre_sorted);

        skip(in, line, true);
        if (in.fail())
        {
            errorMessage("не найдено закрытие ILUT end");
            return false;
        }
    }

    parameters.type = linear_math::PreconditionerType::ILUT;
    parameters.p_max = p_max;
    parameters.threshold = tol;
    parameters.fastType = fast;
    parameters.isAbsoluteThreshold = absolute;
    parameters.fastTypeQueue = pre_sorted;
    parameters.shift_diagonal = shift;
    parameters.coeff_diagonal = coeff;

    return true;
}

InputReader::InputReader(std::istream &in)
{
    precision = 10;
    time_merge = 0.;
    normaDensity=1.;
    normaEnergy=1.;
    normaMagneticField=1.;
    maxwellIntegratePoints = 50;
    std::string line;
    work = true;
    error_message = "";
    numberLine = 0;
    
    {
        threads = 1;
        limit_lin_iterations=200;
        zero_epsilon = 0.;
        parameters.coeff_diagonal = 1.;
        parameters.shift_diagonal = 0.;
        parameters.type = PreconditionerType::NONE;
    }

    bool findMesh = false;
    bool findInitial = false;
    bool findSource = false;
    bool findFStep = false;
    while (getline(in, line, true))
    {
        if (line.empty() || isComment(line))
            continue;

        StringReader::getUnsignedParameter(line, "precision=", precision);
        StringReader::getDoubleParameter(line, "time-merge=", time_merge);
        if (time_merge < 0.)
            time_merge = 0.;

        StringReader::getUnsignedParameter(line, "maxwell-integrate-points", maxwellIntegratePoints);
        if (maxwellIntegratePoints == 0)
            maxwellIntegratePoints = 1;

        if (sources.empty())
        {
            StringReader::getDoubleParameter(line, "normaN=", normaDensity);
            StringReader::getDoubleParameter(line, "normaB=", normaMagneticField);
            StringReader::getDoubleParameter(line, "normaE=", normaEnergy);
            if (normaDensity <= 0)
                normaDensity = 1.;
            if (normaMagneticField <= 0)
                normaMagneticField = 1.;
            if (normaEnergy <= 0)
                normaEnergy = 1.;
        }

        if (!findInitial && !findSource && line.find("ion") != std::string::npos && line.find("ion=") == std::string::npos)
        {
            work = readIon(in);
            if (!work)
                return;
        }
        else if (!findMesh && line.find("mesh") != std::string::npos && line.find("mesh end") == std::string::npos) 
        {
            findMesh = true;
            work = readMesh(in);
            if (!work)
                return;
            mesh.generateCells();
        }
        else if (findMesh && !findInitial && !ions.empty() && line.find("initial-plasma") != std::string::npos && line.find("initial-plasma end") == std::string::npos)
        {
            findInitial = true;
            work = readInitial(in);
            if (!work)
                return;
        }
        else if (!findSource && !findFStep && line.find("time ") != std::string::npos) 
        {
            work = readTime(line);
            if (!work)
                return;
        }
        else if (findMesh && !ions.empty() && !time.empty() && line.find("source") != std::string::npos)
        {
            if (!findSource) {
                findSource = true;
            }
            work = readSource(in);
            if (!work)
                return;
        }
        else if (!time.empty() && line.find("f-step") != std::string::npos && line.find("f-step end") == std::string::npos)
        {
            findFStep = true;
            work = readFStep(in);
            if (!work)
                return;
        }
    }

    if (time.empty())
    {
        work = false;
        errorMessage("не указано разбиение time");
    }

    if (ions.size() == 0)
    {
        work = false;
        errorMessage("не указаны типы ионов");
    }

    if (!findMesh) 
    {
        work = false;
        errorMessage("не указан mesh");
    }

    if (!findInitial)
    {
        work = false;
        errorMessage("не указана initial-plasma");
    }

    if (!findSource)
    {
    }

    if (!findFStep)
        f_step_epsilon.resize(time.size(), 1e-10);

    if (!verify())
    {
        work = false;
        return;
    }
    
}

InputReader::~InputReader()
{
    for (auto source : sources)
        delete source;
}
