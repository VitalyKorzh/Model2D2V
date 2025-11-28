#include "Mesh.h"
#include <iostream>
#include <cmath>

bool Mesh::isMeshLimit(double v0, const darray &v) const
{
    if (v.size() == 0 || v0 - v[v.size()-1] > epsilon)
        return true;
    return false;
}

void Mesh::add(double v0, double v1, uint split, double ratio, bool exp, darray &v)
{
    if (split == 0)
        split = 1;

    v.reserve(v.size()+split);

    if (isMeshLimit(v0, v))
        v.push_back(v0);

    if (split == 1)
        return;

    double deltaV = exp ? (v1 - v0)*(pow(ratio, 1./(split-1))-1.) / (pow(ratio, split/(split-1.))-1.) : 2. / split * (v1 - v0) / (1. + ratio);
    double q = exp ? pow(ratio, 1./(split-1.)) : (ratio - 1.) / (split - 1) * deltaV;

    double current = v0;
    for (uint i = 1; i < split; i++)
    {
        current += deltaV;
        if (isMeshLimit(current, v))
            v.push_back(current); 
        deltaV = exp ? deltaV * q : deltaV + q;
    }

    if (isMeshLimit(v1, v))
        v.push_back(v1);
}

void Mesh::addNeighborhoodsToList(std::list<std::pair<double, uint>> &list, uint k0, uint k1, double xj_0, double xj_1, uint index)
{
    for (uint k = k0; k + 1 < k1; k++) 
    {
        index++;
        const double xk_0 = cells[k];
        const double xk_1 = cells[k+1];

        const double overloop = std::min(xk_1, xj_1) - std::max(xk_0, xj_0);

        if (xk_1 <= xj_0)
            continue;
        
        if (xk_0 >= xj_1)
            break;
        
        list.emplace_back(overloop, index-1);
    }
}

void Mesh::generateNeighborhoods()
{
    numer_neighborhood = 0;
    neighborhood_up.reserve(number_cells);
    neighborhood_down.reserve(number_cells);
    uint index = 0;
    for (uint i = 0; i < getNumberBand(); i++) 
    {
        for (uint j = cells_index[i]; j+1 < cells_index[i+1]; j++) 
        {
            const double xj_0 = cells[j];
            const double xj_1 = cells[j+1];
            std::list <std::pair <double, uint>> list_up;
            std::list <std::pair <double, uint>> list_down;
            if (i < getNumberBand()-1)
                addNeighborhoodsToList(list_up, cells_index[i+1], cells_index[i+2], xj_0, xj_1, index+number_cells_in_line[i]);
            if (i > 0) 
                addNeighborhoodsToList(list_down, cells_index[i-1], cells_index[i], xj_0, xj_1, index-number_cells_in_line[i-1]);
            numer_neighborhood += list_up.size() + list_down.size() + 2;
            neighborhood_up.push_back(list_up);
            neighborhood_down.push_back(list_down);
        }
        index += number_cells_in_line[i];
    }
}

void Mesh::addBandMesh(double min, double max, uint split, double ratio, bool exp)
{
    add(min, max, split, ratio, exp, band);
}

void Mesh::addBandMesh(const darray &band)
{
    this->band.reserve(this->band.size()+band.size());
    for (double v : band) {
        if (this->band.size() == 0 || std::abs(this->band[this->band.size()-1] - v))
            this->band.push_back(v);
    }
}

void Mesh::addCells(quad q)
{
    isGenerateCells=false;
    cells2D.resize(getNumberBand());
    if (q.top.index >= getNumberBand())
        q.top.index = getNumberBand()-1;
    for (uint i = q.bottom.index; i <= q.top.index; i++) 
    {
        double min = q.bottom.index != q.top.index ? q.bottom.min + (q.top.min - q.bottom.min) * (i - q.bottom.index) / (q.top.index - q.bottom.index) : q.bottom.min;
        double max = q.bottom.index != q.top.index ? q.bottom.max + (q.top.max - q.bottom.max) * (i - q.bottom.index) / (q.top.index - q.bottom.index) : q.bottom.max;
        size_t size_prev = cells2D[i].size();
        add(min, max, q.split, q.ratio, q.exp, cells2D[i]);
        number_cells += cells2D[i].size() - size_prev;
        if (size_prev == 0)
            number_cells--;
    }
}

void Mesh::addCells(uint index, const darray &cells_i)
{
    isGenerateCells=false;
    if (index >= getNumberBand())
        index = getNumberBand()-1;
    cells2D.resize(getNumberBand());
    cells2D[index].reserve(cells_i.size()+cells2D[index].size());
    for (size_t i = 0; i < cells_i.size(); i++) 
    {
        if (isMeshLimit(cells_i[i], cells2D[i])) 
        {
            if (cells2D[index].size() != 0)
                number_cells++;
            cells2D[index].push_back(cells_i[i]);
        }
    }
}

uint Mesh::findIndex(double bandValue)
{
    for (uint i = 0; i < getNumberBand(); i++)
    {
        if (bandValue >= band[i] && bandValue < band[i+1])
            return i;
    }
    return getNumberBand()-1;
}

uint Mesh::findCellIndex(double x, double y)
{
    uint cell_index = 0;

    uint i;
    for (i = 0; i < getNumberBand(); i++)
    {
        if (y >= band[i] && y < band[i+1])
            break;
        cell_index += number_cells_in_line[i];
    }

    for (uint j = cells_index[i]; j+1 < cells_index[i+1]; j++) 
    {
        if (x >= cells[j] && x < cells[j+1])
            break;
        cell_index++;
    }

    return cell_index;
}

void Mesh::generateCells()
{
    if (isGenerateCells)
        return;
    isGenerateCells = true;
    cells.reserve(number_cells+getNumberBand());
    number_cells_in_line.reserve(getNumberBand());
    cells_index.reserve(getNumberBand()+1);
    cells_index.push_back(0);

    uint index = 0;
    for (uint i = 0; i < getNumberBand(); i++) 
    {
        for (size_t j = 0; j < cells2D[i].size(); j++)
        {
            cells.push_back(cells2D[i][j]);
            index++;
        }
        number_cells_in_line.push_back(cells2D[i].size()-1);
        cells_index.push_back(index);
    }

    generateNeighborhoods();

    cells2D.clear();
    cells2D.shrink_to_fit();
}

void Mesh::printMesh(const uint precision, std::ostream &os) const
{
    os << "# Work mesh:\n";
    os << "number_band=" << getNumberBand() << " number_cells=" << number_cells << "\n";
    uint cell_index = 0;
    os.precision(precision);
    os << std::scientific;
    for (uint i = 0; i < getNumberBand(); i++) 
    {
        os << number_cells_in_line[i] << " ";
        os << band[i] << " " << band[i+1] << " ";
        os << cells[cells_index[i]] << " ";
        for (uint j = cells_index[i]+1; j < cells_index[i+1]; j++) 
        {
            os << cells[j] << "0 ";
            cell_index++;
        }
        os << "\n";
    }
    os << "#\n";
}

Mesh::~Mesh()
{
}
