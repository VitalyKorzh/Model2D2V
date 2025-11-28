#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include <list>
#include <utility>
#include <ostream>
#include <iostream>
typedef unsigned uint;
typedef std::vector <double> darray;
typedef std::vector <darray> vector2D;

struct base
{
    uint index;
    double min;
    double max;

    base() : index(0), min(0.), max(0.) {}
    base(uint index, double min, double max) : 
    index(index), min(min), max(max) {}   
};

struct quad
{
    uint split;
    double ratio;
    bool exp;
    base bottom;
    base top;
    quad(base bottom, base top, uint split, double ratio=1., bool exp=false) : split(split), ratio(ratio), exp(exp), bottom(bottom), top(top) { }
};

class Mesh
{
private:
    double epsilon;
    double refB;
    darray band;
    vector2D cells2D;
    uint number_cells;

    std::vector <uint> number_cells_in_line;
    std::vector <uint> cells_index;
    darray cells;
    //darray cells_volume;

    uint numer_neighborhood;
    std::vector <std::list <std::pair <double, uint>>> neighborhood_up;
    std::vector <std::list <std::pair <double, uint>>> neighborhood_down;

    bool isMeshLimit(double v0, const darray &v) const;
    void add(double v0, double v1, uint split, double ratio, bool exp, darray &v);
    void addNeighborhoodsToList(std::list <std::pair <double, uint>> &list, uint k0, uint k1, double xj_0, double xj_1, uint index);
    void generateNeighborhoods();
    //void countVolume();

    bool isGenerateCells;

public:
    Mesh(double epsilon=0., double refB=1.) : epsilon(epsilon), refB(refB), number_cells(0), numer_neighborhood(0), isGenerateCells(false) {}
    void addBandMesh(double min, double max, uint split, double ratio=1., bool exp=false);
    void addBandMesh(const darray &band);
    void addCells(quad quad);
    void addCells(uint index, const darray &cells_i);
    void addCells(double bandValue, const darray &cells_i) { addCells(findIndex(bandValue), cells_i); } 
    uint findIndex(double bandValue);
    uint findCellIndex(double x, double y);
    const darray & getBandMesh() const { return band; }
    double getBandMax() const { return band.back(); }
    double getBandMin() const { return band.front(); }
    uint getNumberCells() const { return number_cells; }
    const darray &getCells() const { return cells; }
    //double getCellsVolume(uint index) const { return cells_volume[index]; }
    uint getNumberCellsInLine(uint index) const { return number_cells_in_line[index]; }
    uint getCellsIndex(uint index) const { return cells_index[index]; }
    void generateCells();
    const std::list <std::pair <double, uint>> &getNeighborhoodUp(uint index) const { return neighborhood_up[index]; };
    const std::list <std::pair <double, uint>> &getNeighborhoodDown(uint index) const { return neighborhood_down[index]; };
    uint getNumberBand() const { return band.size()-1; }
    void printMesh(const uint precision=10, std::ostream &os=std::cout) const;
    void setEpsilon(double epsilon) { this->epsilon = epsilon; }
    void setRefB(double refB) { this->refB = refB; } 
    bool isGenerateCellReady() const { return isGenerateCells; }

    uint getNumberNeighborhoods() const { return numer_neighborhood; }
    double getRefB() const { return refB; }
    ~Mesh();
};


#endif