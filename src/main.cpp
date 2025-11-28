#include <iostream>
#include "Solver.h"
#include "help.h"

int command(int argc, char**argv) 
{
    for (int i = 1; i < argc; i++) 
    {
        std::string par = argv[i];
        if (par == "--help")
        {
            help();
            return 0;
        }
        else if (par == "--help-initial") 
        {
            helpInitial();
            return 0;
        }
        else if (par == "--help-source") 
        {
            helpSource();
            return 0;
        }
        else if (par == "--help-f-step")
        {
            helpFStep();
            return 0;
        }
        else if (par == "--help-mesh")
        {
            helpMesh();
            return 0;
        }
        else if (par == "--help-time")
        {
            helpTime();
            return 0;
        }
        else
        {
            std::cerr << "команда " + par + " не найдена" << "\n";
            std::cerr << "--help\n--help-time\n--help-mesh\n--help-f-step\n--help-initial\n--help-source\n";
            return 1;
        }   
    }

    return -1;
}


int main(int argc, char**argv)
{
    int c = command(argc, argv);
    if (c >= 0)
        return c;

    std::ifstream fin("../test.in");
    std::ofstream fout("../test.txt");

    if (fin.is_open() && fout.is_open())
    {
        Solver solver(fin, fout);
        if (!solver.isReadSuccess()) {
            std::cerr << solver.getReader().getError();
            fin.close();
            fout.close();
            return 1;
        }
        solver.printStartInfo();
        solver.count();
    }
    fin.close();
    fout.close();
}