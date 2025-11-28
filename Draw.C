#include <iostream>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TColor.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TPolyLine.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <TBox.h>
#include <TLatex.h>

#define EV_ERG 1.6e-12
#define MP 1.672e-24

typedef unsigned uint;
typedef std::vector <double> darray;


namespace StringReader {

bool findParameterInLine(const std::string &line, const std::string &parameterName)
{
    if (line.empty())
        return false;

    size_t index = line.find(parameterName);
    if (index != std::string::npos && (index == 0 || line[index-1] == ' ')) 
    {
        return true;
    }
    
    return false;
}

bool getDoubleParameter(std::string line, std::string parameterName, double &val)
{
    if (findParameterInLine(line, parameterName)) 
        try {
            val = std::stod(line.substr(line.find(parameterName) + parameterName.size()));
        }catch (const std::invalid_argument &) { 
            val = 0.;                                     
        } catch (const std::out_of_range &) {   
            val = 0.;                                       
        }        
    else
        return false;
    return true;
}

bool getUnsignedParameter(std::string line, std::string parameterName, unsigned &val) 
{
    if (findParameterInLine(line, parameterName)) 
        try {
            val = 0;
            unsigned long int valL = 0;
            valL = std::stoul(line.substr(line.find(parameterName) + parameterName.size()));
            val = static_cast <unsigned> (valL);
        }catch (const std::invalid_argument &) { 
            val = 0;                                     
        } catch (const std::out_of_range &) {   
            val = 0;                                       
        }     
    else
        return false;
    return true;
}

}

TCanvas* makeCanvas(
	const char* const cName)
{
	TCanvas* c;
	TString cTitle(cName);
	TObject* const o = gROOT->FindObject(cName);
	if( o && o->InheritsFrom(TCanvas::Class()) )
	{
		//std::cout << "canvas clear!";
		c = (TCanvas*)o;
		c->Clear();
		c->GetListOfPrimitives()->Delete();
		//c->Close();
		c->SetTitle(cTitle);
	}
	else
		c=new TCanvas(cName,cTitle,1,1,950, 900);
	c->SetBit(kCanDelete);
	c->SetGrid();
	c->cd();
	return c;
}

class DrawMesh {
private:
    uint number_cells;
    uint number_band;
    darray f_plus;
    darray f_minus;
    std::vector <uint> number_cells_in_line;
    std::vector <uint> cells_index;

    std::vector <double> zArray;
    std::vector <double> rArray;
 
    darray volume;
    darray band;
    darray cells;

    double t;
    double tau;

    double M;
    double Z;

    double coeff_to_f;

    uint numberCellsUntil;

    double normaN;
    double normaB;
    double normaE;

    Int_t getColor(double value, double minValue, double maxValue) const
    {
        double f = (value-minValue)/maxValue;
        if (colorMult > 0)
            f = colorMult * (log(f) + colorAdd);
        f = (f < 0) || std::isnan(f) ? 0 : f;
        const TArrayI& colors = TColor::GetPalette();
        int index = int( f * colors.GetSize() + 0.5 );
        if( index >= colors.GetSize() )
            index = colors.GetSize()-1;
        else if( index < 0 )
            index = 0;
        return colors[index];
    }

    bool error;

    double colorMult;
    double colorAdd;
    double axisMin;

    void drawColorBar(double xLo, double xHi, double yLo, double yHi, double min=0, double max=1.) 
    {
        TString colorPadName = "colorBarName";
        delete gROOT->FindObject(colorPadName);
        TPad *pad = new TPad(colorPadName, "", xLo, 0, xHi, 1.);
        pad->SetBit(kCanDelete);
        pad->Draw();
        pad->cd();

        {
			const TArrayI& colors = TColor::GetPalette();
			const Int_t colorn = colors.GetSize();
			Double_t x[4] = { 0, 0.5, 0.5, 0 }, y[4];
			for(Int_t i = 0; i != colorn; ++i )
			{
				y[0] = y[1] = yLo +   i   * ( yHi - yLo ) / colorn;
				y[2] = y[3] = yLo + (i+1) * ( yHi - yLo ) / colorn;
				TPolyLine* const box = new TPolyLine( 4, x, y );
				box->SetFillColor( colors[i] );
				box->SetLineStyle(0);
				box->Draw("f");
				box->SetBit(kCanDelete);
			}
		}


		{
			TGaxis* const a = new TGaxis( 0.5, yLo, 0.5, yHi, axisMin, 1, 50510, colorMult < 0 ? "+L" : "G+L" );
            a->SetBit(kCanDelete);
			a->SetLabelSize(0.25);
			a->SetGridLength(10);
			a->Draw();
			a->SetBit(kCanDelete);
		}

    {
        TLatex* texMin = new TLatex(0.5, 0.08, Form("%.2e", min));
        texMin->SetBit(kCanDelete);
        texMin->SetTextSize(0.25);
        texMin->SetTextAlign(22);  // Центрирование по горизонтали и вертикали
        texMin->Draw();

        TLatex* texMax = new TLatex(0.5, 0.92, Form("%.2e", max));
        texMin->SetBit(kCanDelete);
        texMax->SetTextSize(0.25);
        texMax->SetTextAlign(22);
        texMax->Draw();
    }

        pad->SetEditable(false);
    }

public:

    double countNumberParticle() const
    {
        if (error)
            return 1e100;
        double n = 0;
        for (uint i = 0; i < number_cells; i++) 
        {
            n += (f_plus[i] + f_minus[i])*volume[i];
        }
        return n*normaN;
    }

    double count3VFull() const
    {
        if (error)
            return 0;

        double V = 0;
        for (uint i = 0; i < number_cells; i++)
            V += volume[i];

        return V;
    }

    double count3RV(uint iz, uint ir) const {
        if (error)
            return 0;

        return 2.*M_PI*(zArray[iz+1]-zArray[iz])*(rArray[ir+1]-rArray[ir])*(rArray[ir+1]+rArray[ir]) / 2.;
    }

    static bool readUntil(std::ifstream &fin, const std::string str, std::string *save = nullptr) /*прочесть до заданной строки*/
    {
        std::string line;
        while(std::getline(fin, line)) 
        {
            if (save != nullptr)
                *save = line;
            if(line.size() > str.size())
                line.erase(str.size());
            if (line == str)
                return true;
        }
        return false;
    }

    void readMesh(std::ifstream &fin) 
    {
        cells.reserve(number_cells+number_band);
        number_cells_in_line.reserve(number_band);
        f_plus.resize(number_cells, 0);
        f_minus.resize(number_cells, 0);
        band.reserve(number_band+1);
        volume.reserve(number_cells);
        cells_index.reserve(number_band+1);
        cells_index.push_back(0);

        uint cell_index_0 = 0;
        const double coeff = normaE*1e-3*Z;
        for (uint i = 0; i < number_band; i++) 
        {
            uint index;
            double value;
            fin >> index;
            number_cells_in_line.push_back(index);
            fin >> value;
            band.push_back(value*coeff);
            fin >> value;
            if (i+1 == number_band)
                band.push_back(value*coeff);
            
            fin >> value;
            cells.push_back(value*coeff);

            for (uint j = 0; j < number_cells_in_line.back(); j++) 
            {
                fin >> value;
                cells.push_back(value*coeff);
                fin >> value;
                volume.push_back(value*pow(Z/M, 1.5));
                cell_index_0++;
            }
            cell_index_0++;
            cells_index.push_back(cell_index_0);
        }
    }

    DrawMesh(std::string fileName,  uint it=0, uint iz=0, uint ir=0, uint ionType=0) : colorMult(-1.), colorAdd(0.)
    {
        std::ifstream fin;
        fin.open(fileName);

        if (fin.is_open()) 
        {
            std::string line;


            if (!readUntil(fin, "# normaN=", &line)) {
                std::cerr << "не удалось найти # normaN=\n";
                error = true;
                fin.close();
                return;
            }

            if (!StringReader::getDoubleParameter(line, "normaN=", normaN))
            {
                std::cerr << "не удалось прочитать normaN\n";
                error = true;
                fin.close();
                return;
            }
            
            std::getline(fin, line);

            if (!StringReader::getDoubleParameter(line, "normaE=", normaE))
            {
                std::cerr << "не удалось прочитать normaE\n";
                error = true;
                fin.close();
                return;
            }
            
            coeff_to_f = pow(MP/normaE/EV_ERG, 1.5)*normaN;

            std::getline(fin, line);

            if (!StringReader::getDoubleParameter(line, "normaB=", normaB))
            {
                std::cerr << "не удалось прочитать normaB\n";
                error = true;
                fin.close();
                return;
            }


            if (!readUntil(fin, "# ion " + std::to_string(ionType)))
            {
                std::cerr << "не удалось найти ион " << ionType << "\n";
                error = true;
                fin.close();
                return;
            }

            std::getline(fin, line);

            if (!StringReader::getDoubleParameter(line, "# \tZ=", Z))
            {
                std::cerr << "не удалось найти заряд иона\n";
                error = true;
                fin.close();
                return;
            }
            std::getline(fin, line);

            if (!StringReader::getDoubleParameter(line, "# \tM=", M))
            {
                std::cerr << "не удалось найти массу иона\n";
                error = true;
                fin.close();
                return;
            }


            if (!readUntil(fin, "# time-points")) 
            {
                std::cerr << "# time-points не найдено\n";
                error = true;
                fin.close();
                return;
            }

            getline(fin, line);
            uint nt = 0;
            if(!StringReader::getUnsignedParameter(line, "#    n=", nt)) 
            {
                std::cerr << "не удалось прочитать nt!\n";
                error = true;
                fin.close();
                return;   
            }

            if (it >= nt)
            {
                std::cerr << "it вышел за границу!\n";
                error = true;
                fin.close();
                return;
            }

            if (!readUntil(fin, "# z-axis")) 
            {
                std::cerr << "# z-axis не найдено!\n"; 
                error = true;
                fin.close();
                return;
            }
            
            getline(fin, line);
            uint nz = 0;
            if (!StringReader::getUnsignedParameter(line, "# 	n ", nz))
            {
                std::cerr << "не удалось прочитать nz!\n";
                error = true;
                fin.close();
                return;
            }

            if (iz >= nz)
            {
                std::cerr << "iz вышел за границу!\n";
                error = true;
                fin.close();
                return;
            }
            zArray.reserve(nz+1);

            for (uint i = 0; i < nz+1; i++) 
            {
                char symbol;
                double val;
                fin >> symbol >> val;
                zArray.push_back(val);
            }
            
            if (!readUntil(fin, "# r-axis"))
            {
                std::cerr << "# r-axis не найдено!\n"; 
                error = true;
                fin.close();
                return;
            }
            
            getline(fin, line);
            uint nr = 0;
            
            if (!StringReader::getUnsignedParameter(line, "# 	n ", nr))
            {
                std::cerr << "не удалось прочитать nr!\n";
                error = true;
                fin.close();
                return;
            }


            if (ir >= nr)
            {
                std::cerr << "ir вышел за границу\n";
                error = true;
                fin.close();
                return;
            }

            rArray.reserve(nr+1);

            for (uint i = 0; i < nr+1; i++)
            {
                char symbol;
                double val;
                fin >> symbol >> val;
                rArray.push_back(val);
            }

            if (fin.fail())
            {
                std::cerr << "ошибка чтения axis\n";
                error = true;
                fin.close();
                return;
            }

            std::string timeName = "# it=" + std::to_string(it);

            if (!readUntil(fin, timeName, &line)) {
                std::cerr << timeName << " не найден!\n";
                error = true;
                fin.close();
                return;
            }

            error = false;

            StringReader::getDoubleParameter(line.substr(line.find("it=")+3), "t=", t);
            StringReader::getDoubleParameter(line.substr(line.find("it=")+3), "dt=", tau);

            if (!readUntil(fin, "# Mesh for draw:"))
            {
                std::cerr << "не удалось найти # Mesh for draw:\n";
                error = true;
                fin.close();
                return;
            }

            numberCellsUntil = 0;
            for (uint i = 0; i <= ir+nr*iz; i++)
            {
                if (!readUntil(fin, "number_band=", &line))
                {
                    std::cerr << "не найдена сетка\n";
                    error = true;
                    fin.close();
                    return;
                }

                if (i != ir+nr*iz)
                {
                    uint ncells = 0;
                    if (!StringReader::getUnsignedParameter(line, "number_cells=", ncells)) 
                    {
                        std::cerr << "number cells не найдено!\n";
                        error = true;
                        fin.close();
                        return;
                    }
                    numberCellsUntil += ncells;
                }

            }
            if (!StringReader::getUnsignedParameter(line, "number_band=", number_band)) 
            {
                std::cerr << "number band не найдено!\n";
                error = true;
                fin.close();
                return;
            }

            if (!StringReader::getUnsignedParameter(line, "number_cells=", number_cells)) 
            {
                std::cerr << "number cells не найдено!\n";
                error = true;
                fin.close();
                return;
            }  
            
            std::cout << number_cells << " " << number_band << "\n";

            readMesh(fin);


            if (!readUntil(fin, "# f-count")) {
                std::cerr << "не найда # f-count\n";
                error = true;
                fin.close();
                return;
            }

            if (!readUntil(fin, "ion " +std::to_string(ionType))) {
                std::cerr << "не найден тип иона " << ionType << "\n";
                error = true;
                fin.close();
                return;
            }

            for (uint i = 0; i < numberCellsUntil*2; i++)
            {
                double value;
                fin >> value;
            }

            for (uint i = 0; i < number_cells; i++) {
                fin >> f_plus[i] >> f_minus[i];
            }

            if (fin.fail())
            {
                std::cerr << "не удалось прочиать функцию распределения\n";
                error = true;
                fin.close();
                return;
            }

        }
        else {
            std::cerr << fileName << " не найден!\n";
            error = true;
        }

        fin.close();
    }

    void drawMesh(bool drawGrid=true, bool zeroMinLimit=false, bool coloBar=true, EColorPalette colorMapName=EColorPalette::kRainBow, int drawPlus=0b11) 
    {
        if (error)
            return;
        std::cout << "time: " << t << " tau: " << tau << "\n";
        TColor::SetPalette(colorMapName, 0);
        TString canvasName = "graph_f";
        TCanvas *canvas = makeCanvas(canvasName);
        canvas->SetBit(kCanDelete);

        TString mg_name = "m_graph";
        delete gROOT->FindObject(mg_name);
        TMultiGraph *mg = new TMultiGraph(mg_name, "");

        mg->SetBit(kCanDelete);
        const uint SIZE = 6;
        double x_plus[SIZE];
        double y[SIZE];
        double x_minus[SIZE];
        uint n = 0;
        double max_plus = *std::max_element(f_plus.begin(), f_plus.end());
        double min_plus = *std::min_element(f_plus.begin(), f_plus.end());
        double max_minus = *std::max_element(f_minus.begin(), f_minus.end());
        double min_minus = *std::min_element(f_minus.begin(), f_minus.end());

        double max = std::max(max_minus, max_plus);
        double min = std::min(min_minus, min_plus);

        if (drawPlus >> 0 & 1 && ! drawPlus >> 1 & 1)
        {
            max = max_plus;
            min = min_plus;
        }
        else if (! drawPlus >> 0 & 1 && drawPlus >> 1 & 1) 
        {
            max = max_minus;
            min = min_minus;
        }

        min = zeroMinLimit ? std::min(0., min) : min;
        std::cout << "min: " << min << " max: " << max << "\n";

        for (uint i = 0; i < number_band; i++) 
        {
            for (uint j = cells_index[i]; j+1 < cells_index[i+1]; j++) {

                double x0 = cells[j];
                double x1 = cells[j+1];
                double y0 = band[i];
                double y1 = band[i+1];

                uint size;

                if (x0 <= y0 && x0 <= y1 && x1 <= y0 && x1 <= y1) 
                {
                    size = 5;
                    x_plus[0] = x0, x_plus[1] = x0, x_plus[2] = x1, x_plus[3] = x1, x_plus[4] = x0;
                    y[0] = y0, y[1] = y1, y[2] = y1, y[3] = y0, y[4] = y0;
                }
                else if (x0 > y0 && x0 <= y1 && x1 > y0 && x1 > y1) 
                {
                    size = 4;
                    x_plus[0] = x0, x_plus[1] = x0, x_plus[2] = y1, x_plus[3] = x0;
                    y[0] = x0, y[1] = y1, y[2] = y1, y[3] = x0;
                }
                else if (x0 > y0 && x0 <= y1 && x1 > y0 && x1 <= y1) 
                {
                    size = 5;
                    x_plus[0] = x0, x_plus[1] = x0, x_plus[2] = x1, x_plus[3] = x1, x_plus[4] = x0;
                    y[0] = x0, y[1] = y1, y[2] = y1, y[3] = x1, y[4] = x0;
                }
                else if (x0 <= y0 && x0 <= y1 && x1 > y0 && x1 <= y1) 
                {
                    size = 5;
                    x_plus[0] = x0, x_plus[1] = x0, x_plus[2] = x1, x_plus[3] = x1, x_plus[4] = y0;
                    y[0] = y0, y[1] = y1, y[2] = y1, y[3] = x1, y[4] = y0;
                }
                else if (x0 <= y0 && x0 <= y1 && x1 > y0 && x1 > y1) 
                {
                    size = 5;
                    x_plus[0] = x0, x_plus[1] = x0, x_plus[2] = y1, x_plus[3] = y0, x_plus[4] = x0;
                    y[0] = y0, y[1] = y1, y[2] = y1, y[3] = y0, y[4] = y0;
                }
                else
                {
                    n++;
                    continue;
                }

                for (uint l = 0; l < size; l++)
                    x_minus[l] = -x_plus[l];

                if (drawPlus >> 0 & 1)
                {
                    TGraph *g = new TGraph(size, x_plus, y); 
                    g->SetLineWidth(2);
                    g->SetEditable(false);
                    g->SetBit(kCanDelete);
                    g->SetEditable(kFALSE);
                    auto color = getColor(f_plus[n], min, max);
                    g->SetFillColor(color);
                    g->SetLineColor(drawGrid ? 1 : color);
                    mg->Add(g);
                }

                if (drawPlus >> 1 & 1)
                {
                    TGraph *g = new TGraph(size, x_minus, y); 
                    g->SetLineWidth(2);
                    g->SetEditable(false);
                    g->SetBit(kCanDelete);
                    g->SetEditable(kFALSE);
                    auto color = getColor(f_minus[n], min, max);
                    g->SetFillColor(color);
                    g->SetLineColor(drawGrid ? 1 : color);
                    mg->Add(g);
                }

                n++;
            }
        }

        mg->SetTitle(";#sigma(v_{||})#muB, KeV;K, KeV");

        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();

        mg->Draw("ALF");
        if (coloBar)
           drawColorBar(0.91, 0.99, 0.1, 0.9, min, max);
    }


    void drawFE() {
        if (error)
            return;

        TString canvas_name = "graph_f(E)";
        TCanvas *c = makeCanvas(canvas_name);
        c->SetBit(kCanDelete);
        c->SetLogy();

        TH1D *hist = new TH1D("", ";E, KeV;f(E), s^{3}/cm^{6}", number_band, band.front(), band.back());
        hist->SetLineWidth(2);
        hist->SetBit(kCanDelete);


        darray fE(number_band, 0.);

        uint cell = 0;
        for (uint i = 0; i < number_band; i++)
        {
            double V = 0;
            double f = 0;

            for (uint j = cells_index[i]; j+1 < cells_index[i+1]; j++)
            {
                V += volume[cell];
                f += (f_plus[cell]+f_minus[cell])*volume[cell];
                cell++;
            }
            V *= 2.;
            fE[i] = f/V*coeff_to_f;
            hist->SetBinContent(i+1, fE[i]);
        }

        countTiNi(fE, 0, 50);
        hist->SetStats(kFALSE);
        hist->GetXaxis()->CenterTitle();
        hist->GetYaxis()->CenterTitle();
        hist->Draw();


    }


    void countTiNi(const darray &fE, uint start_points=0, uint step_points=1) {
        double Y = 0;
        double X = 0;
        double Y2 = 0;
        double X2 = 0;
        double XY = 0;
        double N = 0;
        
        for (uint i = start_points; i < std::min(start_points+step_points, number_band); i++)
        {
            double y = log(fE[i]);
            double x = (band[i]+band[i+1]) / 2.; 
            Y += y;
            Y2 += y*y;
            X += x;
            X2 += x*x;
            XY += x*y;
            N += 1;
        }

        double betta = (N*XY - X*Y) / (N*X2-X*X);
        double alpha = (Y - betta*X)/N;
        double Ti = -1./betta;
        double ni = exp(alpha)*pow(MP*M/2./M_PI/Ti/normaE/EV_ERG, -1.5);

        std::cout << "Ti=" << Ti << " ni=" << ni << "\n";
    }

    void Log(double logmin) 
    {
        if (logmin <= 0. || logmin >= 1.) 
        {
            colorMult = -1.;
            colorAdd = 0.;
            axisMin = 0.;
        }
        else 
        {
            const double l = -log(logmin);
            colorAdd = l;
            colorMult = 1./l;
            axisMin = logmin;
        }
    }

    bool isError() { return error; }
};

void Draw(
    std::string fileName, uint it, uint iz, uint ir, uint ionType=0, double logmin=-1., 
    bool drawGrid=false, bool isMinZero=false, bool drawColorBar=true, 
    EColorPalette colorMapName=EColorPalette::kBlueRedYellow, int drawType=0b11, int drawPlus=0b11
) 
{
    if (drawType == 0)
        return;
    DrawMesh ps(fileName, it, iz, ir, ionType);
    ps.Log(logmin);
    if (!ps.isError()) 
    {
        if (drawType >> 0 & 1)
            ps.drawMesh(drawGrid, isMinZero, drawColorBar, colorMapName, drawPlus);
        if (drawType >> 1 & 1)
            ps.drawFE();

        double nakl = ps.countNumberParticle();
        double V3R = ps.count3RV(iz, ir);
        std::cout << "nakl: " << nakl << "\n";
        std::cout << "Vv3: " << 2.*ps.count3VFull() << "\n"; 
        std::cout << "Vr3: " << V3R << "\n";
        std::cout << "Nakl: " << nakl*V3R << "\n";
    }
}