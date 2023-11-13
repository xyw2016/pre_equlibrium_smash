#ifndef GRID_H
#define GRID_H

#include "parameter.h"
#include "reader.h"
#include "eos/eos.h"



#define idx(I,mn) (I)*4+mn
#define GAMMAMAX 60

struct Cell{
    double ed;
    double nb;
    double u0;
    double ux;
    double uy;
    double uz;

    double Temperature;
    double pressure;
    double muB;
    double cs2;
    double Tmnnu[16];
    double Jbmu[4];
    //std::vector<double> Tmnnu(16); //T^{\mu\nu}
    //std::vector<double> Tmnnu(16, 0.0); 

};

typedef std::vector<std::vector<std::vector<Cell>> > Grid3D;

class Grid{

    public:

    //Grid(EventsofTEPaticlelist& Nptclist);
    Grid(int eos_type);

    ~Grid();

    void Laudumatching(double* Tmn, double& ed, double& u0, double& ux, double& uy, double& uz);
    void rootFinding_newton(double& K0,double& M,double& J0, double& ed_find, double& nb_find);

    void GridClear(Grid3D& grid0,int NZ0);
    void GausssmearingTZ(EventsofTEPaticlelist& Nptclist);
    void GausssmearingTauEta(TEPaticlelist& Nptclist);
    void GausssmearingTauEta2(TEPaticlelist& Nptclist);
    

    
    private:

    Grid3D Tgrid;
    Grid3D Taugrid;
    
    EOS eos;
    //const EventsofTEPaticlelist NTptclist;

    const int NX;
    const int NY;
    const int NZ;
    const int NETA;
    

    const double DT;
    const double DX;
    const double DY;
    const double DZ;
    const double DETA;

    

    const double SIGR;
    const double SIGZ;
    const double SIGETA;
    const double TAU0;

    int nmaxtime;
   
    
};


#endif
