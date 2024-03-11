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

    void GridClear(Grid3D& grid0);
    void hydro_ini(Paticlelist_event& ptcl_event);
    // void GausssmearingTZ(EventsofTEPaticlelist& Nptclist);
    // void GausssmearingTauEta(TEPaticlelist& Nptclist);
    // void GausssmearingTauEta2(TEPaticlelist& Nptclist);
    void preequlibirum(Paticlelist_event& ptcl_event);
    void smearing_kernel(Particlelist& ptclist);

    void perform_laudu(double coutevent);

    
    private:

    Grid3D Tgrid;
    int CORES;

    EOS eos;
    //const EventsofTEPaticlelist NTptclist;

    const int NX;
    const int NY;
    const int NETA;
    

    const double DX;
    const double DY;
    const double DETA;

    

    const double SIGR;
    const double SIGETA;
    const double TAU0;
    const double TAU00;
    const int NTAU;


    double one_o_2sigr2; 
    double one_o_2sigz2;

    double w1; 
    double w2; 

   
    
};


#endif
