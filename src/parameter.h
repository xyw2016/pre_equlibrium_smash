#ifndef PARAMETER_H
#define PARAMETER_H
#ifndef HBARC
#define HBARC 0.197326979
#endif
#ifndef small_eps
#define small_eps 1e-16
#endif  
#ifndef M_PI_F
#define M_PI_F 3.1415926
#endif
#include <string>
#include <vector>
#include "INIReader.h"
#include <vector>



class INIReader;

struct Particle{
    double t;
    double x;
    double y;
    double z;
    double mass;
    double e;
    double p0;
    double px;
    double py;
    double pz;
    int pid;
    int charged;
    double tau;
    double etas;
    int ncoll;
    double baryon_number;
};

typedef std::vector<Particle> Particlelist;// Fixed-time Particle Spectrum
typedef std::vector<Particlelist> Paticlelist_tau; //Particle list at different time
typedef std::vector<Paticlelist_tau> Paticlelist_event; //Particle list at different event




namespace Parameter {


extern std::string PATHIN;
extern std::string PATHOUT;
extern int NX;
extern int NY;
extern int NETA;


extern double DX;
extern double DY;
extern double DETA;


extern double SIGR;
extern double SIGETA;

extern double TAU00;
extern double TAU0;
extern int NTAU;

extern double TCUT;

extern int EOS_ID;
extern int turn_on_rhob;
extern int turn_on_shear;
extern int turn_on_bulk;
extern int turn_on_diff;

extern int read_binary;



void Setup(INIReader &reader);
}
#endif
