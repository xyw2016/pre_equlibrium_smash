#ifndef PARAMETER_CPP
#define PARAMETER_CPP
#include "parameter.h"
#include <iostream>

#include "./inih/INIReader.h"

namespace Parameter {

int NX=67;
int NY=67;
int NZ=67;
int NETA=67;

double DX=0.3;
double DY=0.3;
double DZ=0.3;
double DETA=0.3;


double SIGR=0.6;
double SIGZ=0.6;
double SIGETA=0.6;
double TAU0 = 0.4;
std::string PATHIN="./";
std::string PATHOUT="output";

void Setup(INIReader &reader) {
  PATHIN = reader.Get("Input", "PATHIN", PATHIN);
  PATHOUT = reader.Get("Input", "PATHOUT", PATHOUT);
  NX   = reader.GetInteger("Input", "NX", NX);
  NY   = reader.GetInteger("Input", "NY", NY);
  NZ   = reader.GetInteger("Input", "NZ", NZ);
  NETA = reader.GetInteger("Input", "NETA", NETA);


  DX   = reader.GetReal("Input", "DX", DX);
  DY   = reader.GetReal("Input", "DY", DY);
  DZ   = reader.GetReal("Input", "DZ", DZ);
  DETA = reader.GetReal("Input", "DETA", DETA);

  SIGR     = reader.GetReal("Input", "SIGR", SIGR);
  SIGZ     = reader.GetReal("Input", "SIGZ", SIGZ);
  SIGETA   = reader.GetReal("Input", "SIGETA", SIGETA);
  TAU0 = reader.GetReal("Input", "TAU0", TAU0);



  std::cout << "\n ******************************************" ;
  std::cout << "\n The parameter sets" 
            << "  NX     = " << NX << "\n"
            << "  NY     = " << NY << "\n"
            << "  NETA   = " << NETA << "\n"
            << "  DX     = " << DX << "\n"
            << "  DY     = " << DY << "\n"
            << "  DETA   = " << DETA<< "\n"
            << "  SIGR   = " << SIGR << "\n"
            << "  SIGZ   = " << SIGZ << "\n"
            << "  SIGETA = " << SIGETA << "\n"
            << "  TAU0   = " << TAU0 <<  std::endl;
  std::cout << "\n ******************************************" ;
}
}


#endif
