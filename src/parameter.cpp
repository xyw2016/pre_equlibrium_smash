#ifndef PARAMETER_CPP
#define PARAMETER_CPP
#include "parameter.h"
#include <iostream>

#include "INIReader.h"

namespace Parameter {

int NX=67;
int NY=67;
int NETA=67;


double DX=0.3;
double DY=0.3;
double DETA=0.3;


double SIGR=0.6;
double SIGETA=0.6;

double TAU00 = 0.0;
double TAU0  = 0.4;
int NTAU = 10;

double TCUT  = 0.1;

int EOS_ID = 14;

int turn_on_rhob  = 1;
int turn_on_shear = 0;
int turn_on_bulk  = 0;
int turn_on_diff  = 0;

int read_binary = 0;

std::string PATHIN="./";
std::string PATHOUT="output";

void Setup(INIReader &reader) {
  PATHIN = reader.Get("Input", "PATHIN", PATHIN);

  PATHOUT = reader.Get("Input", "PATHOUT", PATHOUT);
  NX   = reader.GetInteger("Input", "NX", NX);
  NY   = reader.GetInteger("Input", "NY", NY);
  NETA = reader.GetInteger("Input", "NETA", NETA);

  DX   = reader.GetReal("Input", "DX", DX);
  DY   = reader.GetReal("Input", "DY", DY);
  DETA = reader.GetReal("Input", "DETA", DETA);

  SIGR     = reader.GetReal("Input", "SIGR", SIGR);
  SIGETA   = reader.GetReal("Input", "SIGETA", SIGETA);
  
  TAU00 = reader.GetReal("Input", "TAU00", TAU00);
  TAU0 = reader.GetReal("Input", "TAU0", TAU0);
  NTAU = reader.GetInteger("Input", "NTAU", NTAU);

  TCUT = reader.GetReal("Input", "TCUT", TCUT);

  EOS_ID   = reader.GetInteger("Input", "EOS_ID", EOS_ID);
  turn_on_rhob   = reader.GetInteger("Input", "turn_on_rhob", turn_on_rhob);
  turn_on_shear   = reader.GetInteger("Input", "turn_on_shear", turn_on_shear);
  turn_on_bulk   = reader.GetInteger("Input", "turn_on_bulk", turn_on_bulk);
  turn_on_diff = reader.GetInteger("Input", "turn_on_diff", turn_on_diff);
  read_binary = reader.GetInteger("Input", "read_binary", turn_on_diff);



  std::cout << "\n ****************************************** \n " ;
  std::cout << "\n The parameter sets\n" 
            << "  NX     = " << NX << "\n"
            << "  NY     = " << NY << "\n"
            << "  NETA   = " << NETA << "\n"
            << "  DX     = " << DX << "\n"
            << "  DY     = " << DY << "\n"
            << "  DETA   = " << DETA<< "\n"
            << "  SIGR   = " << SIGR << "\n"
            << "  SIGETA = " << SIGETA << "\n"
            << "  TAU00  = " << TAU00 <<"\n"
            << "  TAU0   = " << TAU0 <<"\n"
            << "  NTAU   = " << NTAU <<"\n"
            << "  TCUT   = " << TCUT <<"\n"
            << "  turn_on_rhob    = " << turn_on_rhob << "\n"
            << "  turn_on_shear   = " << turn_on_shear << "\n"
            << "  turn_on_bulk    = " << turn_on_bulk << "\n"
            << "  turn_on_diff    = " << turn_on_diff << "\n" 
            << "  read_binary     = " << read_binary << "\n"
            << "  import file     = " << PATHIN << "\n"
            << "  output path     = " << PATHOUT << "\n" <<  std::endl;
  std::cout << "\n ****************************************** \n" ;
}
}


#endif
