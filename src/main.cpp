#include <iostream>
#include <cstring>
//#include <omp.h>

//input file parser
#include "INIReader.h"

#include "parameter.h"
#include "reader.h"
#include "grid.h"
#ifndef _OPENMP
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif

int main(int argc, char **argv) {
  // Try to open input file
  if (argc < 2) {
    std::cerr << "Error: ** " << argv[0] << " ** no input filename " << std::endl;
    std::cerr << "USAGE:" << std::endl;
    std::cerr << "    " << argv[0] << " setup.ini" << std::endl;
    exit(EXIT_FAILURE);
  }

  INIReader reader(argv[1]) ;
  if (reader.ParseError()){
    std::cerr << "Error: ** " << argv[0] << " ** failed to open " << argv[1] << std::endl;
    exit(EXIT_FAILURE);
  }


  Parameter::Setup(reader);



  Events smash_events;
  smash_events.Read_ini();
  
  int eos_type = Parameter::EOS_ID;
  // //Grid gird0(smash_events.allevents);
  Grid grid0(eos_type);
  grid0.preequlibirum(smash_events.ptcl_event);
  grid0.hydro_ini(smash_events.ptcl_event);

  // //gird0.GausssmearingTZ(smash_events.allevents);
  // gird0.GausssmearingTauEta2(smash_events.NIsotauptc);



}

