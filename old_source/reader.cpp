#ifndef READER_CPP
#define READER_CPP

#include "reader.h"

Events::Events(){
  NIsotauptc.clear();
  allevents.clear();

}

Events::~Events(){  
}


void Events::Readprestage(){
  
  std::stringstream fnamefinal;
  fnamefinal<<Parameter::PATHIN<<"/particle_lists.oscar";
  std::ifstream fin(fnamefinal.str());
  Particle particle0;
  std::string dump;
  std::string buf;
  double eps=1e-16;
  int ievent=0;
  Particlelist ptclist;
  TEPaticlelist Nptclist;
  Nptclist.clear();
  ptclist.clear();


  if(fin.is_open() )
  {
      getline(fin, buf);
      getline(fin, buf);
      getline(fin, buf);
      
      std::istringstream iss;
      int nline=0;
      while(fin.good()){
        int nptc=0;
        getline(fin,buf);
        iss.str(buf);
        iss >> dump >> dump >> dump >>dump >> nptc;
        iss.clear();
        if (nptc > 0){
          ptclist.clear();
          for(int iptc=0 ; iptc< nptc; iptc++)
          {
            fin >>  particle0.t    >> particle0.x >> particle0.y>>particle0.z 
                >>  particle0.mass >> particle0.e >> particle0.px >> particle0.py
                >>  particle0.pz   >> particle0.pid >> dump >> dump >> particle0.ncoll
                >> dump>> dump >> dump >> dump >> dump>> dump>> dump>> particle0.baryon_number;
            if (particle0.t < 0 ) continue;
          

            particle0.p0  = sqrt(particle0.px*particle0.px + particle0.py*particle0.py + + particle0.pz*particle0.pz );
            if (fabs(particle0.t) > fabs(particle0.z)){
            particle0.tau = sqrt(particle0.t*particle0.t-particle0.z*particle0.z);
            particle0.etas = 0.5*(log(particle0.t+particle0.z+eps) - log(particle0.t-particle0.z+eps)); 
            
          }
          
          ptclist.push_back(particle0);
          }
          if (ptclist.size()==0) continue;
          Nptclist.push_back(ptclist);
          ptclist.clear();
        }

        if (fin.eof()) break;
        if(dump=="end")
        {
          allevents.push_back(Nptclist);
          Nptclist.clear();
        }

 
    }
	fin.close();
  }
  else{

        std::cerr<<"#Can't open "<< fnamefinal.str()<<"!!!!!! \n" ;
        exit(0);
  }

  std::cout<<" Finish read pre-equlibrium stage data :) "<<std::endl;

}


void Events::Readfinal(){
  
  std::stringstream fnamefinal;
  fnamefinal<<Parameter::PATHIN<<"/SMASH_IC.oscar";
  std::ifstream fin(fnamefinal.str());
  Particle particle0;
  std::string dump;
  std::string buf;

  double eps=1e-16;
  int ievent=0;
  Particlelist Isotauptc;
  if(fin.is_open() )
  {
      getline(fin, buf);
      getline(fin, buf);
      getline(fin, buf);
      
      std::istringstream iss;
      int nline=0;
      
      while(fin.good()){
      
        if (fin.eof()) break;

        while(getline(fin,buf))
        {
           iss.str(buf);
           iss >> dump >> dump >> dump >>dump;
           if (dump == "in"){
            iss.clear();
            Isotauptc.clear();
            continue;
           }
           else if(dump == "end"){
            iss.clear();
            std::cout<<"Finish read event " << ievent <<std::endl;
            ievent++;
            break;
           }
           else{
            iss.seekg(0);
           }

           

            iss >>  particle0.t    >> particle0.x >> particle0.y>>particle0.z 
                >>  particle0.mass >> particle0.e >> particle0.px >> particle0.py
                >>  particle0.pz   >> particle0.pid >> dump >> dump >> particle0.ncoll
                >> dump>> dump >> dump >> dump >> dump>> dump>> dump>> particle0.baryon_number;
            iss.clear();
            particle0.p0  = sqrt(particle0.px*particle0.px + particle0.py*particle0.py + + particle0.pz*particle0.pz );
            if (fabs(particle0.t) > fabs(particle0.z)){
            particle0.tau = sqrt(particle0.t*particle0.t-particle0.z*particle0.z);
            particle0.etas = 0.5*(log(particle0.t+particle0.z+eps) - log(particle0.t-particle0.z+eps));
            }
            
            Isotauptc.push_back(particle0);
            
        }
        if (fin.eof()) break;
        NIsotauptc.push_back(Isotauptc);
    }
	fin.close();
  }
  else{

        std::cerr<<"#Can't open "<< fnamefinal.str()<<"!!!!!! \n" ;
        exit(0);
  }


  std::cout<<" Finish read initial state data :) "<<std::endl;
  
}




#endif
