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
  int read_binary = Parameter::read_binary;
  std::stringstream fnamefinal;
  
  //fnamefinal<<Parameter::PATHIN<<"/particle_lists.oscar";
  fnamefinal<<Parameter::PATHIN;
  
  std::ifstream fin;
  if(read_binary){
      fin.open(fnamefinal.str().c_str(),std::ios::binary);
  }
  else{
      fin.open(fnamefinal.str().c_str());
  }
  //std::ifstream fin(fnamefinal.str());
  Particle particle0;
  std::string dump;
  std::string buf;
  double eps=1e-16;
  int ievent=0;
  Particlelist ptclist;
  TEPaticlelist Nptclist;
  Nptclist.clear();
  ptclist.clear();

  if(!fin.is_open()){
     std::cerr<<"#Can't open "<< fnamefinal.str()<<"!!!!!! \n" ;
      exit(0);
  }

  std::istringstream iss;
  if(read_binary)
  {
    while(fin.good()){
      float fnptc=0;
      fin.read(reinterpret_cast<char*>(&fnptc), sizeof(float));
      int nptc0 = static_cast<int> (fnptc);
      
       if(nptc0>0 && nptc0 != 999999)
          { 
            ptclist.clear();
            float array[13];
            for(int iptc = 0; iptc < nptc0; iptc++ ){

              for (int i = 0; i < 13; i++) {
                float temp = 0.;
                fin.read(reinterpret_cast<char*>(&temp), sizeof(float));
                array[i] = temp;
              }
              particle0.t = array[0];
              particle0.x = array[1];
              particle0.y = array[2];
              particle0.z  = array[3];
              particle0.mass = array[4];
              particle0.e = array[5];
              particle0.px = array[6];
              particle0.py = array[7];
              particle0.pz = array[8];
              particle0.pid = static_cast<int> (array[9]);
              particle0.ncoll = static_cast<int> (array[11]);
              particle0.baryon_number = static_cast<int> (array[12]);

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
        if(nptc0==999999)
        {
          allevents.push_back(Nptclist);
          Nptclist.clear();
        }
    }
    
    


  }
  else{

    
    getline(fin, buf);
    getline(fin, buf);
    getline(fin, buf);
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

  }

  fin.close();
  std::cout<<" Finish read pre-equlibrium stage data :) "<<std::endl;

}


void Events::Readfinal(){
  
  int read_binary = Parameter::read_binary;
  std::stringstream fnamefinal;
  //fnamefinal<<Parameter::PATHIN<<"/SMASH_IC.oscar";
  
  fnamefinal<<Parameter::PATHIN2;
  std::ifstream fin;
  if(read_binary){
      fin.open(fnamefinal.str().c_str(),std::ios::binary);
  }
  else{
      fin.open(fnamefinal.str().c_str());
  }

  Particle particle0;
  std::string dump;
  std::string buf;

  double eps=1e-16;
  int ievent=0;
  Particlelist Isotauptc;
  if(!fin.is_open() ){
    std::cerr<<"#Can't open "<< fnamefinal.str()<<"!!!!!! \n" ;
    exit(0);

  }


  if(read_binary){
     while(fin.good()){

      if (fin.eof()) break;
      
      float fnptc=0;


      while(fin.read(reinterpret_cast<char*>(&fnptc), sizeof(float))){
        int nptc0 = static_cast<int> (fnptc);
       
        if(nptc0==888888){
          Isotauptc.clear();
          std::cout<<"start read event " << ievent <<std::endl;
          continue;
        }
        if(nptc0==999999){
          std::cout<<"Finish read event " << ievent <<std::endl;
          ievent++;
          break;
        }

        float array[12];
        for (int i = 0; i < 12; i++) {
            float temp = 0.;
            fin.read(reinterpret_cast<char*>(&temp), sizeof(float));
            array[i] = temp;
        }
        
        particle0.t = fnptc;
        particle0.x = array[0];
        particle0.y = array[1];
        particle0.z  = array[2];
        particle0.mass = array[3];
        particle0.e = array[4];
        particle0.px = array[5];
        particle0.py = array[6];
        particle0.pz = array[7];
        particle0.pid = static_cast<int> (array[8]);
        particle0.ncoll = static_cast<int> (array[10]);
        particle0.baryon_number = static_cast<int> (array[11]);

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
    
  }
  
  else{

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

  }
  
	fin.close();
  
 


  std::cout<<" Finish read initial state data :) "<<std::endl;
  
}




#endif
