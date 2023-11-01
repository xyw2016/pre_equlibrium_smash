#ifndef GRID_CPP
#define GRID_CPP

#include "grid.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>



Grid::Grid(int use_eos):NX(Parameter::NX), NY(Parameter::NY),NZ(Parameter::NZ),
            DX(Parameter::DX), DY(Parameter::DY),DZ(Parameter::DZ),
            DETA(Parameter::DETA), NETA(Parameter::NETA),
            SIGR(Parameter::SIGR), SIGZ(Parameter::SIGZ),
            SIGETA(Parameter::SIGETA), TAU0(Parameter::TAU0),eos(use_eos){
  
    Cell cellzero;

    std::cout<<" "<<eos.get_temperature(0.4/HBARC,0.0)*HBARC<<std::endl;

    cellzero.ed = 0;
    cellzero.nb = 0;
    cellzero.u0 = 1;
    cellzero.ux = 0;
    cellzero.uy = 0;
    cellzero.uz = 0;
    cellzero.Temperature = 0;
    cellzero.pressure = 0;
    for(int id = 0 ; id < 16 ; id ++)
    {
      cellzero.Tmnnu[id] = 0;
    }  

    for(int id = 0 ; id < 4 ; id ++)
    {
      cellzero.Jbmu[id] = 0;
    }       


    Tgrid.resize(NX);
    Taugrid.resize(NX);
    for (int i = 0; i < NX; i++) {
        Tgrid[i].resize(NY);
        Taugrid[i].resize(NY);
        for (int j = 0; j < NY; j++) {
            Tgrid[i][j].resize(NZ);
            Taugrid[i][j].resize(NETA);
            for(int k = 0; k < NZ; k++){
              Tgrid[i][j][k] = cellzero;
            }

            for(int k = 0; k < NETA; k++){
              Taugrid[i][j][k] = cellzero;
            }

        }
        
    }

    nmaxtime = 0;
    



}

void Grid::GridClear(Grid3D& grid0,int NZ0){

  for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
          for (int k =0 ; k < NZ0 ;k++){
            grid0[i][j][k].ed = 0;
            grid0[i][j][k].nb = 0;
            grid0[i][j][k].u0 = 1;
            grid0[i][j][k].ux = 0;
            grid0[i][j][k].uy = 0;
            grid0[i][j][k].uz = 0;
            grid0[i][j][k].Temperature = 0;
            grid0[i][j][k].pressure = 0;
            for(int id = 0 ; id < 16 ; id ++)
            {
               grid0[i][j][k].Tmnnu[id] = 0;
            }   

            for(int id = 0 ; id < 4 ; id ++)
            {
               grid0[i][j][k].Jbmu[id] = 0;
            }   


          }
        
        }
    }
}


void Grid::Laudumatching(double* Tmn, double& ed, double& u0, double& ux, double& uy, double& uz){

  gsl_matrix * Tmunu;
  Tmunu = gsl_matrix_alloc(4,4); //T^{\mu}_{\nu}

  gsl_matrix_complex *eigen_vectors;
  eigen_vectors = gsl_matrix_complex_alloc(4,4);
  gsl_vector_complex *eigen_values;
  eigen_values = gsl_vector_complex_alloc(4);

  // //set the values of the energy momentum tensor
  gsl_matrix_set(Tmunu, 0, 0,  Tmn[idx(0,0)] ); 
  gsl_matrix_set(Tmunu, 0, 1, -Tmn[idx(0,1)] );
  gsl_matrix_set(Tmunu, 0, 2, -Tmn[idx(0,2)] ); 
  gsl_matrix_set(Tmunu, 0, 3, -Tmn[idx(0,3)] ); 
  
  gsl_matrix_set(Tmunu, 1, 0,  Tmn[idx(1,0)] ); 
  gsl_matrix_set(Tmunu, 1, 1, -Tmn[idx(1,1)] ); 
  gsl_matrix_set(Tmunu, 1, 2, -Tmn[idx(1,2)] ); 
  gsl_matrix_set(Tmunu, 1, 3, -Tmn[idx(1,3)] ); 

  gsl_matrix_set(Tmunu, 2, 0,  Tmn[idx(2,0)] ); 
  gsl_matrix_set(Tmunu, 2, 1, -Tmn[idx(2,1)] ); 
  gsl_matrix_set(Tmunu, 2, 2, -Tmn[idx(2,2)] ); 
  gsl_matrix_set(Tmunu, 2, 3, -Tmn[idx(2,3)] ); 

  gsl_matrix_set(Tmunu, 3, 0,  Tmn[idx(3,0)] ); 
  gsl_matrix_set(Tmunu, 3, 1, -Tmn[idx(3,1)] ); 
  gsl_matrix_set(Tmunu, 3, 2, -Tmn[idx(3,2)] ); 
  gsl_matrix_set(Tmunu, 3, 3, -Tmn[idx(3,3)] ); 
  
  

  
  gsl_eigen_nonsymmv_workspace *eigen_workspace;
  eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
  gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
  gsl_eigen_nonsymmv_free(eigen_workspace);

  int eigenvalue_exists = 0;
  for (int i = 0; i < 4; i++)
  {
    gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);
            
    if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
    {
      gsl_complex v0 = gsl_matrix_complex_get(eigen_vectors, 0 , i);
      gsl_complex v1 = gsl_matrix_complex_get(eigen_vectors, 1 , i);
      gsl_complex v2 = gsl_matrix_complex_get(eigen_vectors, 2 , i);
      gsl_complex v3 = gsl_matrix_complex_get(eigen_vectors, 3 , i);

      double minkowskiLength = GSL_REAL(v0)*GSL_REAL(v0) - (GSL_REAL(v1)*GSL_REAL(v1) + GSL_REAL(v2)*GSL_REAL(v2) + GSL_REAL(v3)*GSL_REAL(v3));

      if (GSL_IMAG(v0) == 0 && minkowskiLength > 0) //choose timelike eigenvector
      {
        
        double flowfactor = 1.0 / sqrt(minkowskiLength);
                    
        if (GSL_REAL(v0) < 0) flowfactor=-flowfactor;
                    
                    //ignore eigenvectors with gamma >~ 60
        if ( (GSL_REAL(v0) * flowfactor) < GAMMAMAX)
        {
            eigenvalue_exists = 1;
            ed = GSL_REAL(eigenvalue);
            u0 = GSL_REAL(v0) * flowfactor;
            ux = GSL_REAL(v1) * flowfactor;
            uy = GSL_REAL(v2) * flowfactor;
            uz = GSL_REAL(v3) * flowfactor;
            //std::cout<< ed <<" "<<u0*u0 - ux*ux - uy*uy-uz*uz<<" "<<ux<<" "<<uy<<" "<<uz<<std::endl;
        }
                    
      } // if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
    } // if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
  } //for (int i = 0; i < 4; ...)







  gsl_vector_complex_free(eigen_values);
	gsl_matrix_complex_free(eigen_vectors);
  gsl_matrix_free(Tmunu);
        

}

void Grid::rootFinding_newton(double& T00,double& M,double& J0, double& ed_find, double& nb_find){

    double vl = 0.0;
    double vh = 1.0;
    double v = 0.5*(vl+vh);
    double ed = T00 - M*v;
    double nb = (J0 )*sqrt(1-v*v);
    double pr = eos.get_pressure(ed/HBARC, nb)*HBARC;
    //real dpe = eos_CS2(ed,nb,eos_table);
    double dpe = 1.0/3.0;
    double f = (T00 + pr)*v - M;
    double df = (T00+pr) - M*v*dpe;
    double dvold = vh-vl;
    double dv = dvold;
    int i = 0;
    while ( true ) {
        if ((f + df * (vh - v)) * (f + df * (vl - v)) > 0.0 ||
            fabs(2. * f) > fabs(dvold * df)) {  // bisection
          dvold = dv;
          dv = 0.5 * (vh - vl);
          v = vl + dv;
        } else {  // Newton
          dvold = dv;
          dv = f / df;
          v -= dv;
        }
        v = std::max(0.0,std::min(v,1.0));
        i ++;
        if ( fabs(dv) < 0.00001 || i > 100 ) break;

        ed = T00 - M*v;
        nb = (J0)*sqrt(1-v*v);
        pr = eos.get_pressure(ed/HBARC, nb)*HBARC;
        f = (T00 + pr)*v - M;
        //dpe = eos_CS2(ed,nb,eos_table);
        dpe = 1.0/3.0;
        df = (T00+pr) - M*v*dpe;
        if ( f > 0.0f ) {
             vh = v;
        } else { 
             vl = v;
        }
    }
    // accelerate the speed of iteration by combining with bisection(low slope) and newton method(high slope)
    ed_find = T00 - M*v;
    nb_find = (J0 )*sqrt(1-v*v);


}


void Grid::GausssmearingTZ(EventsofTEPaticlelist& Nptclist){
   double one_o_2sigr2 = 1.0/(2.0*SIGR*SIGR); 
   double one_o_2sigz2 = 1.0/(2.0*SIGZ*SIGZ);

   double w1 = one_o_2sigr2/M_PI_F;
   double w2 = sqrt(one_o_2sigz2/M_PI_F);
   
   for(int i = 0 ; i < Nptclist.size(); i++){
      int ntime = Nptclist[i].size();
      nmaxtime = std::max(nmaxtime,ntime);
   }
   

  
  for(int itime =0 ; itime < nmaxtime; itime++)
  {
    int coutevent = 0;
    

    GridClear(Tgrid,NZ);
    bool flag_out_grid = 0;
    int out = 0;
    for(int ievent=0; ievent< Nptclist.size(); ievent++) 
    { 
      
      
      if(itime >= Nptclist[ievent].size() ) {
        out += 1;
        continue;
        }

      std::cout<<ievent<<" "<<itime<<" "<<Nptclist[ievent][itime].size()<<std::endl;
      
      double ptc_xmin = Parameter::NX/2* Parameter::DX; 
      double ptc_ymin = Parameter::NY/2* Parameter::DY; 
      double ptc_zmin = Parameter::NZ/2* Parameter::DZ; 
      for(int iptc=0;iptc < Nptclist[ievent][itime].size(); iptc++){
              
            ptc_xmin = std::min(ptc_xmin,fabs(Nptclist[ievent][itime][iptc].x));
            ptc_ymin = std::min(ptc_ymin,fabs(Nptclist[ievent][itime][iptc].y));
            ptc_zmin = std::min(ptc_zmin,fabs(Nptclist[ievent][itime][iptc].z));
      }

       if (fabs(ptc_xmin) >= Parameter::NX/2* Parameter::DX ||
           fabs(ptc_ymin) >= Parameter::NY/2* Parameter::DY ||
           fabs(ptc_zmin) >= Parameter::NZ/2* Parameter::DZ  ) {
            flag_out_grid=1;
            out +=1;
            continue;
           }
      
      std::cout<<ievent<<" "<<itime<<" "<<ptc_zmin<<std::endl;
      
      

      for(int gridi = 0; gridi < NX; gridi++ ){
      double xi= (gridi-NX/2)*DX;
      for(int gridj = 0; gridj < NY; gridj++ ){
          double yj = (gridj-NY/2)*DY;
          for(int gridk = 0; gridk < NZ; gridk++ ){
            double zk = (gridk-NZ/2)*DZ;
            int nncoll=0;
            for(int iptc=0;iptc < Nptclist[ievent][itime].size(); iptc++){
              double ptc_t = Nptclist[ievent][itime][iptc].t;
              double ptc_x = Nptclist[ievent][itime][iptc].x;
              double ptc_y = Nptclist[ievent][itime][iptc].y;
              double ptc_z = Nptclist[ievent][itime][iptc].z;
              double ptc_e = Nptclist[ievent][itime][iptc].e;
              double ptc_px = Nptclist[ievent][itime][iptc].px;
              double ptc_py = Nptclist[ievent][itime][iptc].py;
              double ptc_pz = Nptclist[ievent][itime][iptc].pz;
              double ptc_baryon = Nptclist[ievent][itime][iptc].baryon_number;
              double ptc_p4[4] = {ptc_e,ptc_px,ptc_py,ptc_pz};
              int ncoll = Nptclist[ievent][itime][iptc].ncoll;

            

              double ptc_tausq = ptc_t*ptc_t - ptc_z*ptc_z;
              if(ncoll == 0 || ptc_tausq < 0) continue;
              nncoll++;
              double dxx = ptc_x - xi;
              double dyy = ptc_y - yj;
              double dzz = ptc_z - zk;
              
              double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);

              if ( fabs(dzz) < 10*SIGZ && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
        
                double delta = exp(- distance_sqr)*w1*w2;
                for(int mu = 0 ; mu < 4 ; mu ++){
                  double factorP_bmu = ptc_baryon*ptc_p4[mu]/(ptc_p4[0]+1e-7);
                  Tgrid[gridi][gridj][gridk].Jbmu[mu] += factorP_bmu*delta;
                  for(int nu =0; nu < 4 ; nu ++){
                    double factorPP = ptc_p4[mu]*ptc_p4[nu]/(ptc_p4[0]+1e-7);
                    
                    Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)] += factorPP*delta;

                  }
                }
              }     

            }

            
      }
      }
      }

      coutevent+=1;


    }

    std::cout<<" Out grid event: "<< out<<std::endl;

    if(out == Nptclist.size()){
      break;
    }

   
    std::cout<< coutevent<<std::endl;
    double maxvale2=0.0;

    for(int gridi = 0; gridi < NX; gridi++ ){
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NZ; gridk++ ){
          double maxvale=0.0;
          for(int mu = 0 ; mu < 4 ; mu ++){
                  Tgrid[gridi][gridj][gridk].Jbmu[mu]  /= coutevent;
                  for(int nu =0; nu < 4 ; nu ++){
                    Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]  /= coutevent;
                    
                    maxvale= std::max(maxvale,Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);
                    maxvale2 = std::max(maxvale2,Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);

                  }
                }
          if( maxvale < 1e-4) continue;
          
          //std::cout<<maxvale<<std::endl;
          Laudumatching(Tgrid[gridi][gridj][gridk].Tmnnu,Tgrid[gridi][gridj][gridk].ed,
                        Tgrid[gridi][gridj][gridk].u0,Tgrid[gridi][gridj][gridk].ux,
                        Tgrid[gridi][gridj][gridk].uy,Tgrid[gridi][gridj][gridk].uz);

          Tgrid[gridi][gridj][gridk].nb = Tgrid[gridi][gridj][gridk].Jbmu[0]*Tgrid[gridi][gridj][gridk].u0
                                        - Tgrid[gridi][gridj][gridk].Jbmu[1]*Tgrid[gridi][gridj][gridk].ux
                                        - Tgrid[gridi][gridj][gridk].Jbmu[2]*Tgrid[gridi][gridj][gridk].uy
                                        - Tgrid[gridi][gridj][gridk].Jbmu[3]*Tgrid[gridi][gridj][gridk].uz;

          Tgrid[gridi][gridj][gridk].Temperature = eos.get_temperature(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb)*HBARC;             
        
        }
      }

    }
    
    if(maxvale2 > 0){

    
    std::string filename = Parameter::PATHOUT + std::to_string(itime) + ".dat";   
    std::cout<<" "<<filename<<std::endl; 
    std::ofstream fout(filename);


    for(int gridi = 0; gridi < NX; gridi++ ){
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NZ; gridk++ ){
          fout<< Tgrid[gridi][gridj][gridk].ed
              <<" "<< Tgrid[gridi][gridj][gridk].nb
              <<" "<< Tgrid[gridi][gridj][gridk].Temperature
              <<" "<< Tgrid[gridi][gridj][gridk].u0
              <<" "<< Tgrid[gridi][gridj][gridk].ux
              <<" "<< Tgrid[gridi][gridj][gridk].uy
              <<" "<< Tgrid[gridi][gridj][gridk].uz<<std::endl;     
        }
      }

    }

    fout.close();
    }
 

    
    std::cout<<maxvale2<<std::endl;

    

  }
  

}


void Grid::GausssmearingTauEta(TEPaticlelist& Nptclist){
   double one_o_2sigr2 = 1.0/(2.0*SIGR*SIGR); 
   double one_o_2sigz2 = 1.0/(2.0*SIGETA*SIGETA);

   double w1 = one_o_2sigr2/M_PI_F;
   double w2 = sqrt(one_o_2sigz2/M_PI_F)/TAU0;
  

    int coutevent = 0;


    GridClear(Taugrid,NETA);
    std::cout<<Nptclist.size()<<std::endl;
    for(int ievent=0; ievent< Nptclist.size(); ievent++) 
    //for(int ievent=0; ievent< 2; ievent++) 
    { 
      
      

      std::cout<<ievent<<" "<<Nptclist[ievent].size()<<std::endl;
      double vv = 0;
      for(int gridi = 0; gridi < NX; gridi++ ){
      double xi= (gridi-NX/2)*DX;
      for(int gridj = 0; gridj < NY; gridj++ ){
          double yj = (gridj-NY/2)*DY;
          for(int gridk = 0; gridk < NETA; gridk++ ){
            double zk = (gridk-NETA/2)*DETA;
            int nncoll=0;
           
            for(int iptc=0;iptc < Nptclist[ievent].size(); iptc++){
              double ptc_t = Nptclist[ievent][iptc].t;
              double ptc_x = Nptclist[ievent][iptc].x;
              double ptc_y = Nptclist[ievent][iptc].y;
              double ptc_z = Nptclist[ievent][iptc].z;
              double ptc_e = Nptclist[ievent][iptc].e;
              double ptc_px = Nptclist[ievent][iptc].px;
              double ptc_py = Nptclist[ievent][iptc].py;
              double ptc_pz = Nptclist[ievent][iptc].pz;
              double ptc_mass = Nptclist[ievent][iptc].mass;
              double ptc_baryon = Nptclist[ievent][iptc].baryon_number;
              double ptc_momentum[4] = {ptc_e,ptc_px,ptc_py,ptc_pz};
              double ptc_position[4] = {ptc_t,ptc_x,ptc_y,ptc_z};
              int ncoll = Nptclist[ievent][iptc].ncoll;
             
              double ptc_tausq = ptc_t*ptc_t - ptc_z*ptc_z;
              //std::cout << sqrt(ptc_tausq)<<std::endl;
              if(ncoll == 0 || ptc_tausq < 0) continue;
              nncoll++;
              
              
              double etasi = 0.5f * (log(std::max(ptc_position[0]+ptc_position[3], small_eps)) \
                           - log(std::max(ptc_position[0]-ptc_position[3], small_eps)));
          
              double dxx = ptc_x - xi;
              double dyy = ptc_y - yj;
              double dzz = etasi - zk;
             

              double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);


              if ( fabs(dzz) < 10*SIGETA && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
                 
              
                double delta = exp(- distance_sqr)*w1*w2;
                
                double mt = sqrt(ptc_momentum[1] * ptc_momentum[1] + ptc_momentum[2] * ptc_momentum[2]+ptc_mass*ptc_mass);
                //std::cout << mt<<" "<<delta<<" "<<ptc_baryon <<std::endl;
                double Yi = 0.5f * (log(std::max(ptc_momentum[0]+ptc_momentum[3], small_eps)) \
                           - log(std::max(ptc_momentum[0]-ptc_momentum[3], small_eps)));
                double momentum_miln[4] = { mt*cosh(Yi-zk), ptc_momentum[1],
                                            ptc_momentum[2], mt*sinh(Yi-zk)};
                
                //if(gridi == 70 && gridj == 7 && gridk == 41 && iptc==878)
               


                double factorP_bmu = ptc_baryon;
                Taugrid[gridi][gridj][gridk].Jbmu[0] += ptc_baryon*delta;

                for(int mu = 0 ; mu < 4 ; mu ++){
                 
                  Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)] += momentum_miln[mu]*delta;
                  vv = std::max(vv,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]);
                }

    
              }     

            }

            

      }
      }
      }
      std::cout <<" max value "<< vv <<std::endl;

      coutevent+=1;


    }
    std::cout<<"NEVENT: "<< coutevent<<std::endl;
    double maxvale2=0.0;
    
     
    for(int gridi = 0; gridi < NX; gridi++ ){
      //std::cout<< gridi <<std::endl;
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NETA; gridk++ ){
          double maxvale=0.0;
          Taugrid[gridi][gridj][gridk].Jbmu[0]  /= coutevent;
          for(int mu = 0 ; mu < 4 ; mu ++){
                    Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]  /= coutevent;
                    
                    maxvale= std::max(maxvale,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]);
                    maxvale2 = std::max(maxvale2,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]);

                }
          if( maxvale < 1e-4) continue;
          
         // std::cout<<maxvale<<std::endl;

         

          double K2 = Taugrid[gridi][gridj][gridk].Tmnnu[idx(1,0)]*Taugrid[gridi][gridj][gridk].Tmnnu[idx(1,0)]
                    + Taugrid[gridi][gridj][gridk].Tmnnu[idx(2,0)]*Taugrid[gridi][gridj][gridk].Tmnnu[idx(2,0)]
                    + Taugrid[gridi][gridj][gridk].Tmnnu[idx(3,0)]*Taugrid[gridi][gridj][gridk].Tmnnu[idx(3,0)];
          double K0 = Taugrid[gridi][gridj][gridk].Tmnnu[idx(0,0)];
          double M = sqrt(K2);

          double J0 = Taugrid[gridi][gridj][gridk].Jbmu[0];

          rootFinding_newton(K0, M, J0,Taugrid[gridi][gridj][gridk].ed, Taugrid[gridi][gridj][gridk].nb); 
          
          
          
          double EPV = std::max(small_eps, K0+eos.get_pressure(Taugrid[gridi][gridj][gridk].ed/HBARC,  Taugrid[gridi][gridj][gridk].nb)*HBARC); 
          double vx = Taugrid[gridi][gridj][gridk].Tmnnu[idx(1,0)]/EPV;
          double vy = Taugrid[gridi][gridj][gridk].Tmnnu[idx(2,0)]/EPV;
          double vz = Taugrid[gridi][gridj][gridk].Tmnnu[idx(3,0)]/EPV;
          
          double gamma = 1.0/sqrt(std::max(1.0-vx*vx-vy*vy-vz*vz, small_eps));
          Taugrid[gridi][gridj][gridk].u0 = gamma;
          Taugrid[gridi][gridj][gridk].ux = gamma*vx;
          Taugrid[gridi][gridj][gridk].uy = gamma*vy;
          Taugrid[gridi][gridj][gridk].uz = gamma*vz;

          if(gridi == 95 && gridj == 80 && gridk == 100)
          {
            std::cout<<Taugrid[gridi][gridj][gridk].Tmnnu[idx(0,0)]<<" "
                     <<Taugrid[gridi][gridj][gridk].Tmnnu[idx(1,0)]<<" "
                     <<Taugrid[gridi][gridj][gridk].Tmnnu[idx(2,0)]<<" "
                     <<Taugrid[gridi][gridj][gridk].Tmnnu[idx(3,0)]<<" "
                     <<Taugrid[gridi][gridj][gridk].Jbmu[0]<<std::endl;
          }

        
        }
      }

    }
    
    if(maxvale2 > 0){

    
    std::string filename = Parameter::PATHOUT  + "SMASH_ini.dat";   
    std::cout<<" "<<filename<<std::endl; 
    std::ofstream fout(filename);


    for(int gridi = 0; gridi < NX; gridi++ ){
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NETA; gridk++ ){
          fout<< Taugrid[gridi][gridj][gridk].ed //GeV/fm^3
              <<" "<< Taugrid[gridi][gridj][gridk].nb //1/fm^3
              <<" "<< Taugrid[gridi][gridj][gridk].u0
              <<" "<< Taugrid[gridi][gridj][gridk].ux
              <<" "<< Taugrid[gridi][gridj][gridk].uy
              <<" "<< Taugrid[gridi][gridj][gridk].uz/TAU0<<std::endl;     
        }
      }

    }

    fout.close();
    }
 

    
    std::cout<<maxvale2<<" www"<<std::endl;

    

  
  

}

void Grid::GausssmearingTauEta2(TEPaticlelist& Nptclist){
   double one_o_2sigr2 = 1.0/(2.0*SIGR*SIGR); 
   double one_o_2sigz2 = 1.0/(2.0*SIGETA*SIGETA);

   double w1 = one_o_2sigr2/M_PI_F;
   double w2 = sqrt(one_o_2sigz2/M_PI_F)/TAU0;
  

    int coutevent = 0;


    GridClear(Taugrid,NETA);
    std::cout<<Nptclist.size()<<std::endl;
    for(int ievent=0; ievent< Nptclist.size(); ievent++) 
    //for(int ievent=0; ievent< 2; ievent++) 
    { 
      
      

      std::cout<<ievent<<" "<<Nptclist[ievent].size()<<std::endl;
      double vv = 0;
      
      for(int gridi = 0; gridi < NX; gridi++ ){
      double xi= (gridi-NX/2)*DX;
      for(int gridj = 0; gridj < NY; gridj++ ){
          double yj = (gridj-NY/2)*DY;
          for(int gridk = 0; gridk < NETA; gridk++ ){
            double zk = (gridk-NETA/2)*DETA;
            int nncoll=0;
            double check_nb0=0;
            for(int iptc=0;iptc < Nptclist[ievent].size(); iptc++){
              double ptc_t = Nptclist[ievent][iptc].t;
              double ptc_x = Nptclist[ievent][iptc].x;
              double ptc_y = Nptclist[ievent][iptc].y;
              double ptc_z = Nptclist[ievent][iptc].z;
              double ptc_e = Nptclist[ievent][iptc].e;
              double ptc_px = Nptclist[ievent][iptc].px;
              double ptc_py = Nptclist[ievent][iptc].py;
              double ptc_pz = Nptclist[ievent][iptc].pz;
              double ptc_mass = Nptclist[ievent][iptc].mass;
              double ptc_baryon = Nptclist[ievent][iptc].baryon_number;
              double ptc_momentum[4] = {ptc_e,ptc_px,ptc_py,ptc_pz};
              double ptc_position[4] = {ptc_t,ptc_x,ptc_y,ptc_z};
              int ncoll = Nptclist[ievent][iptc].ncoll;
              
              double ptc_tausq = ptc_t*ptc_t - ptc_z*ptc_z;
              //std::cout << sqrt(ptc_tausq)<<std::endl;
              if(ncoll == 0 || ptc_tausq < 0) continue;
              nncoll++;
              
              
              double etasi = 0.5f * (log(std::max(ptc_position[0]+ptc_position[3], small_eps)) \
                           - log(std::max(ptc_position[0]-ptc_position[3], small_eps)));
          
              double dxx = ptc_x - xi;
              double dyy = ptc_y - yj;
              double dzz = etasi - zk;
             

              double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);


              if ( fabs(dzz) < 10*SIGETA && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
                 
              
                double delta = exp(- distance_sqr)*w1*w2;
                
                double mt = sqrt(ptc_momentum[1] * ptc_momentum[1] + ptc_momentum[2] * ptc_momentum[2]+ptc_mass*ptc_mass);
                double Yi = 0.5f * (log(std::max(ptc_momentum[0]+ptc_momentum[3], small_eps)) \
                           - log(std::max(ptc_momentum[0]-ptc_momentum[3], small_eps)));
                double momentum_miln[4] = { mt*cosh(Yi-zk), ptc_momentum[1],
                                            ptc_momentum[2], mt*sinh(Yi-zk)};
                               


                double factorP_bmu = ptc_baryon;
                
                
                for(int mu = 0 ; mu < 4 ; mu ++){
                  Taugrid[gridi][gridj][gridk].Jbmu[mu] += momentum_miln[mu]*ptc_baryon*delta/(momentum_miln[0]+1e-7);
                  for(int nu = 0; nu <4 ; nu++){
                    double factorPP = momentum_miln[mu]*momentum_miln[nu]/(momentum_miln[0]+1e-7);

                    Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)] += factorPP*delta;
                    vv = std::max(vv,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]);

                  }
                  
                }

    
              }     

            }

      }
      }
      }
      std::cout <<" max value "<< vv <<std::endl;

      coutevent+=1;


    }
    std::cout<<"NEVENT: "<< coutevent<<std::endl;
    double maxvale2=0.0;
    
    double check_nb = 0.0;
    for(int gridi = 0; gridi < NX; gridi++ ){
      //std::cout<< gridi <<std::endl;
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NETA; gridk++ ){
          double maxvale=0.0;
         
          for(int mu = 0 ; mu < 4 ; mu ++){
            Taugrid[gridi][gridj][gridk].Jbmu[mu]  /= coutevent;
            if (mu==0) check_nb+= Taugrid[gridi][gridj][gridk].Jbmu[mu];
             for(int nu =0; nu < 4 ; nu ++){

                    Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]  /= coutevent;
                    
                    maxvale= std::max(maxvale,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);
                    maxvale2 = std::max(maxvale2,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);

             }

                    

          }
          if( maxvale < 1e-4) continue;
          
         // std::cout<<maxvale<<std::endl;
          Laudumatching(Taugrid[gridi][gridj][gridk].Tmnnu,Taugrid[gridi][gridj][gridk].ed,
                        Taugrid[gridi][gridj][gridk].u0,Taugrid[gridi][gridj][gridk].ux,
                        Taugrid[gridi][gridj][gridk].uy,Taugrid[gridi][gridj][gridk].uz);

          Taugrid[gridi][gridj][gridk].nb = Taugrid[gridi][gridj][gridk].Jbmu[0]*Taugrid[gridi][gridj][gridk].u0
                                        - Taugrid[gridi][gridj][gridk].Jbmu[1]*Taugrid[gridi][gridj][gridk].ux
                                        - Taugrid[gridi][gridj][gridk].Jbmu[2]*Taugrid[gridi][gridj][gridk].uy
                                        - Taugrid[gridi][gridj][gridk].Jbmu[3]*Taugrid[gridi][gridj][gridk].uz;
          
          Taugrid[gridi][gridj][gridk].Temperature = eos.get_temperature(Taugrid[gridi][gridj][gridk].ed/HBARC,Taugrid[gridi][gridj][gridk].nb)*HBARC;    

          if(gridi == 95 && gridj == 80 && gridk == 100)
          {
            for(int mu=0;mu<4;mu++)
            {
              for(int nu=0;nu<4;nu++){
                 std::cout<<Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]<<" ";
                
              }
            }
            std::cout<<std::endl;
            std::cout<<Taugrid[gridi][gridj][gridk].ed<<" "
                     <<Taugrid[gridi][gridj][gridk].nb<<" "
                     <<Taugrid[gridi][gridj][gridk].ux<<" "
                     <<Taugrid[gridi][gridj][gridk].uy<<" "
                     <<Taugrid[gridi][gridj][gridk].uz<<" "<<check_nb*0.15*0.15*0.15*3.2<<std::endl;
          }

        
        }
      }

    }
    std::cout<<check_nb*0.15*0.15*0.15*3.2<<std::endl;
    
    if(maxvale2 > 0){

    
    std::string filename = Parameter::PATHOUT  + "SMASH_ini.dat";   
    std::cout<<" "<<filename<<std::endl; 
    std::ofstream fout(filename);


    for(int gridi = 0; gridi < NX; gridi++ ){
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NETA; gridk++ ){
          fout<< Taugrid[gridi][gridj][gridk].ed //GeV/fm^3
              <<" "<< Taugrid[gridi][gridj][gridk].nb //1/fm^3
              <<" "<< Taugrid[gridi][gridj][gridk].u0
              <<" "<< Taugrid[gridi][gridj][gridk].ux
              <<" "<< Taugrid[gridi][gridj][gridk].uy
              <<" "<< Taugrid[gridi][gridj][gridk].uz/TAU0<<std::endl;     
        }
      }

    }

    fout.close();
    }
 

    
    std::cout<<maxvale2<<" www"<<std::endl;

    

  
  

}

Grid::~Grid(){
  
}





#endif
