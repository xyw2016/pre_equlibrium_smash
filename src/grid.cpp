#ifndef GRID_CPP
#define GRID_CPP

#include "grid.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#ifndef _OPENMP
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif


Grid::Grid(int use_eos):NX(Parameter::NX), NY(Parameter::NY), NETA(Parameter::NETA),
                        DX(Parameter::DX), DY(Parameter::DY),DETA(Parameter::DETA),
                        SIGR(Parameter::SIGR), SIGETA(Parameter::SIGETA), TAU0(Parameter::TAU0),
                        TAU00(Parameter::TAU00),NTAU(Parameter::NTAU),eos(use_eos){
  
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
    cellzero.muB = 0;
    cellzero.cs2 = 0;
    
    for(int id = 0 ; id < 16 ; id ++)
    {
      cellzero.Tmnnu[id] = 0;
    }  

    for(int id = 0 ; id < 4 ; id ++)
    {
      cellzero.Jbmu[id] = 0;
    }       


    Tgrid.resize(NX);
   
  
    for (int i = 0; i < NX; i++) {
        Tgrid[i].resize(NY);     
        for (int j = 0; j < NY; j++) {
            Tgrid[i][j].resize(NETA);
            for(int k = 0; k < NETA; k++){
              Tgrid[i][j][k] = cellzero;
            }
        }
        
    }

    CORES = 1;

    #ifdef _OPENMP
      CORES = omp_get_max_threads();
    #endif

std::cout<<" use the number of threads: " <<CORES<<std::endl;


}

void Grid::GridClear(Grid3D& grid0){

  for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
          for (int k =0 ; k < NETA ;k++){
            grid0[i][j][k].ed = 0;
            grid0[i][j][k].nb = 0;
            grid0[i][j][k].u0 = 1;
            grid0[i][j][k].ux = 0;
            grid0[i][j][k].uy = 0;
            grid0[i][j][k].uz = 0;
            grid0[i][j][k].Temperature = 0;
            grid0[i][j][k].pressure = 0;
            grid0[i][j][k].muB = 0;
            grid0[i][j][k].cs2 = 0;
            
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

void Grid::smearing_kernel(Particlelist& ptclist){
        
        #pragma omp parallel for
        for(int iptc=0;iptc < ptclist.size(); iptc++){
        
        double norm_facor = 0.0;
        double ptc_t = ptclist[iptc].t;
        double ptc_x = ptclist[iptc].x;
        double ptc_y = ptclist[iptc].y;
        double ptc_z = ptclist[iptc].z;
        double ptc_e = ptclist[iptc].e;
        double ptc_px = ptclist[iptc].px;
        double ptc_py = ptclist[iptc].py;
        double ptc_pz = ptclist[iptc].pz;
        double ptc_mass = ptclist[iptc].mass;
        double ptc_baryon = ptclist[iptc].baryon_number;

        double ptc_momentum[4] = {ptc_e,ptc_px,ptc_py,ptc_pz};
        double ptc_position[4] = {ptc_t,ptc_x,ptc_y,ptc_z};
        int ncoll = ptclist[iptc].ncoll;

        double ptc_tausq = ptc_t*ptc_t - ptc_z*ptc_z;
        if(ncoll == 0 || ptc_tausq < 0) continue;

        one_o_2sigr2 = 1.0/(2.0*SIGR*SIGR); 
        one_o_2sigz2 = 1.0/(2.0*SIGETA*SIGETA);

        w1 = one_o_2sigr2/M_PI_F;
        w2 = sqrt(one_o_2sigz2/M_PI_F)/sqrt(ptc_tausq);
        
        double etasi = 0.5f * (log(std::max(ptc_position[0]+ptc_position[3], small_eps)) \
                           - log(std::max(ptc_position[0]-ptc_position[3], small_eps)));
        
        #pragma omp parallel for collapse(3) reduction(+:norm_facor)
         for(int gridi = 0; gridi < NX; gridi++ ){
          for(int gridj = 0; gridj < NY; gridj++ ){
              for(int gridk = 0; gridk < NETA; gridk++ ){
                double xi= (gridi-NX/2)*DX;
                double yj = (gridj-NY/2)*DY;
                double zk = (gridk-NETA/2)*DETA;         
                double dxx = ptc_x - xi;
                double dyy = ptc_y - yj;
                double dzz = etasi - zk;
             

                double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);


                if ( fabs(dzz) < 10*SIGETA && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
                  double delta = exp(- distance_sqr)*w1*w2;
                  norm_facor += delta;
                }
              }
          }
         }

        norm_facor =  norm_facor*DX*DY*DETA*TAU0;
        

     

        #pragma omp parallel for collapse(3) 
        for(int gridi = 0; gridi < NX; gridi++ ){
          for(int gridj = 0; gridj < NY; gridj++ ){
              for(int gridk = 0; gridk < NETA; gridk++ ){
                double xi= (gridi-NX/2)*DX;
                double yj = (gridj-NY/2)*DY;
                double zk = (gridk-NETA/2)*DETA;         
                double dxx = ptc_x - xi;
                double dyy = ptc_y - yj;
                double dzz = etasi - zk;
             

                double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);


                if ( fabs(dzz) < 10*SIGETA && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
                 
              
                  double delta = exp(- distance_sqr)*w1*w2/norm_facor;
                  double mt = sqrt(ptc_momentum[1] * ptc_momentum[1] + ptc_momentum[2] * ptc_momentum[2]+ptc_mass*ptc_mass);
                  double Yi = 0.5f * (log(std::max(ptc_momentum[0]+ptc_momentum[3], small_eps)) \
                      - log(std::max(ptc_momentum[0]-ptc_momentum[3], small_eps)));
                  double momentum_miln[4] = { mt*cosh(Yi-zk), ptc_momentum[1],
                                            ptc_momentum[2], mt*sinh(Yi-zk)};
                               


                  double factorP_bmu = ptc_baryon;
                
                
                  for(int mu = 0 ; mu < 4 ; mu ++){
                    double jmu_temp = momentum_miln[mu]*ptc_baryon*delta/(momentum_miln[0]+1e-7);
                    #pragma omp atomic
                    Tgrid[gridi][gridj][gridk].Jbmu[mu] += jmu_temp;
                     for(int nu = 0; nu <4 ; nu++){
                        double factorPP = momentum_miln[mu]*momentum_miln[nu]/(momentum_miln[0]+1e-7);
                        double Tmn_temp = factorPP*delta;
                        #pragma omp atomic
                        Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)] += Tmn_temp;
                  

                    }
                  
                  }

    
                }     

            

              }
          }
        }


      } //end ptc loops in one event



}

void Grid::perform_laudu(double coutevent){
    double maxvale2=0.0;
    double check_nb = 0.0;
    for(int gridi = 0; gridi < NX; gridi++ ){
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NETA; gridk++ ){
          double maxvale=0.0;
         
          for(int mu = 0 ; mu < 4 ; mu ++){
            Tgrid[gridi][gridj][gridk].Jbmu[mu]  /= coutevent;
            if (mu==0) check_nb+= Tgrid[gridi][gridj][gridk].Jbmu[mu];
             for(int nu =0; nu < 4 ; nu ++){

                    Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]  /= coutevent;
                    
                    maxvale= std::max(maxvale,Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);
                    maxvale2 = std::max(maxvale2,Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);

             }

                    

          }
          if( maxvale < 1e-4) continue;
          
         // std::cout<<maxvale<<std::endl;
          Laudumatching(Tgrid[gridi][gridj][gridk].Tmnnu,Tgrid[gridi][gridj][gridk].ed,
                        Tgrid[gridi][gridj][gridk].u0,Tgrid[gridi][gridj][gridk].ux,
                        Tgrid[gridi][gridj][gridk].uy,Tgrid[gridi][gridj][gridk].uz);

          Tgrid[gridi][gridj][gridk].nb = Tgrid[gridi][gridj][gridk].Jbmu[0]*Tgrid[gridi][gridj][gridk].u0
                                        - Tgrid[gridi][gridj][gridk].Jbmu[1]*Tgrid[gridi][gridj][gridk].ux
                                        - Tgrid[gridi][gridj][gridk].Jbmu[2]*Tgrid[gridi][gridj][gridk].uy
                                        - Tgrid[gridi][gridj][gridk].Jbmu[3]*Tgrid[gridi][gridj][gridk].uz;
          
          Tgrid[gridi][gridj][gridk].Temperature = eos.get_temperature(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb)*HBARC;  

          Tgrid[gridi][gridj][gridk].pressure = eos.get_pressure(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb)*HBARC;  

          Tgrid[gridi][gridj][gridk].muB = eos.get_muB(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb)*HBARC;

          Tgrid[gridi][gridj][gridk].cs2 = eos.get_cs2(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb);   

        
        }
      }

    }
    //std::cout<<check_nb*DX*DY*DETA*TAU0<<std::endl;
    std::cout<<maxvale2<<" www"<<std::endl;
}
void Grid::preequlibirum(Paticlelist_event& ptcl_event){


    
   double dtau = (TAU0-TAU00)/NTAU;
   
  std::string out_open_mode;

  std::string filename = Parameter::PATHOUT  + "pre_equlibrium_evo_txyz.dat";   
  std::cout<<" "<<filename<<std::endl; 
 
  for(int itime = 0 ; itime < NTAU; itime+=1)
  {

    std::string out_open_mode;

    if (itime == 0) {
        out_open_mode = "wb";
    } else {
        out_open_mode = "ab";
    }


    FILE *out_file_xyz;
    out_file_xyz = fopen(filename.c_str(), out_open_mode.c_str());
    


    int coutevent = 0;
    GridClear(Tgrid);
    bool flag_out_grid = 0;
    int out = 0;
    double total_energy_ptc = 0.0;
    double total_px_ptc = 0.0;
    double total_py_ptc = 0.0;
    double total_pz_ptc = 0.0;
    double total_baryon_ptc = 0.0;
    std::cout<<"t = "<<(itime+1)*dtau<<" fm"<<std::endl;
    
    for(int ievent=0; ievent< ptcl_event.size(); ievent++) 
    { 
      
      // std::cout<<"t = "<<(itime+1)*dtau<<" fm"<<" "
      //          <<"EVENT ID: "<<ievent<<" "
      //          <<"Number of particle: "<<ptcl_event[ievent][itime].size()<< std::endl;

      smearing_kernel(ptcl_event[ievent][itime]);
          
      double total_energy_grid = 0.0;
      double total_px_grid = 0.0;
      double total_py_grid = 0.0;
      double total_pz_grid = 0.0;
      
      double total_baryon_grid = 0.0;
      
      for(int gridi = 0; gridi < NX; gridi++ ){

      for(int gridj = 0; gridj < NY; gridj++ ){
        
          for(int gridk = 0; gridk < NETA; gridk++ ){
            total_energy_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,0)];
            total_px_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,1)];
            total_py_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,2)];
            total_pz_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,3)];
            total_baryon_grid +=Tgrid[gridi][gridj][gridk].Jbmu[0];       
            }     
      }
      }

      
      // std::cout <<"Total energy of grid: "<< total_energy_grid*DX*DY*DETA*(itime+1)*dtau<<std::endl;
      // std::cout <<"Total baryon of grid: "<< total_baryon_grid*DX*DY*DETA*(itime+1)*dtau<<std::endl;

      coutevent+=1;


    }

   
    // std::cout<< coutevent<<std::endl;
    perform_laudu(coutevent);
    
    const int nVar_per_cell = (11 + Parameter::turn_on_rhob *2 + Parameter::turn_on_shear*5
                                  + Parameter::turn_on_bulk*1 + Parameter::turn_on_diff*3);

    const int output_nx        = NX;
    const int output_ny        = NY;
    const int output_nz        = NETA;
    const double output_dx     = DX;
    const double output_dy     = DY;
    const double output_dz     = DETA;
    const double output_xmin   = - (NX-1)*DX/2.;
    const double output_ymin   = - (NY-1)*DY/2.;
    const double output_zmin   = - (NETA-1)*DETA/2.;

    if (itime == 0) {
        float header[] = {
            static_cast<float>(0), static_cast<float>(dtau),
            static_cast<float>(output_nx), static_cast<float>(output_dx),
            static_cast<float>(output_xmin),
            static_cast<float>(output_ny), static_cast<float>(output_dy),
            static_cast<float>(output_ymin),
            static_cast<float>(output_nz), static_cast<float>(output_dz),
            static_cast<float>(output_zmin),
            static_cast<float>(Parameter::turn_on_rhob),
            static_cast<float>(Parameter::turn_on_shear),
            static_cast<float>(Parameter::turn_on_bulk),
            static_cast<float>(Parameter::turn_on_diff),
            static_cast<float>(nVar_per_cell)};
        fwrite(header, sizeof(float), 16, out_file_xyz);
    }

     for(int gridi = 0; gridi < NX; gridi++ ){
      for(int gridj = 0; gridj < NY; gridj++ ){
        for(int gridk = 0; gridk < NETA; gridk++ ){
          if(Tgrid[gridi][gridj][gridk].ed > Parameter::TCUT)
          {
          
          
          // Tgrid[gridi][gridj][gridk].Temperature=0.2;
          // Tgrid[gridi][gridj][gridk].muB=0.0;
          // Tgrid[gridi][gridj][gridk].nb=0.0;
          // Tgrid[gridi][gridj][gridk].ux=0.0;
          // Tgrid[gridi][gridj][gridk].uy=0.0;
          // Tgrid[gridi][gridj][gridk].uz=0.0;
          
          float ideal[] = {static_cast<float>(itime+1),
                            static_cast<float>(gridi),
                            static_cast<float>(gridj),
                            static_cast<float>(gridk),
                            static_cast<float>(Tgrid[gridi][gridj][gridk].ed),
                            static_cast<float>(Tgrid[gridi][gridj][gridk].pressure),
                            static_cast<float>(Tgrid[gridi][gridj][gridk].Temperature),
                            static_cast<float>(Tgrid[gridi][gridj][gridk].cs2),
                            static_cast<float>(Tgrid[gridi][gridj][gridk].ux),
                            static_cast<float>(Tgrid[gridi][gridj][gridk].uy),
                            static_cast<float>(Tgrid[gridi][gridj][gridk].uz)};
          
          fwrite(ideal, sizeof(float), 11, out_file_xyz);
          
          if (Parameter::turn_on_rhob == 1) {
              float ouput_mu[] = {static_cast<float>(Tgrid[gridi][gridj][gridk].nb),
                                  static_cast<float>(Tgrid[gridi][gridj][gridk].muB)};
              fwrite(ouput_mu, sizeof(float), 2, out_file_xyz);
          }
        }

        }
      }

    }



            

    fclose(out_file_xyz);

 

    
   

    

  }
   



}

void Grid::hydro_ini(Paticlelist_event& ptcl_event){
    int select_ptc = NTAU - 1;
    GridClear(Tgrid);
    int coutevent = 0;
    std::cout<<"Total events: "<< ptcl_event.size()<<std::endl;
    
    for(int ievent=0; ievent< ptcl_event.size(); ievent++) 
    { 
      
      std::cout<<" start to smearing, event "<< ievent
               <<" ptc: "<<ptcl_event[ievent][select_ptc].size()<<std::endl;
      smearing_kernel(ptcl_event[ievent][select_ptc]);
      coutevent+=1;
    }
    std::cout<<"Laudu matching.... "<< coutevent<<std::endl;

    perform_laudu(coutevent);

    
    // std::string filename = Parameter::PATHOUT  + "SMASH_ini.dat";   
    // // std::ofstream fout(filename,"wb");

    // FILE *fout;
    // fout = fopen(filename.c_str(), "w");


    // for(int gridi = 0; gridi < NX; gridi++ ){
    //   for(int gridj = 0; gridj < NY; gridj++ ){
    //     for(int gridk = 0; gridk < NETA; gridk++ ){


    //       float ideal[] = {static_cast<float>(Taugrid[gridi][gridj][gridk].ed),
    //                         static_cast<float>(Taugrid[gridi][gridj][gridk].nb),
    //                         static_cast<float>(Taugrid[gridi][gridj][gridk].u0),
    //                         static_cast<float>(Taugrid[gridi][gridj][gridk].ux),
    //                         static_cast<float>(Taugrid[gridi][gridj][gridk].uy),
    //                         static_cast<float>(Taugrid[gridi][gridj][gridk].uz/TAU0)};
    //       fwrite(ideal, sizeof(float), 6, fout);


    //     }
    //   }

    // }

    // fclose(fout);




    std::string filename = Parameter::PATHOUT + "SMASH_ini.dat";
    FILE *fout;
    fout = fopen(filename.c_str(), "w");

    if(fout == nullptr) {
    // Handle file open error
    }

    for(int gridi = 0; gridi < NX; gridi++) {
      for(int gridj = 0; gridj < NY; gridj++) {
        for(int gridk = 0; gridk < NETA; gridk++) {

            // Format the output to write in text format
            fprintf(fout, "%f %f %f %f %f %f\n",
                    static_cast<float>(Tgrid[gridi][gridj][gridk].ed),
                    static_cast<float>(Tgrid[gridi][gridj][gridk].nb),
                    static_cast<float>(Tgrid[gridi][gridj][gridk].u0),
                    static_cast<float>(Tgrid[gridi][gridj][gridk].ux),
                    static_cast<float>(Tgrid[gridi][gridj][gridk].uy),
                    static_cast<float>(Tgrid[gridi][gridj][gridk].uz/TAU0));
        }
      }
    }

    fclose(fout);

}


// void Grid::GausssmearingTZ(EventsofTEPaticlelist& Nptclist){
//    double one_o_2sigr2 = 1.0/(2.0*SIGR*SIGR); 
//    double one_o_2sigz2 = 1.0/(2.0*SIGZ*SIGZ);

//    double w1 = one_o_2sigr2/M_PI_F;
//    double w2 = sqrt(one_o_2sigz2/M_PI_F);
   
//    for(int i = 0 ; i < Nptclist.size(); i++){
//       int ntime = Nptclist[i].size();
//       nmaxtime = std::max(nmaxtime,ntime);
//    }
   
//   std::string out_open_mode;

//   std::string filename = Parameter::PATHOUT  + "pre_equlibrium_evo_txyz.dat";   
//   std::cout<<" "<<filename<<std::endl; 
 
//   int itstart = 0;
//   int itend = 100;
//   //for(int itime =0 ; itime < nmaxtime; itime++)
//   for(int itime =itstart ; itime < itend; itime+=1)
//   {

//     std::string out_open_mode;

//     if (itime == itstart) {
//         out_open_mode = "wb";
//     } else {
//         out_open_mode = "ab";
//     }


//     FILE *out_file_xyz;
//     out_file_xyz = fopen(filename.c_str(), out_open_mode.c_str());
    


//     int coutevent = 0;
//     GridClear(Tgrid,NZ);
//     bool flag_out_grid = 0;
//     int out = 0;
//     double total_energy_ptc = 0.0;
//     double total_px_ptc = 0.0;
//     double total_py_ptc = 0.0;
//     double total_pz_ptc = 0.0;
//     double total_baryon_ptc = 0.0;
//     for(int ievent=0; ievent< Nptclist.size(); ievent++) 
//     { 
      
      
//       if(itime >= Nptclist[ievent].size() ) {
//         out += 1;
//         continue;
//         }

   
      
//       double ptc_xmin = Parameter::NX/2* Parameter::DX; 
//       double ptc_ymin = Parameter::NY/2* Parameter::DY; 
//       double ptc_zmin = Parameter::NZ/2* Parameter::DZ; 
//       for(int iptc=0;iptc < Nptclist[ievent][itime].size(); iptc++){
              
//             ptc_xmin = std::min(ptc_xmin,fabs(Nptclist[ievent][itime][iptc].x));
//             ptc_ymin = std::min(ptc_ymin,fabs(Nptclist[ievent][itime][iptc].y));
//             ptc_zmin = std::min(ptc_zmin,fabs(Nptclist[ievent][itime][iptc].z));
//       }

//        if (fabs(ptc_xmin) >= Parameter::NX/2* Parameter::DX ||
//            fabs(ptc_ymin) >= Parameter::NY/2* Parameter::DY ||
//            fabs(ptc_zmin) >= Parameter::NZ/2* Parameter::DZ  ) {
//             flag_out_grid=1;
//             out +=1;
//             continue;
//            }
      

//       std::cout<<"t = "<<itime*DT<<" fm"<<" "
//                <<"EVENT ID: "<<ievent<<" "
//                <<"Zmin: "<<ptc_zmin<<" "
//                <<"Number of particle: "<<Nptclist[ievent][itime].size()<< std::endl;


     
//       for(int iptc=0;iptc < Nptclist[ievent][itime].size(); iptc++){
      
//       int nncoll=0;
//       double norm_facor = 0.0;
//       double ptc_t = Nptclist[ievent][itime][iptc].t;
//       double ptc_x = Nptclist[ievent][itime][iptc].x;
//       double ptc_y = Nptclist[ievent][itime][iptc].y;
//       double ptc_z = Nptclist[ievent][itime][iptc].z;
//       double ptc_e = Nptclist[ievent][itime][iptc].e;
//       double ptc_px = Nptclist[ievent][itime][iptc].px;
//       double ptc_py = Nptclist[ievent][itime][iptc].py;
//       double ptc_pz = Nptclist[ievent][itime][iptc].pz;
//       double ptc_baryon = Nptclist[ievent][itime][iptc].baryon_number;
//       double ptc_p4[4] = {ptc_e,ptc_px,ptc_py,ptc_pz};
//       int ncoll = Nptclist[ievent][itime][iptc].ncoll;
//       double ptc_tausq = ptc_t*ptc_t - ptc_z*ptc_z;
//       if(ncoll == 0 || ptc_tausq < 0) continue;
      
//       double local_time = itime*DT;
//       total_energy_ptc = total_energy_ptc + ptc_e;
//       total_px_ptc = total_px_ptc + ptc_px;
//       total_py_ptc = total_py_ptc + ptc_py;
//       total_pz_ptc = total_pz_ptc + ptc_pz;
//       total_baryon_ptc = total_baryon_ptc + ptc_baryon;

//       for(int gridi = 0; gridi < NX; gridi++ ){
//       double xi= (gridi-NX/2)*DX;
//       for(int gridj = 0; gridj < NY; gridj++ ){
//           double yj = (gridj-NY/2)*DY;
//           for(int gridk = 0; gridk < NZ; gridk++ ){
//             double zk = (gridk-NZ/2)*DZ;
//               if (fabs(zk) > local_time ) continue;
//               double dxx = ptc_x - xi;
//               double dyy = ptc_y - yj;
//               double dzz = ptc_z - zk;
              
//               double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);

//               if ( fabs(dzz) < 10*SIGZ && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
            
//                 double delta = exp(- distance_sqr)*w1*w2;
//                 norm_facor += delta;
//               }     

//             }

            
//       }
//       }
//       norm_facor =  norm_facor*DX*DY*DZ;
//       //std::cout<<" nrom " << norm_facor<<std::endl;


//       for(int gridi = 0; gridi < NX; gridi++ ){
//       double xi= (gridi-NX/2)*DX;
//       for(int gridj = 0; gridj < NY; gridj++ ){
//           double yj = (gridj-NY/2)*DY;
//           for(int gridk = 0; gridk < NZ; gridk++ ){
//             double zk = (gridk-NZ/2)*DZ;
//               if (fabs(zk) > local_time ) continue;
//               double dxx = ptc_x - xi;
//               double dyy = ptc_y - yj;
//               double dzz = ptc_z - zk;
              
//               double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);

//               if ( fabs(dzz) < 10*SIGZ && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
            
//                 double delta = exp(- distance_sqr)*w1*w2/norm_facor;
//                 for(int mu = 0 ; mu < 4 ; mu ++){
//                   double factorP_bmu = ptc_baryon*ptc_p4[mu]/(ptc_p4[0]+1e-7);
//                   Tgrid[gridi][gridj][gridk].Jbmu[mu] += factorP_bmu*delta;
//                   for(int nu =0; nu < 4 ; nu ++){
//                     double factorPP = ptc_p4[mu]*ptc_p4[nu]/(ptc_p4[0]+1e-7);
                    
//                     Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)] += factorPP*delta;

//                   }
//                 }
//               }     

//             }

            
//       }
//       }

//       }
      
//       double total_energy_grid = 0.0;
//       double total_px_grid = 0.0;
//       double total_py_grid = 0.0;
//       double total_pz_grid = 0.0;
      
//       double total_baryon_grid = 0.0;
      
//       for(int gridi = 0; gridi < NX; gridi++ ){

//       for(int gridj = 0; gridj < NY; gridj++ ){
        
//           for(int gridk = 0; gridk < NZ; gridk++ ){
//             total_energy_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,0)];
//             total_px_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,1)];
//             total_py_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,2)];
//             total_pz_grid +=Tgrid[gridi][gridj][gridk].Tmnnu[idx(0,3)];
//             total_baryon_grid +=Tgrid[gridi][gridj][gridk].Jbmu[0];       
//             }     
//       }
//       }

//       std::cout <<"Total energy of ptc: " << total_energy_ptc<<"\t";
//       std::cout <<"Total energy of grid: "<< total_energy_grid*DX*DY*DZ<<std::endl;

//       // std::cout <<"Total px of ptc: " << total_px_ptc<<"\t";
//       // std::cout <<"Total px of grid: "<< total_px_grid*DX*DY*DZ<<std::endl;

//       // std::cout <<"Total py of ptc: " << total_py_ptc<<"\t";
//       // std::cout <<"Total py of grid: "<< total_py_grid*DX*DY*DZ<<std::endl;

//       // std::cout <<"Total pz of ptc: " << total_pz_ptc<<"\t";
//       // std::cout <<"Total pz of grid: "<< total_pz_grid*DX*DY*DZ<<std::endl;

//       std::cout <<"Total baryon of ptc: " << total_baryon_ptc<<"\t";
//       std::cout <<"Total baryon of grid: "<< total_baryon_grid*DX*DY*DZ<<std::endl;




//       coutevent+=1;


//     }

//     std::cout<<" Out grid event: "<< out<<std::endl;

//     if(out == Nptclist.size()){
//       fclose(out_file_xyz);
//       break;
//     }

   
//     std::cout<< coutevent<<std::endl;
//     double maxed=0.0;
//     double maxnb=0.0;
//     double maxT=0.0;

//     for(int gridi = 0; gridi < NX; gridi++ ){
//       for(int gridj = 0; gridj < NY; gridj++ ){
//         for(int gridk = 0; gridk < NZ; gridk++ ){
//           double maxvale=0.0;
//           for(int mu = 0 ; mu < 4 ; mu ++){
//                   Tgrid[gridi][gridj][gridk].Jbmu[mu]  /= coutevent;
//                   for(int nu =0; nu < 4 ; nu ++){
//                     Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]  /= coutevent;
                    
//                     maxvale= std::max(maxvale,Tgrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);
                    

//                   }
//                 }
//           if( maxvale < 1e-4) continue;
          
          
//           Laudumatching(Tgrid[gridi][gridj][gridk].Tmnnu,Tgrid[gridi][gridj][gridk].ed,
//                         Tgrid[gridi][gridj][gridk].u0,Tgrid[gridi][gridj][gridk].ux,
//                         Tgrid[gridi][gridj][gridk].uy,Tgrid[gridi][gridj][gridk].uz);

//           Tgrid[gridi][gridj][gridk].nb = Tgrid[gridi][gridj][gridk].Jbmu[0]*Tgrid[gridi][gridj][gridk].u0
//                                         - Tgrid[gridi][gridj][gridk].Jbmu[1]*Tgrid[gridi][gridj][gridk].ux
//                                         - Tgrid[gridi][gridj][gridk].Jbmu[2]*Tgrid[gridi][gridj][gridk].uy
//                                         - Tgrid[gridi][gridj][gridk].Jbmu[3]*Tgrid[gridi][gridj][gridk].uz;
         
//           Tgrid[gridi][gridj][gridk].Temperature = eos.get_temperature(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb)*HBARC;

//           Tgrid[gridi][gridj][gridk].pressure = eos.get_pressure(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb)*HBARC;  

//           Tgrid[gridi][gridj][gridk].muB = eos.get_muB(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb)*HBARC;

//           Tgrid[gridi][gridj][gridk].cs2 = eos.get_cs2(Tgrid[gridi][gridj][gridk].ed/HBARC,Tgrid[gridi][gridj][gridk].nb);  

//           maxed = std::max(maxed,Tgrid[gridi][gridj][gridk].ed);
//           maxnb = std::max(maxnb,Tgrid[gridi][gridj][gridk].nb);
//           maxT = std::max(maxnb,Tgrid[gridi][gridj][gridk].Temperature);

        
//         }
//       }

//     }
    
   
    
    

//     const int nVar_per_cell = (11 + Parameter::turn_on_rhob *2 + Parameter::turn_on_shear*5
//                                   + Parameter::turn_on_bulk*1 + Parameter::turn_on_diff*3);

//     const int output_nx        = NX;
//     const int output_ny        = NY;
//     const int output_nz        = NZ;
//     const double output_dx     = DX;
//     const double output_dy     = DY;
//     const double output_dz     = DZ;
//     const double output_xmin   = - (NX-1)*DX/2.;
//     const double output_ymin   = - (NY-1)*DY/2.;
//     const double output_zmin   = - (NZ-1)*DZ/2.;

//     if (itime == itstart) {
//         float header[] = {
//             static_cast<float>(0), static_cast<float>(DT),
//             static_cast<float>(output_nx), static_cast<float>(output_dx),
//             static_cast<float>(output_xmin),
//             static_cast<float>(output_ny), static_cast<float>(output_dy),
//             static_cast<float>(output_ymin),
//             static_cast<float>(output_nz), static_cast<float>(output_dz),
//             static_cast<float>(output_zmin),
//             static_cast<float>(Parameter::turn_on_rhob),
//             static_cast<float>(Parameter::turn_on_shear),
//             static_cast<float>(Parameter::turn_on_bulk),
//             static_cast<float>(Parameter::turn_on_diff),
//             static_cast<float>(nVar_per_cell)};
//         fwrite(header, sizeof(float), 16, out_file_xyz);
//     }

//      for(int gridi = 0; gridi < NX; gridi++ ){
//       for(int gridj = 0; gridj < NY; gridj++ ){
//         for(int gridk = 0; gridk < NZ; gridk++ ){
//           if(Tgrid[gridi][gridj][gridk].ed > Parameter::TCUT)
//           {
          
          
//           // Tgrid[gridi][gridj][gridk].Temperature=0.2;
//           // Tgrid[gridi][gridj][gridk].muB=0.0;
//           // Tgrid[gridi][gridj][gridk].nb=0.0;
//           // Tgrid[gridi][gridj][gridk].ux=0.0;
//           // Tgrid[gridi][gridj][gridk].uy=0.0;
//           // Tgrid[gridi][gridj][gridk].uz=0.0;
          
//           float ideal[] = {static_cast<float>(itime),
//                             static_cast<float>(gridi),
//                             static_cast<float>(gridj),
//                             static_cast<float>(gridk),
//                             static_cast<float>(Tgrid[gridi][gridj][gridk].ed),
//                             static_cast<float>(Tgrid[gridi][gridj][gridk].pressure),
//                             static_cast<float>(Tgrid[gridi][gridj][gridk].Temperature),
//                             static_cast<float>(Tgrid[gridi][gridj][gridk].cs2),
//                             static_cast<float>(Tgrid[gridi][gridj][gridk].ux),
//                             static_cast<float>(Tgrid[gridi][gridj][gridk].uy),
//                             static_cast<float>(Tgrid[gridi][gridj][gridk].uz)};
          
//           fwrite(ideal, sizeof(float), 11, out_file_xyz);
          
//           if (Parameter::turn_on_rhob == 1) {
//               float ouput_mu[] = {static_cast<float>(Tgrid[gridi][gridj][gridk].nb),
//                                   static_cast<float>(Tgrid[gridi][gridj][gridk].muB)};
//               fwrite(ouput_mu, sizeof(float), 2, out_file_xyz);
//           }
//         }

//         }
//       }

//     }



            

//     fclose(out_file_xyz);

 

    
//     std::cout<<"Edmax: "<< maxed <<std::endl;
//     std::cout<<"Nbmax: "<< maxnb <<std::endl;

    

//   }
  

// }


// void Grid::GausssmearingTauEta(TEPaticlelist& Nptclist){
//    double one_o_2sigr2 = 1.0/(2.0*SIGR*SIGR); 
//    double one_o_2sigz2 = 1.0/(2.0*SIGETA*SIGETA);

//    double w1 = one_o_2sigr2/M_PI_F;
//    double w2 = sqrt(one_o_2sigz2/M_PI_F)/TAU0;
  

//     int coutevent = 0;


//     GridClear(Taugrid,NETA);
//     std::cout<<Nptclist.size()<<std::endl;
//     for(int ievent=0; ievent< Nptclist.size(); ievent++) 
//     //for(int ievent=0; ievent< 2; ievent++) 
//     { 
      
      

//       std::cout<<ievent<<" "<<Nptclist[ievent].size()<<std::endl;
//       double vv = 0;
//       for(int gridi = 0; gridi < NX; gridi++ ){
//       double xi= (gridi-NX/2)*DX;
//       for(int gridj = 0; gridj < NY; gridj++ ){
//           double yj = (gridj-NY/2)*DY;
//           for(int gridk = 0; gridk < NETA; gridk++ ){
//             double zk = (gridk-NETA/2)*DETA;
//             int nncoll=0;
           
//             for(int iptc=0;iptc < Nptclist[ievent].size(); iptc++){
//               double ptc_t = Nptclist[ievent][iptc].t;
//               double ptc_x = Nptclist[ievent][iptc].x;
//               double ptc_y = Nptclist[ievent][iptc].y;
//               double ptc_z = Nptclist[ievent][iptc].z;
//               double ptc_e = Nptclist[ievent][iptc].e;
//               double ptc_px = Nptclist[ievent][iptc].px;
//               double ptc_py = Nptclist[ievent][iptc].py;
//               double ptc_pz = Nptclist[ievent][iptc].pz;
//               double ptc_mass = Nptclist[ievent][iptc].mass;
//               double ptc_baryon = Nptclist[ievent][iptc].baryon_number;
//               double ptc_momentum[4] = {ptc_e,ptc_px,ptc_py,ptc_pz};
//               double ptc_position[4] = {ptc_t,ptc_x,ptc_y,ptc_z};
//               int ncoll = Nptclist[ievent][iptc].ncoll;
             
//               double ptc_tausq = ptc_t*ptc_t - ptc_z*ptc_z;
//               //std::cout << sqrt(ptc_tausq)<<std::endl;
//               if(ncoll == 0 || ptc_tausq < 0) continue;
//               nncoll++;
              
              
//               double etasi = 0.5f * (log(std::max(ptc_position[0]+ptc_position[3], small_eps)) \
//                            - log(std::max(ptc_position[0]-ptc_position[3], small_eps)));
          
//               double dxx = ptc_x - xi;
//               double dyy = ptc_y - yj;
//               double dzz = etasi - zk;
             

//               double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);


//               if ( fabs(dzz) < 10*SIGETA && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
                 
              
//                 double delta = exp(- distance_sqr)*w1*w2;
                
//                 double mt = sqrt(ptc_momentum[1] * ptc_momentum[1] + ptc_momentum[2] * ptc_momentum[2]+ptc_mass*ptc_mass);
//                 //std::cout << mt<<" "<<delta<<" "<<ptc_baryon <<std::endl;
//                 double Yi = 0.5f * (log(std::max(ptc_momentum[0]+ptc_momentum[3], small_eps)) \
//                            - log(std::max(ptc_momentum[0]-ptc_momentum[3], small_eps)));
//                 double momentum_miln[4] = { mt*cosh(Yi-zk), ptc_momentum[1],
//                                             ptc_momentum[2], mt*sinh(Yi-zk)};
                
//                 //if(gridi == 70 && gridj == 7 && gridk == 41 && iptc==878)
               


//                 double factorP_bmu = ptc_baryon;
//                 Taugrid[gridi][gridj][gridk].Jbmu[0] += ptc_baryon*delta;

//                 for(int mu = 0 ; mu < 4 ; mu ++){
                 
//                   Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)] += momentum_miln[mu]*delta;
//                   vv = std::max(vv,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]);
//                 }

    
//               }     

//             }

            

//       }
//       }
//       }
//       std::cout <<" max value "<< vv <<std::endl;

//       coutevent+=1;


//     }
//     std::cout<<"NEVENT: "<< coutevent<<std::endl;
//     double maxvale2=0.0;
    
     
//     for(int gridi = 0; gridi < NX; gridi++ ){
//       //std::cout<< gridi <<std::endl;
//       for(int gridj = 0; gridj < NY; gridj++ ){
//         for(int gridk = 0; gridk < NETA; gridk++ ){
//           double maxvale=0.0;
//           Taugrid[gridi][gridj][gridk].Jbmu[0]  /= coutevent;
//           for(int mu = 0 ; mu < 4 ; mu ++){
//                     Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]  /= coutevent;
                    
//                     maxvale= std::max(maxvale,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]);
//                     maxvale2 = std::max(maxvale2,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,0)]);

//                 }
//           if( maxvale < 1e-4) continue;
          
//          // std::cout<<maxvale<<std::endl;

         

//           double K2 = Taugrid[gridi][gridj][gridk].Tmnnu[idx(1,0)]*Taugrid[gridi][gridj][gridk].Tmnnu[idx(1,0)]
//                     + Taugrid[gridi][gridj][gridk].Tmnnu[idx(2,0)]*Taugrid[gridi][gridj][gridk].Tmnnu[idx(2,0)]
//                     + Taugrid[gridi][gridj][gridk].Tmnnu[idx(3,0)]*Taugrid[gridi][gridj][gridk].Tmnnu[idx(3,0)];
//           double K0 = Taugrid[gridi][gridj][gridk].Tmnnu[idx(0,0)];
//           double M = sqrt(K2);

//           double J0 = Taugrid[gridi][gridj][gridk].Jbmu[0];

//           rootFinding_newton(K0, M, J0,Taugrid[gridi][gridj][gridk].ed, Taugrid[gridi][gridj][gridk].nb); 
          
          
          
//           double EPV = std::max(small_eps, K0+eos.get_pressure(Taugrid[gridi][gridj][gridk].ed/HBARC,  Taugrid[gridi][gridj][gridk].nb)*HBARC); 
//           double vx = Taugrid[gridi][gridj][gridk].Tmnnu[idx(1,0)]/EPV;
//           double vy = Taugrid[gridi][gridj][gridk].Tmnnu[idx(2,0)]/EPV;
//           double vz = Taugrid[gridi][gridj][gridk].Tmnnu[idx(3,0)]/EPV;
          
//           double gamma = 1.0/sqrt(std::max(1.0-vx*vx-vy*vy-vz*vz, small_eps));
//           Taugrid[gridi][gridj][gridk].u0 = gamma;
//           Taugrid[gridi][gridj][gridk].ux = gamma*vx;
//           Taugrid[gridi][gridj][gridk].uy = gamma*vy;
//           Taugrid[gridi][gridj][gridk].uz = gamma*vz;

       
        
//         }
//       }

//     }
    
//     if(maxvale2 > 0){

    
//     std::string filename = Parameter::PATHOUT  + "SMASH_ini.dat";   
//     std::cout<<" "<<filename<<std::endl; 
//     std::ofstream fout(filename);


//     for(int gridi = 0; gridi < NX; gridi++ ){
//       for(int gridj = 0; gridj < NY; gridj++ ){
//         for(int gridk = 0; gridk < NETA; gridk++ ){
//           fout<< Taugrid[gridi][gridj][gridk].ed //GeV/fm^3
//               <<" "<< Taugrid[gridi][gridj][gridk].nb //1/fm^3
//               <<" "<< Taugrid[gridi][gridj][gridk].u0
//               <<" "<< Taugrid[gridi][gridj][gridk].ux
//               <<" "<< Taugrid[gridi][gridj][gridk].uy
//               <<" "<< Taugrid[gridi][gridj][gridk].uz/TAU0<<std::endl;     
//         }
//       }

//     }

//     fout.close();
//     }
//     std::cout<<maxvale2<<" www"<<std::endl;
// }

// void Grid::GausssmearingTauEta2(TEPaticlelist& Nptclist){
//    double one_o_2sigr2 = 1.0/(2.0*SIGR*SIGR); 
//    double one_o_2sigz2 = 1.0/(2.0*SIGETA*SIGETA);

//    double w1 = one_o_2sigr2/M_PI_F;
//    double w2 = sqrt(one_o_2sigz2/M_PI_F)/TAU0;
  

//     int coutevent = 0;


//     GridClear(Taugrid,NETA);
//     std::cout<<Nptclist.size()<<std::endl;
//     for(int ievent=0; ievent< Nptclist.size(); ievent++) 
//     { 
      
//       std::cout<<ievent<<" "<<Nptclist[ievent].size()<<std::endl;
     
//       for(int iptc=0;iptc < Nptclist[ievent].size(); iptc++){
        
//         double norm_facor = 0.0;
//         double ptc_t = Nptclist[ievent][iptc].t;
//         double ptc_x = Nptclist[ievent][iptc].x;
//         double ptc_y = Nptclist[ievent][iptc].y;
//         double ptc_z = Nptclist[ievent][iptc].z;
//         double ptc_e = Nptclist[ievent][iptc].e;
//         double ptc_px = Nptclist[ievent][iptc].px;
//         double ptc_py = Nptclist[ievent][iptc].py;
//         double ptc_pz = Nptclist[ievent][iptc].pz;
//         double ptc_mass = Nptclist[ievent][iptc].mass;
//         double ptc_baryon = Nptclist[ievent][iptc].baryon_number;

//         double ptc_momentum[4] = {ptc_e,ptc_px,ptc_py,ptc_pz};
//         double ptc_position[4] = {ptc_t,ptc_x,ptc_y,ptc_z};
//         int ncoll = Nptclist[ievent][iptc].ncoll;

//         double ptc_tausq = ptc_t*ptc_t - ptc_z*ptc_z;
//         if(ncoll == 0 || ptc_tausq < 0) continue;
        
//         double etasi = 0.5f * (log(std::max(ptc_position[0]+ptc_position[3], small_eps)) \
//                            - log(std::max(ptc_position[0]-ptc_position[3], small_eps)));

//          for(int gridi = 0; gridi < NX; gridi++ ){
//           double xi= (gridi-NX/2)*DX;
//           for(int gridj = 0; gridj < NY; gridj++ ){
//               double yj = (gridj-NY/2)*DY;
//               for(int gridk = 0; gridk < NETA; gridk++ ){
//                 double zk = (gridk-NETA/2)*DETA;         
//                 double dxx = ptc_x - xi;
//                 double dyy = ptc_y - yj;
//                 double dzz = etasi - zk;
             

//                 double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);


//                 if ( fabs(dzz) < 10*SIGETA && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
//                   double delta = exp(- distance_sqr)*w1*w2;
//                   norm_facor += delta;
//                 }
//               }
//           }
//          }

//         norm_facor =  norm_facor*DX*DY*DETA*TAU0;
                  
      
//         for(int gridi = 0; gridi < NX; gridi++ ){
//           double xi= (gridi-NX/2)*DX;
//           for(int gridj = 0; gridj < NY; gridj++ ){
//               double yj = (gridj-NY/2)*DY;
//               for(int gridk = 0; gridk < NETA; gridk++ ){
//                 double zk = (gridk-NETA/2)*DETA;         
//                 double dxx = ptc_x - xi;
//                 double dyy = ptc_y - yj;
//                 double dzz = etasi - zk;
             

//                 double distance_sqr = one_o_2sigr2*(dxx*dxx + dyy*dyy) + one_o_2sigz2*(dzz*dzz);


//                 if ( fabs(dzz) < 10*SIGETA && fabs(dyy) < 10*SIGR && fabs(dxx)< 10*SIGR ){
                 
              
//                   double delta = exp(- distance_sqr)*w1*w2/norm_facor;
//                   double mt = sqrt(ptc_momentum[1] * ptc_momentum[1] + ptc_momentum[2] * ptc_momentum[2]+ptc_mass*ptc_mass);
//                   double Yi = 0.5f * (log(std::max(ptc_momentum[0]+ptc_momentum[3], small_eps)) \
//                       - log(std::max(ptc_momentum[0]-ptc_momentum[3], small_eps)));
//                   double momentum_miln[4] = { mt*cosh(Yi-zk), ptc_momentum[1],
//                                             ptc_momentum[2], mt*sinh(Yi-zk)};
                               


//                   double factorP_bmu = ptc_baryon;
                
                
//                   for(int mu = 0 ; mu < 4 ; mu ++){
//                      Taugrid[gridi][gridj][gridk].Jbmu[mu] += momentum_miln[mu]*ptc_baryon*delta/(momentum_miln[0]+1e-7);
//                      for(int nu = 0; nu <4 ; nu++){
//                         double factorPP = momentum_miln[mu]*momentum_miln[nu]/(momentum_miln[0]+1e-7);

//                         Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)] += factorPP*delta;
                  

//                     }
                  
//                   }

    
//                 }     

            

//               }
//           }
//         }
//       } //end ptc loops in one event

//       coutevent+=1;


//     }
//     std::cout<<"NEVENT: "<< coutevent<<std::endl;
//     double maxvale2=0.0;
    
//     double check_nb = 0.0;
//     for(int gridi = 0; gridi < NX; gridi++ ){
//       for(int gridj = 0; gridj < NY; gridj++ ){
//         for(int gridk = 0; gridk < NETA; gridk++ ){
//           double maxvale=0.0;
         
//           for(int mu = 0 ; mu < 4 ; mu ++){
//             Taugrid[gridi][gridj][gridk].Jbmu[mu]  /= coutevent;
//             if (mu==0) check_nb+= Taugrid[gridi][gridj][gridk].Jbmu[mu];
//              for(int nu =0; nu < 4 ; nu ++){

//                     Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]  /= coutevent;
                    
//                     maxvale= std::max(maxvale,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);
//                     maxvale2 = std::max(maxvale2,Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]);

//              }

                    

//           }
//           if( maxvale < 1e-4) continue;
          
//          // std::cout<<maxvale<<std::endl;
//           Laudumatching(Taugrid[gridi][gridj][gridk].Tmnnu,Taugrid[gridi][gridj][gridk].ed,
//                         Taugrid[gridi][gridj][gridk].u0,Taugrid[gridi][gridj][gridk].ux,
//                         Taugrid[gridi][gridj][gridk].uy,Taugrid[gridi][gridj][gridk].uz);

//           Taugrid[gridi][gridj][gridk].nb = Taugrid[gridi][gridj][gridk].Jbmu[0]*Taugrid[gridi][gridj][gridk].u0
//                                         - Taugrid[gridi][gridj][gridk].Jbmu[1]*Taugrid[gridi][gridj][gridk].ux
//                                         - Taugrid[gridi][gridj][gridk].Jbmu[2]*Taugrid[gridi][gridj][gridk].uy
//                                         - Taugrid[gridi][gridj][gridk].Jbmu[3]*Taugrid[gridi][gridj][gridk].uz;
          
//           Taugrid[gridi][gridj][gridk].Temperature = eos.get_temperature(Taugrid[gridi][gridj][gridk].ed/HBARC,Taugrid[gridi][gridj][gridk].nb)*HBARC;    

//           // if(gridi == 95 && gridj == 80 && gridk == 100)
//           // {
//           //   // for(int mu=0;mu<4;mu++)
//           //   // {
//           //   //   for(int nu=0;nu<4;nu++){
//           //   //      std::cout<<Taugrid[gridi][gridj][gridk].Tmnnu[idx(mu,nu)]<<" ";
                
//           //   //   }
//           //   // }
//           //   std::cout<<std::endl;
//           //   std::cout<<Taugrid[gridi][gridj][gridk].ed<<" "
//           //            <<Taugrid[gridi][gridj][gridk].nb<<" "
//           //            <<Taugrid[gridi][gridj][gridk].ux<<" "
//           //            <<Taugrid[gridi][gridj][gridk].uy<<" "
//           //            <<Taugrid[gridi][gridj][gridk].uz<<" "<<check_nb*0.15*0.15*0.15*3.2<<std::endl;
//           // }

        
//         }
//       }

//     }
//     std::cout<<check_nb*DX*DY*DETA*TAU0<<std::endl;
    
//     if(maxvale2 > 0){

    
//     // std::string filename = Parameter::PATHOUT  + "SMASH_ini.dat";   
//     // // std::ofstream fout(filename,"wb");

//     // FILE *fout;
//     // fout = fopen(filename.c_str(), "w");


//     // for(int gridi = 0; gridi < NX; gridi++ ){
//     //   for(int gridj = 0; gridj < NY; gridj++ ){
//     //     for(int gridk = 0; gridk < NETA; gridk++ ){


//     //       float ideal[] = {static_cast<float>(Taugrid[gridi][gridj][gridk].ed),
//     //                         static_cast<float>(Taugrid[gridi][gridj][gridk].nb),
//     //                         static_cast<float>(Taugrid[gridi][gridj][gridk].u0),
//     //                         static_cast<float>(Taugrid[gridi][gridj][gridk].ux),
//     //                         static_cast<float>(Taugrid[gridi][gridj][gridk].uy),
//     //                         static_cast<float>(Taugrid[gridi][gridj][gridk].uz/TAU0)};
//     //       fwrite(ideal, sizeof(float), 6, fout);


//     //     }
//     //   }

//     // }

//     // fclose(fout);




//     std::string filename = Parameter::PATHOUT + "SMASH_ini.dat";
//     FILE *fout;
//     fout = fopen(filename.c_str(), "w");

//     if(fout == nullptr) {
//     // Handle file open error
//     }

//     for(int gridi = 0; gridi < NX; gridi++) {
//       for(int gridj = 0; gridj < NY; gridj++) {
//         for(int gridk = 0; gridk < NETA; gridk++) {

//             // Format the output to write in text format
//             fprintf(fout, "%f %f %f %f %f %f\n",
//                     static_cast<float>(Taugrid[gridi][gridj][gridk].ed),
//                     static_cast<float>(Taugrid[gridi][gridj][gridk].nb),
//                     static_cast<float>(Taugrid[gridi][gridj][gridk].u0),
//                     static_cast<float>(Taugrid[gridi][gridj][gridk].ux),
//                     static_cast<float>(Taugrid[gridi][gridj][gridk].uy),
//                     static_cast<float>(Taugrid[gridi][gridj][gridk].uz/TAU0));
//         }
//       }
//     }

// fclose(fout);




//     }
 

    
//     std::cout<<maxvale2<<" www"<<std::endl;

    

  
  

// }

Grid::~Grid(){
  
}





#endif
