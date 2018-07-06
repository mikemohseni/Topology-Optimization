
#ifndef __organize_volume_h_h
#define __organize_volume_h_h

#include "prelim.h"

using namespace dealii;

template <int dim>
void elastic<dim>::org_vol() 
{
  Vector<double>  HSP ;

  HSP.reinit(phase_number);

  if (PN ==0 ) 
  {
      mu.reinit(phase_number);
      landa.reinit(phase_number);
      vol_bound.reinit(phase_number);
 
      mu(0)=0.0000000001; 
      mu(1)=1.0;
      mu(2)=2.0;
      //mu(3)=0.0000000001;
      //mu(4)=0.0000000001;
      landa(0)=0.0000000001;
      landa(1)=1.0;
      landa(2)=2.0;
      //landa(3)=0.0000000001;
      //landa(4)=0.0000000001;
      vol_bound(0)= 0.4 ;
      vol_bound(1)= 0.25 ;
      vol_bound(2)= 0.35 ;
      //vol_bound(3)= 0.4 ;
      //vol_bound(4)= 0.4 ;

      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),endc = dof_handler.end();

      for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)

          cell_dist.push_back(vol_bound) ; 

      w.reinit(triangulation.n_active_cells());
      wnew.reinit(triangulation.n_active_cells());
      fixed_mu.reinit(triangulation.n_active_cells());
      fixed_landa.reinit(triangulation.n_active_cells());
      fixed_vol.reinit(triangulation.n_active_cells());

      change_ok = false ;


  }

  if (PN != 0) 
  {
     typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                    endc = dof_handler.end();
     for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
     {
	 cell_dist[icell](PN-1) = 1.0 - w(icell) ;
         
         HSP = cell_dist[icell] ;

         double uUV = 0. ;

         double P_V = 0. ; 
     
         for (unsigned int f = 0 ; f < PN ; ++f)    
             uUV += HSP(f);

         for (unsigned int f = PN ; f < phase_number ; ++f)    
             P_V += HSP(f);

         for (unsigned int f = PN ; f < phase_number ; ++f)  
         {
             double HL = ((1.-uUV)/P_V) * HSP(f) ; 
             HSP(f) = HL ;
         }     

         for (unsigned int f = PN ; f < phase_number ; ++f)

             cell_dist[icell](f) = HSP(f) ;
     }  
  }

  vol_frac = 0. ;

  for (unsigned int f = PN+1 ; f < phase_number ; ++f)

      vol_frac += vol_bound(f) ;
 
  w = wnew =fixed_mu = fixed_landa = fixed_vol = 0.0 ;

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
  {  
      for (unsigned int f = PN+1 ; f < phase_number ; ++f)

	  w(icell) += cell_dist[icell](f) ;

      wnew(icell) = w(icell) ;

      if (PN != 0 )

         for (unsigned int f =0 ; f < PN ; ++f)
         {
	     fixed_mu(icell) += mu(f) * cell_dist[icell](f) ;
	     fixed_landa(icell) += mu(f) * cell_dist[icell](f) ;
      	     fixed_vol(icell) += cell_dist[icell](f) ;
         }
  }
	 //unsigned int n_des_var = triangulation.n_active_cells();
	 //if(PN==1&&iter==0)for (unsigned int k = 0; k < n_des_var ; k++)printf("%f   ",w(k));
}

#endif
