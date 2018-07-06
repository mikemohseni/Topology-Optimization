
#ifndef __opt_h_h
#define __opt_h_h

#include "prelim.h"

template <int dim>
void elastic<dim>::opt()

{
     double l1=0, l2=100000, move=0.2, lmid, sum;//, constraint;

	 unsigned int n_des_var = triangulation.n_active_cells();

	 while (l2-l1 > 0.00000001) 
         {
              lmid= 0.5*(l2+l1);

	      for (unsigned int k = 0; k < n_des_var ; ++k) { 

		   double step = w(k)*sqrt(-gradian(k)/lmid);

		   wnew(k) = step;

		   if ( w(k)+move < step )			wnew(k) = w(k)+move;

		   if ( (1.0-fixed_vol(k)) < wnew(k)  )	        wnew(k)= (1.0-fixed_vol(k)) ;
			   
		   if (w(k)-move > wnew(k)  )			wnew(k) = w(k)-move ;
 
		   if ( 0.0000001 > wnew(k) )			wnew(k) = 0.0000001 ;
	   }
	   
	   sum=0.;  
	   for (unsigned int k = 0; k < n_des_var ; k++)	sum+= (wnew(k) * cell_volume(k));   
	   
	   if ( (sum - vol_frac * volume_total) > 0)	l1 = lmid;	
	   else	l2 = lmid;
	 }

	 for (unsigned int k = 0; k < n_des_var ; k++)	w(k) = wnew(k);
//if(iter==iter_number-1 )for (unsigned int k = 0; k < n_des_var ; ++k)printf("%f        ",w(k));

}

#endif
