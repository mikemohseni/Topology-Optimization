
#ifndef __output_h_h
#define __output_h_h

#include "prelim.h"

template <int dim>
  void elastic<dim>::output ()
 /* {
    std::string filename = "solution-";
    filename += ('0' + iter);
    
    filename += ".plt";
    std::ofstream output (filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);

    data_out.add_data_vector (w, "w");
    data_out.build_patches ();
    data_out.write_gmv (output);
    
  }*/

    {   



    /*Vector<double>  toprint ;
    toprint.reinit(triangulation.n_active_cells());
    unsigned int n_des_var = triangulation.n_active_cells();
    for (unsigned int k = 0; k < n_des_var ; ++k)
        toprint(k) = cell_dist[k](1);*/



     std::string filename = "solution-";

  if(iter<10 ) filename += '0';
  if(iter>=10 && iter < 20) filename += '1';
  if(iter>=20 && iter < 30) filename += '2';
  if(iter>=30 && iter < 40) filename += '3';
  if(iter>=40 && iter < 50) filename += '4';
  if(iter>=50 && iter < 60) filename += '5';
  if(iter>=60 && iter < 70) filename += '6';
  if(iter>=70 && iter < 80) filename += '7';
  if(iter>=80 && iter < 90) filename += '8';
  if(iter>=90 && iter < 100) filename += '9';

  unsigned int filenumber = iter % 10;
    
     filename += ".dat"; 
     std::ofstream output (filename.c_str());

     DataOut<dim> data_out;
 
     data_out.attach_dof_handler (dof_handler);
     data_out.add_data_vector (color, "w");
     data_out.build_patches ();
     data_out.write_gnuplot (output);
   }

#endif

