#ifndef __run_h_h
#define __run_h_h

#include "prelim.h"


template <int dim>
void elastic<dim>::run () 
   
  {

        if ( iter == 0 && PN == 0)

           grid_gen (0);

/*if (iter == 0)  penal=2 ;
if (iter == 50) penal=3 ;
if (iter == 100)penal=4 ;
if (iter == 150)penal=5 ;*/
   /* if ( iter == 0)
      {
       std::vector<unsigned int> subdivisions (dim, 1);
       subdivisions[0] = 3;
       subdivisions[1] = 2;

       const Point<dim> bottom_left = Point<dim>(0,0);
       const Point<dim> top_right   = Point<dim>(3,2) ;

       GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                  subdivisions,
                                                  bottom_left,
                                                  top_right);

      // GridGenerator::hyper_cube (triangulation, 0, 1);
       triangulation.refine_global (golobal_refine_number);
      
       for (typename Triangulation<dim>::active_cell_iterator
                         cell = triangulation.begin_active();
                         cell != triangulation.end(); ++cell)
       for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
           if (cell->face(f) -> at_boundary())   
           {
               const Point<dim> face_center = cell->face(f)->center();
               if (face_center(1) < 0.05)
                 if (face_center(0) < 0.05 || face_center(0) > 2.95)
 
                   cell->face(f)->set_all_boundary_indicators(1);
              }
   }*/

	setup_system () ;

        if(iter==0) org_vol() ;

       	assemble_direct () ;

	solve_direct () ;

	grad () ;

	grad_filter ();

        //grad_smoothing ();         

	opt () ;
 
        change_ok = true ;

        grad ();

        unsigned int n_des_var = triangulation.n_active_cells();

        double  energy ;

        for (unsigned int k = 0; k < n_des_var ; k++)  energy += (gradian(k)/n_des_var) ;
        printf("%E  \n", energy*1000000);
        
        change_ok = false ;

/////////////interpolate last phase fractions //////////////////

        if ( PN == phase_number - 2 && iter == (iter_number -1) )           
        {
           typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                          endc = dof_handler.end();
           for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
           {
	       cell_dist[icell](phase_number -1 ) = w(icell) ;
	       cell_dist[icell](phase_number -2 ) = 1. - w(icell) - fixed_vol(icell) ;
               if (cell_dist[icell](phase_number -2) < 0.) 
                   cell_dist[icell](phase_number -2 ) = 0.;
           }
        }

        if ( PN == phase_number - 2 && iter == (iter_number -1) ) //////coloring all cells and output/////////////
        {    
           typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                          endc = dof_handler.end();
           for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
           {
                unsigned int x  = 0 ;

                for (unsigned int f = 1 ; f < phase_number ; ++f)  
                {   
		    if (cell_dist[icell](f) > 0.5)
                    x = f ;
                }
                
                color(icell) = x ; 
            }

            output () ;
        }
	// unsigned int n_des_var = triangulation.n_active_cells();
	 //if(PN==0&&iter==0)for (unsigned int k = 0; k < n_des_var ; k++)printf("%E   ",gradian(k));
//////////////////////////////////////////////////////////////////////////////

  }

#endif

