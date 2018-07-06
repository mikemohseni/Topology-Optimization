#ifndef __grad_filter_h_h
#define __grad_filter_h_h

#include "prelim.h"

using namespace dealii ;

template <int dim>
void elastic<dim>::grad_filter () 
{ 
 /*for(int ifilter=0; ifilter<1; ifilter++) {


 filter_gradient = 0.;

 typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						                        endc = dof_handler.end();

 
 std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
          
 active_neighbors.reserve (GeometryInfo<dim>::faces_per_cell *
                          GeometryInfo<dim>::max_children_per_face);


  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
  {

	  active_neighbors.clear ();
              
	  for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
                
		  if (! cell->at_boundary(n))      
		  {
			  
			  const typename DoFHandler<dim>::cell_iterator 
				  neighbor = cell->neighbor(n);

			  if (neighbor->active())
                      
				  active_neighbors.push_back (neighbor);
                    
			  else {

				  for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                            
					  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                             
						  if (neighbor->child(c)->neighbor(f) == cell)    
						  {  
							  active_neighbors.push_back (neighbor->child(c));
                                  
							  break;     
						  }
			  }
		  }*/

 //for(unsigned int i = 0; i<triangulation.n_active_cells(); ++i) 
   //              if( PN==0 && iter==99 ) printf ("  %f    " ,gradian(i) );


  filter_gradient = gradian ;

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell) 
      if ( ! cell->at_boundary() )
      {
         double sum1 = 1.2 * w(icell) * gradian(icell);
         double sum2 = 1.2;

       for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
            int celdex = cell -> neighbor_index(face_no) ;
	    sum1 += (0.2 * w(celdex) * gradian(celdex));
	    sum2 += 0.2;
          }

         filter_gradient(icell) = sum1 / (sum2 * w(icell));
     } 

  for (unsigned int i = 0; i<triangulation.n_active_cells(); ++i)   gradian(i) = filter_gradient(i);  
}

#endif
