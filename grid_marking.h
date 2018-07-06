#ifndef __grid_marking_h_h
#define __grid_marking_h_h

#include "prelim.h"

template <int dim>
void elastic<dim>::grid_gen (int gn )
{
   if (gn==0)////////////////////////bridge////////////////////////////
   {
      std::vector<unsigned int> subdivisions (dim, 1);
      subdivisions[0] = 2 ;
      //subdivisions[1] = 2 ;

      const Point<dim> bottom_left = Point<dim>(0,0);
      const Point<dim> top_right   = Point<dim>(2,1) ;

      GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                 subdivisions,
                                                 bottom_left,
                                                 top_right);

      triangulation.refine_global (golobal_refine_number);  

      for (typename Triangulation<dim>::active_cell_iterator
                        cell = triangulation.begin_active();
                        cell != triangulation.end(); ++cell)
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->face(f) -> at_boundary())   
              {
                 const Point<dim> face_center = cell->face(f)->center();

                    if (face_center(1) < 0.05)
                    {
		       if (face_center(0) < 0.05 ) 
			  cell->face(f)->set_all_boundary_indicators(1);
		       if (face_center(0) > 1.95 ) 
			  cell->face(f)->set_all_boundary_indicators(2);
                       if (face_center(0) >0.45 && face_center(0) < 0.55)  
			  cell->face(f)->set_all_boundary_indicators(10);
                       if (face_center(0) >0.95 && face_center(0) < 1.05)  
			  cell->face(f)->set_all_boundary_indicators(4);
                       if (face_center(0) >1.45 && face_center(0) < 1.55)  
			  cell->face(f)->set_all_boundary_indicators(10);
		    }
              }
   }///////////******************************///////////////////////

   if (gn==1)//////////////////cantilever beam /////////////////////
   {
      std::vector<unsigned int> subdivisions (dim, 1);
      subdivisions[0] = 2 ;
      //subdivisions[1] = 2 ;

      const Point<dim> bottom_left = Point<dim>(0,0);
      const Point<dim> top_right   = Point<dim>(2,1) ;

      GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                 subdivisions,
                                                 bottom_left,
                                                 top_right);

      triangulation.refine_global (golobal_refine_number);  

      typename Triangulation<dim>:: cell_iterator cell = triangulation.begin_active (),
                                      endc = triangulation.end(); 
			
		for (; cell!=endc; ++cell)
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              //if (cell->face(f) -> at_boundary())   
              {
                 const Point<dim> face_center = cell->face(f)->center();
                 if (face_center(0) == 0.0)	
	            cell->face(f)->set_boundary_indicator(1);
                 if (face_center[1] < 0.05 )
                    if (face_center[0] > 1.95)
                       cell->face(f)->set_all_boundary_indicators(4); 
              }

   }

   if (gn == 2 )///////////////////////////////////////
   {
      GridGenerator::hyper_cube(triangulation,0,1) ;
      triangulation.refine_global (golobal_refine_number); 

   }///////////******************************///////////////////////

   if (gn == 3 )///////////////tensin test//////////////////////////////////
   {
      std::vector<unsigned int> subdivisions (dim, 1);
      subdivisions[0] = 6 ;

      const Point<dim> bottom_left = Point<dim>(0,0);
      const Point<dim> top_right   = Point<dim>(6,1) ;

      GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                 subdivisions,
                                                 bottom_left,
                                                 top_right);
      triangulation.refine_global (golobal_refine_number);  

      for (typename Triangulation<dim>::active_cell_iterator
                        cell = triangulation.begin_active();
                        cell != triangulation.end(); ++cell)
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->face(f) -> at_boundary())   
              {
                 const Point<dim> face_center = cell->face(f)->center();
                 if (face_center[0] < 0.05 )
		    cell->face(f)->set_all_boundary_indicators(6);
                 if (face_center[0] > 2.95)
		    cell->face(f)->set_all_boundary_indicators(7);
              }
   }///////////******************************///////////////////////

   if (gn == 4)////////////MBB Beam/////////////////////////////////
   {
            std::vector<unsigned int> subdivisions (dim, 1);
      subdivisions[0] = 4 ;
      //subdivisions[1] = 2 ;

      const Point<dim> bottom_left = Point<dim>(0,0);
      const Point<dim> top_right   = Point<dim>(4,1) ;

      GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                 subdivisions,
                                                 bottom_left,
                                                 top_right);

      triangulation.refine_global (golobal_refine_number);  

      for (typename Triangulation<dim>::active_cell_iterator
                        cell = triangulation.begin_active();
                        cell != triangulation.end(); ++cell)
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->face(f) -> at_boundary())   
              {
                 const Point<dim> face_center = cell->face(f)->center();
                 if (face_center[1] < 0.05)	
		 {
                    if (face_center[0] < 0.05 )
                       cell->face(f)->set_all_boundary_indicators(1);
                    if (face_center[0] > 3.95)
                       cell->face(f)->set_all_boundary_indicators(2);
		 }

		 if (face_center[1] > 0.95)
		    if (face_center[0] > 1.95 && face_center[0] < 2.05)
		       cell->face(f)->set_all_boundary_indicators(4);
              }
   }///////////******************************///////////////////////

}

#endif
