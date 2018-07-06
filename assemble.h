
#ifndef __assemble_h_h
#define __assemble_h_h

#include "prelim.h"
/*
template <int dim>
  class RightHandSide :  public Function<dim>
  {
    public:
      RightHandSide ();

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &values) const;

      virtual void vector_value_list (const std::vector<Point<dim> > &points,
				      std::vector<Vector<double> >   &value_list) const;
  };


  template <int dim>
  RightHandSide<dim>::RightHandSide ()
		  :
		  Function<dim> (dim)
  {}


  template <int dim>
  inline
  void RightHandSide<dim>::vector_value (const Point<dim> &p,
					 Vector<double>   &values) const
  {
     if (p(0) > 1.45 && p(0) < 1.55) {
        if (p(1) < 0.05 ){
            values(0) = 0.0;
            values(1) = -100.0 ;}
      /* else if ( p(1)<0.05)
          {
           values(0) = -100.0;
           values(1) = 0.0;
          }  
     }
     else 
     {   
         values(0) = 0.0;
         values(1) = 0.0;  
     }
   // if (p.square() < 0.2*0.2)

  //  else
    //  values(1) = 20.0;

  }

  template <int dim>
  void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
					      std::vector<Vector<double> >   &value_list) const
  {
    Assert (value_list.size() == points.size(),
	    ExcDimensionMismatch (value_list.size(), points.size()));

    const unsigned int n_points = points.size();


    for (unsigned int p=0; p<n_points; ++p)
      RightHandSide<dim>::vector_value (points[p],value_list[p]);
	
  }  

///////////////////////boundary condition///////////////////////////
 template <int dim>
  class DirichletBoundaryValues : public Function<dim>
  {
    public:
      DirichletBoundaryValues() : Function<dim> (2) {};

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &values) const;

      virtual void vector_value_list (const std::vector<Point<dim> > &points,
				      std::vector<Vector<double> >   &value_list) const;
  };


  template <int dim>
  inline
  void DirichletBoundaryValues<dim>::vector_value (const Point<dim> &p,
						   Vector<double>   &values) const
  {
    if (p(0) < 0.05 ) {
       values(0) = 0.;
       values(1) = 0.;
    }

    if (p(0) > 2.95 )
       values(1) = 0. ;
  }


  template <int dim>
  void DirichletBoundaryValues<dim>::vector_value_list (const std::vector<Point<dim> > &points,
							std::vector<Vector<double> >   &value_list) const
  {
    Assert (value_list.size() == points.size(),
	    ExcDimensionMismatch (value_list.size(), points.size()));

    for (unsigned int p=0; p<points.size(); ++p)
      DirichletBoundaryValues<dim>::vector_value (points[p], value_list[p]);
  }

*/

template <int dim>
void elastic<dim>::assemble_direct () 
{  
    QGauss<dim>  quadrature_formula(2);
    QGauss<dim-1> face_quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
                                     UpdateFlags(update_values  |
                                     update_q_points       |
                                     update_normal_vectors |
                                     update_JxW_values));


    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();
    const unsigned int   n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    //RightHandSide<dim>      right_hand_side;
    //std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim));

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();

    for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
    {
        cell_matrix = 0;

        cell_rhs = 0. ;

        fe_values.reinit (cell) ; 

     // right_hand_side.vector_value_list (fe_values.get_quadrature_points(),rhs_values);

        double  cell_mu = (pow(w(icell) , penal)*(mu(phase_number -1) - mu(PN)))+mu(PN) + fixed_mu(icell) ;
        double  cell_landa = (pow(w(icell) , penal)*(landa(phase_number-1) - landa(PN)))+landa(PN) + fixed_landa(icell) ;

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            const unsigned int   component_i = fe.system_to_component_index(i).first;

            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
                const unsigned int component_j = fe.system_to_component_index(j).first;

                for (unsigned int q_point=0; q_point<n_q_points; ++q_point) 
                {
                    cell_matrix(i,j)+=(

                        (fe_values.shape_grad(i,q_point)[component_i] *
                         fe_values.shape_grad(j,q_point)[component_j] *
                         cell_landa) +
                        (fe_values.shape_grad(i,q_point)[component_j] *
                         fe_values.shape_grad(j,q_point)[component_i] *
                         cell_mu)    +
                        ((component_i == component_j) ?
                         (fe_values.shape_grad(i,q_point) *
                          fe_values.shape_grad(j,q_point) *
                          cell_mu)  :
                         0.)
                      )
                      *
                      fe_values.JxW(q_point);
                  }
              }
          }

       
/* cell->get_dof_indices (local_dof_indices);
       
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int
              component_i = fe.system_to_component_index(i).first;

            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_rhs(i) += fe_values.shape_value(i,q_point) *
                             rhs_values[q_point](component_i) *
                             fe_values.JxW(q_point);
          }  */   


	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face) 
          {
	      unsigned int mark = cell->face(face)->boundary_indicator();
	      fe_face_values.reinit (cell, face);

	      for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		  const unsigned int component_i = fe.system_to_component_index(i).first;

		  for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point) 
		  {
		      if (mark == 1) 
		      { 					
			 cell_rhs(i) = 0.; 
			 cell_matrix(i,i) += 1.e30;
		      }

		      if (mark == 2) // y=0 x-free 
			 if (component_i ==1) 
			 {
			    cell_rhs(i) = 0.;
			    cell_matrix(i,i) += 1.e30;
		         }

		      if (mark == 3) // x=0 y-free 
			 if (component_i ==0) 
			 {
			    cell_rhs(i) = 0.;
			    cell_matrix(i,i) += 1.e30;
			 }

	              if (mark == 4) 
		      {
			 cell_rhs(i) += ( (component_i == 0 ? 0.0:-100.0) *     
			 fe_values.shape_value(i,q_point) *	
			 fe_values.JxW(q_point));
		      }

	              if (mark == 5) 
			 if (component_i ==1)
			 {
			    cell_rhs(i) += ( 1.0 *     
			    fe_face_values.shape_value(i,q_point) *	
			    fe_face_values.JxW(q_point));
			 }

	              if (mark == 6) 
			 if (component_i ==0)
			 {
			    cell_rhs(i) += ( -1.0 *     
			    fe_face_values.shape_value(i,q_point) *	
			    fe_face_values.JxW(q_point));
			 }

	              if (mark == 7) 
			 if (component_i ==0)
			 {
			    cell_rhs(i) += ( 1.0 *     
			    fe_face_values.shape_value(i,q_point) *	
			    fe_face_values.JxW(q_point));
			 }
	              if (mark == 10) 
		      {
			 cell_rhs(i) += ( (component_i == 0 ? 0.0:-50.0) *     
			 fe_values.shape_value(i,q_point) *	
			 fe_values.JxW(q_point));
		      }

		  }
	      }
	  }   

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            for (unsigned int j=0; j<dofs_per_cell; ++j)

              system_matrix.add (local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i,j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    hanging_node_constraints.condense (system_matrix);
    hanging_node_constraints.condense (system_rhs);

   /* std::map<unsigned int,double>    boundary_values;

    VectorTools::interpolate_boundary_values (dof_handler,       
					      1,
					      DirichletBoundaryValues<dim>(),
					      boundary_values);

    MatrixTools::apply_boundary_values (boundary_values,
				        system_matrix,
				        solution_d,
				        system_rhs);*/

}

/*template <int dim>
  void elastic<dim>::assemble_adjoint () 
  {  

    QGauss<dim>  quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);
    Vector<double>       direct_ans(dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    RightHandSide<dim>      right_hand_side;
    std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim));

    //system_matrix.reinit (sparsity_pattern);
    //system_rhs.reinit (dof_handler.n_dofs());

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();

    for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);
    
        double  mu = (pow(w(icell) , penal)*(mu_1 - mu_2))+mu_2 ;
        double  landa = (pow(w(icell) , penal)*(landa_1 - landa_2))+landa_2 ;

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int
              component_i = fe.system_to_component_index(i).first;

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const unsigned int component_j = fe.system_to_component_index(j).first;

                for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
                    cell_matrix(i,j)-=((

                        (fe_values.shape_grad(i,q_point)[component_i] *
                         fe_values.shape_grad(j,q_point)[component_j] *
                         landa)
                        -
                        (fe_values.shape_grad(i,q_point)[component_j] *
                         fe_values.shape_grad(j,q_point)[component_i] *
                         mu)
                        -
                        ((component_i == component_j) ?
                         (fe_values.shape_grad(i,q_point) *
                          fe_values.shape_grad(j,q_point) *
                          mu)  :
                         0)
                      )
                      *
                      fe_values.JxW(q_point));

                  }
              }
          }

       for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += cell_matrix(i,j)* solution_d(icell) ; 

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              system_matrix_2.add (local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i,j));

            system_rhs_2(local_dof_indices[i]) += cell_rhs(i);
          }
      }

    hanging_node_constraints.condense (system_matrix);
    hanging_node_constraints.condense (system_rhs);

  std::map<unsigned int,double> boundary_values;


 if(0) { VectorTools::interpolate_boundary_values (dof_handler,
					    1,
					    DirichletBoundaryValues<dim>(),
					    boundary_values);

  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution_a,
				      system_rhs);
       }
  }*/

#endif
