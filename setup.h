
#ifndef __setup_h_h
#define __setup_h_h

#include "prelim.h"

using namespace dealii ;

template <int dim>

  elastic<dim>::elastic(): dof_handler (triangulation), fe (FE_Q<dim>(1),dim)
               
  {}

template <int dim>

  elastic<dim>::~elastic ()
  {
    dof_handler.clear ();
    
  }


template <int dim>
void elastic<dim>::calc_volume () 
{  
  
  QGauss<dim>  quadrature_formula(2);
				   
  FEValues<dim> fe_values (fe, quadrature_formula, 
			               UpdateFlags(update_values    |
				           update_gradients |
				           update_q_points  |
				           update_JxW_values));
			   
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),endc = dof_handler.end();
  
  volume_total = 0.0;
  cell_volume = 0.0;

  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
  {

      fe_values.reinit (cell);

	  cell->get_dof_indices (local_dof_indices); 
	  
	  double temp  = 0.0;

	  for (unsigned int i=0; i<dofs_per_cell; ++i)

		  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)

			  temp  += fe_values.shape_value(i,q_point) * fe_values.JxW(q_point);

		  cell_volume (icell) = temp;

		  volume_total += temp;
  }

}

/*template <int dim>
   void elastic<dim>::refine_grid ()
   {
     Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
 
     typename FunctionMap<dim>::type neumann_boundary;
     KellyErrorEstimator<dim>::estimate (dof_handler,QGauss<dim-1>(2),
                                         neumann_boundary,
                                         solution_d,
                                         estimated_error_per_cell);
 
     GridRefinement::refine_and_coarsen_fixed_number (triangulation,estimated_error_per_cell, 0.3, 0.03);
 
     triangulation.execute_coarsening_and_refinement ();
   }*/



template <int dim>
  void elastic<dim>::setup_system ()
{
   // else
   //   refine_grid ();

    unsigned int n_des_var= triangulation.n_active_cells();

    dof_handler.distribute_dofs (fe);
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
					     hanging_node_constraints);
    hanging_node_constraints.close ();
    sparsity_pattern.reinit (dof_handler.n_dofs(),
			     dof_handler.n_dofs(),
			     dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    hanging_node_constraints.condense (sparsity_pattern);

    sparsity_pattern.compress();

    system_matrix.reinit (sparsity_pattern);
    system_matrix_2.reinit (sparsity_pattern);

    solution_d.reinit(dof_handler.n_dofs());
    solution_a.reinit(dof_handler.n_dofs());
    gradian.reinit(triangulation.n_active_cells());
    color.reinit(triangulation.n_active_cells());
   // dist.reinit(phase_number);
    
    system_rhs.reinit (dof_handler.n_dofs());
    system_rhs_2.reinit (dof_handler.n_dofs());


    cell_volume.reinit (n_des_var);
    cell_index.reinit (n_des_var);

    calc_volume () ;

}

#endif

