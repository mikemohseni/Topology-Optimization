#ifndef __grad_smoothing_h_h
#define __grad_smoothing_h_h

#include "prelim.h"
//#include <deal.II/lac/solver_cg.h>

using namespace dealii ;

Vector<double>    solution_h , gradian_h ;



template <int dim>
class to_hilbert
{
  public:
    to_hilbert ();
    ~to_hilbert ();
    void run ();

  private:
    void setup_system ();
    void assemble_system ();

    Triangulation<dim>     triangulation;
    FE_Q<dim>              fe_h;
    DoFHandler<dim>        dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    ConstraintMatrix     hanging_node_constraints;
    Vector<double>       solution;
    Vector<double>       system_rhs;
};

  template <int dim>
  to_hilbert<dim>::to_hilbert (): fe_h (1) , 
                             dof_handler (triangulation) {}

template <int dim>
to_hilbert<dim>::~to_hilbert ()
   {
     dof_handler.clear ();
   }

template <int dim>
void to_hilbert<dim>::run ()
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

      // GridGenerator::hyper_cube (triangulation, 0, 1);

  triangulation.refine_global (golobal_refine_number);

  dof_handler.distribute_dofs (fe_h);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());


  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);

  hanging_node_constraints.close ();

  CompressedSparsityPattern  c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);

  hanging_node_constraints.condense (c_sparsity);

  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit (sparsity_pattern);

  const QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe_h, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe_h.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);
  Vector<double>       grad (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),endc = dof_handler.end();
  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += ((0.01*fe_values.shape_grad(i,q_point) *
                                        fe_values.shape_grad(j,q_point)*
                                        fe_values.JxW(q_point))
                                       +
                                       (fe_values.shape_value(i,q_point) *
                                       fe_values.shape_value(j,q_point) *
                                      fe_values.JxW(q_point)));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            gradian_h(icell) *
                            fe_values.JxW(q_point));
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

    SparseDirectUMFPACK umfpack;
    umfpack.factorize (system_matrix);
    umfpack.solve (system_rhs);
    solution = system_rhs;
    hanging_node_constraints.distribute (solution);

       

      cell = dof_handler.begin_active();
      endc = dof_handler.end();
  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
         {

	  fe_values.reinit (cell);
	  
	  cell->get_dof_indices (local_dof_indices);

	  for (unsigned int i=0; i<dofs_per_cell; ++i) 
		  grad(i) = solution(local_dof_indices[i]);
     

	   double tmp =0. ;
           for (unsigned int i=0; i<dofs_per_cell; ++i)
           for (unsigned int q_point=0; q_point<n_q_points; ++q_point)      
              tmp += (fe_values.shape_value(i,q_point) *
                       grad(i) );

 
               solution_h(icell) = tmp/n_q_points ;
         }
}

template <int dim>
void elastic<dim>::grad_smoothing ()
{
  //solution_h = gradian ;
  //gradian_h =gradian ;

  gradian_h.reinit(triangulation.n_active_cells());
  solution_h.reinit(triangulation.n_active_cells());

  unsigned  int n_des_var = triangulation.n_active_cells();
  for (unsigned int k = 0; k < n_des_var ; ++k)   gradian_h(k) = gradian(k) ;

  to_hilbert<dim>     Hilbertproblem;

  Hilbertproblem.run() ;

  typename DoFHandler<dim>::active_cell_iterator  cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)
     gradian(icell)=solution_h(icell);
if(iter==4)for (unsigned int k = 0; k < n_des_var ; ++k) printf("%E    ",gradian(k));
}

#endif
