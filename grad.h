#ifndef __grad_h_h
#define __grad_h_h
#include "prelim.h"

using namespace dealii;

template <int dim>
void elastic<dim>::grad () 
{ 
  QGauss<dim>  quadrature_formula(2);
  const unsigned int   n_q_points    = quadrature_formula.size();
				   
  FEValues<dim> fe_values (fe, quadrature_formula, 
			                   UpdateFlags(update_values |
				           update_gradients |
				           update_q_points  |
				           update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double>   d_cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  gradian =0;

  FEValuesExtractors::Vector displacement(0);
  std::vector<Tensor<2,dim> >  solution_grads (n_q_points);
  std::vector<SymmetricTensor<2,dim> > strain (n_q_points);
  std::vector<SymmetricTensor<2,dim> > stress (n_q_points);


  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int icell = 0; cell!=endc; ++cell, ++icell)   
  {
      double  d_mu = penal * pow(w(icell) ,(penal-1))*(mu(phase_number-1.) - mu(PN)) ;
      double  d_landa = penal * pow(w(icell) ,(penal-1))*(landa(phase_number-1.) - landa(PN)) ;
      double  grad = 0. ;

      if (change_ok)/// if calculating compliance theses are added////////
      {
	 d_mu = w(icell)*(mu(phase_number-1.) - mu(PN)) + fixed_mu(icell) ;
         d_landa = w(icell)*(landa(phase_number-1.) - landa(PN)) + fixed_landa(icell) ;
      }

      fe_values.reinit (cell);
      fe_values[displacement].get_function_gradients (solution_d, solution_grads);

      for (unsigned int q=0; q<n_q_points ; ++q) 
      {
          strain[q] = symmetrize(solution_grads[q]);
   
          stress[q] =  (2*d_mu*strain[q] + d_landa*trace(strain[q])*
                       unit_symmetric_tensor<dim>());
          grad += stress[q]*strain[q]*fe_values.JxW(q);
      }
  
      gradian(icell) = -0.5*grad ;
  }
}

#endif
