#ifndef __solver_h_h
#define __solver_h_h

#include "prelim.h"
#include <deal.II/lac/sparse_direct.h>

template <int dim>
   void elastic<dim>::solve_direct ()

   {

    SparseDirectUMFPACK umfpack;
    umfpack.factorize (system_matrix);
    umfpack.solve (system_rhs);
    solution_d = system_rhs;
    hanging_node_constraints.distribute (solution_d);

   }

template <int dim>
   void elastic<dim>::solve_adjoint ()
   {
    SparseDirectUMFPACK umfpack;
    umfpack.factorize (system_matrix_2);
    umfpack.solve (system_rhs_2);
    solution_a = system_rhs;
    hanging_node_constraints.distribute (solution_a);
   }

#endif
