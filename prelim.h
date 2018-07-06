
#ifndef __prelim_h_h
#define __prelim_h_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
//#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>


#include <deal.II/fe/fe_system.h>

#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iostream>

using namespace dealii;
 
template <int dim>
class elastic
{
  public:

	elastic();
	 ~elastic() ;
	
        void org_vol();
        void  run() ;

  private:
	void grid_gen (int gn);
        void setup_system ();
	void assemble_direct ();
	void solve_direct ();
	void assemble_adjoint ();
	void solve_adjoint();
	void grad ();
	void calc_volume ();
	void grad_filter ();
	void opt();
	void output (); 
        void refine_grid ();
        void grad_smoothing ();

	Triangulation<dim>    triangulation ; 

	DoFHandler<dim>      dof_handler;
	 
        FESystem<dim>        fe;
	ConstraintMatrix     hanging_node_constraints;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix , system_matrix_2;

	Vector<double> 		      solution_a , solution_d , gradian,filter_gradient;
	Vector<double>		      system_rhs , system_rhs_2;
	Vector<double>		      cell_grad ;
	Vector<double>		      cell_volume;
	Vector<unsigned int> 	      cell_index;

        Vector<double> 		      dist ;  /////////////phase distribution in a cell//
        Vector<double> 		      mu ;    ////////////contain all mu values//////////
        Vector<double> 		      landa ; ////////////contain all landa values///////       
        Vector<double> 		      vol_bound ;
        Vector<double> 		      color ;
        Vector<double>                w ,wnew , fixed_mu ,fixed_landa ,fixed_vol ;
 
        std::vector<Vector<double> >  cell_dist;

	double		              volume_total , vol_frac ;

        //Vector<double>                mu_HS , landa_HS ;

};

bool             change_ok ;  //////if change be permited mu & landa change  (in case of energy caculation)//////

const double     golobal_refine_number = 6.0;

const double     phase_number = 3.; ////////total number of diffrent phases //////

const double     penal = 3. ;

const double     teta = 0.5 ;       //////////wheight for upper HS bound///////////

int	 	 iter;

const int        iter_number = 400; 

unsigned int     PN ;               //////////phase counter used in main///////////

double           energy ; 


#include "grid_marking.h"
#include "setup.h"
#include "organize_volume.h"
//#include "calculate_HS.h"
#include "assemble.h"
#include "solver.h"
#include "grad.h"
#include "grad_filter.h"
#include "grad_smoothing.h"
#include "opt.h"
#include "output.h"
#include "run.h"

#endif


