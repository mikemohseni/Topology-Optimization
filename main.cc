#include "prelim.h"

using namespace dealii;

int main()

{
  PN = 0 ;

  dealii::deallog.depth_console (0);

  elastic<2>  elast_opt;

  for(PN = 0 ; PN < phase_number-1 ; ++PN) 
  {

     for ( iter=0; iter< iter_number ; ++iter)

            elast_opt.run();
  }

  return 0;
}

