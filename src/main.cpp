#include <iostream>

#include "cell.h"
#include "grid2D.h"
#include "isentropicVortex.h"
#include "solver2D.h"
#include "timeIntegrator.h"


int main()
{
   grid2D grid(50, 50, 1 / 50.0, 1 / 50.0);
   std::cout << "hi" << std::endl;
   isentropicVortex vort(1.4, 0.0, 0.0, 0.5, 0.5, 3.2);

   std::cout << "hi" << std::endl;
   vort.initialize(grid);
}
