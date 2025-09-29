#include <iostream>

#include "cell.h"
#include "grid2D.h"
#include "isentropicVortex.h"
#include "solver2D.h"
#include "timeIntegrator.h"
#include "rungeKutta3.h"

int main()
{
   // +4 for WENO5 ghost cells

   isentropicVortex vort(1.4, 0.0, 0.0, 0.5, 0.5, 10.2);

   rungeKutta3 integrator;

   solver2D solver(20, 20, 1 / 20.0, 1 / 20.0, 3, vort, integrator, 1.0);
   solver.run();
}
