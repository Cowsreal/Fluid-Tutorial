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

   isentropicVortex vort(1.4, 5.0, 5.0, 0.5, 0.5, 5.0);

   rungeKutta3 integrator;

   solver2D solver(300, 300, 1 / 300.0, 1 / 300.0, 3, vort, integrator, 1.0);
   solver.run();
}
