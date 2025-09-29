#include "rungeKutta3.h"
#include "timeIntegrator.h"

void rungeKutta3::step(grid2D& grid, physics& physics, double dt, int nStep, double t)
{

   // RUNGE KUTTA 3, SHU OSHER SCHEME
   int ghostCells = grid.getGhostCells();

   // Stage 1
   grid2D U0 = grid;
   grid2D F0 = grid;
   grid2D U1 = grid;
   physics.computeFluxes(F0, nStep, t);
   
   // Time step
   
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = ghostCells; i < grid.getNx() - ghostCells; i++)
      {
         cell& c = U0.getCell(i, j);
         cell& cNew = U1.getCell(i, j);
         std::vector<double> flux0 = F0.getCell(i, j).varState(5);

         for(int k = 0; k < 4; k++)
         {
            cNew.var(k) = c.var(k) + dt * flux0[k]; 
         }
      }
   }


   physics.applyBC(U1);

   // Stage 2

   grid2D F1 = U1;
   grid2D U2 = U1;
   
   physics.computeFluxes(F1, nStep, t);
   
   // Time step
   
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = ghostCells; i < grid.getNx() - ghostCells; i++)
      {
         cell& c0 = U0.getCell(i, j);
         cell& c1= U1.getCell(i, j);
         cell& cNew = U2.getCell(i, j);
         std::vector<double> flux0 = F0.getCell(i, j).varState(5);
         std::vector<double> flux1 = F1.getCell(i, j).varState(5);

         for(int k = 0; k < 4; k++)
         {
            cNew.var(k) = c0.var(k) + 1.0 / 4.0 * dt * flux0[k] + 1.0 / 4.0 * dt * flux1[k]; 
         }
      }
   }

   // Reset boundary conditions

   physics.applyBC(U1);

   // Final Stage

   grid2D F2 = U2;

   physics.computeFluxes(F2, nStep, t);
   
   // Time step
   
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = ghostCells; i < grid.getNx() - ghostCells; i++)
      {
         cell& c0 = U0.getCell(i, j);
         cell& c1 = U1.getCell(i, j);
         cell& c2 = U2.getCell(i, j);
         cell& cNew = grid.getCell(i, j);
         std::vector<double> flux0 = F0.getCell(i, j).varState(5);
         std::vector<double> flux1 = F1.getCell(i, j).varState(5);
         std::vector<double> flux2 = F2.getCell(i, j).varState(5);

         for(int k = 0; k < 4; k++)
         {
            cNew.var(k) = c0.var(k) + 1.0 / 6.0 * dt * flux0[k] + 1.0 / 6.0 * dt * flux1[k] + 2.0 / 3.0 * dt * flux2[k]; 
         }
      }
   }

   // Reset boundary conditions

   physics.applyBC(U1);
}
