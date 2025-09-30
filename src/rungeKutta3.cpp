#include "rungeKutta3.h"
#include "timeIntegrator.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>


void outputData(grid2D& grid, int step, double time, int field)
{
   std::ostringstream filename;
   switch(field)
   {
      case(0):
         filename << "data/rho/step" << std::setw(4) << std::setfill('0') << step << ".csv";
         break;
      case(1):
         filename << "data/u/step" << std::setw(4) << std::setfill('0') << step << ".csv";
         break;
      case(2):
         filename << "data/v/step" << std::setw(4) << std::setfill('0') << step << ".csv";
         break;
      case(3):
         filename << "data/pressure/step" << std::setw(4) << std::setfill('0') << step << ".csv";
         break;
   }

   std::ofstream file(filename.str());

   int nx = grid.getNx();
   int ny = grid.getNy();

   file << "# time = " << time << "\n";

   for (int j = 0; j < ny; j++)
   {
      for (int i = 0; i < nx; i++)
      {
         double rho = grid.getCell(i, j).var(field);
         file << rho;
         if (i < nx - 1) file << ",";
      }
      file << "\n";
   }
}

void rungeKutta3::step(grid2D& grid, physics& physics, double dt, int nStep, double t)
{

   physics.applyBC(grid);
   // RUNGE KUTTA 3, SHU OSHER SCHEME
   int ghostCells = grid.getGhostCells();

   for(int j = 0; j < grid.getNy(); j++)
   {
      for(int i = 0; i < grid.getNx(); i++)
      {
         std::vector<double>& s = grid.getCell(i, j).varState(0);
         physics.primToCons(s);
      }
   }

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

   physics.applyBC(F0);

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

   physics.applyBC(U2);

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

   physics.applyBC(grid);

   for(int j = 0; j < grid.getNy(); j++)
   {
      for(int i = 0; i < grid.getNx(); i++)
      {
         std::vector<double>& s = grid.getCell(i, j).varState(0);
         physics.consToPrim(s);
      }
   }

   for(int i = 0; i < 4; i++)
   {
      outputData(grid, nStep, t, i);
   }
}
