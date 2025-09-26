#include "isentropicVortex.h"
#include <fstream>

#define SMALL_VAL 1e-6

isentropicVortex::isentropicVortex(double gamma, double vel1, double vel2, double x0, double y0, double epsilon)
   : m_gamma(gamma), m_vel1(vel1), m_vel2(vel2), m_x0(x0), m_y0(y0), m_epsilon(epsilon) {}

void isentropicVortex::initialize(grid2D& grid)
{
   for(int i = 2; i < grid.getNx() - 2; i++)
   {
      double x = (i - 2) * grid.getDx() + grid.getDx() * 0.5;
      std::cout << i << ", "<< grid.getDx()<< ", " << std::endl;
      for(int j = 2; j < grid.getNy() - 2; j++)
      {
         cell& currCell = grid.getCell(i, j);
         // get cell center
         double y = (j - 2) * grid.getDy() + grid.getDy() * 0.5;

         currCell.var(1) = m_vel1 - (x - m_x0) * m_epsilon / (2 * M_PI) * exp(1 - (x - m_x0) * (x - m_x0) - (y - m_y0) * (y - m_y0)); 
         currCell.var(2) = m_vel1 + (y - m_y0) * m_epsilon / (2 * M_PI) * exp(1 - (x - m_x0) * (x - m_x0) - (y - m_y0) * (y - m_y0)); 

         // temperature
         double temp = 1 - (m_gamma - 1) * m_epsilon * m_epsilon / (8 * m_gamma * M_PI * M_PI);
         temp *= exp(1 - (x - m_x0) * (x - m_x0) - (y - m_y0) * (y - m_y0));

         currCell.var(0) = pow(temp, 1 / (m_gamma - 1));
         currCell.var(3) = currCell.var(0) * temp;
      }
   }

   applyBC(grid);

   std::ofstream f("grid.csv");

   for(int i = 0; i < grid.getNx(); i++)
   {
      for(int j = 0; j < grid.getNy(); j++)
      {
         cell& currCell = grid.getCell(i, j);
         f << currCell.var(0);
         if(j != grid.getNy() - 1)
         {
            f << ", ";
         }
      }
      f << "\n";
   }
   f.close();
}

void isentropicVortex::computeFluxes(grid2D& grid)
{
   int ghostCells = grid.getGhostCells();
   // TRANSFORM TO CCONSERVATIVE
   
   primToCons(grid);

   // STEP 1: WENO RECONSTRUCTION

   // COLS 
   for(int i = ghostCells; i < grid.getNx() - ghostCells; i++)
   {
      for(int j = ghostCells; j < grid.getNx() - ghostCells; j++)
      {
         cell& cm2 = grid.getCell(i, j - 2);
         cell& cm1 = grid.getCell(i, j - 1);
         cell& c = grid.getCell(i, j);
         cell& cp1 = grid.getCell(i, j + 1);
         cell& cp2 = grid.getCell(i, j + 2);
         for(int k = 0; k < 4; k++)
         {
            // Calculate i + 1/2 value
            double poly0 = 1.0 / 3.0 * c.var(k) + 5.0 / 6.0 * cp1.var(k) - 1.0 / 6.0 * cp2.var(k);  
            double poly1 = -1.0 / 6.0 * cm1.var(k) + 5.0 / 6.0 * c.var(k) + 1.0 / 3.0 * cp1.var(k);  
            double poly2 = 1.0 / 3.0 * cm2.var(k) - 7.0 / 6.0 * cm1.var(k) + 11.0 / 6.0 * c.var(k);  

            double d0 = 3.0 / 10.0;
            double d1 = 3.0 / 5.0;
            double d2 = 1.0 / 10.0;

            double beta0 = 13.0 / 12.0 * (c.var(k) - 2 * cp1.var(k) + cp2.var(k)) * (c.var(k) - 2 * cp1.var(k) + cp2.var(k)) + 1.0 / 4.0 * (3 * c.var(k) - 4 * cp1.var(k) + cp2.var(k)) * (3 * c.var(k) - 4 * cp1.var(k) + cp2.var(k));
            double beta1 = 13.0 / 12.0 * (cm1.var(k) - 2 * c.var(k) + cp1.var(k)) * (cm1.var(k) - 2 * c.var(k) + cp1.var(k)) + 1.0 / 4.0 * (cm1.var(k) - cp1.var(k)) * (cm1.var(k) - cp1.var(k));
            double beta2 = 13.0 / 12.0 * (cm2.var(k) - 2 * cm1.var(k) + c.var(k)) * (cm2.var(k) - 2 * cm1.var(k) + c.var(k)) + 1.0 / 4.0 * (cm2.var(k) - 4 * cm1.var(k) + 3 * c.var(k)) * (cm2.var(k) - 4 * cm1.var(k) + 3 * c.var(k));

            double alpha0 = d0 / ((SMALL_VAL + beta0) * (SMALL_VAL + beta0));
            double alpha1 = d1 / ((SMALL_VAL + beta1) * (SMALL_VAL + beta1));
            double alpha2 = d2 / ((SMALL_VAL + beta2) * (SMALL_VAL + beta2));

            double alphaSum = alpha0 + alpha1 + alpha2;
            
            double w0 = alpha0 / alphaSum;
            double w1 = alpha1 / alphaSum;
            double w2 = alpha2 / alphaSum;
            c.varU(k) = w0 * poly0 + w1 * poly1 + w2 * poly2;

            // Calculate i - 1/2 value
            


         }
      }
   }



};
void isentropicVortex::applyBC(grid2D& grid)
{
   int ghostCells = grid.getGhostCells();
   // LEFT
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = 0; i < ghostCells; i++)
      {
         cell& targetCell = grid.getCell(grid.getNx() - 2 * ghostCells + i, j);
         cell& ghostCell = grid.getCell(i, j);
         ghostCell = targetCell;
      }
   }

   // RIGHT 
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = grid.getNx() - ghostCells; i < grid.getNx(); i++)
      {
         cell& targetCell = grid.getCell(i - (grid.getNx() - 2 * ghostCells), j);
         cell& ghostCell = grid.getCell(i, j);
         ghostCell = targetCell;
      }
   }

   // TOP 
   for(int j = ghostCells; j < grid.getNx() - ghostCells; j++)
   {
      for(int i = grid.getNy() - ghostCells; i < grid.getNy(); i++)
      {
         cell& targetCell = grid.getCell(j, i - (grid.getNy() - 2 * ghostCells));
         cell& ghostCell = grid.getCell(j, i);
         ghostCell = targetCell;
      }
   }

   // BOTTOM 
   for(int j = ghostCells; j < grid.getNx() - ghostCells; j++)
   {
      for(int i = 0; i < ghostCells; i++)
      {
         cell& targetCell = grid.getCell(j, grid.getNy() - 2 * ghostCells + i);
         cell& ghostCell = grid.getCell(j, i);
         ghostCell = targetCell;
      }
   }
};

void isentropicVortex::consToPrim(grid2D& grid)
{
   for(int i = 2; i < grid.getNx() - 2; i++)
   {
      for(int j = 2; j < grid.getNy() - 2; j++)
      {
         cell& currCell = grid.getCell(i, j);

         double rho = currCell.var(0);
         double energy = currCell.var(3);

         currCell.var(1) /= rho;
         currCell.var(2) /= rho;
         currCell.var(3) = (m_gamma - 1) * (energy  - 0.5 * rho * (currCell.var(1) * currCell.var(1) + currCell.var(2) * currCell.var(2)));
      }
   }
};

void isentropicVortex::primToCons(grid2D& grid)
{
   for(int i = 2; i < grid.getNx() - 2; i++)
   {
      for(int j = 2; j < grid.getNy() - 2; j++)
      {
         cell& currCell = grid.getCell(i, j);

         double rho = currCell.var(0);
         double pres = currCell.var(3);
         double energy = pres / (m_gamma - 1) + 0.5 * rho * (currCell.var(1) * currCell.var(1) + currCell.var(2) * currCell.var(2));

         currCell.var(1) *= rho;
         currCell.var(2) *= rho;
         currCell.var(3) = (m_gamma - 1) * (energy  - 0.5 * rho * (currCell.var(1) * currCell.var(1) + currCell.var(2) * currCell.var(2)));
      }
   }
};

