#include "isentropicVortex.h"
#include <fstream>

isentropicVortex::isentropicVortex(double gamma, double vel1, double vel2, double x0, double y0, double epsilon)
   : m_gamma(gamma), m_vel1(vel1), m_vel2(vel2), m_x0(x0), m_y0(y0), m_epsilon(epsilon) {}

void isentropicVortex::initialize(grid2D& grid)
{
   for(int i = 0; i < grid.getNx(); i++)
   {
      double x = i * grid.getDx() + grid.getDx() * 0.5;
      for(int j = 0; j < grid.getNy(); j++)
      {
         cell& currCell = grid.getCell(i, j);
         // get cell center
         double y = i * grid.getDy() + grid.getDy() * 0.5;

         currCell.var(1) = m_vel1 - (x - m_x0) * m_epsilon / (2 * M_PI) * exp(1 - (x - m_x0) * (x - m_x0) - (y - m_y0) * (y - m_y0)); 
         currCell.var(2) = m_vel1 + (y - m_y0) * m_epsilon / (2 * M_PI) * exp(1 - (x - m_x0) * (x - m_x0) - (y - m_y0) * (y - m_y0)); 

         std::cout << i << ", " << j << std::endl;
         // temperature
         double temp = 1 - (m_gamma - 1) * m_epsilon * m_epsilon / (8 * m_gamma * M_PI * M_PI);
         temp *= exp(1 - (x - m_x0) * (x - m_x0) - (y - m_y0) * (y - m_y0));

         currCell.var(0) = pow(temp, 1 / (m_gamma - 1));
         currCell.var(3) = currCell.var(0) * temp;
      }
   }

   std::ofstream f("grid.csv");
   for(int i = 0; i < grid.getNx(); i++)
   {
      for(int j = 0; j < grid.getNy(); j++)
      {
         cell& currCell = grid.getCell(i, j);
         f << currCell.var(0);

         if(i == grid.getNx() - 1)
         {
            f << "\n";
         }
      }
   }
   f.close();
}

void isentropicVortex::computeFluxes(grid2D& grid)
{

};
void isentropicVortex::applyBC(grid2D& grid){

};

void isentropicVortex::consToPrim(grid2D& grid)
{

};
void isentropicVortex::primToCons(grid2D& grid)
{

};

