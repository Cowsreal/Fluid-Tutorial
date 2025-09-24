#include "isentropicVortex.h"

isentropicVortex::isentropicVortex(double gamma, double vel1, double vel2, double x0, double y0, double epsilon, )
   : m_gamma(gamma), m_vel1(vel1), m_vel2(vel2), m_x0(x0), m_y0(y0), m_epsilon(epsilon) {}

void isentropicVortex::initialize(grid2D& grid)
{
   for(int i = 0; i < grid.getNx(); i++)
   {
      for(int j = 0; j < grid.getNy(); j++)
      {
         cell& currCell = grid.getCell(i, j);
         double x = i * grid.getDx() + grid.getDx() * 0.5;
         double y = i * grid.getDy() + grid.getDy() * 0.5;


      }
   }
}
