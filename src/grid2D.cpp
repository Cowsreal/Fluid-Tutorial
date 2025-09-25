#include "grid2D.h"


grid2D::grid2D(int nx, int ny, double dx, double dy)
   : m_nx(nx), m_ny(ny), m_dx(dx), m_dy(dy)
{
   m_grid.resize(nx);
   for(auto& row : m_grid)
   {
      row.resize(ny);
   }
}

cell& grid2D::getCell(int i, int j)
{
   return m_grid[i][j];
}
