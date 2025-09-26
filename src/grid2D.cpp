#include "grid2D.h"


grid2D::grid2D(int nx, int ny, double dx, double dy, int ghostCells)
   : m_nx(nx), m_ny(ny), m_dx(dx), m_dy(dy), m_ghostCells(ghostCells)
{
   m_grid.resize(nx + ghostCells * 2);
   for(auto& row : m_grid)
   {
      row.resize(ny + ghostCells * 2);
   }
}

cell& grid2D::getCell(int i, int j)
{
   return m_grid[i][j];
}
