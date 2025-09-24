#include "grid2D.h"


grid2D::grid2D(int nx, int ny, int dx, int dy)
   : m_nx(nx), m_ny(ny), m_dx(dx), m_dy(dy)
{
   m_grid.resize(nx, std::vector<cell>(ny, 4));
}

cell& grid2D::getCell(int i, int j)
{
   return m_grid[i][j];
}
