#pragma once

#include "cell.h"
#include <vector>

class grid2D
{
   public:
      // constructor
      grid2D(int nx, int ny, double dx, double dy);

      // return pointer to cell at i, j 
      cell& getCell(int i, int j);

      int getNx(){return m_nx;}
      int getNy(){return m_ny;}
      double getDx(){return m_dx;}
      double getDy(){return m_dy;}

   private:
      int m_nx, m_ny;
      double m_dx, m_dy;
      std::vector<std::vector<cell>> m_grid;

};
