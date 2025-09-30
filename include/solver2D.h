#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "cell.h"
#include "grid2D.h"
#include "physics.h"
#include "timeIntegrator.h"

class solver2D
{
   public: 
      solver2D(int nx, int ny, double dx, double dy, int ghostCells, physics& phys, timeIntegrator& timeIntegrator, double tStop);
      void run();
      double cfl(grid2D& grid, physics& phys);
      
      void outputData(int step, double time);
   private:
      grid2D m_grid;
      physics& m_physics;
      timeIntegrator& m_timeIntegrator;
      double m_dt;
      double m_tStop;
      double m_c = 0.2;
};
