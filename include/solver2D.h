#pragma once

#include "grid2D.h"
#include "physics.h"
#include "timeIntegrator.h"

class solver2D
{
   public: 
      solver2D(int nx, int ny, double dx, double dy, physics& phys, timeIntegrator& timeIntegrator, double dt, double tStop);
      void run();
      
   private:
      grid2D m_grid;
      physics& m_physics;
      timeIntegrator& m_timeIntegrator;
      double m_dt;
      double m_tStop;
};
