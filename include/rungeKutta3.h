#pragma once

#include "timeIntegrator.h"


class rungeKutta3 : public timeIntegrator
{
   public: 
      rungeKutta3(){};
      void step(grid2D& grid, physics& physics, double dt, int nStep, double t) override; 
};

