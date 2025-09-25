#pragma once
#include "grid2D.h"
#include "physics.h"

class timeIntegrator
{
   public:
      virtual void step(grid2D& grid, physics& physics, double dt) = 0;
      virtual ~timeIntegrator() {}
};
