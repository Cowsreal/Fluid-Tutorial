#pragma once
#include "grid2D.h"

class physics
{
   public: 
      virtual void initialize(grid2D& grid) = 0;
      virtual void computeFluxes(grid2D& grid) = 0;
      virtual void applyBC(grid2D& grid) = 0;

      virtual void consToPrim(grid2D& grid) = 0;
      virtual void primToCons(grid2D& grid) = 0;

      virtual ~physics(){}
}
