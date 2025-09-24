#pragma once

#include "physics.h"

#include <cmath>

class isentropicVortex : public physics
{
   public:
      isentropicVortex();

      void initialize(grid2D& grid) override;
      void computeFluxes(grid2D& grid) override;
      void applyBC(grid2D& grid) override;

      void consToPrim(grid2D grid);
      void primToCons(grid2D grid);

   private:
      double m_gamma;
      double m_epsilon;
      double m_vel1;
      double m_vel2;
      double m_x0, m_y0;
}
