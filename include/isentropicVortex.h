#pragma once
#define _USE_MATH_DEFINES

#include "physics.h"

#include <cmath>
#include <algorithm>
#include <fstream>
#include <vector>

class isentropicVortex : public physics
{
   public:
      isentropicVortex(double gamma, double vel1, double vel2, double x0, double y0, double epsilon);

      void initialize(grid2D& grid) override;
      void computeFluxes(grid2D& grid, int nStep, double t) override;
      void applyBC(grid2D& grid) override;
      std::vector<double> HLLC(std::vector<double>& L, std::vector<double>& R, int dir);

      double getGamma() override { return m_gamma; } 

      void consToPrim(std::vector<double>& s) override;
      void primToCons(std::vector<double>& s) override;

   private:
      double m_gamma;
      double m_epsilon;
      double m_vel1;
      double m_vel2;
      double m_x0, m_y0;
};
