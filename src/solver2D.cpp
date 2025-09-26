#include "solver2D.h"

solver2D::solver2D(int nx, int ny, double dx, double dy, int ghostCells, physics& phys, timeIntegrator& timeIntegrator, double dt, double tStop)
   : m_grid(nx, ny, dx, dy, ghostCells), m_physics(phys), m_timeIntegrator(timeIntegrator), m_dt(dt), m_tStop(tStop)
{

}


void solver2D::run()
{
   double t = 0.0;
   int nStep = 0;

   while(t < m_tStop)
   {
      m_timeIntegrator.step(m_grid, m_physics, m_dt);
      t += m_dt;
      nStep++;
      if(nStep % 10 == 0)
      {
         std::cout << "Step " << nStep << ", t = " << t << std::endl;
      }
   }
}
