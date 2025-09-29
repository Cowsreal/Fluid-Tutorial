#include "solver2D.h"
/*
void outputData(grid2D& grid, int step, double time)
{
    std::ostringstream filename;
    filename << "data/density_step" << std::setw(4) << std::setfill('0') << step << ".csv";

    std::ofstream file(filename.str());

    int nx = grid.getNx();
    int ny = grid.getNy();

    file << "# time = " << time << "\n";

    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            double rho = grid.getCell(i, j).var(0);
            file << rho;
            if (i < nx - 1) file << ",";
        }
        file << "\n";
    }
}
*/
solver2D::solver2D(int nx, int ny, double dx, double dy, int ghostCells, physics& phys, timeIntegrator& timeIntegrator, double tStop)
   : m_grid(nx, ny, dx, dy, ghostCells), m_physics(phys), m_timeIntegrator(timeIntegrator), m_tStop(tStop)
{
   m_physics.initialize(m_grid);
}


void solver2D::run()
{
   double t = 0.0;
   int nStep = 0;
   m_dt = cfl(m_grid, m_physics);

   while(t < m_tStop)
   {
      m_timeIntegrator.step(m_grid, m_physics, m_dt, nStep, t);
      t += m_dt;
      nStep++;
      if(nStep % 10 == 0)
      {
         std::cout << "Step " << nStep << ", t = " << t << std::endl;
      }
//      outputData(nStep, t);
   } 
}


double solver2D::cfl(grid2D& grid, physics& phys)
{
   int ghostCells = grid.getGhostCells();

   double minVal = 1e9; 
   double dx = grid.getDx();
   double dy = grid.getDy();
   double gamma = phys.getGamma();

   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = ghostCells; i < grid.getNx() - ghostCells; i++)
      {
         cell& c = grid.getCell(i , j);
         double rho = c.var(0);
         double p = c.var(3);
         double u = std::abs(c.var(1));
         double v = std::abs(c.var(2));
         double a = sqrt(gamma * p / rho);
         minVal = std::min(minVal, std::min(dx / (u + a), dy / (v + a)));
      }
   }

   return m_c * minVal;
}

