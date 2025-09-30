#include "isentropicVortex.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#define SMALL_VAL 1e-6

/*
void outputData(grid2D& grid, int step, double time)
{
    std::ostringstream filename;
    filename << "data/density_step" << std::setw(4) << std::setfill('0') << step << ".csv";

   std::cout << "saving data at step " << step << std::endl;
    std::ofstream file(filename.str());

    int nx = grid.getNx();
    int ny = grid.getNy();

    file << "# time = " << time << "\n";

    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            double rho = grid.getCell(i, j).varState(5)[0];
            file << rho;
            if (i < nx - 1) file << ",";
        }
        file << "\n";
    }
}
*/
isentropicVortex::isentropicVortex(double gamma, double vel1, double vel2, double x0, double y0, double epsilon)
   : m_gamma(gamma), m_vel1(vel1), m_vel2(vel2), m_x0(x0), m_y0(y0), m_epsilon(epsilon) {}

void isentropicVortex::initialize(grid2D& grid)
{
   double beta = 0.2;
   for(int i = 2; i < grid.getNx() - 2; i++)
   {
      double x = (i - 3) * grid.getDx() + grid.getDx() * 0.5;
      for(int j = 2; j < grid.getNy() - 2; j++)
      {
         // get cell center
         double y = (j - 3) * grid.getDy() + grid.getDy() * 0.5;
         double r2 = ((x - m_x0) * (x - m_x0) + (y - m_y0) * (y - m_y0)) / (beta * beta);
         double expTerm = exp(1.0 - r2);

         cell& c = grid.getCell(i,j);

         // velocity
         c.var(1) = m_vel1 - (y - m_y0) * m_epsilon / (2.0 * M_PI) * expTerm;
         c.var(2) = m_vel2 + (x - m_x0) * m_epsilon / (2.0 * M_PI) * expTerm;

         double T = 1.0 - (m_gamma - 1) * m_epsilon * m_epsilon / (8.0 * m_gamma * M_PI * M_PI) * expTerm;
         c.var(0) = pow(T, 1.0 / (m_gamma - 1.0));   // density
         c.var(3) = pow(c.var(0), m_gamma);       
      }
   }

   applyBC(grid);

   std::ofstream f("grid.csv");

   for(int i = 0; i < grid.getNx(); i++)
   {
      for(int j = 0; j < grid.getNy(); j++)
      {
         cell& currCell = grid.getCell(i, j);
         f << currCell.var(0);
         if(j != grid.getNy() - 1)
         {
            f << ", ";
         }
      }
      f << "\n";
   }
   f.close();
}

void isentropicVortex::computeFluxes(grid2D& grid, int nStep, double t)
{
   int ghostCells = grid.getGhostCells();
   // TRANSFORM TO CCONSERVATIVE
/*   
   for(int j = 0; j < grid.getNy(); j++)
   {
      for(int i = 0; i < grid.getNx(); i++)
      {
         std::vector<double>& s = grid.getCell(i, j).varState(0);
         primToCons(s);
      }
   }
*/
   // STEP 1: WENO RECONSTRUCTION

   // COLS 
   for(int j = ghostCells - 1; j < grid.getNy() - ghostCells + 1; j++)
   {
      for(int i = ghostCells - 1; i < grid.getNx() - ghostCells + 1; i++)
      {
         cell& cm2 = grid.getCell(i, j - 2);
         cell& cm1 = grid.getCell(i, j - 1);
         cell& c = grid.getCell(i, j);
         cell& cp1 = grid.getCell(i, j + 1);
         cell& cp2 = grid.getCell(i, j + 2);
         for(int k = 0; k < 4; k++)
         {
            // Calculate i + 1/2 value
            double poly0 = 1.0 / 3.0 * c.var(k) + 5.0 / 6.0 * cp1.var(k) - 1.0 / 6.0 * cp2.var(k);  
            double poly1 = -1.0 / 6.0 * cm1.var(k) + 5.0 / 6.0 * c.var(k) + 1.0 / 3.0 * cp1.var(k);  
            double poly2 = 1.0 / 3.0 * cm2.var(k) - 7.0 / 6.0 * cm1.var(k) + 11.0 / 6.0 * c.var(k);  

            double d0 = 3.0 / 10.0;
            double d1 = 3.0 / 5.0;
            double d2 = 1.0 / 10.0;

            double beta0 = 13.0 / 12.0 * (c.var(k) - 2 * cp1.var(k) + cp2.var(k)) * (c.var(k) - 2 * cp1.var(k) + cp2.var(k)) + 1.0 / 4.0 * (3 * c.var(k) - 4 * cp1.var(k) + cp2.var(k)) * (3 * c.var(k) - 4 * cp1.var(k) + cp2.var(k));
            double beta1 = 13.0 / 12.0 * (cm1.var(k) - 2 * c.var(k) + cp1.var(k)) * (cm1.var(k) - 2 * c.var(k) + cp1.var(k)) + 1.0 / 4.0 * (cm1.var(k) - cp1.var(k)) * (cm1.var(k) - cp1.var(k));
            double beta2 = 13.0 / 12.0 * (cm2.var(k) - 2 * cm1.var(k) + c.var(k)) * (cm2.var(k) - 2 * cm1.var(k) + c.var(k)) + 1.0 / 4.0 * (cm2.var(k) - 4 * cm1.var(k) + 3 * c.var(k)) * (cm2.var(k) - 4 * cm1.var(k) + 3 * c.var(k));

            double alpha0 = d0 / ((SMALL_VAL + beta0) * (SMALL_VAL + beta0));
            double alpha1 = d1 / ((SMALL_VAL + beta1) * (SMALL_VAL + beta1));
            double alpha2 = d2 / ((SMALL_VAL + beta2) * (SMALL_VAL + beta2));

            double alphaSum = alpha0 + alpha1 + alpha2;
            
            double w0 = alpha0 / alphaSum;
            double w1 = alpha1 / alphaSum;
            double w2 = alpha2 / alphaSum;
            c.varU(k) = w0 * poly0 + w1 * poly1 + w2 * poly2;

            // Calculate i - 1/2 value
            
            poly0 = 11.0 / 6.0 * c.var(k) - 7.0 / 6.0 * cp1.var(k) + 1.0 / 3.0 * cp2.var(k);  
            poly1 = 1.0 / 3.0 * cm1.var(k) + 5.0 / 6.0 * c.var(k) - 1.0 / 6.0 * cp1.var(k);  
            poly2 = - 1.0 / 6.0 * cm2.var(k) + 5.0 / 6.0 * cm1.var(k) + 1.0 / 3.0 * c.var(k);  

            d0 = 1.0 / 10.0;
            d1 = 3.0 / 5.0;
            d2 = 3.0 / 10.0;

            alpha0 = d0 / ((SMALL_VAL + beta0) * (SMALL_VAL + beta0));
            alpha1 = d1 / ((SMALL_VAL + beta1) * (SMALL_VAL + beta1));
            alpha2 = d2 / ((SMALL_VAL + beta2) * (SMALL_VAL + beta2));

            alphaSum = alpha0 + alpha1 + alpha2;
            
            w0 = alpha0 / alphaSum;
            w1 = alpha1 / alphaSum;
            w2 = alpha2 / alphaSum;
            c.varD(k) = w0 * poly0 + w1 * poly1 + w2 * poly2;
         }
      }
   }
   

   // ROWS 
   for(int j = ghostCells - 1; j < grid.getNy() - ghostCells + 1; j++)
   {
      for(int i = ghostCells - 1; i < grid.getNx() - ghostCells + 1; i++)
      {
         cell& cm2 = grid.getCell(i - 2, j);
         cell& cm1 = grid.getCell(i - 1, j);
         cell& c = grid.getCell(i, j);
         cell& cp1 = grid.getCell(i + 1, j);
         cell& cp2 = grid.getCell(i + 2, j);
         for(int k = 0; k < 4; k++)
         {
            // Calculate i + 1/2 value
            double poly0 = 1.0 / 3.0 * c.var(k) + 5.0 / 6.0 * cp1.var(k) - 1.0 / 6.0 * cp2.var(k);  
            double poly1 = -1.0 / 6.0 * cm1.var(k) + 5.0 / 6.0 * c.var(k) + 1.0 / 3.0 * cp1.var(k);  
            double poly2 = 1.0 / 3.0 * cm2.var(k) - 7.0 / 6.0 * cm1.var(k) + 11.0 / 6.0 * c.var(k);  

            double d0 = 3.0 / 10.0;
            double d1 = 3.0 / 5.0;
            double d2 = 1.0 / 10.0;

            double beta0 = 13.0 / 12.0 * (c.var(k) - 2 * cp1.var(k) + cp2.var(k)) * (c.var(k) - 2 * cp1.var(k) + cp2.var(k)) + 1.0 / 4.0 * (3 * c.var(k) - 4 * cp1.var(k) + cp2.var(k)) * (3 * c.var(k) - 4 * cp1.var(k) + cp2.var(k));
            double beta1 = 13.0 / 12.0 * (cm1.var(k) - 2 * c.var(k) + cp1.var(k)) * (cm1.var(k) - 2 * c.var(k) + cp1.var(k)) + 1.0 / 4.0 * (cm1.var(k) - cp1.var(k)) * (cm1.var(k) - cp1.var(k));
            double beta2 = 13.0 / 12.0 * (cm2.var(k) - 2 * cm1.var(k) + c.var(k)) * (cm2.var(k) - 2 * cm1.var(k) + c.var(k)) + 1.0 / 4.0 * (cm2.var(k) - 4 * cm1.var(k) + 3 * c.var(k)) * (cm2.var(k) - 4 * cm1.var(k) + 3 * c.var(k));

            double alpha0 = d0 / ((SMALL_VAL + beta0) * (SMALL_VAL + beta0));
            double alpha1 = d1 / ((SMALL_VAL + beta1) * (SMALL_VAL + beta1));
            double alpha2 = d2 / ((SMALL_VAL + beta2) * (SMALL_VAL + beta2));

            double alphaSum = alpha0 + alpha1 + alpha2;
            
            double w0 = alpha0 / alphaSum;
            double w1 = alpha1 / alphaSum;
            double w2 = alpha2 / alphaSum;
            c.varR(k) = w0 * poly0 + w1 * poly1 + w2 * poly2;

            // Calculate i - 1/2 value
            
            poly0 = 11.0 / 6.0 * c.var(k) - 7.0 / 6.0 * cp1.var(k) + 1.0 / 3.0 * cp2.var(k);  
            poly1 = 1.0 / 3.0 * cm1.var(k) + 5.0 / 6.0 * c.var(k) - 1.0 / 6.0 * cp1.var(k);  
            poly2 = - 1.0 / 6.0 * cm2.var(k) + 5.0 / 6.0 * cm1.var(k) + 1.0 / 3.0 * c.var(k);  

            d0 = 1.0 / 10.0;
            d1 = 3.0 / 5.0;
            d2 = 3.0 / 10.0;

            alpha0 = d0 / ((SMALL_VAL + beta0) * (SMALL_VAL + beta0));
            alpha1 = d1 / ((SMALL_VAL + beta1) * (SMALL_VAL + beta1));
            alpha2 = d2 / ((SMALL_VAL + beta2) * (SMALL_VAL + beta2));

            alphaSum = alpha0 + alpha1 + alpha2;
            
            w0 = alpha0 / alphaSum;
            w1 = alpha1 / alphaSum;
            w2 = alpha2 / alphaSum;
            c.varL(k) = w0 * poly0 + w1 * poly1 + w2 * poly2;
         }
      }
   }

   /*
   for(int j = ghostCells - 1; j < grid.getNy() - ghostCells + 1; j++)
   {
      for(int i = ghostCells - 1; i < grid.getNx() - ghostCells + 1; i++)
      {
         for(int k = 0; k < 5; k++)
         {
            auto& s = grid.getCell(i, j).varState(k);
            consToPrim(s);
         }
      }
   }
   */

   // HLLC SOLVER
   
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = ghostCells - 1; i < grid.getNx() - ghostCells; i++)
      {
         cell& cl = grid.getCell(i, j);
         cell& cr = grid.getCell(i + 1, j);
         
         // RIEMANN PROBLEM RIGHT EDGE
         cl.varState(4) = HLLC(cl.varState(4), cr.varState(3), 0);
         cr.varState(3) = cl.varState(4);
      }
   }

   for(int i = ghostCells; i < grid.getNx() - ghostCells; i++)
   {
      for(int j = ghostCells - 1; j < grid.getNy() - ghostCells; j++)
      {
         cell& cd = grid.getCell(i, j);
         cell& cu = grid.getCell(i, j + 1);
         
         // RIEMANN PROBLEM TOP EDGE
         cd.varState(1) = HLLC(cd.varState(1), cu.varState(2), 1);
         cu.varState(2) = cd.varState(1);
      }
   }
      
   // FLUX average estimates

   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = ghostCells; i < grid.getNx() - ghostCells; i++)
      {
         cell& c = grid.getCell(i, j);

         std::vector<double> hoizontalDifference(4);
         std::vector<double> verticalDifference(4);

         std::transform(
            c.varState(3).begin(), c.varState(3).end(),
            c.varState(4).begin(),
            hoizontalDifference.begin(),
            [&grid](double a, double b) { return -(b - a) / grid.getDx(); }
         );
         std::transform(
            c.varState(2).begin(), c.varState(2).end(),
            c.varState(1).begin(),
            verticalDifference.begin(),
            [&grid](double a, double b) { return -(b - a) / grid.getDy(); }
         );        
         std::transform(
            hoizontalDifference.begin(), hoizontalDifference.end(), 
            verticalDifference.begin(), 
            c.varState(5).begin(),
            [](double a, double b) { return a + b; }
         );
      }
   }

   //outputData(grid, nStep, t);
};


std::vector<double> isentropicVortex::HLLC(std::vector<double>& L, std::vector<double>& R, int dir)
{
   double rhoL = L[0];
   double rhoR = R[0];
   double uL;
   double uR;
   double vL;
   double vR;
   // Total energy
   double eL = L[3];
   double eR = R[3];
   switch(dir)
   {
      case(0):   
         uL = L[1] / rhoL;
         uR = R[1] / rhoR;
         vL = L[2] / rhoL;
         vR = R[2] / rhoR;
         break;
      case(1):   
         uL = L[2] / rhoL;
         uR = R[2] / rhoR;
         vL = L[1] / rhoL;
         vR = R[1] / rhoR;
         break;
   } 

   double pL = (m_gamma - 1) * (eL - 0.5 * rhoL * (uL * uL + vL * vL));
   double pR = (m_gamma - 1) * (eL - 0.5 * rhoL * (uL * uL + vL * vL));

   double aL = sqrt(m_gamma * pL / rhoL);
   double aR = sqrt(m_gamma * pR / rhoR);

   double pvrs = 0.5 * (pL + pR) - 0.5 * (uR - uL) * 0.5 * (rhoL + rhoR) * 0.5 * (aL + aR); 
   double p = std::max(pvrs, 0.0);

   double qL;
   double qR;

   if(p <= pL)
   {
      qL = 1;
   }
   else
   {
      qL = sqrt(1 + (m_gamma + 1) / (2 * m_gamma) * (p / pL - 1));
   }
   if(p <= pR)
   {
      qR = 1;
   }
   else
   {
      qR = sqrt(1 + (m_gamma + 1) / (2 * m_gamma) * (p / pR - 1));
   }

   // Wave speed estimates

   double SL = uL - aL * qL;
   double SR = uR + aR * qR;
   double S = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR * (SR - uR));

   /*
   // Total energy
   double eL = pL / (m_gamma - 1) + 0.5 * rhoL * (uL * uL + vL * vL);
   double eR = pR / (m_gamma - 1) + 0.5 * rhoR * (uR * uR + vR * vR);
   */

   // UL and UR in the normal directions
   std::vector<double> UL = {rhoL, rhoL * uL, rhoL * vL, eL};
   std::vector<double> UR = {rhoR, rhoR * uR, rhoR * vR, eR};

   // FL and FR in the normal directions
   std::vector<double> FL = {rhoL * uL, rhoL * uL * uL + pL, rhoL * uL * vL, uL * (eL + pL)};
   std::vector<double> FR = {rhoR * uR, rhoR * uR * uR + pR, rhoR * uR * vR, uR * (eR + pR)};

   std::vector<double> HLLCFlux(4);
   std::vector<double> D = {0, 1, 0, S};

   if(0 <= SL)
   {
      HLLCFlux = FL;
   }
   else if(SL < 0 && 0 <= S)
   {
      
      double scale = rhoL * (SL - uL) / (SL - S);
      // Difference between center (star region) and state
      std::vector<double> UStarL = {
         scale * 1.0,
         scale * S,
         scale * vL,
         scale * eL / rhoL + (S - uL) * (S + (pL) / (rhoL * (SL - uL)))
      }; 
      for(int i = 0; i < 4; i++)
      {
         HLLCFlux[i] = FL[i] + SL * (UStarL[i] - UL[i]);
      }
      
      // Difference between center (star region) and state
      /*
      for(int i = 0; i < 4; i++)
      {
         double add = SL * (pL + rhoL * (SL - uL) * (S - uL) * D[i]);
         HLLCFlux[i] = (S * (SL * UL[i] - FL[i]) + add) / (SL - S);
      }
      */
   }
   else if(S < 0 && 0 <= SR)
   {
      
      double scale = rhoR * (SR - uR) / (SR - S);
      // Difference between center (star region) and state
      std::vector<double> UStarR = {
         scale * 1.0,
         scale * S,
         scale * vR,
         scale * eR / rhoR + (S - uR) * (S + (pR) / (rhoR * (SR - uR)))
      }; 
      for(int i = 0; i < 4; i++)
      {
         HLLCFlux[i] = FR[i] + SR * (UStarR[i] - UR[i]);
      }
      
      /*
      for(int i = 0; i < 4; i++)
      {
         double add = SR * (pR + rhoR * (SR - uR) * (S - uR) * D[i]);
         HLLCFlux[i] = (S * (SR * UR[i] - FR[i]) + add) / (SR - S);
      }
      */
   }
   else if(0 >= SR)
   {
      HLLCFlux = FR;
   }

   // Map flux depending on vertical or horizontal calculation
   if(dir == 0)
   {
      return HLLCFlux;
   }
   else
   {
      return {HLLCFlux[0], HLLCFlux[2], HLLCFlux[1], HLLCFlux[3]};
   }
}

void isentropicVortex::applyBC(grid2D& grid)
{
   int ghostCells = grid.getGhostCells();
   // LEFT
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = 0; i < ghostCells; i++)
      {
         cell& targetCell = grid.getCell(grid.getNx() - 2 * ghostCells + i, j);
         cell& ghostCell = grid.getCell(i, j);
         ghostCell = targetCell;
      }
   }

   // RIGHT 
   for(int j = ghostCells; j < grid.getNy() - ghostCells; j++)
   {
      for(int i = grid.getNx() - ghostCells; i < grid.getNx(); i++)
      {
         cell& targetCell = grid.getCell(i - (grid.getNx() - 2 * ghostCells), j);
         cell& ghostCell = grid.getCell(i, j);
         ghostCell = targetCell;
      }
   }

   // TOP 
   for(int j = ghostCells; j < grid.getNx() - ghostCells; j++)
   {
      for(int i = grid.getNy() - ghostCells; i < grid.getNy(); i++)
      {
         cell& targetCell = grid.getCell(j, i - (grid.getNy() - 2 * ghostCells));
         cell& ghostCell = grid.getCell(j, i);
         ghostCell = targetCell;
      }
   }

   // BOTTOM 
   for(int j = ghostCells; j < grid.getNx() - ghostCells; j++)
   {
      for(int i = 0; i < ghostCells; i++)
      {
         cell& targetCell = grid.getCell(j, grid.getNy() - 2 * ghostCells + i);
         cell& ghostCell = grid.getCell(j, i);
         ghostCell = targetCell;
      }
   }
};

void isentropicVortex::consToPrim(std::vector<double>& state)
{

   double rho = state[0];
   double energy = state[3];

   state[1] /= rho;
   state[2] /= rho;
   state[3] = (m_gamma - 1) * (energy - 0.5 * rho * (state[1] * state[1] + state[2] * state[2]));
};

void isentropicVortex::primToCons(std::vector<double>& state)
{
   double rho = state[0];
   double pres = state[3];
   double energy = pres / (m_gamma - 1) + 0.5 * rho * (state[1] * state[1] + state[2] * state[2]);

   state[1] *= rho;
   state[2] *= rho;
   state[3] = energy;
};

