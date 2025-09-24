#pragma once
#include <iostream>
#include <vector>
#include <math>

class Cell
{
   public:
      Cell(){}

      double& var(int i) {return m_vars[i];}

   private:
      std::vector<double> m_vars;
}
