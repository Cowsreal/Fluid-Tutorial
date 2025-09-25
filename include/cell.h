#pragma once
#include <iostream>
#include <vector>
#include <cmath>

class cell
{
   public:
      cell(){m_vars.resize(4);}

      double& var(int i) {return m_vars[i];}
      // 
      // 

   private:
      std::vector<double> m_vars;
};
