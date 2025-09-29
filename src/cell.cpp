#include "cell.h"

cell::cell()
{
   m_vars.resize(4);
   m_varsR.resize(4);
   m_varsL.resize(4);
   m_varsU.resize(4);
   m_varsD.resize(4);
   m_flux.resize(4);
}

std::vector<double>& cell::varState(int i)
{
   switch (i)
   {
      case 0:
         return m_vars;
      case 1:
         return m_varsU;
      case 2:
         return m_varsD;
      case 3:
         return m_varsL;
      case 4:
         return m_varsR;
      case 5:
         return m_flux;
      default:
      throw std::out_of_range("Invalid index for varState.");
   }
}
