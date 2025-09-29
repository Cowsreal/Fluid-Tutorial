#pragma once
#include <iostream>
#include <vector>
#include <cmath>

class cell
{
   public:
      cell();

      // copy constructor
      cell(const cell& other)
      {
         m_vars = other.m_vars;
         m_varsR = other.m_varsR;
         m_varsL = other.m_varsL;
         m_varsU = other.m_varsU;
         m_varsD = other.m_varsD;
         m_flux = other.m_flux;
      }

      // Assignment operator
      cell& operator=(const cell& other)
      {
         if(this != &other)
         {
            m_vars = other.m_vars;
            m_varsR = other.m_varsR;
            m_varsL = other.m_varsL;
            m_varsU = other.m_varsU;
            m_varsD = other.m_varsD;
            m_flux = other.m_flux;
         }

         return *this;
      }

      double& var(int i) {return m_vars[i];}
      double& varR(int i) {return m_varsR[i];}
      double& varL(int i) {return m_varsL[i];}
      double& varU(int i) {return m_varsU[i];}
      double& varD(int i) {return m_varsD[i];}

      std::vector<double>& varState(int i);

   private:
      std::vector<double> m_vars;
      std::vector<double> m_varsR;
      std::vector<double> m_varsL;
      std::vector<double> m_varsU;
      std::vector<double> m_varsD;

      std::vector<double> m_flux;
};
