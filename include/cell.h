#pragma once
#include <iostream>
#include <vector>
#include <cmath>

class cell
{
   public:
      cell(){m_vars.resize(4);}

      // copy constructor
      cell(const cell& other)
      {
         m_vars = other.m_vars;
      }

      // Assignment operator
      cell& operator=(const cell& other)
      {
         if(this != &other)
         {
            m_vars = other.m_vars;
         }

         return *this;
      }

      double& var(int i) {return m_vars[i];}
      double& varR(int i) {return m_varsR[i];}
      double& varL(int i) {return m_varsL[i];}
      double& varU(int i) {return m_varsU[i];}
      double& varD(int i) {return m_varsD[i];}
      // 
      // 

   private:
      std::vector<double> m_vars;
      std::vector<double> m_varsR;
      std::vector<double> m_varsL;
      std::vector<double> m_varsU;
      std::vector<double> m_varsD;

};
