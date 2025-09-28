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

      std::vector<double>& varState(int i)
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
            default:
            throw std::out_of_range("Invalid index for varState.");
         }
      }

   private:
      std::vector<double> m_vars;
      std::vector<double> m_varsR;
      std::vector<double> m_varsL;
      std::vector<double> m_varsU;
      std::vector<double> m_varsD;

};
