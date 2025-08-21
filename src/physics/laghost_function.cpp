#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <cmath>
#include "laghost_function.hpp"

namespace mfem
{
   double rho0(const Vector &x)
   {
      return 2800.0; // This case in initialized in main().
   }

   double gamma_func(const Vector &x)
   {
      return 1.0; // This case in initialized in main().
   }

   void v0(const Vector &x, Vector &v)
   {
      // const double atn = dim!=1 ? pow((x(0)*(1.0-x(0))*4*x(1)*(1.0-x(1))*4.0),W
      //                               0.4) : 0.0;
      // const double s = 0.1/64.;
      // v(0) = 0.0;
      // v(1) = 1.0e-16;

      v =0.0;
      // if(x(0) <= 1)
      // {
      //    v(0) = -1*0.003/86400/365.25;
      // } 

      // if(x(1) >= 35e3)
      // {
      //    v(1) = -1e-20;
      // } 
   }

   double e0(const Vector &x)
   {
      return 0.0; // This case in initialized in main().
   }

   double p0(const Vector &x)
   {
      double r = sqrt(pow((x(0)-50e3), 2) + pow((x(1)-2e3), 2));

      if(r <= 1.0e3)
      {
         return 0.5;
      }

      return 0.0;
   }

   double depth0(const Vector &x)
   {
      double r = sqrt(pow((x(0)-200e3/2), 2) + pow((x(1)-100e3/2), 2));

      if(r <= 10.0e3)
      {
         return 1;
      }

      return 0.0;
   }

   void surf_norm(const Vector &x, Vector &v)
   {
      if(x.Size() == 2)
      {
         v(0) = x(0);
         v(1) = x(1);
      } 
      else
      {
         v(0) = x(0);
         v(1) = x(1);
         v(2) = x(2);
      }
   }

   double zero_func(const Vector &x)
   {
      return 0.0;
   }

   void xyz0(const Vector &x, Vector &v)
   {
      if(x.Size() == 2)
      {
         v(0) = x(0);
         v(1) = x(1);
      } 
      else
      {
         v(0) = x(0);
         v(1) = x(1);
         v(2) = x(2);
      }
   }

   double x_l2(const Vector &x)
   {
      return x(0);
   }

   double y_l2(const Vector &x)
   {
      return x(1);
   }

   double z_l2(const Vector &x)
   {
      return x.Size() >= 3 ? x(2) : 0.0;
   }

}
