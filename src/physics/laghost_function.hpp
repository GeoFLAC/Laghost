#include "mfem.hpp"
#include <ctime>
#include <cstdlib>

namespace mfem
{
   double e0(const Vector &);
   double p0(const Vector &);
   double depth0(const Vector &);
   double rho0(const Vector &);
   double gamma_func(const Vector &);
   double zero_func(const Vector &);
   void v0(const Vector &, Vector &);
   void xyz0(const Vector &, Vector &);
   double x_l2(const Vector &);
   double y_l2(const Vector &);
   double z_l2(const Vector &);

   static void principal_stresses2(const double* s, double p[2],
                                double& cos2t, double& sin2t)
   {
      /* 's' is a flattened stress vector, with the components {XX, ZZ, XZ}.
      * Returns the eigenvalues 'p', and the direction cosine of the
      * eigenvectors in the X-Z plane.
      * The eigenvalues are ordered such that p[0] <= p[1].
      */

      // center and radius of Mohr circle
      double s0 = 0.5 * (s[0] + s[1]);
      double rad =  std::sqrt(0.25*(s[0] -  s[1])*(s[0] -  s[1]) + s[2]*s[2]);

      // principal stresses in the X-Z plane
      p[0] = s0 - rad;
      p[1] = s0 + rad;

      {
         // direction cosine and sine of 2*theta
         const double eps = 1e-15;
         double a = 0.5 * (s[0] - s[1]);
         double b = - rad; // always negative
         if (b < -eps) {
               cos2t = a / b;
               sin2t = s[2] / b;
         }
         else {
               cos2t = 1;
               sin2t = 0;
         }
      }
   }

   class PlasticCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &xyz;
      int dim;
      Vector location;
      double rad, ini_pls;

   public:
      PlasticCoefficient (int &_dim, ParGridFunction &_xyz, Vector &_location, double &_rad, double &_ini_pls)
         : VectorCoefficient(1), xyz(_xyz)
         {
            dim=_dim; location = _location; rad = _rad; ini_pls = _ini_pls;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(1);
         double r = 0.0;
         double xc = xyz.GetValue(T, ip,1);
         double yc = xyz.GetValue(T, ip,2);
         double zc = 0.0;
         double randomNumber;
         
         double xc_0 = 0.0;
         double yc_0 = 0.0;
         double zc_0 = 0.0;
         double dist = 0.0; 
         if(dim == 3){zc = xyz.GetValue(T, ip,3);}
         
         if(dim == 2){ r = sqrt(pow((xc-location[0]), 2) + pow((yc-location[1]), 2));}
         else if(dim == 3){r = sqrt(pow((xc-location[0]), 2) + pow((yc-location[1]), 2) + pow((zc-location[2]), 2));}

         if(r <= rad)
         {
            K(0) = ini_pls;
         }
         else
         {
            K(0) = 0.0;
         }

         // // Define the line's angle (in radians)
         // double angle_rad = M_PI / 3; // 60 degrees in radians

         // // Define the slope and intercept of the line
         // double m = tan(angle_rad);
         // double b = 0; // Assuming the line passes through the origin

         // // Define the width of the band
         // double width = 3000;
         // xc_0 = xc-location[0]; yc_0 = yc-location[1];

         // dist = fabs((m * xc_0 - yc_0 + b) / sqrt(m * m + 1)); // distance between points and line

         // if ((dist <= width / 2) & (yc <= 5000)) 
         // {
         //    K(0) = ini_pls;
         // }
         // else
         // {
         //    K(0) = 0.0;
         // }

         // if((xc >= 49375.0) & (xc <= 50625.0) & (yc <= 1250.0))
         // {
         //    K(0) = ini_pls;
         // }
         // else
         // {
         //    K(0) = 0.0;
         // }

         // srand(time(0));
         // randomNumber =  (rand() % 101)/100.0;
         // if((xc >= 97.5e3) & (xc <= 102.5e3) & (yc <= 5e3))
         // // if((xc >= 80e3) & (xc <= 120e3))
         // // if(xc <= 100e3)
         // {
         //    // K(0) = ini_pls*randomNumber;
         //    K(0) = ini_pls;
         //    // std::cout << randomNumber << std::endl;
         // }
         // else
         // {
         //    K(0) = 0.0;
         // }

         // K(0) = ini_pls*randomNumber;
         
      }
      virtual ~PlasticCoefficient() { }
   };

   class LithostaticCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &xyz, &rho;
      int dim;
      double gravity, thickness;

   public:
      LithostaticCoefficient (int &_dim, ParGridFunction &_xyz, ParGridFunction &_rho, double &_gravity, double &_thickness)
         : VectorCoefficient(3*(_dim-1)), xyz(_xyz), rho(_rho)  
         {
            dim=_dim; gravity = _gravity; thickness = _thickness; 
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(3*(dim-1));

         double xc = xyz.GetValue(T, ip, 1);
         double zc = xyz.GetValue(T, ip, dim);
         double denc = rho.GetValue(T, ip);
         double atm = -101325; // 1 atm in Pa

         K = 0.0;
         for (int i = 0; i < dim; i++)
            K(i) = -1.0*fabs(thickness - zc)*denc*gravity + atm;
      }
      virtual ~LithostaticCoefficient() { }
   };

   class ATMCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &xyz, &rho;
      int dim;
      double gravity, thickness;

   public:
      ATMCoefficient (int &_dim, ParGridFunction &_xyz, ParGridFunction &_rho, double &_gravity, double &_thickness)
         : VectorCoefficient(3*(_dim-1)), xyz(_xyz), rho(_rho)  
         {
            dim=_dim; gravity = _gravity; thickness = _thickness; 
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(3*(dim-1));
         // K.SetSize(dim);

         double atm = -101325; // 1 atm in Pa
         // atm =0.0;
         double zc = xyz.GetValue(T, ip, dim);
         // 
         K = 0.0;
         for (int i = 0; i < dim; i++)
            K(i) = atm;
      }
      virtual ~ATMCoefficient() { }
   };

   class CompoCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &mat;
      int mat_num;
      int comp;

   public:
      CompoCoefficient (int &_mat_num, ParGridFunction &_mat)
         : VectorCoefficient(_mat_num), mat(_mat)  
         {
           mat_num=_mat_num; comp = 0;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(mat_num);
         K = 0.0;
         comp = mat.GetValue(T, ip);
         for (int i = 0; i < mat_num; i++)
         {
            if(comp == i)
            {
               K(i) = 1.0;
            }
         }

      }
      virtual ~CompoCoefficient() { }
   };
   
   class CompMassCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &comp_ref_gf, vol_ini_gf, quality;
      int mat_num;

   public:
      CompMassCoefficient (int &_mat_num, ParGridFunction &_comp_ref_gf, ParGridFunction &_vol_ini_gf, ParGridFunction &_quality)
         : VectorCoefficient(_mat_num), comp_ref_gf(_comp_ref_gf), vol_ini_gf(_vol_ini_gf), quality(_quality)
         {
           mat_num=_mat_num;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(mat_num);
         K = 0.0;
         for (int i = 0; i < mat_num; i++)
         {
            K(i) = comp_ref_gf.GetValue(T, ip, i+1) * vol_ini_gf.GetValue(T, ip) / quality.GetValue(T, ip, 1);
         }
      }
      virtual ~CompMassCoefficient() { }
   };

   class Temp_steessCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &s_gf;
      int component;

   public:
      Temp_steessCoefficient (int &_component, ParGridFunction &_s_gf)
         : VectorCoefficient(_component), s_gf(_s_gf)
         {
           component=_component;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         
         K(0) = s_gf.GetValue(T, ip, component);
      }
      virtual ~Temp_steessCoefficient() { }
   };

   class Temp_compCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &comp_gf;
      int component;

   public:
      Temp_compCoefficient (int &_component, ParGridFunction &_comp_gf)
         : VectorCoefficient(_component), comp_gf(_comp_gf)
         {
           component=_component;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K(0) = comp_gf.GetValue(T, ip, component+1);
      }
      virtual ~Temp_compCoefficient() { }
   };

   class StressMappingCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &temp_gf;
      int dim;
      
   public:
      StressMappingCoefficient (int &_dim, ParGridFunction &_temp_gf)
         : VectorCoefficient(3*(_dim-1)), temp_gf(_temp_gf)
         {
           dim = _dim; 
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(3*(dim-1));
         if(dim == 2)
         {
            K(0) = temp_gf.GetValue(T, ip, 1);
            K(1) = temp_gf.GetValue(T, ip, 2);
            K(2) = temp_gf.GetValue(T, ip, 3);
         }
         else
         {
            K(0) = temp_gf.GetValue(T, ip, 1);
            K(1) = temp_gf.GetValue(T, ip, 2);
            K(2) = temp_gf.GetValue(T, ip, 3);
            K(3) = temp_gf.GetValue(T, ip, 4);
            K(4) = temp_gf.GetValue(T, ip, 5);
            K(5) = temp_gf.GetValue(T, ip, 6);
         }
         
      }
      virtual ~StressMappingCoefficient() { }
   };

   class PlasticityMappingCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &temp_gf;
      int dim;
      
   public:
      PlasticityMappingCoefficient (int &_dim, ParGridFunction &_temp_gf)
         : VectorCoefficient(1), temp_gf(_temp_gf)
         {
           dim = _dim; 
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(1);
         if(dim == 2)
         {
            K(0) = temp_gf.GetValue(T, ip, 4);
         }
         else
         {
            K(0) = temp_gf.GetValue(T, ip, 7);
         }
         
      }
      virtual ~PlasticityMappingCoefficient() { }
   };

   class ReturnMapping3DCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &comp_gf, &s_gf, &s_old_gf, &p_gf, &mat_gf;
      int dim;
      double h_min, dt_old;
      bool viscoplastic;
      Vector rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity;

   public:
      ReturnMapping3DCoefficient (int &_dim, double &_h_min, double &_dt_old, bool &_viscoplastic, ParGridFunction &_comp_gf, ParGridFunction &_s_gf, ParGridFunction &_s_old_gf, ParGridFunction &_p_gf, ParGridFunction &_mat_gf, Vector &_rho, Vector &_lambda, Vector &_mu, Vector &_tension_cutoff, Vector &_cohesion0, Vector &_cohesion1, Vector &_pls0, Vector &_pls1, Vector &_friction_angle0, Vector &_friction_angle1, Vector &_dilation_angle0, Vector &_dilation_angle1, Vector &_plastic_viscosity)
         : VectorCoefficient(_dim), comp_gf(_comp_gf), s_gf(_s_gf), s_old_gf(_s_old_gf), p_gf(_p_gf), mat_gf(_mat_gf)
         {
           dim = _dim; h_min = _h_min; dt_old = _dt_old; viscoplastic = _viscoplastic;
           rho = _rho; lambda = _lambda; mu = _mu; tension_cutoff = _tension_cutoff; cohesion0 = _cohesion0; cohesion1 = _cohesion1; pls0 = _pls0; pls1 = _pls1; friction_angle0 = _friction_angle0; friction_angle1 = _friction_angle1; dilation_angle0 = _dilation_angle0; dilation_angle1 = _dilation_angle1; plastic_viscosity = _plastic_viscosity;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(3*(dim-1)+1);
         K = 0.0;
         DenseMatrix esig(3);
         DenseMatrix esig_old(3);
         DenseMatrix esig_inc(3);
         DenseMatrix plastic_sig(3);
         DenseMatrix plastic_str(3);
         esig=0.0; plastic_sig=0.0; plastic_str=0.0;
         double eig_sig_var[3], eig_sig_vec[9];

         double sig1{0.0};
         double sig3{0.0};
         double  syy{0.0}; // Syy is non-zero value in plane strain condition
         double msig{0.0}; // mean stress
         double evol{0.0}; // volumetric strain
         double DEG2RAD{M_PI/180.0};
         double depls{0.0}; // 2nd invariant of plastic strain

         double fs{0.0};
         double ft{0.0};
         double fh{0.0};
         double N_phi{0.0};
         double st_N_phi{0.0};
         double N_psi{0.0};
         double beta{0.0};
         double viscosity{0.0};
         double relax{0.0};
         // double relax_limit{1.0};
         double dt_scaled{0.0};
         double numerator   = {0.0};
         double denominator = {0.0};
         double pls_old = {0.0}; // cumulative 2nd invariant of plastic strain
         double p_slope = {0.0}; 
         double fri_str = {0.0}; // strain_dependent friction angle
         double dil_str = {0.0}; // strain_dependent dilation angle
         double coh_str = {0.0}; // strain_dependent cohesion
         double ten_cut = {0.0};
         int mat{0};
         int nsize{mat_gf.Size()};
         int mat_num{lambda.Size()};
      
         double pls0_c = {0.0}; 
         double pls1_c = {0.0}; 
         double rho_c = {0.0}; 
         double lambda_c = {0.0}; 
         double mu_c = {0.0}; 
         double tension_cutoff_c = {0.0}; 
         double cohesion0_c = {0.0}; 
         double cohesion1_c = {0.0}; 
         double friction_angle0_c ={0.0}; 
         double friction_angle1_c ={0.0}; 
         double dilation_angle0_c ={0.0}; 
         double dilation_angle1_c ={0.0}; 
         double plastic_viscosity_c ={0.0};
         double pwave_speed ={0.0};
         double time_scale ={1.0};
         
         // element-wise loop
         esig=0.0; plastic_sig=0.0; plastic_str=0.0;
         esig_old=0.0; esig_inc=0.0;
         pls_old = p_gf.GetValue(T, ip);
         if(pls_old < 0.0){pls_old =0.0;}

         pls0_c =0.0; pls1_c =0.0; rho_c = 0.0; lambda_c = 0.0; mu_c = 0.0; time_scale = 1.0;
         tension_cutoff_c = 0.0; cohesion0_c = 0.0; cohesion1_c = 0.0; friction_angle0_c = 0.0; friction_angle1_c = 0.0;
         dilation_angle0_c = 0.0; dilation_angle1_c = 0.0; plastic_viscosity_c = 0.0;

         for( int i = 0; i < mat_num; i++ )
         {
            pls0_c = pls0_c + comp_gf.GetValue(T, ip, i+1)*pls0[i];
            pls1_c = pls1_c + comp_gf.GetValue(T, ip, i+1)*pls1[i];
            rho_c = rho_c + comp_gf.GetValue(T, ip, i+1)*rho[i];
            lambda_c = lambda_c + comp_gf.GetValue(T, ip, i+1)*lambda[i];
            mu_c = mu_c + comp_gf.GetValue(T, ip, i+1)*mu[i];
            tension_cutoff_c = tension_cutoff_c + comp_gf.GetValue(T, ip, i+1)*tension_cutoff[i];
            cohesion0_c = cohesion0_c + comp_gf.GetValue(T, ip, i+1)*cohesion0[i];
            cohesion1_c = cohesion1_c + comp_gf.GetValue(T, ip, i+1)*cohesion1[i];
            friction_angle0_c = friction_angle0_c + comp_gf.GetValue(T, ip, i+1)*friction_angle0[i];
            friction_angle1_c = friction_angle1_c + comp_gf.GetValue(T, ip, i+1)*friction_angle1[i];
            dilation_angle0_c = dilation_angle0_c + comp_gf.GetValue(T, ip, i+1)*dilation_angle0[i];
            dilation_angle1_c = dilation_angle1_c + comp_gf.GetValue(T, ip, i+1)*dilation_angle1[i];
         }

         // linear weakening
         p_slope = (pls_old - pls0_c)/(pls1_c - pls0_c);
         pwave_speed = sqrt((lambda_c + 2*mu_c)/rho_c);
         if(h_min  > 0){time_scale = h_min / pwave_speed;}
         plastic_viscosity_c = time_scale * mu_c;
         
         if(dim ==2) 
         {
            // 2D in plane strain condition
            // sxx, syy, szz, sxz are non zeros.
            // sxy, syz are zero.
            msig = (s_gf.GetValue(T, ip, 1) + s_gf.GetValue(T, ip, 2))*0.5;
            evol = msig/(lambda_c+mu_c);
            syy  = evol * lambda_c;
            esig(0,0) = s_gf.GetValue(T, ip, 1); esig(0,1) = 0.0; esig(0,2) = s_gf.GetValue(T, ip, 3); 
            esig(1,0) =                     0,0; esig(1,1) = syy; esig(1,2) =                     0.0;
            esig(2,0) = s_gf.GetValue(T, ip, 3); esig(2,1) = 0.0; esig(2,2) = s_gf.GetValue(T, ip, 2);

            // Caushy stress at previous time step
            msig = (s_old_gf.GetValue(T, ip, 1) + s_old_gf.GetValue(T, ip, 2))*0.5;
            evol = msig/(lambda_c+mu_c);
            syy  = evol * lambda_c;
            esig_old(0,0) = s_old_gf.GetValue(T, ip, 1); esig_old(0,1) = 0.0; esig_old(0,2) = s_old_gf.GetValue(T, ip, 3); 
            esig_old(1,0) =                         0,0; esig_old(1,1) = syy; esig_old(1,2) =                         0.0;
            esig_old(2,0) = s_old_gf.GetValue(T, ip, 3); esig_old(2,1) = 0.0; esig_old(2,2) = s_old_gf.GetValue(T, ip, 2);
         }
         else
         {
            esig(0,0) = s_gf.GetValue(T, ip, 1); esig(0,1) = s_gf.GetValue(T, ip, 4); esig(0,2) = s_gf.GetValue(T, ip, 5); 
            esig(1,0) = s_gf.GetValue(T, ip, 4); esig(1,1) = s_gf.GetValue(T, ip, 2); esig(1,2) = s_gf.GetValue(T, ip, 6);
            esig(2,0) = s_gf.GetValue(T, ip, 5); esig(2,1) = s_gf.GetValue(T, ip, 6); esig(2,2) = s_gf.GetValue(T, ip, 3);

            // Caushy stress at previous time step
            esig_old(0,0) = s_gf.GetValue(T, ip, 1); esig_old(0,1) = s_gf.GetValue(T, ip, 4); esig_old(0,2) = s_gf.GetValue(T, ip, 5); 
            esig_old(1,0) = s_gf.GetValue(T, ip, 4); esig_old(1,1) = s_gf.GetValue(T, ip, 2); esig_old(1,2) = s_gf.GetValue(T, ip, 6);
            esig_old(2,0) = s_gf.GetValue(T, ip, 5); esig_old(2,1) = s_gf.GetValue(T, ip, 6); esig_old(2,2) = s_gf.GetValue(T, ip, 3);
         }

         // Elastic stress increment   
         esig_inc(0,0) = esig(0,0) - esig_old(0,0); esig_inc(0,1) = esig(0,1) - esig_old(0,1); esig_inc(0,2) = esig(0,2) - esig_old(0,2); 
         esig_inc(1,0) = esig(1,0) - esig_old(1,0); esig_inc(1,1) = esig(1,1) - esig_old(1,1); esig_inc(1,2) = esig(1,2) - esig_old(1,2);
         esig_inc(2,0) = esig(2,0) - esig_old(2,0); esig_inc(2,1) = esig(2,1) - esig_old(2,1); esig_inc(2,2) = esig(2,2) - esig_old(2,2);
         
         esig.CalcEigenvalues(eig_sig_var, eig_sig_vec); 

         Vector sig_var(eig_sig_var, 3);
         Vector sig_dir(eig_sig_vec, 3);

         auto max_it = std::max_element(sig_var.begin(), sig_var.end()); // find iterator to max element
         auto min_it = std::min_element(sig_var.begin(), sig_var.end()); // find iterator to min element
         
         int max_index = std::distance(sig_var.begin(), max_it); // calculate index of max element
         int min_index = std::distance(sig_var.begin(), min_it); // calculate index of min element
         
         int itm_index = 0; // calculate index of intermediate element
         if (max_index + min_index == 1) {itm_index = 2;}
         else if(max_index + min_index == 2) {itm_index = 1;}
         else {itm_index = 0;}

         sig1 = sig_var[min_index]; // most compressive pincipal stress
         sig3 = sig_var[max_index]; // least compressive pincipal stress

         // linear strain weaking on cohesion, friction and dilation angles.
         coh_str = cohesion0_c; fri_str = friction_angle0_c; dil_str = dilation_angle0_c;

         if (pls_old < pls0_c) {
            // no weakening yet
            coh_str = cohesion0_c;
            fri_str = friction_angle0_c;
            dil_str = dilation_angle0_c;
         }
         else if (pls_old < pls1_c) {
            // linear weakening
            coh_str = cohesion0_c + p_slope * (cohesion1_c - cohesion0_c);
            fri_str = friction_angle0_c + p_slope * (friction_angle1_c - friction_angle0_c);
            dil_str = dilation_angle0_c + p_slope * (dilation_angle1_c - dilation_angle0_c);
         }
         else {
            // saturated weakening
            coh_str = cohesion1_c; fri_str = friction_angle1_c; dil_str = dilation_angle1_c;
         }

         N_phi = (1+sin(DEG2RAD*fri_str))/(1-sin(DEG2RAD*fri_str));
         st_N_phi = cos(DEG2RAD*fri_str)/(1-sin(DEG2RAD*fri_str));
         N_psi = -1*(1+sin(DEG2RAD*dil_str))/(1-sin(DEG2RAD*dil_str)); // partial_g/partial_sig3

         if(tension_cutoff_c == 0)
         {
            ten_cut = coh_str/tan(DEG2RAD*fri_str);
         }
         else{ten_cut = tension_cutoff_c;}

         // shear failure function
         fs = sig1 - N_phi*sig3 + 2*coh_str*st_N_phi;
         // tension failure function
         ft = sig3 - ten_cut;
         // bisects the obtuse angle made by two yield function
         fh = sig3 - ten_cut + (sqrt(N_phi*N_phi + 1.0)+ N_phi)*(sig1 - N_phi*ten_cut + 2*coh_str*st_N_phi);

         depls = 0.0;
         if((fs < 0 & fh < 0) | (ft > 0 & fh > 0))
         { 
            if(fs < 0 & fh < 0) // stress correction at shear failure
            {
               // Equations 28 and 30 from Choi et al. (2013; DynEarthSol2D: An efficient unstructured finite element method to study long-term tectonic deformation). 
               beta = fs;
               beta = beta / (((lambda_c+2*mu_c)*1 - N_phi*lambda_c*1) + (lambda_c*N_psi - N_phi*(lambda_c+2*mu_c)*N_psi));
               
               plastic_str(0,0) = (lambda_c + 2*mu_c + lambda_c*N_psi) * beta; 
               plastic_str(1,1) = (lambda_c + lambda_c*N_psi) * beta;
               plastic_str(2,2) = (lambda_c + (lambda_c+2*mu_c)*N_psi) * beta;
               // reduced form of 2nd invariant
               if(dim ==2)
               {
                  depls = std::fabs(beta) * std::sqrt((3 - 2*N_psi + 3*N_psi*N_psi) / 8); 
               }
               else
               {
                  depls = std::fabs(beta) * std::sqrt((7 - 4*N_psi + 7*N_psi*N_psi) / 18);
               }
               
            }
            else if (ft > 0 & fh > 0) // stress correction at tension failure
            {
               beta = ft;
               beta = beta/(lambda_c+2*mu_c);

               plastic_str(0,0) = lambda_c * beta * 1;
               plastic_str(1,1) = lambda_c * beta * 1;
               plastic_str(2,2) = (lambda_c+2*mu_c) * beta * 1;

               // reduced form of 2nd invariant
               if(dim ==2)
               {
                  depls = std::fabs(beta) * std::sqrt(3. / 8);
               }
               else
               {
                  depls = std::fabs(beta) * std::sqrt(7. / 18);
               }       
            }

            // Rotating Principal axis to XYZ axis
            plastic_sig(0,0) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[0+min_index*3]*sig_dir[0+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[0+itm_index*3]*sig_dir[0+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[0+max_index*3]*sig_dir[0+max_index*3]);
            plastic_sig(0,1) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[0+min_index*3]*sig_dir[1+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[0+itm_index*3]*sig_dir[1+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[0+max_index*3]*sig_dir[1+max_index*3]);
            plastic_sig(0,2) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[0+min_index*3]*sig_dir[2+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[0+itm_index*3]*sig_dir[2+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[0+max_index*3]*sig_dir[2+max_index*3]);
            plastic_sig(1,0) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[1+min_index*3]*sig_dir[0+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[1+itm_index*3]*sig_dir[0+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[1+max_index*3]*sig_dir[0+max_index*3]);
            plastic_sig(1,1) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[1+min_index*3]*sig_dir[1+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[1+itm_index*3]*sig_dir[1+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[1+max_index*3]*sig_dir[1+max_index*3]);
            plastic_sig(1,2) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[1+min_index*3]*sig_dir[2+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[1+itm_index*3]*sig_dir[2+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[1+max_index*3]*sig_dir[2+max_index*3]);
            plastic_sig(2,0) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[2+min_index*3]*sig_dir[0+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[2+itm_index*3]*sig_dir[0+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[2+max_index*3]*sig_dir[0+max_index*3]);
            plastic_sig(2,1) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[2+min_index*3]*sig_dir[1+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[2+itm_index*3]*sig_dir[1+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[2+max_index*3]*sig_dir[1+max_index*3]);
            plastic_sig(2,2) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[2+min_index*3]*sig_dir[2+min_index*3] + (sig_var[itm_index]-plastic_str(1,1))*sig_dir[2+itm_index*3]*sig_dir[2+itm_index*3] + (sig_var[max_index]-plastic_str(2,2))*sig_dir[2+max_index*3]*sig_dir[2+max_index*3]);

            // Updating new stress to grid function
            viscosity = plastic_viscosity_c;
            relax = exp(-dt_old/viscosity);
            numerator  = fabs(1 - relax);
            denominator= fabs(dt_old/viscosity);
            dt_scaled = dt_old/viscosity;

            // Based on L'Hôpital's rule, 0/0 become 1 
            // if((numerator > 1e-15) & (denominator > 1e-15)){relax_limit = (1 - relax)/(dt_old/viscosity);}

         
            if(dim ==2)
            {
               if(viscoplastic)
               {
                  // Closed-form algorithm for viscoplasticity from Computational Inelasticity on p. 218 (Simo and Hughes, 1998)
                  // s_gf[i+nsize*0]=esig_old(0,0)*relax + (1.0 - relax)*plastic_sig(0,0) + relax_limit*esig_inc(0,0); 
                  // s_gf[i+nsize*2]=esig_old(0,2)*relax + (1.0 - relax)*plastic_sig(0,2) + relax_limit*esig_inc(0,2); 
                  // s_gf[i+nsize*2]=esig_old(2,0)*relax + (1.0 - relax)*plastic_sig(2,0) + relax_limit*esig_inc(2,0); 
                  // s_gf[i+nsize*1]=esig_old(2,2)*relax + (1.0 - relax)*plastic_sig(2,2) + relax_limit*esig_inc(2,2);

                  // Implicit backward Euler algorithm
                  K(0)=((esig_old(0,0) + esig_inc(0,0)) + dt_scaled*plastic_sig(0,0))/(1.0+dt_scaled);
                  K(1)=((esig_old(2,2) + esig_inc(2,2)) + dt_scaled*plastic_sig(2,2))/(1.0+dt_scaled);
                  K(2)=((esig_old(2,0) + esig_inc(2,0)) + dt_scaled*plastic_sig(2,0))/(1.0+dt_scaled);
                  depls = dt_scaled*depls/(1.0+dt_scaled); K(3) = depls;
               }
               else
               {
                  K(0)=plastic_sig(0,0); 
                  K(1)=plastic_sig(2,2);
                  K(2)=plastic_sig(0,2);
                  K(3)=depls;
               }
            }
            else
            {
               if(viscoplastic)
               {
                  // // Closed-form algorithm for viscoplasticity from Computational Inelasticity on p. 218 (Simo and Hughes, 1998)
                  // s_gf[i+nsize*0]=esig_old(0,0)*relax + (1.0 - relax)*plastic_sig(0,0) + relax_limit*esig_inc(0,0); 
                  // s_gf[i+nsize*3]=esig_old(0,1)*relax + (1.0 - relax)*plastic_sig(0,1) + relax_limit*esig_inc(0,1); 
                  // s_gf[i+nsize*4]=esig_old(0,2)*relax + (1.0 - relax)*plastic_sig(0,2) + relax_limit*esig_inc(0,2); 
                  // s_gf[i+nsize*3]=esig_old(1,0)*relax + (1.0 - relax)*plastic_sig(1,0) + relax_limit*esig_inc(1,0); 
                  // s_gf[i+nsize*1]=esig_old(1,1)*relax + (1.0 - relax)*plastic_sig(1,1) + relax_limit*esig_inc(1,1); 
                  // s_gf[i+nsize*5]=esig_old(1,2)*relax + (1.0 - relax)*plastic_sig(1,2) + relax_limit*esig_inc(1,2); 
                  // s_gf[i+nsize*4]=esig_old(2,0)*relax + (1.0 - relax)*plastic_sig(2,0) + relax_limit*esig_inc(2,0); 
                  // s_gf[i+nsize*5]=esig_old(2,1)*relax + (1.0 - relax)*plastic_sig(2,1) + relax_limit*esig_inc(2,1); 
                  // s_gf[i+nsize*2]=esig_old(2,2)*relax + (1.0 - relax)*plastic_sig(2,2) + relax_limit*esig_inc(2,2);

                  // Implicit backward Euler algorithm
                  K(0)=((esig_old(0,0) + esig_inc(0,0)) + dt_scaled*plastic_sig(0,0))/(1.0+dt_scaled);
                  K(1)=((esig_old(1,1) + esig_inc(1,1)) + dt_scaled*plastic_sig(1,1))/(1.0+dt_scaled);
                  K(2)=((esig_old(2,2) + esig_inc(2,2)) + dt_scaled*plastic_sig(2,2))/(1.0+dt_scaled);
                  K(3)=((esig_old(0,1) + esig_inc(0,1)) + dt_scaled*plastic_sig(0,1))/(1.0+dt_scaled);
                  K(4)=((esig_old(0,2) + esig_inc(0,2)) + dt_scaled*plastic_sig(0,2))/(1.0+dt_scaled);
                  K(5)=((esig_old(1,2) + esig_inc(1,2)) + dt_scaled*plastic_sig(1,2))/(1.0+dt_scaled);
                  depls = dt_scaled*depls/(1.0+dt_scaled); K(6) = depls;
               }
               else
               {
                  K(0)=plastic_sig(0,0); 
                  K(1)=plastic_sig(1,1); 
                  K(2)=plastic_sig(2,2);
                  K(3)=plastic_sig(0,1); 
                  K(4)=plastic_sig(0,2); 
                  K(5)=plastic_sig(1,2);
                  K(6)=depls;
               }
            }
         }
         else
         {
            if(dim ==2)
            {
               K(0)=esig(0,0); 
               K(1)=esig(2,2);
               K(2)=esig(0,2);
               K(3)=0.0;
            }
            else
            {
               K(0)=esig(0,0); 
               K(1)=esig(1,1); 
               K(2)=esig(2,2);
               K(3)=esig(0,1); 
               K(4)=esig(0,2); 
               K(5)=esig(1,2);
               K(6)=0.0;
            }
         }

      }
      virtual ~ReturnMapping3DCoefficient() { }
   };

   class ReturnMapping2DCoefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &comp_gf, &s_gf, &s_old_gf, &p_gf, &mat_gf;
      int dim;
      double h_min, dt_old;
      bool viscoplastic;
      Vector rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity;

   public:
      ReturnMapping2DCoefficient (int &_dim, double &_h_min, double &_dt_old, bool &_viscoplastic, ParGridFunction &_comp_gf, ParGridFunction &_s_gf, ParGridFunction &_s_old_gf, ParGridFunction &_p_gf, ParGridFunction &_mat_gf, Vector &_rho, Vector &_lambda, Vector &_mu, Vector &_tension_cutoff, Vector &_cohesion0, Vector &_cohesion1, Vector &_pls0, Vector &_pls1, Vector &_friction_angle0, Vector &_friction_angle1, Vector &_dilation_angle0, Vector &_dilation_angle1, Vector &_plastic_viscosity)
         : VectorCoefficient(_dim), comp_gf(_comp_gf), s_gf(_s_gf), s_old_gf(_s_old_gf), p_gf(_p_gf), mat_gf(_mat_gf)
         {
           dim = _dim; h_min = _h_min; dt_old = _dt_old; viscoplastic = _viscoplastic;
           rho = _rho; lambda = _lambda; mu = _mu; tension_cutoff = _tension_cutoff; cohesion0 = _cohesion0; cohesion1 = _cohesion1; pls0 = _pls0; pls1 = _pls1; friction_angle0 = _friction_angle0; friction_angle1 = _friction_angle1; dilation_angle0 = _dilation_angle0; dilation_angle1 = _dilation_angle1; plastic_viscosity = _plastic_viscosity;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(3*(dim-1)+1);
         K = 0.0;
         DenseMatrix esig(2);
         DenseMatrix esig_old(2);
         DenseMatrix esig_inc(2);
         DenseMatrix plastic_sig(2);
         DenseMatrix plastic_str(2);
         esig=0.0; plastic_sig=0.0; plastic_str=0.0;
         double eig_sig_var[2], eig_sig_vec[4];

         double sig1{0.0};
         double sig3{0.0};
         double  syy{0.0}; // Syy is non-zero value in plane strain condition
         double  syy_old{0.0}; // Syy is non-zero value in plane strain condition
         double msig{0.0}; // mean stress
         double evol{0.0}; // volumetric strain
         double DEG2RAD{M_PI/180.0};
         double depls{0.0}; // 2nd invariant of plastic strain

         double fs{0.0};
         double ft{0.0};
         double fh{0.0};
         double N_phi{0.0};
         double st_N_phi{0.0};
         double N_psi{0.0};
         double beta{0.0};
         double viscosity{0.0};
         double relax{0.0};
         // double relax_limit{1.0};
         double dt_scaled{0.0};
         double numerator   = {0.0};
         double denominator = {0.0};
         double pls_old = {0.0}; // cumulative 2nd invariant of plastic strain
         double p_slope = {0.0}; 
         double fri_str = {0.0}; // strain_dependent friction angle
         double dil_str = {0.0}; // strain_dependent dilation angle
         double coh_str = {0.0}; // strain_dependent cohesion
         double ten_cut = {0.0};
         int mat{0};
         int nsize{mat_gf.Size()};
         int mat_num{lambda.Size()};
      
         double pls0_c = {0.0}; 
         double pls1_c = {0.0}; 
         double rho_c = {0.0}; 
         double lambda_c = {0.0}; 
         double mu_c = {0.0}; 
         double tension_cutoff_c = {0.0}; 
         double cohesion0_c = {0.0}; 
         double cohesion1_c = {0.0}; 
         double friction_angle0_c ={0.0}; 
         double friction_angle1_c ={0.0}; 
         double dilation_angle0_c ={0.0}; 
         double dilation_angle1_c ={0.0}; 
         double plastic_viscosity_c ={0.0};
         double pwave_speed ={0.0};
         double time_scale ={1.0};

         // element-wise loop
         esig=0.0; plastic_sig=0.0; plastic_str=0.0;
         esig_old=0.0; esig_inc=0.0;
         pls_old = p_gf.GetValue(T, ip);
         if(pls_old < 0.0){pls_old =0.0;}

         pls0_c =0.0; pls1_c =0.0; rho_c = 0.0; lambda_c = 0.0; mu_c = 0.0; time_scale = 1.0;
         tension_cutoff_c = 0.0; cohesion0_c = 0.0; cohesion1_c = 0.0; friction_angle0_c = 0.0; friction_angle1_c = 0.0;
         dilation_angle0_c = 0.0; dilation_angle1_c = 0.0; plastic_viscosity_c = 0.0;

         for( int i = 0; i < mat_num; i++ )
         {
            pls0_c = pls0_c + comp_gf.GetValue(T, ip, i+1)*pls0[i];
            pls1_c = pls1_c + comp_gf.GetValue(T, ip, i+1)*pls1[i];
            rho_c = rho_c + comp_gf.GetValue(T, ip, i+1)*rho[i];
            lambda_c = lambda_c + comp_gf.GetValue(T, ip, i+1)*lambda[i];
            mu_c = mu_c + comp_gf.GetValue(T, ip, i+1)*mu[i];
            tension_cutoff_c = tension_cutoff_c + comp_gf.GetValue(T, ip, i+1)*tension_cutoff[i];
            cohesion0_c = cohesion0_c + comp_gf.GetValue(T, ip, i+1)*cohesion0[i];
            cohesion1_c = cohesion1_c + comp_gf.GetValue(T, ip, i+1)*cohesion1[i];
            friction_angle0_c = friction_angle0_c + comp_gf.GetValue(T, ip, i+1)*friction_angle0[i];
            friction_angle1_c = friction_angle1_c + comp_gf.GetValue(T, ip, i+1)*friction_angle1[i];
            dilation_angle0_c = dilation_angle0_c + comp_gf.GetValue(T, ip, i+1)*dilation_angle0[i];
            dilation_angle1_c = dilation_angle1_c + comp_gf.GetValue(T, ip, i+1)*dilation_angle1[i];
         }
         

         // linear weakening
         p_slope = (pls_old - pls0_c)/(pls1_c - pls0_c);
         pwave_speed = sqrt((lambda_c + 2*mu_c)/rho_c);
         if(h_min  > 0){time_scale = h_min / pwave_speed;}
         plastic_viscosity_c = time_scale * mu_c;

         msig = (s_gf.GetValue(T, ip, 1) + s_gf.GetValue(T, ip, 2))*0.5;
         evol = msig/(lambda_c+mu_c);
         syy  = evol * lambda_c;
         esig(0,0) = s_gf.GetValue(T, ip, 1); esig(0,1) = s_gf.GetValue(T, ip, 3);  
         esig(1,0) = s_gf.GetValue(T, ip, 3); esig(1,1) = s_gf.GetValue(T, ip, 2); 


         // std::cout << mat_gf.GetValue(T, ip) << " " << s_gf.GetValue(T, ip, 1) << " " << s_gf.GetValue(T, ip, 2) << " " << s_gf.GetValue(T, ip, 3) << std::endl;

         // Caushy stress at previous time step
         msig = (s_old_gf.GetValue(T, ip, 1) + s_old_gf.GetValue(T, ip, 2))*0.5;
         evol = msig/(lambda_c+mu_c);
         syy_old  = evol * lambda_c;
         esig_old(0,0) = s_old_gf.GetValue(T, ip, 1); esig_old(0,1) = s_old_gf.GetValue(T, ip, 3);  
         esig_old(1,0) = s_old_gf.GetValue(T, ip, 3); esig_old(1,1) = s_old_gf.GetValue(T, ip, 2); 

         // Elastic stress increment
         esig_inc(0,0) = esig(0,0) - esig_old(0,0); esig_inc(0,1) = esig(0,1) - esig_old(0,1); 
         esig_inc(1,0) = esig(1,0) - esig_old(1,0); esig_inc(1,1) = esig(1,1) - esig_old(1,1); 
         esig.CalcEigenvalues(eig_sig_var, eig_sig_vec); 

         Vector sig_var(eig_sig_var, 2);
         Vector sig_dir(eig_sig_vec, 2);

         auto max_it = std::max_element(sig_var.begin(), sig_var.end()); // find iterator to max element
         auto min_it = std::min_element(sig_var.begin(), sig_var.end()); // find iterator to min element
         
         int max_index = std::distance(sig_var.begin(), max_it); // calculate index of max element
         int min_index = std::distance(sig_var.begin(), min_it); // calculate index of min element
         
         sig1 = sig_var[min_index]; // most compressive pincipal stress
         sig3 = sig_var[max_index]; // least compressive pincipal stress

         // linear strain weaking on cohesion, friction and dilation angles.
         coh_str = cohesion0_c; fri_str = friction_angle0_c; dil_str = dilation_angle0_c;

         if (pls_old < pls0_c) {
            // no weakening yet
            coh_str = cohesion0_c;
            fri_str = friction_angle0_c;
            dil_str = dilation_angle0_c;
         }
         else if (pls_old < pls1_c) {
            // linear weakening
            coh_str = cohesion0_c + p_slope * (cohesion1_c - cohesion0_c);
            fri_str = friction_angle0_c + p_slope * (friction_angle1_c - friction_angle0_c);
            dil_str = dilation_angle0_c + p_slope * (dilation_angle1_c - dilation_angle0_c);
         }
         else {
            // saturated weakening
            coh_str = cohesion1_c; fri_str = friction_angle1_c; dil_str = dilation_angle1_c;
         }

         // std::cout << coh_str <<"," << p_slope << "," << fri_str << std::endl;
         N_phi = (1+sin(DEG2RAD*fri_str))/(1-sin(DEG2RAD*fri_str));
         st_N_phi = cos(DEG2RAD*fri_str)/(1-sin(DEG2RAD*fri_str));
         N_psi = -1*(1+sin(DEG2RAD*dil_str))/(1-sin(DEG2RAD*dil_str)); // partial_g/partial_sig3

         if(tension_cutoff_c == 0)
         {
            ten_cut = coh_str/tan(DEG2RAD*fri_str);
         }
         else{ten_cut = tension_cutoff_c;}

         // shear failure function
         fs = sig1 - N_phi*sig3 + 2*coh_str*st_N_phi;
         // tension failure function
         ft = sig3 - ten_cut;
         // bisects the obtuse angle made by two yield function
         fh = sig3 - ten_cut + (sqrt(N_phi*N_phi + 1.0)+ N_phi)*(sig1 - N_phi*ten_cut + 2*coh_str*st_N_phi);
         depls = 0.0;

         if((fs < 0 & fh < 0) | (ft > 0 & fh > 0))
         {  
            // std::cout <<"failure occurs" << std::endl;

            if(fs < 0 & fh < 0) // stress correction at shear failure
            {
               // Equations 28 and 30 from Choi et al. (2013; DynEarthSol2D: An efficient unstructured finite element method to study long-term tectonic deformation). 
               beta = fs;
               beta = beta / (((lambda_c+2*mu_c)*1 - N_phi*lambda_c*1) + (lambda_c*N_psi - N_phi*(lambda_c+2*mu_c)*N_psi));
               
               plastic_str(0,0) = (lambda_c + 2*mu_c + lambda_c*N_psi) * beta; 
               syy -= (lambda_c + lambda_c*N_psi) * beta;
               plastic_str(1,1) = (lambda_c + (lambda_c+2*mu_c)*N_psi) * beta;
               // reduced form of 2nd invariant
               depls = std::fabs(beta) * std::sqrt((3 - 2*N_psi + 3*N_psi*N_psi) / 8); 
               
            }
            else if (ft > 0 & fh > 0) // stress correction at tension failure
            {
               beta = ft;
               beta = beta/(lambda_c+2*mu_c);

               plastic_str(0,0) = lambda_c * beta * 1;
               syy -= lambda_c * beta * 1;
               plastic_str(1,1) = (lambda_c+2*mu_c) * beta * 1;

               // reduced form of 2nd invariant
               depls = std::fabs(beta) * std::sqrt(7. / 18);
            }

            // Rotating Principal axis to XYZ axis
            plastic_sig(0,0) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[0+min_index*2]*sig_dir[0+min_index*2]  + (sig_var[max_index]-plastic_str(1,1))*sig_dir[0+max_index*2]*sig_dir[0+max_index*2]);
            plastic_sig(0,1) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[0+min_index*2]*sig_dir[1+min_index*2]  + (sig_var[max_index]-plastic_str(1,1))*sig_dir[0+max_index*2]*sig_dir[1+max_index*2]);
            plastic_sig(1,0) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[1+min_index*2]*sig_dir[0+min_index*2]  + (sig_var[max_index]-plastic_str(1,1))*sig_dir[1+max_index*2]*sig_dir[0+max_index*2]);
            plastic_sig(1,1) = ((sig_var[min_index]-plastic_str(0,0))*sig_dir[1+min_index*2]*sig_dir[1+min_index*2]  + (sig_var[max_index]-plastic_str(1,1))*sig_dir[1+max_index*2]*sig_dir[1+max_index*2]);

            // Updating new stress to grid function
            viscosity = plastic_viscosity_c;
            relax = exp(-dt_old/viscosity);
            numerator  = fabs(1 - relax);
            denominator= fabs(dt_old/viscosity);
            dt_scaled = dt_old/viscosity;

            // Based on L'Hôpital's rule, 0/0 become 1 
            // if((numerator > 1e-15) & (denominator > 1e-15)){relax_limit = (1 - relax)/(dt_old/viscosity);}

         
            if(viscoplastic)
            {
               // Implicit backward Euler algorithm
               K(0)=((esig_old(0,0) + esig_inc(0,0)) + dt_scaled*plastic_sig(0,0))/(1.0+dt_scaled);
               K(1)=((esig_old(1,1) + esig_inc(1,1)) + dt_scaled*plastic_sig(1,1))/(1.0+dt_scaled);
               K(2)=((esig_old(0,1) + esig_inc(0,1)) + dt_scaled*plastic_sig(0,1))/(1.0+dt_scaled);
               depls = dt_scaled*depls/(1.0+dt_scaled); 
               K(3)=depls;
            }
            else
            {
               K(0)=plastic_sig(0,0);
               K(1)=plastic_sig(1,1);
               K(2)=plastic_sig(0,1);
               K(3)=depls;
            }
         }
         else
         {
            K(0)=esig(0,0);
            K(1)=esig(1,1);
            K(2)=esig(0,1);
            K(3)=0.0;
         }

      }
      virtual ~ReturnMapping2DCoefficient() { }
   };

   class ReturnMapping2D_simple_Coefficient : public VectorCoefficient
   {
   private:
      ParGridFunction &comp_gf, &s_gf, &s_old_gf, &p_gf, &mat_gf;
      int dim;
      double h_min, dt_old;
      bool viscoplastic;
      Vector rho, lambda, mu, tension_cutoff, cohesion0, cohesion1, pls0, pls1, friction_angle0, friction_angle1, dilation_angle0, dilation_angle1, plastic_viscosity;

   public:
      ReturnMapping2D_simple_Coefficient (int &_dim, double &_h_min, double &_dt_old, bool &_viscoplastic, ParGridFunction &_comp_gf, ParGridFunction &_s_gf, ParGridFunction &_s_old_gf, ParGridFunction &_p_gf, ParGridFunction &_mat_gf, Vector &_rho, Vector &_lambda, Vector &_mu, Vector &_tension_cutoff, Vector &_cohesion0, Vector &_cohesion1, Vector &_pls0, Vector &_pls1, Vector &_friction_angle0, Vector &_friction_angle1, Vector &_dilation_angle0, Vector &_dilation_angle1, Vector &_plastic_viscosity)
         : VectorCoefficient(_dim), comp_gf(_comp_gf), s_gf(_s_gf), s_old_gf(_s_old_gf), p_gf(_p_gf), mat_gf(_mat_gf)
         {
           dim = _dim; h_min = _h_min; dt_old = _dt_old; viscoplastic = _viscoplastic;
           rho = _rho; lambda = _lambda; mu = _mu; tension_cutoff = _tension_cutoff; cohesion0 = _cohesion0; cohesion1 = _cohesion1; pls0 = _pls0; pls1 = _pls1; friction_angle0 = _friction_angle0; friction_angle1 = _friction_angle1; dilation_angle0 = _dilation_angle0; dilation_angle1 = _dilation_angle1; plastic_viscosity = _plastic_viscosity;
         }
      virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
      {
         K.SetSize(3*(dim-1)+1);
         K = 0.0;
         DenseMatrix esig(2);
         DenseMatrix esig_old(2);
         DenseMatrix esig_inc(2);
         DenseMatrix plastic_sig(2);
         DenseMatrix plastic_str(2);
         esig=0.0; plastic_sig=0.0; plastic_str=0.0;
         double eig_sig_var[2], eig_sig_vec[4];

         double sig1{0.0};
         double sig3{0.0};
         double  syy{0.0}; // Syy is non-zero value in plane strain condition
         double  syy_old{0.0}; // Syy is non-zero value in plane strain condition
         double msig{0.0}; // mean stress
         double evol{0.0}; // volumetric strain
         double DEG2RAD{M_PI/180.0};
         double depls{0.0}; // 2nd invariant of plastic strain

         double fs{0.0};
         double ft{0.0};
         double fh{0.0};
         double N_phi{0.0};
         double st_N_phi{0.0};
         double N_psi{0.0};
         double beta{0.0};
         double viscosity{0.0};
         double relax{0.0};
         // double relax_limit{1.0};
         double dt_scaled{0.0};
         double numerator   = {0.0};
         double denominator = {0.0};
         double pls_old = {0.0}; // cumulative 2nd invariant of plastic strain
         double p_slope = {0.0}; 
         double fri_str = {0.0}; // strain_dependent friction angle
         double dil_str = {0.0}; // strain_dependent dilation angle
         double coh_str = {0.0}; // strain_dependent cohesion
         double ten_cut = {0.0};
         int mat{0};
         int nsize{mat_gf.Size()};
         int mat_num{lambda.Size()};
      
         double pls0_c = {0.0}; 
         double pls1_c = {0.0}; 
         double rho_c = {0.0}; 
         double lambda_c = {0.0}; 
         double mu_c = {0.0}; 
         double tension_cutoff_c = {0.0}; 
         double cohesion0_c = {0.0}; 
         double cohesion1_c = {0.0}; 
         double friction_angle0_c ={0.0}; 
         double friction_angle1_c ={0.0}; 
         double dilation_angle0_c ={0.0}; 
         double dilation_angle1_c ={0.0}; 
         double plastic_viscosity_c ={0.0};
         double pwave_speed ={0.0};
         double time_scale ={1.0};
         double rad ={0.0};
         double cos2t ={1.0};
         double sin2t ={0.0};

         // element-wise loop
         esig=0.0; plastic_sig=0.0; plastic_str=0.0;
         esig_old=0.0; esig_inc=0.0;
         pls_old = p_gf.GetValue(T, ip);
         if(pls_old < 0.0){pls_old =0.0;}

         pls0_c =0.0; pls1_c =0.0; rho_c = 0.0; lambda_c = 0.0; mu_c = 0.0; time_scale = 1.0;
         tension_cutoff_c = 0.0; cohesion0_c = 0.0; cohesion1_c = 0.0; friction_angle0_c = 0.0; friction_angle1_c = 0.0;
         dilation_angle0_c = 0.0; dilation_angle1_c = 0.0; plastic_viscosity_c = 0.0;

         for( int i = 0; i < mat_num; i++ )
         {
            pls0_c = pls0_c + comp_gf.GetValue(T, ip, i+1)*pls0[i];
            pls1_c = pls1_c + comp_gf.GetValue(T, ip, i+1)*pls1[i];
            rho_c = rho_c + comp_gf.GetValue(T, ip, i+1)*rho[i];
            lambda_c = lambda_c + comp_gf.GetValue(T, ip, i+1)*lambda[i];
            mu_c = mu_c + comp_gf.GetValue(T, ip, i+1)*mu[i];
            tension_cutoff_c = tension_cutoff_c + comp_gf.GetValue(T, ip, i+1)*tension_cutoff[i];
            cohesion0_c = cohesion0_c + comp_gf.GetValue(T, ip, i+1)*cohesion0[i];
            cohesion1_c = cohesion1_c + comp_gf.GetValue(T, ip, i+1)*cohesion1[i];
            friction_angle0_c = friction_angle0_c + comp_gf.GetValue(T, ip, i+1)*friction_angle0[i];
            friction_angle1_c = friction_angle1_c + comp_gf.GetValue(T, ip, i+1)*friction_angle1[i];
            dilation_angle0_c = dilation_angle0_c + comp_gf.GetValue(T, ip, i+1)*dilation_angle0[i];
            dilation_angle1_c = dilation_angle1_c + comp_gf.GetValue(T, ip, i+1)*dilation_angle1[i];
         }
         

         // linear weakening
         p_slope = (pls_old - pls0_c)/(pls1_c - pls0_c);
         pwave_speed = sqrt((lambda_c + 2*mu_c)/rho_c);
         if(h_min  > 0){time_scale = h_min / pwave_speed;}
         plastic_viscosity_c = time_scale * mu_c;

         // sqrt(pow((xc-location[0]), 2) + pow((yc-location[1]), 2));

         msig = (s_gf.GetValue(T, ip, 1) + s_gf.GetValue(T, ip, 2))*0.5;
         rad  = std::sqrt(0.25*(s_gf.GetValue(T, ip, 1) -  s_gf.GetValue(T, ip, 2))*(s_gf.GetValue(T, ip, 1) -  s_gf.GetValue(T, ip, 2)) + s_gf.GetValue(T, ip, 3)*s_gf.GetValue(T, ip, 3));
         
         Vector sig_var(2);
         Vector sig_dir(2);

         sig_var[0] = msig - rad;
         sig_var[1] = msig + rad;

         {
            // direction cosine and sine of 2*theta
            const double eps = 1e-15;
            double a = 0.5 * (sig_var[0] - sig_var[1]);
            double b = - rad; // always negative
            if (b < -eps) {
                  cos2t = a / b;
                  sin2t = s_gf.GetValue(T, ip, 3) / b;
            }
            else {
                  cos2t = 1.0;
                  sin2t = 0.0;
            }
         }

         sig1 = sig_var[0]; // most compressive pincipal stress
         sig3 = sig_var[1]; // least compressive pincipal stress

         // linear strain weaking on cohesion, friction and dilation angles.
         coh_str = cohesion0_c; fri_str = friction_angle0_c; dil_str = dilation_angle0_c;

         if (pls_old < pls0_c) {
            // no weakening yet
            coh_str = cohesion0_c;
            fri_str = friction_angle0_c;
            dil_str = dilation_angle0_c;
         }
         else if (pls_old < pls1_c) {
            // linear weakening
            coh_str = cohesion0_c + p_slope * (cohesion1_c - cohesion0_c);
            fri_str = friction_angle0_c + p_slope * (friction_angle1_c - friction_angle0_c);
            dil_str = dilation_angle0_c + p_slope * (dilation_angle1_c - dilation_angle0_c);
         }
         else {
            // saturated weakening
            coh_str = cohesion1_c; fri_str = friction_angle1_c; dil_str = dilation_angle1_c;
         }

         // std::cout << coh_str <<"," << p_slope << "," << fri_str << std::endl;
         N_phi = (1+sin(DEG2RAD*fri_str))/(1-sin(DEG2RAD*fri_str));
         st_N_phi = cos(DEG2RAD*fri_str)/(1-sin(DEG2RAD*fri_str));
         N_psi = -1*(1+sin(DEG2RAD*dil_str))/(1-sin(DEG2RAD*dil_str)); // partial_g/partial_sig3

         if(tension_cutoff_c == 0)
         {
            ten_cut = coh_str/tan(DEG2RAD*fri_str);
         }
         else{ten_cut = tension_cutoff_c;}

         // shear failure function
         fs = sig1 - N_phi*sig3 + 2*coh_str*st_N_phi;
         // tension failure function
         ft = sig3 - ten_cut;
         // bisects the obtuse angle made by two yield function
         fh = sig3 - ten_cut + (sqrt(N_phi*N_phi + 1.0)+ N_phi)*(sig1 - N_phi*ten_cut + 2*coh_str*st_N_phi);
         depls = 0.0;

         if((fs < 0 & fh < 0) | (ft > 0 & fh > 0))
         {  
            // std::cout <<"failure occurs" << std::endl;

            if(fs < 0 & fh < 0) // stress correction at shear failure
            {
               // Equations 28 and 30 from Choi et al. (2013; DynEarthSol2D: An efficient unstructured finite element method to study long-term tectonic deformation). 
               beta = fs;
               beta = beta / (((lambda_c+2*mu_c)*1 - N_phi*lambda_c*1) + (lambda_c*N_psi - N_phi*(lambda_c+2*mu_c)*N_psi));
               
               plastic_str(0,0) = (lambda_c + 2*mu_c + lambda_c*N_psi) * beta; 
               // syy -= (lambda_c + lambda_c*N_psi) * beta;
               plastic_str(1,1) = (lambda_c + (lambda_c+2*mu_c)*N_psi) * beta;
               // reduced form of 2nd invariant
               depls = std::fabs(beta) * std::sqrt((3 - 2*N_psi + 3*N_psi*N_psi) / 8); 
               
            }
            else if (ft > 0 & fh > 0) // stress correction at tension failure
            {
               beta = ft;
               beta = beta/(lambda_c+2*mu_c);

               plastic_str(0,0) = lambda_c * beta * 1;
               // syy -= lambda_c * beta * 1;
               plastic_str(1,1) = (lambda_c+2*mu_c) * beta * 1;

               // reduced form of 2nd invariant
               depls = std::fabs(beta) * std::sqrt(7. / 18);
            }

            // Rotating Principal axis to XYZ axis
            double dc2 = ((sig1-plastic_str(0,0)) - (sig3-plastic_str(1,1))) * cos2t;
            double dss =  (sig1-plastic_str(0,0)) + (sig3-plastic_str(1,1));

            // Updating new stress to grid function
            viscosity = plastic_viscosity_c;
            relax = exp(-dt_old/viscosity);
            numerator  = fabs(1 - relax);
            denominator= fabs(dt_old/viscosity);
            dt_scaled = dt_old/viscosity;

            // Based on L'Hôpital's rule, 0/0 become 1 
            // if((numerator > 1e-15) & (denominator > 1e-15)){relax_limit = (1 - relax)/(dt_old/viscosity);}

         
            if(viscoplastic)
            {
               // Caushy stress at previous time step
               esig_old(0,0) = s_old_gf.GetValue(T, ip, 1); esig_old(0,1) = s_old_gf.GetValue(T, ip, 3);  
               esig_old(1,0) = s_old_gf.GetValue(T, ip, 3); esig_old(1,1) = s_old_gf.GetValue(T, ip, 2); 

               // Elastic stress increment
               esig_inc(0,0) = esig(0,0) - esig_old(0,0); esig_inc(0,1) = esig(0,1) - esig_old(0,1); 
               esig_inc(1,0) = esig(1,0) - esig_old(1,0); esig_inc(1,1) = esig(1,1) - esig_old(1,1); 

               // Implicit backward Euler algorithm
               K(0)=((esig_old(0,0) + esig_inc(0,0)) + dt_scaled*0.5 * (dss + dc2))/(1.0+dt_scaled);
               K(1)=((esig_old(1,1) + esig_inc(1,1)) + dt_scaled*0.5 * (dss - dc2))/(1.0+dt_scaled);
               K(2)=((esig_old(0,1) + esig_inc(0,1)) + dt_scaled*0.5 * ((sig1-plastic_str(0,0)) - (sig3-plastic_str(1,1))) * sin2t)/(1.0+dt_scaled);
               depls = dt_scaled*depls/(1.0+dt_scaled); 
               K(3)=depls;
            }
            else
            {
               K(0)=0.5 * (dss + dc2);
               K(1)=0.5 * (dss - dc2);
               K(2)=0.5 * ((sig1-plastic_str(0,0)) - (sig3-plastic_str(1,1))) * sin2t;
               K(3)=depls;
            }
         }
         else
         {
            K(0)=esig(0,0);
            K(1)=esig(1,1);
            K(2)=esig(0,1);
            K(3)=0.0;
         }

      }
      virtual ~ReturnMapping2D_simple_Coefficient() { }
   };

   // class IcVelCoefficient : public VectorCoefficient
   // {
   // private:
   //    ParGridFunction &xyz, &rho;
   //    int dim;
   //    double gravity, thickness;

   // public:
   //    IcVelCoefficient (int &_dim, ParGridFunction &_xyz, Vector &_location, double &_rad, double &_ini_pls)
   //       : VectorCoefficient(_dim), xyz(_xyz)
   //       {
   //          dim=_dim; location = _location; rad = _rad; ini_pls = _ini_pls;
   //       }
   //    virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
   //    {
   //       K.SetSize(dim);
   //       double r = 0.0;
   //       double xc = xyz.GetValue(T, ip,1);
   //       double yc = xyz.GetValue(T, ip,2);
   //       double zc = 0.0;
   //       if(dim == 3){zc = xyz.GetValue(T, ip,3);}
         
   //    }
   //    virtual ~IcVelCoefficient() { }
   // };

   // class PlasticCoefficient : public VectorCoefficient
   // {
   // private:
   //    ParGridFunction &x, &y, &z;
   //    int dim;
   //    Vector location;
   //    double rad, ini_pls;

   // public:
   //    PlasticCoefficient (int &_dim, ParGridFunction &_x, ParGridFunction &_y, ParGridFunction &_z, Vector &_location, double &_rad, double &_ini_pls)
   //       : VectorCoefficient(_dim), x(_x), y(_y), z(_z)  
   //       {
   //          dim=_dim; location = _location; rad = _rad; ini_pls = _ini_pls;
   //       }
   //    virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
   //    {
   //       K.SetSize(1);
   //       double r = 0.0;
   //       double xc = x.GetValue(T, ip);
   //       double yc = y.GetValue(T, ip);
   //       double zc = z.GetValue(T, ip);

   //       if(dim == 2){ r = sqrt(pow((xc-location[0]), 2) + pow((yc-location[1]), 2));}
   //       else if(dim == 3){r = sqrt(pow((xc-location[0]), 2) + pow((yc-location[1]), 2) + pow((zc-location[2]), 2));}

   //       if(r <= rad)
   //       {
   //          K(0) = ini_pls;
   //       }
   //       else
   //       {
   //          K(0) = 0.0;
   //       }
   //    }
   //    virtual ~PlasticCoefficient() { }
   // };

   // class PlasticCoefficient : public VectorCoefficient
   // {
   // private:
   //    ParGridFunction &x, &y, &z;
   //    int dim;
   //    Vector location;
   //    double rad, ini_pls;

   // public:
   //    PlasticCoefficient (int &_dim, ParGridFunction &_x, ParGridFunction &_y, ParGridFunction &_z, Vector &_location, double &_rad, double &_ini_pls)
   //       : VectorCoefficient(_dim), x(_x), y(_y), z(_z)  
   //       {
   //          dim=_dim; location = _location; rad = _rad; ini_pls = _ini_pls;
   //       }
   //    virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
   //    {
   //       K.SetSize(1);
   //       double r = 0.0;
   //       double xc = x.GetValue(T, ip);
   //       double yc = y.GetValue(T, ip);
   //       double zc = z.GetValue(T, ip);

   //       if(dim == 2){ r = sqrt(pow((xc-location[0]), 2) + pow((yc-location[1]), 2));}
   //       else if(dim == 3){r = sqrt(pow((xc-location[0]), 2) + pow((yc-location[1]), 2) + pow((zc-location[2]), 2));}

   //       if(r <= rad)
   //       {
   //          K(0) = ini_pls;
   //       }
   //       else
   //       {
   //          K(0) = 0.0;
   //       }
   //    }
   //    virtual ~PlasticCoefficient() { }
   // };

   // class LithostaticCoefficient : public VectorCoefficient
   // {
   // private:
   //    ParGridFunction &y, &z, &rho;
   //    int dim;
   //    double gravity, thickness;

   // public:
   //    LithostaticCoefficient (int &_dim, ParGridFunction &_y, ParGridFunction &_z, ParGridFunction &_rho, double &_gravity, double &_thickness)
   //       : VectorCoefficient(_dim), y(_y), z(_z), rho(_rho)  
   //       {
   //          dim=_dim; gravity = _gravity; thickness = _thickness; 
   //       }
   //    virtual void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip)
   //    {
   //       K.SetSize(3*(dim-1));

   //       double yc = y.GetValue(T, ip);
   //       double zc = z.GetValue(T, ip);
   //       double denc = rho.GetValue(T, ip);

   //       if(dim == 2)
   //       {
   //          K(0) = -1.0*fabs(thickness - yc)*denc*gravity;
   //          K(1) = -1.0*fabs(thickness - yc)*denc*gravity;
   //       }
   //       else if(dim ==3)
   //       {
   //          K(0) = -1.0*fabs(thickness - zc)*denc*gravity;
   //          K(1) = -1.0*fabs(thickness - zc)*denc*gravity;
   //          K(2) = -1.0*fabs(thickness - zc)*denc*gravity;
   //       }
   //    }
   //    virtual ~LithostaticCoefficient() { }
   // };
}
