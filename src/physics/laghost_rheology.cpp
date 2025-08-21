#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <cmath>
#include "laghost_rheology.hpp"
namespace mfem
{
   void Returnmapping3d (Vector &comp_gf, Vector &s_gf, Vector &s_old_gf, Vector &p_gf, Vector &mat_gf, int &dim, double &h_min, Vector &rho, Vector &lambda, Vector &mu, Vector &tension_cutoff, Vector &cohesion0, Vector &cohesion1, Vector &pls0, Vector &pls1, Vector &friction_angle0, Vector &friction_angle1, Vector &dilation_angle0, Vector &dilation_angle1, Vector &plastic_viscosity, bool &viscoplastic, double &dt_old) 
   {
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
      // bool viscoplastic = true;
      // bool viscoplastic = false;

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

      for( int i = 0; i < nsize; i++ )
      {  
            esig=0.0; plastic_sig=0.0; plastic_str=0.0;
            esig_old=0.0; esig_inc=0.0;
            double eig_sig_var[3], eig_sig_vec[9];

            // mat = mat_gf[i];
            pls_old = p_gf[i];
            if(pls_old < 0.0){pls_old =0.0; p_gf[i] = 0.0;}

            pls0_c =0.0; pls1_c =0.0; rho_c = 0.0; lambda_c = 0.0; mu_c = 0.0; time_scale = 1.0;
            tension_cutoff_c = 0.0; cohesion0_c = 0.0; cohesion1_c = 0.0; friction_angle0_c = 0.0; friction_angle1_c = 0.0;
            dilation_angle0_c = 0.0; dilation_angle1_c = 0.0; plastic_viscosity_c = 0.0;
            for( int ii = 0; ii < mat_num; ii++ )
            {
               pls0_c = pls0_c + comp_gf[i+nsize*ii]*pls0[ii];
               pls1_c = pls1_c + comp_gf[i+nsize*ii]*pls1[ii];
               rho_c = rho_c + comp_gf[i+nsize*ii]*rho[ii];
               lambda_c = lambda_c + comp_gf[i+nsize*ii]*lambda[ii];
               mu_c = mu_c + comp_gf[i+nsize*ii]*mu[ii];
               tension_cutoff_c = tension_cutoff_c + comp_gf[i+nsize*ii]*tension_cutoff[ii];
               cohesion0_c = cohesion0_c + comp_gf[i+nsize*ii]*cohesion0[ii];
               cohesion1_c = cohesion1_c + comp_gf[i+nsize*ii]*cohesion1[ii];
               friction_angle0_c = friction_angle0_c + comp_gf[i+nsize*ii]*friction_angle0[ii];
               friction_angle1_c = friction_angle1_c + comp_gf[i+nsize*ii]*friction_angle1[ii];
               dilation_angle0_c = dilation_angle0_c + comp_gf[i+nsize*ii]*dilation_angle0[ii];
               dilation_angle1_c = dilation_angle1_c + comp_gf[i+nsize*ii]*dilation_angle1[ii];
               // plastic_viscosity_c = plastic_viscosity_c + comp_gf[i+nsize*ii]*plastic_viscosity[ii];
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
               msig = (s_gf[i+nsize*0] + s_gf[i+nsize*1])*0.5;
               evol = msig/(lambda_c+mu_c);
               syy  = evol * lambda_c;
               esig(0,0) = s_gf[i+nsize*0]; esig(0,1) = 0.0; esig(0,2) = s_gf[i+nsize*2]; 
               esig(1,0) =             0,0; esig(1,1) = syy; esig(1,2) =             0.0;
               esig(2,0) = s_gf[i+nsize*2]; esig(2,1) = 0.0; esig(2,2) = s_gf[i+nsize*1];

               // Caushy stress at previous time step
               msig = (s_old_gf[i+nsize*0] + s_old_gf[i+nsize*1])*0.5;
               evol = msig/(lambda_c+mu_c);
               syy  = evol * lambda_c;
               esig_old(0,0) = s_old_gf[i+nsize*0]; esig_old(0,1) = 0.0; esig_old(0,2) = s_old_gf[i+nsize*2]; 
               esig_old(1,0) =                 0,0; esig_old(1,1) = syy; esig_old(1,2) =                 0.0;
               esig_old(2,0) = s_old_gf[i+nsize*2]; esig_old(2,1) = 0.0; esig_old(2,2) = s_old_gf[i+nsize*1];
            }
            else
            {
               esig(0,0) = s_gf[i+nsize*0]; esig(0,1) = s_gf[i+nsize*3]; esig(0,2) = s_gf[i+nsize*4]; 
               esig(1,0) = s_gf[i+nsize*3]; esig(1,1) = s_gf[i+nsize*1]; esig(1,2) = s_gf[i+nsize*5];
               esig(2,0) = s_gf[i+nsize*4]; esig(2,1) = s_gf[i+nsize*5]; esig(2,2) = s_gf[i+nsize*2];

               // Caushy stress at previous time step
               esig_old(0,0) = s_old_gf[i+nsize*0]; esig_old(0,1) = s_old_gf[i+nsize*3]; esig_old(0,2) = s_old_gf[i+nsize*4]; 
               esig_old(1,0) = s_old_gf[i+nsize*3]; esig_old(1,1) = s_old_gf[i+nsize*1]; esig_old(1,2) = s_old_gf[i+nsize*5];
               esig_old(2,0) = s_old_gf[i+nsize*4]; esig_old(2,1) = s_old_gf[i+nsize*5]; esig_old(2,2) = s_old_gf[i+nsize*2];
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
               // std::cout << "tension cutoff " << ten_cut << std::endl;
            }
            else{ten_cut = tension_cutoff_c;}

            // shear failure function
            fs = sig1 - N_phi*sig3 + 2*coh_str*st_N_phi;
            // tension failure function
            ft = sig3 - ten_cut;
            // bisects the obtuse angle made by two yield function
            fh = sig3 - ten_cut + (sqrt(N_phi*N_phi + 1.0)+ N_phi)*(sig1 - N_phi*ten_cut + 2*coh_str*st_N_phi);

            depls = 0.0;
                    
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
                  // depls = std::fabs(alam) * std::sqrt((3 + 2*anpsi + 3*anpsi*anpsi) / 8);
               }
               else
               {
                  depls = std::fabs(beta) * std::sqrt((7 - 4*N_psi + 7*N_psi*N_psi) / 18);
                  // depls = std::fabs(alam) * std::sqrt((7 + 4*anpsi + 7*anpsi*anpsi) / 18);
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
                  // depls = std::fabs(alam) * std::sqrt(3. / 8);
               }
               else
               {
                  depls = std::fabs(beta) * std::sqrt(7. / 18);
                  // depls = std::fabs(alam) * std::sqrt(7. / 18);
               }

               // std::cout << i << ", tensile failure occures,  ft = " << ft << ", plastic_str(0,0) = " << plastic_str(0,0) << ", plastic_str(1,1) = " << \
               // plastic_str(1,1) << "plastic_str(2,2) = " << plastic_str(2,2) << ", depls = " << depls << std::endl;
                  
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
            // viscosity = 10000.0; 
            // viscosity = 0.0;  // rate_independent
            // viscosity = 1.0e+38; //  elastic
            relax = exp(-dt_old/viscosity);
            numerator  = fabs(1 - relax);
            denominator= fabs(dt_old/viscosity);
            dt_scaled = dt_old/viscosity;

            // Based on L'Hôpital's rule, 0/0 become 1 
            // if((numerator > 1e-15) & (denominator > 1e-15)){relax_limit = (1 - relax)/(dt_old/viscosity);}

            // std::cout 
            // if(viscosity > 1e+100){viscoplastic = false;}

            if((fs < 0 & fh < 0) | (ft > 0 & fh > 0))
            {
               // if(i ==0){std::cout << fs <<","<< beta*1e6 << ","<< dt_old << "," << relax <<","<< relax_limit << std::endl;}
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
                     s_gf[i+nsize*0]=((esig_old(0,0) + esig_inc(0,0)) + dt_scaled*plastic_sig(0,0))/(1.0+dt_scaled);
                     s_gf[i+nsize*2]=((esig_old(0,2) + esig_inc(0,2)) + dt_scaled*plastic_sig(0,2))/(1.0+dt_scaled);
                     s_gf[i+nsize*2]=((esig_old(2,0) + esig_inc(2,0)) + dt_scaled*plastic_sig(2,0))/(1.0+dt_scaled);
                     s_gf[i+nsize*1]=((esig_old(2,2) + esig_inc(2,2)) + dt_scaled*plastic_sig(2,2))/(1.0+dt_scaled);
                  }
                  else
                  {
                     s_gf[i+nsize*0]=plastic_sig(0,0); s_gf[i+nsize*2]=plastic_sig(0,2);
                     s_gf[i+nsize*2]=plastic_sig(2,0); s_gf[i+nsize*1]=plastic_sig(2,2);
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
                     s_gf[i+nsize*0]=((esig_old(0,0) + esig_inc(0,0)) + dt_scaled*plastic_sig(0,0))/(1.0+dt_scaled);
                     s_gf[i+nsize*3]=((esig_old(0,1) + esig_inc(0,1)) + dt_scaled*plastic_sig(0,1))/(1.0+dt_scaled);
                     s_gf[i+nsize*4]=((esig_old(0,2) + esig_inc(0,2)) + dt_scaled*plastic_sig(0,2))/(1.0+dt_scaled);
                     s_gf[i+nsize*3]=((esig_old(1,0) + esig_inc(1,0)) + dt_scaled*plastic_sig(1,0))/(1.0+dt_scaled);
                     s_gf[i+nsize*1]=((esig_old(1,1) + esig_inc(1,1)) + dt_scaled*plastic_sig(1,1))/(1.0+dt_scaled);
                     s_gf[i+nsize*5]=((esig_old(1,2) + esig_inc(1,2)) + dt_scaled*plastic_sig(1,2))/(1.0+dt_scaled);
                     s_gf[i+nsize*4]=((esig_old(2,0) + esig_inc(2,0)) + dt_scaled*plastic_sig(2,0))/(1.0+dt_scaled);
                     s_gf[i+nsize*5]=((esig_old(2,1) + esig_inc(2,1)) + dt_scaled*plastic_sig(2,1))/(1.0+dt_scaled);
                     s_gf[i+nsize*2]=((esig_old(2,2) + esig_inc(2,2)) + dt_scaled*plastic_sig(2,2))/(1.0+dt_scaled);

                  }
                  else
                  {
                     s_gf[i+nsize*0]=plastic_sig(0,0); s_gf[i+nsize*3]=plastic_sig(0,1); s_gf[i+nsize*4]=plastic_sig(0,2); 
                     s_gf[i+nsize*3]=plastic_sig(1,0); s_gf[i+nsize*1]=plastic_sig(1,1); s_gf[i+nsize*5]=plastic_sig(1,2);
                     s_gf[i+nsize*4]=plastic_sig(2,0); s_gf[i+nsize*5]=plastic_sig(2,1); s_gf[i+nsize*2]=plastic_sig(2,2);
                  }
               }
            }
            
            // Adding 2nd invariant of plastic strain increment

            if(viscoplastic)
            {
               // depls = (1 - relax)*depls; 
               depls = dt_scaled*depls/(1.0+dt_scaled); 
               p_gf[i] += depls;
            }
            else
            {
               // p_gf[i] += depls;
               p_gf[i] += std::fabs(depls);
            }
            
      }
   }

   void Returnmapping2d (Vector &comp_gf, Vector &s_gf, Vector &s_old_gf, Vector &p_gf, Vector &mat_gf, int &dim, double &h_min, Vector &rho, Vector &lambda, Vector &mu, Vector &tension_cutoff, Vector &cohesion0, Vector &cohesion1, Vector &pls0, Vector &pls1, Vector &friction_angle0, Vector &friction_angle1, Vector &dilation_angle0, Vector &dilation_angle1, Vector &plastic_viscosity, bool &viscoplastic, double &dt_old) 
   {
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
      double viscosity{1.0};
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
      // bool viscoplastic = true;
      // bool viscoplastic = false;

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

      for( int i = 0; i < nsize; i++ )
      {  
            esig=0.0; plastic_sig=0.0; plastic_str=0.0;
            esig_old=0.0; esig_inc=0.0;
            // double eig_sig_var[2], eig_sig_vec[4];
            Vector eig_sig_var(2), eig_sig_vec(4);

            // mat = mat_gf[i];
            pls_old = p_gf[i];
            if(pls_old < 0.0){pls_old =0.0; p_gf[i] = 0.0;}

            pls0_c =0.0; pls1_c =0.0; rho_c = 0.0; lambda_c = 0.0; mu_c = 0.0; time_scale = 1.0;
            tension_cutoff_c = 0.0; cohesion0_c = 0.0; cohesion1_c = 0.0; friction_angle0_c = 0.0; friction_angle1_c = 0.0;
            dilation_angle0_c = 0.0; dilation_angle1_c = 0.0; plastic_viscosity_c = 0.0;
            for( int ii = 0; ii < mat_num; ii++ )
            {
               pls0_c = pls0_c + comp_gf[i+nsize*ii]*pls0[ii];
               pls1_c = pls1_c + comp_gf[i+nsize*ii]*pls1[ii];
               rho_c = rho_c + comp_gf[i+nsize*ii]*rho[ii];
               lambda_c = lambda_c + comp_gf[i+nsize*ii]*lambda[ii];
               mu_c = mu_c + comp_gf[i+nsize*ii]*mu[ii];
               tension_cutoff_c = tension_cutoff_c + comp_gf[i+nsize*ii]*tension_cutoff[ii];
               cohesion0_c = cohesion0_c + comp_gf[i+nsize*ii]*cohesion0[ii];
               cohesion1_c = cohesion1_c + comp_gf[i+nsize*ii]*cohesion1[ii];
               friction_angle0_c = friction_angle0_c + comp_gf[i+nsize*ii]*friction_angle0[ii];
               friction_angle1_c = friction_angle1_c + comp_gf[i+nsize*ii]*friction_angle1[ii];
               dilation_angle0_c = dilation_angle0_c + comp_gf[i+nsize*ii]*dilation_angle0[ii];
               dilation_angle1_c = dilation_angle1_c + comp_gf[i+nsize*ii]*dilation_angle1[ii];
               // plastic_viscosity_c = plastic_viscosity_c + comp_gf[i+nsize*ii]*plastic_viscosity[ii];
            }
            // linear weakening
            p_slope = (pls_old - pls0_c)/(pls1_c - pls0_c);
            pwave_speed = sqrt((lambda_c + 2*mu_c)/rho_c);
            if(h_min  > 0){time_scale = h_min / pwave_speed;}
            plastic_viscosity_c = time_scale * mu_c;

            // std::cout << time_scale <<',' << plastic_viscosity_c << std::endl;
            // 2D in plane strain condition
            // sxx, syy, szz, sxz are non zeros.
            // sxy, syz are zero.
            msig = (s_gf[i+nsize*0] + s_gf[i+nsize*1])*0.5;
            evol = msig/(lambda_c+mu_c);
            syy  = evol * lambda_c;
            esig(0,0) = s_gf[i+nsize*0]; esig(0,1) = s_gf[i+nsize*2];  
            esig(1,0) = s_gf[i+nsize*2]; esig(1,1) = s_gf[i+nsize*1]; 

            // Caushy stress at previous time step
            msig = (s_old_gf[i+nsize*0] + s_old_gf[i+nsize*1])*0.5;
            evol = msig/(lambda_c+mu_c);
            syy_old  = evol * lambda_c;
            esig_old(0,0) = s_old_gf[i+nsize*0]; esig_old(0,1) = s_old_gf[i+nsize*2]; 
            esig_old(1,0) = s_old_gf[i+nsize*2]; esig_old(1,1) = s_old_gf[i+nsize*1]; 

            // Elastic stress increment
            esig_inc(0,0) = esig(0,0) - esig_old(0,0); esig_inc(0,1) = esig(0,1) - esig_old(0,1); 
            esig_inc(1,0) = esig(1,0) - esig_old(1,0); esig_inc(1,1) = esig(1,1) - esig_old(1,1); 
            // esig.CalcEigenvalues(eig_sig_var, eig_sig_vec); 
            esig.CalcEigenvalues(eig_sig_var.GetData(), eig_sig_vec.GetData()); 

            // Vector sig_var(eig_sig_var, 2);
            // Vector eig_sig_dir(eig_sig_vec, 4);

            auto max_it = std::max_element(eig_sig_var.begin(), eig_sig_var.end()); // find iterator to max element
            auto min_it = std::min_element(eig_sig_var.begin(), eig_sig_var.end()); // find iterator to min element
            
            int max_index = std::distance(eig_sig_var.begin(), max_it); // calculate index of max element
            int min_index = std::distance(eig_sig_var.begin(), min_it); // calculate index of min element
            
            // int itm_index = 0; // calculate index of intermediate element
            // if (max_index + min_index == 1) {itm_index = 2;}
            // else if(max_index + min_index == 2) {itm_index = 1;}
            // else {itm_index = 0;}

            sig1 = eig_sig_var[min_index]; // most compressive pincipal stress
            sig3 = eig_sig_var[max_index]; // least compressive pincipal stress

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
            plastic_sig(0,0) = ((eig_sig_var[min_index]-plastic_str(0,0))*eig_sig_vec[0+min_index*2]*eig_sig_vec[0+min_index*2]  + (eig_sig_var[max_index]-plastic_str(1,1))*eig_sig_vec[0+max_index*2]*eig_sig_vec[0+max_index*2]);
            plastic_sig(0,1) = ((eig_sig_var[min_index]-plastic_str(0,0))*eig_sig_vec[0+min_index*2]*eig_sig_vec[1+min_index*2]  + (eig_sig_var[max_index]-plastic_str(1,1))*eig_sig_vec[0+max_index*2]*eig_sig_vec[1+max_index*2]);
            plastic_sig(1,0) = ((eig_sig_var[min_index]-plastic_str(0,0))*eig_sig_vec[1+min_index*2]*eig_sig_vec[0+min_index*2]  + (eig_sig_var[max_index]-plastic_str(1,1))*eig_sig_vec[1+max_index*2]*eig_sig_vec[0+max_index*2]);
            plastic_sig(1,1) = ((eig_sig_var[min_index]-plastic_str(0,0))*eig_sig_vec[1+min_index*2]*eig_sig_vec[1+min_index*2]  + (eig_sig_var[max_index]-plastic_str(1,1))*eig_sig_vec[1+max_index*2]*eig_sig_vec[1+max_index*2]);

            // Updating new stress to grid function
            viscosity = plastic_viscosity_c;
            relax = exp(-dt_old/viscosity);
            numerator  = fabs(1 - relax);
            denominator= fabs(dt_old/viscosity);
            dt_scaled = dt_old/viscosity;

            // Based on L'Hôpital's rule, 0/0 become 1 
            // if((numerator > 1e-15) & (denominator > 1e-15)){relax_limit = (1 - relax)/(dt_old/viscosity);}
            // if(viscosity > 1e+100){viscoplastic = false;}

            if((fs < 0 & fh < 0) | (ft > 0 & fh > 0))
            {
               if(viscoplastic)
               {

                  // Implicit backward Euler algorithm
                  s_gf[i+nsize*0]=((esig_old(0,0) + esig_inc(0,0)) + dt_scaled*plastic_sig(0,0))/(1.0+dt_scaled);
                  s_gf[i+nsize*2]=((esig_old(0,1) + esig_inc(0,1)) + dt_scaled*plastic_sig(0,1))/(1.0+dt_scaled);
                  s_gf[i+nsize*2]=((esig_old(1,0) + esig_inc(1,0)) + dt_scaled*plastic_sig(1,0))/(1.0+dt_scaled);
                  s_gf[i+nsize*1]=((esig_old(1,1) + esig_inc(1,1)) + dt_scaled*plastic_sig(1,1))/(1.0+dt_scaled);
               }
               else
               {
                  s_gf[i+nsize*0]=plastic_sig(0,0); s_gf[i+nsize*2]=plastic_sig(0,1);
                  s_gf[i+nsize*2]=plastic_sig(1,0); s_gf[i+nsize*1]=plastic_sig(1,1);
               }
            }
            
            // Adding 2nd invariant of plastic strain increment

            if(viscoplastic)
            {
               // depls = (1 - relax)*depls; 
               depls = dt_scaled*depls/(1.0+dt_scaled); 
               p_gf[i] += depls;
            }
            else
            {
               p_gf[i] += std::fabs(depls);
            }
            
      }
   }  
}
