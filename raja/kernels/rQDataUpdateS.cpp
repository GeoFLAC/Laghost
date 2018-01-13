// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.
#include "raja.hpp"

// *****************************************************************************
template<const int NUM_DIM,
         const int NUM_QUAD,
         const int NUM_QUAD_1D,
         const int NUM_DOFS_1D>
__global__
void rUpdateQuadratureData2DS(const double GAMMA,
                              const double H0,
                              const double CFL,
                              const bool USE_VISCOSITY,
                              const int numElements,
                              const double* restrict dofToQuad,
                              const double* restrict dofToQuadD,
                              const double* restrict quadWeights,
                              const double* restrict v,
                              const double* restrict e,
                              const double* restrict rho0DetJ0w,
                              const double* restrict invJ0,
                              const double* restrict J,
                              const double* restrict invJ,
                              const double* restrict detJ,
                              double* restrict stressJinvT,
                              double* restrict dtEst) {
  const int NUM_QUAD_2D = NUM_QUAD_1D*NUM_QUAD_1D;
  const int NUM_QUAD_DOFS_1D = (NUM_QUAD_1D * NUM_DOFS_1D);
  const int NUM_MAX_1D = (NUM_QUAD_1D<NUM_DOFS_1D)?NUM_DOFS_1D:NUM_QUAD_1D;

  for (int el = 0; el < numElements; ++el) {
    __shared__ double s_dofToQuad[NUM_QUAD_DOFS_1D];//@dim(NUM_QUAD_1D, NUM_DOFS_1D);
    __shared__ double s_dofToQuadD[NUM_QUAD_DOFS_1D];//@dim(NUM_QUAD_1D, NUM_DOFS_1D);

    __shared__ double s_xy[NUM_DIM * NUM_QUAD_DOFS_1D];//@dim(NUM_DIM, NUM_DOFS_1D, NUM_QUAD_1D);
    __shared__ double s_xDy[NUM_DIM * NUM_QUAD_DOFS_1D];//@dim(NUM_DIM, NUM_DOFS_1D, NUM_QUAD_1D);

    __shared__ double s_gradv[NUM_DIM * NUM_DIM * NUM_QUAD_2D];//@dim(NUM_DIM, NUM_DIM, NUM_QUAD_2D);

    double r_v[NUM_DIM * NUM_DOFS_1D];//@dim(NUM_DIM, NUM_DOFS_1D);

    for (int x = 0; x < NUM_MAX_1D; ++x) {
      for (int id = x; id < NUM_QUAD_DOFS_1D; id += NUM_MAX_1D) {
        s_dofToQuad[id]  = dofToQuad[id];
        s_dofToQuadD[id] = dofToQuadD[id];
      }
    }

    for (int dx = 0; dx < NUM_MAX_1D; ++dx) {
      if (dx < NUM_DOFS_1D) {
        for (int qy = 0; qy < NUM_QUAD_1D; ++qy) {
          for (int vi = 0; vi < NUM_DIM; ++vi) {
            s_xy[ijkNM(vi, dx, qy,NUM_DIM,NUM_DOFS_1D)] = 0;
            s_xDy[ijkNM(vi, dx, qy,NUM_DIM,NUM_DOFS_1D)] = 0;
          }
        }
        for (int dy = 0; dy < NUM_DOFS_1D; ++dy) {
          for (int vi = 0; vi < NUM_DIM; ++vi) {
            r_v[ijN(vi, dy,NUM_DIM)] = v[_ijklNM(vi,dx,dy,el,NUM_DOFS_1D,numElements)];
          }
        }
        for (int qy = 0; qy < NUM_QUAD_1D; ++qy) {
          double xy[NUM_DIM];
          double xDy[NUM_DIM];
          for (int vi = 0; vi < NUM_DIM; ++vi) {
            xy[vi]  = 0;
            xDy[vi] = 0;
          }
          for (int dy = 0; dy < NUM_DOFS_1D; ++dy) {
            for (int vi = 0; vi < NUM_DIM; ++vi) {
              xy[vi]  += r_v[ijN(vi, dy,NUM_DIM)] * s_dofToQuad[ijN(qy,dy,NUM_QUAD_1D)];
              xDy[vi] += r_v[ijN(vi, dy,NUM_DIM)] * s_dofToQuadD[ijN(qy,dy,NUM_QUAD_1D)];
            }
          }
          for (int vi = 0; vi < NUM_DIM; ++vi) {
            s_xy[ijkNM(vi, dx, qy,NUM_DIM,NUM_DOFS_1D)]  = xy[vi];
            s_xDy[ijkNM(vi, dx, qy,NUM_DIM,NUM_DOFS_1D)] = xDy[vi];
          }
        }
      }
    }

    for (int qy = 0; qy < NUM_MAX_1D; ++qy) {
      if (qy < NUM_QUAD_1D) {
        for (int qx = 0; qx < NUM_MAX_1D; ++qx) {
          double gradX[NUM_DIM];
          double gradY[NUM_DIM];
          for (int vi = 0; vi < NUM_DIM; ++vi) {
            gradX[vi] = 0;
            gradY[vi] = 0;
          }
          for (int dx = 0; dx < NUM_DOFS_1D; ++dx) {
            for (int vi = 0; vi < NUM_DIM; ++vi) {
              gradX[vi] += s_xy[ijkNM(vi, dx, qy,NUM_DIM,NUM_DOFS_1D)]  * s_dofToQuadD[ijN(qx, dx,NUM_QUAD_1D)];
              gradY[vi] += s_xDy[ijkNM(vi, dx, qy,NUM_DIM,NUM_DOFS_1D)] * s_dofToQuad[ijN(qx, dx,NUM_QUAD_1D)];
            }
          }
          for (int vi = 0; vi < NUM_DIM; ++vi) {
            s_gradv[ijkN(vi, 0, qx + qy*NUM_QUAD_1D,NUM_DIM)] = gradX[vi];
            s_gradv[ijkN(vi, 1, qx + qy*NUM_QUAD_1D,NUM_DIM)] = gradY[vi];
          }
        }
      }
    }

    for (int qBlock = 0; qBlock < NUM_MAX_1D; ++qBlock) {
      for (int q = qBlock; q < NUM_QUAD; q += NUM_MAX_1D) {
        double q_gradv[NUM_DIM * NUM_DIM];//@dim(NUM_DIM, NUM_DIM);
        double q_stress[NUM_DIM * NUM_DIM];//@dim(NUM_DIM, NUM_DIM);

        const double invJ_00 = invJ[ijklNM(0,0,q,el,NUM_DIM,NUM_QUAD)];
        const double invJ_10 = invJ[ijklNM(1,0,q,el,NUM_DIM,NUM_QUAD)];
        const double invJ_01 = invJ[ijklNM(0,1,q,el,NUM_DIM,NUM_QUAD)];
        const double invJ_11 = invJ[ijklNM(1,1,q,el,NUM_DIM,NUM_QUAD)];

        q_gradv[ijN(0,0,2)] = ((s_gradv[ijkN(0,0,q,2)]*invJ_00) + (s_gradv[ijkN(1,0,q,2)]*invJ_01));
        q_gradv[ijN(1,0,2)] = ((s_gradv[ijkN(0,0,q,2)]*invJ_10) + (s_gradv[ijkN(1,0,q,2)]*invJ_11));
        q_gradv[ijN(0,1,2)] = ((s_gradv[ijkN(0,1,q,2)]*invJ_00) + (s_gradv[ijkN(1,1,q,2)]*invJ_01));
        q_gradv[ijN(1,1,2)] = ((s_gradv[ijkN(0,1,q,2)]*invJ_10) + (s_gradv[ijkN(1,1,q,2)]*invJ_11));

        const double q_Jw = detJ[ijN(q,el,NUM_QUAD)]*quadWeights[q];

        const double q_rho = rho0DetJ0w[ijN(q,el,NUM_QUAD)]/q_Jw;
        const double q_e   = fmax(0.0,e[ijN(q,el,NUM_QUAD)]);

        // TODO: Input OccaVector eos(q,e) -> (stress, soundSpeed)
        const double s = -(GAMMA - 1.0) * q_rho * q_e;
        q_stress[ijN(0,0,2)] = s; q_stress[ijN(1,0,2)] = 0;
        q_stress[ijN(0,1,2)] = 0; q_stress[ijN(1,1,2)] = s;

        const double gradv00 = q_gradv[ijN(0,0,2)];
        const double gradv11 = q_gradv[ijN(1,1,2)];
        const double gradv10 = 0.5 * (q_gradv[ijN(1,0,2)] + q_gradv[ijN(0,1,2)]);
        q_gradv[ijN(1,0,2)] = gradv10;
        q_gradv[ijN(0,1,2)] = gradv10;

        double comprDirX = 1;
        double comprDirY = 0;
        double minEig = 0;
        // linalg/densemat.cpp: Eigensystem2S()
        if (gradv10 == 0) {
          minEig = (gradv00 < gradv11) ? gradv00 : gradv11;
        } else {
          const double zeta  = (gradv11 - gradv00) / (2.0 * gradv10);
          const double azeta = fabs(zeta);
          double t = 1.0 / (azeta + sqrt(1.0 + zeta*zeta));
          if ((t < 0) != (zeta < 0)) {
            t = -t;
          }

          const double c = sqrt(1.0 / (1.0 + t*t));
          const double s = c * t;
          t *= gradv10;

          if ((gradv00 - t) <= (gradv11 + t)) {
            minEig = gradv00 - t;
            comprDirX = c;
            comprDirY = -s;
          } else {
            minEig = gradv11 + t;
            comprDirX = s;
            comprDirY = c;
          }
        }

        // Computes the initial->physical transformation Jacobian.
        const double J_00 = J[ijklNM(0,0,q,el,NUM_DIM,NUM_QUAD)];
        const double J_10 = J[ijklNM(1,0,q,el,NUM_DIM,NUM_QUAD)];
        const double J_01 = J[ijklNM(0,1,q,el,NUM_DIM,NUM_QUAD)];
        const double J_11 = J[ijklNM(1,1,q,el,NUM_DIM,NUM_QUAD)];

        const double invJ0_00 = invJ0[ijklNM(0,0,q,el,NUM_DIM,NUM_QUAD)];
        const double invJ0_10 = invJ0[ijklNM(1,0,q,el,NUM_DIM,NUM_QUAD)];
        const double invJ0_01 = invJ0[ijklNM(0,1,q,el,NUM_DIM,NUM_QUAD)];
        const double invJ0_11 = invJ0[ijklNM(1,1,q,el,NUM_DIM,NUM_QUAD)];

        const double Jpi_00 = ((J_00 * invJ0_00) + (J_10 * invJ0_01));
        const double Jpi_10 = ((J_00 * invJ0_10) + (J_10 * invJ0_11));
        const double Jpi_01 = ((J_01 * invJ0_00) + (J_11 * invJ0_01));
        const double Jpi_11 = ((J_01 * invJ0_10) + (J_11 * invJ0_11));

        const double physDirX = (Jpi_00 * comprDirX) + (Jpi_10 * comprDirY);
        const double physDirY = (Jpi_01 * comprDirX) + (Jpi_11 * comprDirY);

        const double q_h = H0 * sqrt((physDirX * physDirX) + (physDirY * physDirY));

        // TODO: soundSpeed will be an input as well (function call or values per q)
        const double soundSpeed = sqrt(GAMMA * (GAMMA - 1.0) * q_e);
        dtEst[ijN(q, el,NUM_QUAD)] = CFL * q_h / soundSpeed;

        if (USE_VISCOSITY) {
          // TODO: Check how we can extract outside of kernel
          const double mu = minEig;
          double coeff = 2.0 * q_rho * q_h * q_h * fabs(mu);
          if (mu < 0) {
            coeff += 0.5 * q_rho * q_h * soundSpeed;
          }
          for (int y = 0; y < NUM_DIM; ++y) {
            for (int x = 0; x < NUM_DIM; ++x) {
              q_stress[ijN(x,y,2)] += coeff * q_gradv[ijN(x,y,2)];
            }
          }
        }
        const double S00 = q_stress[ijN(0,0,2)]; const double S10 = q_stress[ijN(1,0,2)];
        const double S01 = q_stress[ijN(0,1,2)]; const double S11 = q_stress[ijN(1,1,2)];

        stressJinvT[ijklNM(0,0,q,el,NUM_DIM,NUM_QUAD)] = q_Jw * ((S00 * invJ_00) + (S10 * invJ_01));
        stressJinvT[ijklNM(1,0,q,el,NUM_DIM,NUM_QUAD)] = q_Jw * ((S00 * invJ_10) + (S10 * invJ_11));

        stressJinvT[ijklNM(0,1,q,el,NUM_DIM,NUM_QUAD)] = q_Jw * ((S01 * invJ_00) + (S11 * invJ_01));
        stressJinvT[ijklNM(1,1,q,el,NUM_DIM,NUM_QUAD)] = q_Jw * ((S01 * invJ_10) + (S11 * invJ_11));
      }
    }
  }
}
// *****************************************************************************
typedef void (*fUpdateQuadratureDataS)(const double GAMMA,
                                       const double H0,
                                       const double CFL,
                                       const bool USE_VISCOSITY,
                                       const int numElements,
                                       const double* restrict dofToQuad,
                                       const double* restrict dofToQuadD,
                                       const double* restrict quadWeights,
                                       const double* restrict v,
                                       const double* restrict e,
                                       const double* restrict rho0DetJ0w,
                                       const double* restrict invJ0,
                                       const double* restrict J,
                                       const double* restrict invJ,
                                       const double* restrict detJ,
                                       double* restrict stressJinvT,
                                       double* restrict dtEst);

// *****************************************************************************
void rUpdateQuadratureDataS(const double GAMMA,
                            const double H0,
                            const double CFL,
                            const bool USE_VISCOSITY,
                            const int NUM_DIM,
                            const int NUM_QUAD,
                            const int NUM_QUAD_1D,
                            const int NUM_DOFS_1D,
                            const int nzones,
                            const double* restrict dofToQuad,
                            const double* restrict dofToQuadD,
                            const double* restrict quadWeights,
                            const double* restrict v,
                            const double* restrict e,
                            const double* restrict rho0DetJ0w,
                            const double* restrict invJ0,
                            const double* restrict J,
                            const double* restrict invJ,
                            const double* restrict detJ,
                            double* restrict stressJinvT,
                            double* restrict dtEst){
  const unsigned int id =
    (NUM_DIM<<24)|
    (NUM_QUAD<<16)|
    (NUM_QUAD_1D<<8)|
    (NUM_DOFS_1D);
  assert(LOG2(NUM_DIM)<=8);//printf("NUM_DIM:%d ",(NUM_DIM));
  assert(LOG2(NUM_QUAD)<=8);//printf("NUM_QUAD:%d ",(NUM_QUAD));
  assert(LOG2(NUM_QUAD_1D)<=8);//printf("NUM_QUAD_1D:%d ",(NUM_QUAD_1D));
  assert(LOG2(NUM_DOFS_1D)<=8);//printf("NUM_DOFS_1D:%d ",(NUM_DOFS_1D));
  static std::unordered_map<unsigned int, fUpdateQuadratureDataS> call = {
    // 2D
    {0x2040202,&rUpdateQuadratureData2DS<2,4,2,2>},
    {0x2100403,&rUpdateQuadratureData2DS<2,16,4,3>},
    {0x2100404,&rUpdateQuadratureData2DS<2,16,4,4>},
    
    {0x2190503,&rUpdateQuadratureData2DS<2,25,5,3>},
    {0x2190504,&rUpdateQuadratureData2DS<2,25,5,4>},
    
    {0x2240603,&rUpdateQuadratureData2DS<2,36,6,3>},
    {0x2240604,&rUpdateQuadratureData2DS<2,36,6,4>},
    
    {0x2310704,&rUpdateQuadratureData2DS<2,49,7,4>},
    
    {0x2400805,&rUpdateQuadratureData2DS<2,64,8,5>},
    
    {0x2510905,&rUpdateQuadratureData2DS<2,81,9,5>},
    
    {0x2640A06,&rUpdateQuadratureData2DS<2,100,10,6>},
    {0x2900C07,&rUpdateQuadratureData2DS<2,144,12,7>},
    
    // 3D
    //{0x3100403,&rUpdateQuadratureData3D<3,16,4,3>},
    //{0x3400403,&rUpdateQuadratureData3D<3,64,4,3>},

  };
  if (!call[id]){
    printf("\n[rUpdateQuadratureDataS] id \033[33m0x%X\033[m ",id);
    fflush(stdout);
  }
  assert(call[id]);
  call[id](GAMMA,H0,CFL,USE_VISCOSITY,
           nzones,
           dofToQuad,
           dofToQuadD,
           quadWeights,
           v,
           e,
           rho0DetJ0w,
           invJ0,
           J,
           invJ,
           detJ,
           stressJinvT,
           dtEst);
}
