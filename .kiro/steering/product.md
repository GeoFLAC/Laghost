# Product Overview

**Laghost** (LAGrangian High-Order Solver for Tectonics) is a high-performance computational geodynamics solver that simulates the time-dependent momentum balance of geological media using Lagrangian finite element methods.

## Purpose
- Solves dynamic momentum balance equations for geological materials in a moving Lagrangian frame
- Targets long-term brittle and ductile deformation of rocks coupled with thermal evolution
- Extends capabilities of the MFEM Laghos miniapp for geodynamics applications
- Supports quasi-static solutions with dynamic relaxation and large time steps via mass scaling

## Key Features
- High-order finite element spatial discretization using MFEM
- Explicit high-order time-stepping (Runge-Kutta methods)
- Rock rheologies: elastic, Mohr-Coulomb plasticity, viscoplasticity
- Multi-material tracking and remeshing capabilities
- Parallel MPI execution with optional GPU support
- Adaptive mesh refinement and mesh optimization (TMOP)
- Mass scaling for year-length time steps in geological simulations

## Target Applications
- Long-term geological deformation modeling
- Tectonic processes simulation
- Rock mechanics and geomechanics
- Coupled thermal-mechanical problems in geosciences