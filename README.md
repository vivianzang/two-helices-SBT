# two-helices-SBT

This repository contains MATLAB codes implementing **Slender-Body Theory (SBT)** for low-Reynolds-number flow around one or two helices. 

## Adapted From
The underlying approach is adapted from:
> Bruce (2025). *Helical Swimming Simulator* (https://www.mathworks.com/matlabcentral/fileexchange/39265-helical-swimming-simulator), MATLAB Central File Exchange. Retrieved February 25, 2025.

---

## 1. Overview

- **Core Idea**: We discretize a slender filament (helix) into \(N\) segments and solve for the force density via a linear system \(M \mathbf{f} = \mathbf{U}\), where \(M\) includes both the local diagonal approximation and off-diagonal stokeslets.
- **Functions**:
  - `solveOneHelixKnSBT.m` – Discretizes a **single** helix and solves for the force density.
  - `solveTwoHelixKnSBT.m` – Discretizes **two** helices (offset or phase-shifted).
  - `calculateVelocity.m` – Computes the velocity at user-specified spatial points due to the solved force distribution.

---

## 2. Files

1. **`main_double_helix.m`** (example)  
   Shows how to set parameters, call the double-helix solver, and do a parametric sweep over separation distance or phase shift.

2. **`timeAvgAxialVelDouble.m`**  
   Time-averages the axial velocity over one rotation for a **two-helix** configuration.

3. **`timeAvgAxialVzSingle.m`**  
   Time-averages the axial velocity for a **single** helix.

4. **`solveOneHelixKnSBT.m`**  
   Main solver for a single helix in dimensionless form, solves for the force density, computes total thrust and torque and velocity fields.

5. **`solveTwoHelixKnSBT.m`**  
   Similar to the single-helix solver but sets up two helices (with separation distance, phase shift, etc.) and returns combined forces/torques.

6. **`calculateVelocity.m`**  
   Given force distributions and positions, computes the velocity at query points.

---

## 3. Basic Usage

1. **Clone or Download** the repository.
2. **Open MATLAB** in the downloaded folder.
3. Run one of the main scripts (e.g. `main_double_helix.m`):
   ```matlab
   % Example
   >> main_double_helix
