# two-helices-SBT

This repository contains MATLAB codes implementing **Slender-Body Theory (SBT)** for low-Reynolds-number flow around one or two helices. 

---

## 1. Functions
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

## 3. Acknowledgments

This implementation follows the methods discussed in:

- **Kim & Powers (2004)**: Hydrodynamic interactions between rotating helices. *Phys. Rev. E* **69**, 061910.  
- **Lighthill (1976)**: Flagellar hydrodynamics. *SIAM Rev.* **18**, 161--230.  
- **Rodenborn et al. (2013)**: Propulsion of microorganisms by a helical flagellum. *Proc. Natl. Acad. Sci. USA* **110**, 338--347.

The code is also partially adapted from:

> **Bruce (2025).** *Helical Swimming Simulator* (https://www.mathworks.com/matlabcentral/fileexchange/39265-helical-swimming-simulator),  
> MATLAB Central File Exchange. Retrieved February 25, 2025.
