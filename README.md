# Vectorized Computational Framework for Floating Flexible Bodies Using Component Mode Synthesis

This repository provides the MATLAB implementation and numerical simulation code for the paper "Vectorized Computational Framework for Floating Flexible Bodies Using Component Mode Synthesis" by Younghwa Park and Christoffer Sloth. The framework was developed at The Maersk Mc-Kinney Moller Institute, University of Southern Denmark, and is designed to address the dynamic modeling requirements of fusion reactor blanket maintenance systems, such as the DEMOnstration power plant (DEMO).

## Overview
Real-time control of floating flexible multi-body systems requires dynamic models that balance high-fidelity structural representation with strict computational efficiency. This codebase provides an optimized, fully vectorized computational framework for assembling the state-dependent mass matrix and generalized force vectors. 

The model integrates the floating frame of reference formulation (FFRF) with Craig-Bampton component mode synthesis to reduce elastic coordinates. To overcome the computational bottleneck of sequential node-by-node summation loops, the implementation systematically reformulates the assembly algorithm to utilize bulk matrix outer products and native array operations. This explicitly delegates the arithmetic workload to SIMD-accelerated underlying linear algebra libraries.

## Performance Benchmarks
Based on a 1-second numerical simulation of an 80t flexible payload ($E=200$ GPa, $\nu=0.3$, $\rho=3650~kg/m^3$), the vectorized framework demonstrates significant computational advantages:
* **Execution Speed:** Achieves a Real-Time Factor (RTF) of 1.0989, running thirteen times faster (13.3x speedup) than a standard iterative baseline (RTF of 0.0829).
* **Commercial Solver Superiority:** Executes approximately three to six times faster than Simulink ODE15s (RTF of 0.1887) and vastly outperforms the Simulink DAESSC solver (RTF of 0.0710).
* **Accuracy:** Accurately captures transient modal deflections and maintains properly decoupled macro-dynamics without the constraint deviations characteristic of penalty methods.

## System Requirements
The simulations and benchmarks were executed using the following software and hardware environment:
* **Software:** MathWorks MATLAB R2026a (version 26.1.0.3346908).
* **Required Toolboxes:** Partial Differential Equation Toolbox and Simulink.
* **Hardware Profile:** Apple M5 processor with 16 GB of unified memory.

## License (Enforced Open Source)
This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. 

This is a strong copyleft license that ensures the software remains free and open-source. **If you use, modify, or distribute this software in your own projects, you are strictly required to release your derivative work under the same open-source GPL-3.0 license.** This enforces that any advancements or adaptations built upon this framework remain open and accessible to the research community. 

## Acknowledgments
This work has been carried out within the framework of the EUROfusion Consortium, funded by the European Union via the Euratom Research and Training Programme (Grant Agreement No 101052200 - EUROfusion). 

## Citation
If you utilize this code or framework in your academic work, please cite the original publication:
> Park, Y., & Sloth, C. (2026). *Vectorized Computational Framework for Floating Flexible Bodies Using Component Mode Synthesis*. 

**Repository Link:** [https://github.com/yhpark64/iccma2026](https://github.com/yhpark64/iccma2026)