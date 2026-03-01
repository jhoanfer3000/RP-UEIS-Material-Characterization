# RP-UEIS-Material-Characterization
This repository contains the MATLAB scripts, COMSOL Multiphysics models, and experimental data required to run the Resonant Piezoelectric Spectroscopy coupled with Ultrasound Electrical Impedance Spectroscopy (RP-UEIS) framework. [cite_start]This low-cost, impedance-driven framework combines Finite Element (FE) modeling and inverse optimization to extract the complex electromechanical constants of piezoelectric actuators and the viscoelastic moduli of polymer layers[cite: 10].


## ⚙️ Prerequisites and Software Requirements

To successfully run the scripts in this repository, you need:
* **MATLAB**: Version **R2024a** (or higher). The **Optimization Toolbox** must be installed (required for constrained optimization functions).
* **COMSOL Multiphysics**: Version **6.1a** (or higher).
* **COMSOL LiveLink™ for MATLAB®**: You must launch MATLAB via the COMSOL LiveLink shortcut to establish the connection between the numerical FEM models and the optimization scripts.

## 📂 Repository Structure

The repository is organized into three main directories to cleanly separate the source code, finite element models, and experimental/computed data:

* `Scripts/`: Contains all MATLAB source code, optimization algorithms, and objective cost functions.
* `models/`: Contains the COMSOL Multiphysics (`.mph`) finite element models for the standalone PZT and the coupled material assemblies.
* `data/`: Contains the raw experimental electrical impedance data (`.mat`). During execution, the converged PZT properties (`theta.mat`) and the final viscoelastic properties of the tested materials (`theta_s.mat`) are automatically saved here.

## 🧠 Main Algorithms

The core methodology is divided into three sequential MATLAB scripts located in the `Scripts` folder. [cite_start]This pipeline bridges the gap between raw component evaluation and predictive multi-physics modeling[cite: 43]:

### 1. Algorithm 1 (`Algorithm_1.m`) - PZT Sensitivity Analysis
This script performs a sensitivity analysis on the material parameters of the piezoelectric disk. [cite_start]By analyzing the local response of the cost function to individual perturbations, it quantifies how variations in these parameters propagate to the model outputs[cite: 353]. [cite_start]This crucial step reduces the high-dimensional parameter space to only the most influential constants, significantly improving parameter identifiability for the subsequent inversion[cite: 35].

### 2. Algorithm 2 (`Algorithm_2.m`) - PZT Model Calibration
[cite_start]This script calibrates the standalone PZT disc model by fitting the simulated electrical impedance to the measured spectra[cite: 34]. [cite_start]Using the reduced parameter set identified in Algorithm 1, it runs an inverse parameter estimation procedure to reconstruct the full frequency response of the PZT stack[cite: 203]. The final, converged electromechanical properties of the actuator are saved as `theta.mat` in the `data` folder.

### 3. Algorithm 3 (`Algorithm_3.m`) - Viscoelastic Material Adjustment
This is the main material extraction script. [cite_start]It extends the model to include a polymeric sample acoustically coupled to the actuator via an ultrasonic gel layer[cite: 36]. [cite_start]By comparing the coupled FEM model against experimental data, it runs a constrained optimization to adjust and extract the frequency-dependent complex viscoelastic moduli of the material while accounting for interface damping[cite: 37]. The final tested material parameters (Storage and Loss Bulk/Shear moduli) are saved as `theta_s.mat` in the `data` folder.

## 🚀 How to Run

1. Open **COMSOL Multiphysics 6.1a with MATLAB** (This starts the LiveLink server).
2. In MATLAB, navigate your "Current Folder" to the `Scripts` directory of this repository.
3. Run `Algorithm_1.m` to perform the sensitivity analysis on the base piezoelectric sensor and determine the active optimization parameters.
4. Run `Algorithm_2.m` to calibrate the PZT model. This will generate and save the baseline `theta.mat` file.
5. Run `Algorithm_3.m` to perform the final inversion and extract the viscoelastic properties of your coupled polymer sample, which will be saved as `theta_s.mat`.

