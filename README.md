# MeSA 3D Diffusion and Communication Effectiveness Indicator (CEI) Simulation System

This project is a MATLAB-based 3D convection-diffusion Partial Differential Equation (PDE) solver. It is designed to simulate the spatial diffusion of molecules like Methyl Salicylate (MeSA) under various wind scenarios and transmitter (TX) topologies, and to evaluate their Communication Effectiveness Indicator (CEI) in macroscopic molecular communication.

## ⚙️ Features

* **3D Convection-Diffusion Solver**: Utilizes an Explicit Euler scheme for time stepping, central difference for the diffusion term, and a first-order upwind scheme for the advection (convection) term.
* **Multiple Wind Scenarios**: Supports 4 predefined wind environments, incorporating diurnal variation factors and random noise:
    1.  *Scenario 1*: Baseline open field
    2.  *Scenario 2*: Sheltered / calm
    3.  *Scenario 3*: Windy day
    4.  *Scenario 4*: Gusty / unstable
* **Multiple TX Topologies**: Automatically generates and maps 4 source node physical layouts to the grid: `central`, `uniform`, `corners`, and `perimeter`.
* **Dynamic Release Model**: Uses an experimentally fitted model for the dynamic mass release rate `q(t)` over a 24-hour period.
* **CEI Evaluation System**: Calculates the volume fraction of the spatial coverage at different concentration thresholds, as well as the dynamic evolution of this indicator over time.

---

## 📂 File Structure

### Main Scripts
* **`main_fig6.m`**: Reproduces the diffusion heatmaps (Fig. 6). It simulates the `central` and `corners` topologies across the 4 wind scenarios and outputs 2D concentration heatmaps (Z-axis slice average) at specific hours (1h, 11h, 22h).
* **`main_fig7.m`**: Reproduces the CEI evaluation curves (Fig. 7). It iterates through all 4 TX topologies across the 4 wind scenarios, calculating and plotting **CEI vs. Concentration Threshold** (24h average) and **CEI vs. Time**.

### Core Physics & Environmental Models
* **`step_pde.m`**: The core PDE step solver. It includes 3D Laplacian diffusion, upwind advection calculation, and boundary conditions (absorbing ground at z=0 and zero-gradient Neumann boundaries on the sides/top).
* **`wind_field.m`**: Wind field generation model. Generates the 3D wind velocity field `(vx, vy, vz)` based on the input time `t` and the 3D spatial grid.
* **`release_model.m`**: Dynamic source release model. Returns the step-by-step mass release rate (in g/s) over a total of 24 hours.
* **`build_source.m`**: Source term builder. Allocates the continuous-time source release rate into the corresponding discrete 3D spatial grid matrix.

### Spatial & Layout Tools
* **`tx_positions.m`**: Generates the physical coordinates `(x, y, z)` for the source transmitters based on the specified topology mode.
* **`map_tx_to_grid.m`**: Maps and deduplicates the continuous physical coordinates of the transmitters to discrete 3D PDE grid indices.
* **`run_CEI_for_mode.m`**: A wrapper script that runs a full 24-hour CEI calculation for a single specified TX topology mode and returns the results.

---

## 🚀 Quick Start

### 1. Requirements
This project only requires a basic MATLAB environment (R2020a or newer is recommended). No additional toolboxes (like Deep Learning or PDE toolboxes) are needed.

### 2. Running the Simulations

* **Generate Concentration Heatmaps (Fig. 6)**:
  Run `main_fig6` directly in the MATLAB Command Window or Editor. The script will sequentially calculate and display the concentration evolution for the 4 wind scenarios.
  
* **Generate CEI Evaluation Curves (Fig. 7)**:
  Run `main_fig7`. This script iterates through 16 combinations (4 wind scenarios x 4 layouts), simulating a 24-hour physical process (43,200 time steps) for each. Once completed, it will output the corresponding statistical line charts.

### ⚠️ Performance Note
Because the explicit PDE solver for the 3D grid must satisfy the CFL stability condition, the current settings (grid `dx=5, dy=5, dz=0.5`, time step `dt=2` seconds) are computationally heavy. If you only wish to quickly verify the code logic, you can temporarily change `T_hours = 24` to a smaller value (e.g., `1` or `2` hours) inside `main_fig6.m` and `main_fig7.m`.
