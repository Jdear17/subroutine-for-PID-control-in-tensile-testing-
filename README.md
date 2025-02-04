# PID Temperature Control in Fully Coupled Thermomechanical-Electrical-Structural FEA

## Overview
This project implements a **PID-based temperature control** system within a fully coupled **thermomechanical-electrical-structural** **Finite Element Analysis (FEA)** model in **Abaqus**. The control system is defined using **Abaqus user subroutines** (`UAMP` and `UHARD`), which adjust the applied **heat flux** based on sensor feedback to maintain a target temperature profile.

## Features
- **Custom PID Controller** for precise temperature regulation.
- **Fully coupled physics** integration in Abaqus.
- **Dynamic adaptation of material properties** based on strain rate and temperature.
- **Robust numerical implementation** with safeguards against integral windup and derivative noise.

## Subroutines
### `UAMP`
- Implements the **PID control algorithm**.
- Adjusts the **applied heat flux** based on the temperature sensor feedback.
- Includes safeguards such as **integral windup prevention** and **derivative smoothing**.
- Uses an **exponential function for control output** with clamping to avoid overflow.

### `UHARD`
- Defines **strain rate- and temperature-dependent material hardening**.
- Uses **Johnson-Cook plasticity model** with modifications for strain rate and temperature.
- Ensures stability by enforcing **yield strength lower bounds**.

## Installation
1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-repo-name.git
   cd your-repo-name

## Installation

### Setup Abaqus User Subroutines:
- Place `UAMP` and `UHARD` in the Abaqus working directory.
- Ensure Abaqus is configured to compile Fortran subroutines.

### Run the Abaqus simulation:
```bash
abaqus job=your_simulation_name user=UAMP,UHARD
```

## Usage
- Modify the **PID gains** (`Kp`, `Ki`, `Kd`) in `UAMP` for better control tuning.
- Adjust **material properties** in `UHARD` according to specific material data.
- Modify **temperature setpoints** in `UAMP` to change the control profile.

## Tuning the PID Controller
- Start with **Kd = 0**, then increase **Kp** until stable control is achieved.
- If oscillations occur, increase **Kd** to dampen response.
- If steady-state error remains, increase **Ki**, but beware of integral windup.

## License
This project is licensed under the **MIT License**.

## Contact
For any issues or contributions, please create an **issue** or submit a **pull request** on GitHub.
