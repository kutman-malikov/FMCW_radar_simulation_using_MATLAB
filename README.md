# FMCW Radar Simulation (MATLAB)

## Overview
This repository contains a MATLAB-based simulation of an **FMCW (Frequency-Modulated Continuous Wave) radar system**.  
The project demonstrates the fundamental principles of FMCW radar operation, including chirp generation, signal reflection, beat frequency extraction, and range/velocity estimation.

The simulator is intended for **educational and research purposes**, bridging electromagnetic theory, signal processing, and practical radar algorithms.

---

## Objectives
- Demonstrate the operating principle of FMCW radar
- Simulate range and velocity estimation using beat frequency
- Implement dechirp processing and FFT-based spectral analysis
- Provide clear visualization of radar signals and results
- Serve as a foundation for more advanced radar simulations

---

## FMCW Radar Principle
In FMCW radar, a continuous signal with **linearly modulated frequency (chirp)** is transmitted.  
The received signal is a delayed and Doppler-shifted version of the transmitted signal.

After mixing (dechirping), the resulting **beat frequency** is proportional to the target range:

\[
R = \frac{c f_{beat}}{2 S}
\]

where:
- \( c \) — speed of light  
- \( f_{beat} \) — beat frequency  
- \( S \) — chirp slope  

Using triangular chirps (up and down), the Doppler component can be separated, enabling velocity estimation.

---

## Features
- Baseband FMCW radar simulation
- Linear chirp signal generation
- Target motion modeling (range & velocity)
- Dechirp processing
- FFT-based beat frequency extraction
- Range estimation
- Real-time visualization
- Adjustable radar and target parameters

---

## Project Structure
