# FMCW Radar Simulation (MATLAB)

## Overview
This repository contains a MATLAB simulation of an **FMCW (Frequency-Modulated Continuous Wave) radar**.  
The project demonstrates how an FMCW radar measures target distance and velocity using frequency-modulated signals and digital signal processing.

The main goal is to clearly show the **operating principle of FMCW radar** and how it can be implemented in software.

---

## Project Goals
- Explain the working principle of FMCW radar in a simple and visual way
- Simulate radar signal transmission and reception
- Extract target range and velocity from the received signal
- Connect radar theory with practical signal processing
- Provide a clean and understandable reference project

---

## How the Radar Works
The radar transmits a continuous signal whose frequency changes linearly over time (chirp).  
The signal reflected from a target returns with a delay and a frequency shift caused by motion.

By comparing the transmitted and received signals, the system obtains a beat signal.  
The frequency of this beat signal contains information about the target distance, and its change allows estimation of target velocity.

---

## What the Code Does
1. Generates a frequency-modulated radar signal  
2. Models a moving target with distance and velocity  
3. Applies time delay and Doppler shift to the received signal  
4. Mixes transmitted and received signals (dechirp processing)  
5. Uses FFT to extract beat frequencies  
6. Estimates target range and velocity  
7. Visualizes signals and results in real time

---

## Features
- Baseband FMCW radar simulation
- Adjustable radar and target parameters
- Range and velocity estimation
- FFT-based signal analysis
- Real-time plots and visualization
- Clean and readable MATLAB implementation

---

## Project Structure
├── code/ # MATLAB simulation code

├── report/ # Technical report

└── README.md

---

## Assumptions and Limitations
- Single target only
- No noise or clutter
- Ideal frequency modulation
- Simplified radar channel model

These assumptions were made to keep the focus on understanding the core radar principles.

---

## Future Improvements
- Multi-target support
- Noise and interference modeling
- Range–Doppler processing
- CFAR detection
- Hardware or SDR-based implementation

---

## Author
Kutman Malikov 

Bahcesehir Cyprus University

Electronics Engineering

kutman.engineer@gmail.com

---

## License
This project is intended for educational and research use.
