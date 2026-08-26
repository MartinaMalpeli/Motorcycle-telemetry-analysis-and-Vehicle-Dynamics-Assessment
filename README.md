# 🏍️ Motorcycle Dynamics and Stability Analysis

<div align="center">

![Status](https://img.shields.io/badge/status-complete-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)
![Python](https://img.shields.io/badge/Python-3.8+-blue)
![MATLAB](https://img.shields.io/badge/MATLAB-R2021b+-red)

**Analysis of motorcycle weave and wobble modes through experimental acquisition and numerical modeling**

[📋 Description](#description) • [🎯 Objectives](#objectives) • [📊 Methodology](#methodology) • [📁 Structure](#structure) • [🚀 Usage](#usage)

</div>

---

## 📋 Description

This research project analyzes the dynamic behavior of motorcycles under external excitations, with particular focus on the **weave mode** (low-frequency oscillation 2-4.5 Hz involving coupled roll, yaw, and steering motions).

The study combines:
- ✅ **Experimental acquisition** on road with MEMS sensors
- ✅ **Advanced data processing** (FFT, Hilbert Transform, filtering)
- ✅ **Numerical modal analysis** with 4-DOF model
- ✅ **Parametric simulations** for validation and optimization

**Author:** Martina Malpeli  
**Institution:** Politecnico di Milano - Master in Mechanical Engineering (Sports Engineering)  
**Academic Year:** 2025/2026

---

## 🎯 Objectives

| Objective | Description |
|-----------|-------------|
| 🔍 **Identification** | Isolation and characterization of vibration modes (weave, wobble) |
| 📈 **Quantification** | Measurement of natural frequencies, damping ratios, and decay rates |
| 🔄 **Validation** | Comparison between experimental data and numerical predictions |
| 🛠️ **Optimization** | Identification of critical parameters for stability at various speeds |
🏍️ Data Acquisition
├── Gyroscope on handlebar (steering axis)
├── Gyroscope on rear frame (yaw axis)
├── GPS (position and velocity)
└── IMU 6-DOF

📍 Test Route: 6 segments at increasing speeds
├── Seg 1-2: ~10-40 km/h
├── Seg 3-4: ~90-110 km/h
└── Seg 5-6: ~100-160 km/h
---

## 📊 Methodology

### 1️⃣ **Experimental Phase**
