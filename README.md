# Analysis of Motorcycle Dynamics and Stability 🏍️

**Institution:** Politecnico di Milano[cite: 3]  
**Program:** Master of Science in Mechanical Engineering - Sports Engineering (CC6)[cite: 3]  
**Course:** Dynamics of mechanical systems[cite: 3]  
**Author:** Martina Malpeli (A.Y. 2025/2026)[cite: 3]  

---

## 📖 Project Overview
This repository contains the analytical tools, methodologies, and findings developed to evaluate the dynamic behavior and stability of motorcycles under external excitations[cite: 3]. The objective is to identify setups and operating conditions that maximize rider safety, vehicle control, and stability[cite: 3].

The study focuses primarily on isolating and evaluating the **Weave** mode, a low-frequency (2-4.5 Hz) global eigenmode characterized by coupled oscillations in roll, yaw, and steer involving the entire vehicle chassis[cite: 3]. 

### State of the Art & Scope
Literature identifies three main modes of motorcycle vibration:
1. **Weave:** Low-frequency coupled yaw/roll motion, highly dependent on frame stiffness, tire properties, and aerodynamics[cite: 3].
2. **Wobble:** High-frequency (6-10 Hz) oscillation of the front assembly[cite: 3].
3. **Capsize:** Non-oscillatory eigenmode dominated by roll motion[cite: 3].

This project intentionally deprioritizes Wobble because accurate representation requires complex elastic multibody simulations (EMBS) and tire relaxation properties, focusing instead on Weave as the primary handling limit for global trajectory control[cite: 3]. The mathematical foundation relies on a 4-Degrees of Freedom (4-DOF) model established by R.S. Sharp (1971) considering lateral displacement, yaw, roll, and steer[cite: 3].

---

## 🔬 1. Experimental Analysis

The experimental phase relies on a real road test designed to excite the motorcycle at various racing speeds to collect transient dynamic responses[cite: 3].

### Data Acquisition & Telemetry
The vehicle was equipped with a comprehensive sensor suite:
*   **Handlebar & Tail Triaxle MEMS**[cite: 3]
*   **Front Triaxle MEMS**[cite: 3]
*   **Steering Gyroscope (Voltmeter):** To measure angular speed on the steering axis[cite: 3].
*   **IMU 6DOF + GPS:** To obtain latitude, longitude, and yaw rate[cite: 3].

The test track consisted of 6 distinct segments (Go and Return phases), merging GPS peaks to map the spatial velocity gradient[cite: 3]. Two noticeable stops occurred during the acquisition, which impacted the vehicle's setup and subsequent data[cite: 3].

### Signal Processing
1.  **Filtering:** Gyroscope signals were filtered using a high-pass filter (2 Hz cutoff) to remove low-frequency noise and drift[cite: 3].
2.  **Impulse Detection:** Multiple excitation events (introduced by the rider) were manually isolated via interactive MATLAB windows overlapping steer, yaw rate, and speed[cite: 3].
3.  **Hilbert Transform:** The peaks were approximated using the Hilbert Transform to evaluate the instantaneous amplitude envelope[cite: 3].
4.  **Damping Calculation:** The decay rate ($\sigma$) was extracted via linear regression on the natural logarithm of the envelope[cite: 3]. The modal damping ratio ($\zeta$) was computed as[cite: 3]:
    $$\zeta = \frac{\sigma}{\sqrt{\sigma^2 + \omega_d^2}} * 100$$
    *(where $\omega_d$ is the dominant damped natural frequency obtained via FFT)*[cite: 3].

### Experimental Findings
*   **Frequency vs. Speed:** Frequencies hovered around 2.5 Hz initially but showed a progressive shift toward higher frequencies in the later segments of the test[cite: 3].
*   **Stability (Root Locus):** Experimental Gauss Diagrams mapped complex eigenvalues ($\lambda = -\sigma \pm j\omega_d$)[cite: 3]. As speed increased, the real part ($\sigma$) moved dangerously close to the $0 \text{ rad/s}$ stability limit, indicating reduced damping capacity[cite: 3].
*   **Wobble Identification:** Supplementary band-pass filtering (8-11 Hz) detected wobble peaks primarily in the 60-80 km/h range, confirming theoretical expectations[cite: 3].

---

## 💻 2. Numerical Analysis & Simulation

To validate experimental anomalies (like the frequency shift in later segments), a 4-DOF numerical model was used to simulate how geometric variations and load conditions impact stability[cite: 3]. 

### Stability Drivers Evaluated
*   **Geometric & Mass Parameters:** Center of Mass (CoM), wheel inertia, steering inertia, and rake angle[cite: 3].
*   **Structural Properties:** Front fork and frame torsional stiffness[cite: 3].
*   **Tires & Damping:** Steering damper and tire relaxation length[cite: 3].
*   **Rider Influence:** Rider body damping and suspended mass[cite: 3].

### Sensitivity Testing Results
*   **Steering Damper (+/- 20%):** Proved to be an almost negligible parameter, though reducing it slightly delayed the corrective steering response[cite: 3].
*   **Front Mass Center (+/- 40%):** Lowering the CoM optimized the feedback loop between steering and lean angle, improving mid-speed stability[cite: 3]. A high CoM introduced a destabilizing phase lag[cite: 3].
*   **Wheelbase (+/- 40%):** Increasing the wheelbase led to significant instability at high speeds[cite: 3].
*   **Load Conditions:** 
    *   *Rider Only (75 kg):* Provided aerodynamic damping and grip at high speeds, while slightly raising the CoG at low speeds (managed by rider balance)[cite: 3].
    *   *Extreme Loads:* Adding a passenger and bags created a dangerous region of instability around 120 km/h[cite: 3].

### Time-Domain Analysis & Optimization
Simulations of an impact response at 120 km/h (the critical speed threshold) were performed[cite: 3]. Spectral diagrams of the standard setup confirmed peaks for Weave (2-4 Hz) and Wobble (9.5 Hz)[cite: 3].

Several setups were tested to flatten the weave peaks:
*   **Excessive Inertia Modifications:** Disrupt the balance between kinetic energy and damping capacity[cite: 3]. The phase lag between perturbation and corrective tire forces transforms stabilizing feedback into positive reinforcement, causing divergent, unstable oscillations[cite: 3].
*   **Optimal Setup ("Try 2"):** A combination of modified front trail (+3%), lowered rear/front CoM (-10%), and reduced wheelbase (-10%) proved to be the best compromise[cite: 3]. It significantly reduced weave vibrations and provided a convergent, stable response[cite: 3].

---

## 🎯 Conclusions
1.  **Validation:** The experimental methodology successfully isolated the Weave mode, confirming that it becomes critical at high speeds as damping diminishes[cite: 3].
2.  **Test Anomalies Explained:** The numerical model confirmed that the shifts in frequency and damping observed in the later segments of the experimental test were likely caused by a load-reduction of the bike during the brief track stops[cite: 3].
3.  **Design Impact:** The project highlights that while theoretical setups can flatten frequency spectra, careful balance is required in the time domain to avoid divergent impulse responses[cite: 3].

---
## 📚 Bibliography
*   Sharp, R. S. (1971). *Stability, control and steering responses of motorcycles.*[cite: 3]
*   Sharp, R. S. (2001). *The stability and control of motorcycles.*[cite: 3]
*   Passigato, F. (2015). *Study of motorcycle dynamics for the improvement of the stability of weave and wobble eigenmodes.*[cite: 3]
*   Dębowski, A. (2017). *Analysis of the effect of mass parameters on motorcycle vibration and stability.*[cite: 3]
