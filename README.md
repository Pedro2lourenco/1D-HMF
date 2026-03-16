# Hamiltonian Mean Field (HMF) Model Simulation

This repository contains a **numerical simulation of the Hamiltonian Mean Field (HMF) model**, a paradigmatic system in statistical physics used to study **long-range interacting systems, phase transitions, and quasi-stationary states**.

The system is evolved using a **fourth-order symplectic integrator**, ensuring good energy conservation during long simulations.

---

# Model Description

The **Hamiltonian Mean Field (HMF) model** describes a system of \(N\) classical rotors interacting through a mean-field potential.

The Hamiltonian is

$$
H = \sum_{i=1}^{N} \frac{p_i^2}{2} + \frac{1}{2N}\sum_{i,j}\left[1 - \cos(\theta_i - \theta_j)\right]
$$

Where:

- \( \theta_i \) — angular position of rotor \(i\)
- \( p_i \) — momentum of rotor \(i\)
- \(N\) — number of particles

The potential energy can be written in terms of the **magnetization vector**

$$
M_x = \frac{1}{N}\sum_i \cos(\theta_i)
$$

$$
M_y = \frac{1}{N}\sum_i \sin(\theta_i)
$$

$$
M = \sqrt{M_x^2 + M_y^2}
$$

Thus,

$$
V = \frac{1}{2}(1 - M^2)
$$

This formulation allows the computation of the potential energy with **O(N)** complexity instead of **O(N²)**.

---

# Numerical Method

The simulation uses a **fourth-order symplectic integrator (Forest–Ruth / Yoshida scheme)**.

Symplectic integrators are ideal for Hamiltonian systems because they preserve the **phase-space structure and improve long-term energy conservation**.

The equations of motion are

$$
\dot{q}_i = p_i
$$

$$
\dot{p}_i = -\frac{\partial V}{\partial q_i}
$$

---

# Initial Conditions

The simulation starts from a **water-bag distribution**, commonly used in studies of long-range interacting systems.

Properties of the initialization:

- uniform distribution of angles
- uniform distribution of momenta
- rescaling of momenta to match the desired **initial energy density**
- removal of center-of-mass momentum


# Requirements

The code requires:

- Python 3
- NumPy
- Matplotlib

Install dependencies:

```bash
pip install numpy matplotlib
