from statistics import mean
import numpy as np
import matplotlib.pyplot as plt 
import time

# Start timer to measure total execution time
begin = time.time()

# Number of particles (or rotors) in the system
N = 10000

# Integration time step
h = 0.05

# Initial energy density
E_in = 0.5

# Arrays to store observables during the simulation
Trotter = []   # time evolution
Kinet = []     # kinetic energy
Pot = []       # potential energy
erro = []      # energy conservation error

Mag = []       # magnetization
E_Mag = []     # magnetization energy

# Time array
T = np.arange(0,100,h)

#### Initial condition ####

# Position (angles of rotors)
q = np.zeros(N)

# Momentum of rotors
p = np.zeros(N)


# ---------------------------------------------------
# Function that evolves positions using momenta
# q(t+dt) = q(t) + p*dt
# ---------------------------------------------------
def evolveQ(Q,P,dt):
    Q += P*dt
    return Q


# ---------------------------------------------------
# Function that evolves momenta using the mean-field force
# The force depends on the global magnetization
# ---------------------------------------------------
def evolveP(Q,P,dt):

    # Magnetization components
    Mx = np.sum(np.cos(Q))/float(N)
    My = np.sum(np.sin(Q))/float(N)   

    # Mean-field force acting on each rotor
    F = np.cos(Q)*My - np.sin(Q)*Mx
    
    # Momentum update
    P += F*dt

    return P


# ---------------------------------------------------
# Kinetic energy per particle
# K = (1/2N) Σ p_i²
# ---------------------------------------------------
def K(P):
    return 0.5*np.sum(np.square(P))/float(N)


# ---------------------------------------------------
# Potential energy per particle
# V = (1/2)(1 - M²)
# where M² = Mx² + My²
# ---------------------------------------------------
def V(Q):

    Mx = np.sum(np.cos(Q))/float(N)     
    My = np.sum(np.sin(Q))/float(N)

    M2 = Mx*Mx + My*My 
    
    return(0.5*(1.- M2))


# ---------------------------------------------------
# Total energy per particle
# H = K + V
# ---------------------------------------------------
def H(K,V):
    return(K+V)


# ---------------------------------------------------
# Water-bag initial condition
# This initializes positions and momenta uniformly
# inside a finite region of phase space.
# The routine rescales the momenta to ensure the
# desired initial energy density.
# ---------------------------------------------------
def water_bag(Q,P):
    
    deltaq = 0.5*np.pi
    
    # Initial magnetization corresponding to deltaq
    M0 = (2./deltaq)*( np.sin(deltaq/2.) ) 
    
    # Momentum width obtained from the energy constraint
    deltap = np.sqrt( 24.*E_in - 12.*( 1. - M0*M0 ) )

    deltap = 0.5*deltap
    deltaq = 0.5*deltaq
    
    while True:
        
        # Uniform distribution of angles
        Q = deltaq*(2.*np.random.rand(N) - 1.)

        # Uniform distribution of momenta
        P = deltap*(2.*np.random.rand(N) - 1.)
        
        # Raw kinetic energy
        k_w = 0.5*np.sum(P*P)
        
        # Remove center-of-mass momentum
        p_cm = np.sum(P)/float(N)
        P -= np.ones(N)*p_cm  

        # Kinetic energy after removing momentum drift
        k_s = 0.5*np.sum(P*P)
        
        # Rescaling factor
        A = np.sqrt(k_w/k_s)
        
        # Rescale momenta
        P = A*P
        
        # Compute total energy
        bac = H(K(P),V(Q))
        
        # Accept configuration if energy is correct
        if( np.abs(E_in - bac) < 1e-3 ):
            break

    return(Q,P)


# Generate initial condition
q,p = water_bag(q,p)

# Initial energy
e0 = H(K(p),V(q))

print(e0)

# Initial total momentum (should be approximately zero)
p0 = np.sum(p)/N


# Initial magnetization
Mx = np.sum(np.cos(q))/float(N)     
My = np.sum(np.sin(q))/float(N)

M_in = np.sqrt(Mx*Mx + My*My) 


# ---------------------------------------------------
# Fourth-order symplectic integrator coefficients
# (Forest-Ruth / Yoshida scheme)
# ---------------------------------------------------
s = 1/(4-4**(1/3))

delta1 = h*s
delta2 = 0.5*(h*s)
delta3 = 0.5*h*(1-3*s)
delta4 = (1-4*s)*h


# ---------------------------------------------------
# Time evolution loop
# ---------------------------------------------------
for i in range(len(T)):

    print(i)

    # Symplectic integration steps
    evolveQ(q,p,delta2)
    evolveP(q,p,delta1)

    evolveQ(q,p,delta1)
    evolveP(q,p,delta1)

    evolveQ(q,p,delta3)
    evolveP(q,p,delta4)

    evolveQ(q,p,delta3)
    evolveP(q,p,delta1)

    evolveQ(q,p,delta1)
    evolveP(q,p,delta1)

    evolveQ(q,p,delta2)
    
    # Compute magnetization
    Mx = np.sum(np.cos(q))/float(N)     
    My = np.sum(np.sin(q))/float(N)

    M2 = np.sqrt(Mx*Mx + My*My)