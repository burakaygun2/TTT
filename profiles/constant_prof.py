import numpy as np
from scipy.optimize import fsolve, root

rho_ice   = 930.0
rho_ocean = 1100.0
rho_hpice = 1250.0
rho_core  = 3300.0

R_ice     = 2575.0e3
R_mantle  = 2000.0e3

Mass = 1.345e23
MOI  = 0.3414

comp   = [0.5, 1.0]
bice   = [167]
# docean = np.linspace(1.0, 50.0, 96)
docean = np.linspace(1.0, 20.0, 48)
eta    = [2.0]
print(eta[0:10])
print(eta[10:21])

def equations(vars):
    x1, x2 = vars
    AA = Mass*3.0/(4.0*np.pi) - (rho_ice*(R_ice**3 - R_ocean**3) + rho_ocean*(R_ocean**3 - R_hpice**3) + rho_hpice*(R_hpice**3 - R_mantle**3))
    eq1 = x1 * (R_mantle**3 - x2**3) + rho_core * (x2**3) - AA
    BB = MOI*Mass*15.0*R_ice**2/(8.0*np.pi) -(rho_ice*(R_ice**5 - R_ocean**5) + rho_ocean*(R_ocean**5 - R_hpice**5) + rho_hpice*(R_hpice**5 - R_mantle**5))
    eq2 = x1 * (R_mantle**5 - x2**5) + rho_core * (x2**5) - BB
    return [eq1, eq2]
    # x1, x2 = vars
    # AA = Mass*3.0/(4.0*np.pi) - (rho_ice*(R_ice**3 - R_ocean**3) + rho_ocean*(R_ocean**3 - R_hpice**3) + rho_hpice*R_hpice**3 + rho_core*R_core**3)
    # eq1 = -rho_hpice * x2**3 + x1 * (x2**3 - R_core**3) - AA
    # BB = MOI*Mass*15.0*R_ice**2/(8.0*np.pi) - (rho_ice*(R_ice**5 - R_ocean**5) + rho_ocean*(R_ocean**5 - R_hpice**5) + rho_hpice*R_hpice**5 + rho_core*R_core**5)
    # eq2 = -rho_hpice * x2**5 + x1 * (x2**5 - R_core**5) - BB
    # return [eq1, eq2]

def chebyshev_nodes(a, b, N, rr):
        """
        Generate N Chebyshev nodes in the interval [a, b], sorted from b (outer) to a (inner).
        """
        k = np.arange(N)
        # x = np.cos(np.pi * (2*k + 1) / (2*N))
        for ir in range(N):
            if ir == 0:
                x = 1.0
            elif ir == N-1:
                x = -1.0
            else:
                x = np.cos( (2.0*(ir-1)+1.0)*np.pi/(2.0*(N-1)) )
            rr[ir] = 0.5 * (a + b) + 0.5 * (b - a) * x
        # Map from [-1, 1] to [a, b]
        return rr[::-1]

def grid(R0, R1, n, type):
    if type == 'h':
        return np.linspace(R0, R1, n)
    elif type == 'n':
        return chebyshev_nodes(R0, R1, n, np.zeros(n))

for c in comp:
    for b in bice:
        R_ocean = R_ice - b*1e3
        for d in docean:
            R_hpice = R_ocean - d*1e3
            initial_guess = [2000.0e0, 500.0e3]
            rho_mantle, R_core = fsolve(equations, initial_guess)
            print(rho_mantle, R_core)
            total_mass = (4.0/3.0)*np.pi*(rho_ice*(R_ice**3 - R_ocean**3) + rho_ocean*(R_ocean**3 - R_hpice**3) + rho_hpice*(R_hpice**3 - R_mantle**3) + rho_mantle*(R_mantle**3 - R_core**3) + rho_core*R_core**3)
            tota_moi   = (8.0/15.0)*np.pi*(rho_ice*(R_ice**5 - R_ocean**5) + rho_ocean*(R_ocean**5 - R_hpice**5) + rho_hpice*(R_hpice**5 - R_mantle**5) + rho_mantle*(R_mantle**5 - R_core**5) + rho_core*R_core**5)
            print(total_mass, Mass, tota_moi/(total_mass*R_ice**2), MOI)
            for vis in eta:
                file = open(f"./prof_set1/Titan_{b}km_{d:.2f}km_{vis:.2f}Pas.dat", "w")
                file.write(f"{'#Radius (m)':>30}{'Density (kg/m^3)':>31}{'Bulk modulus (Pa)':>31}{'Shear modulus (Pa)':>31}{'Viscosity (Pa.s)':>31}\n")
                Nr = 70
                r = grid(R_ice, R_ice-20., Nr, 'h')
                for ir in range(Nr):
                    file.write(f"{r[ir]:30.6f} {rho_ice:30.6f} {1.0e30:30.6e} {3.5e9:30.6e} {24.0:30.6f}\n")
                Nr = 130
                r = grid(R_ice-20.0, R_ocean, Nr, 'h')
                for ir in range(Nr):
                    file.write(f"{r[ir]:30.6f} {rho_ice:30.6f} {1.0e30:30.6e} {3.5e9:30.6e} {15.0:30.6f}\n")
                r = grid(R_ocean, R_hpice, Nr, 'n')
                for ir in range(Nr):
                    file.write(f"{r[ir]:30.6f} {rho_ocean:30.6f} {1.0e30:30.6e} {0.0:30.6e} {vis:30.6f}\n")
                r = grid(R_hpice, R_mantle, Nr, 'h')
                for ir in range(Nr):
                    file.write(f"{r[ir]:30.6f} {rho_hpice:30.6f} {1.0e30:30.6e} {6.0e9:30.6e} {15.0:30.6f}\n")
                r = grid(R_mantle, R_core, Nr, 'h')
                for ir in range(Nr):
                    file.write(f"{r[ir]:30.6f} {rho_mantle:30.6f} {1.0e30:30.6e} {30.0e9:30.6e} {17.0:30.6f}\n")
                r = grid(R_core, 10.0, Nr, 'h')
                for ir in range(Nr):
                    file.write(f"{r[ir]:30.6f} {rho_core:30.6f} {1.0e30:30.6e} {50.0e9:30.6e} {22.0:30.6f}\n")