"""
1D Finite Element Solver for a Stepped Bar
Static analysis (displacements, stresses, reactions) + free-vibration (modal) analysis
Units: SI (N, m, Pa, kg)
"""

import numpy as np
from scipy.linalg import eigh

# ---------------- Input data ----------------
E = 207e9                                   # Young's modulus, Pa
rho = 7850.0                                # steel density, kg/m^3
A = np.array([5e-4, 3e-4, 1e-4])            # areas, m^2 (5, 3, 1 cm^2)
L = np.array([0.20, 0.30, 0.40])            # lengths, m
P = 2000.0                                  # axial load at node 3, N

n_nodes = 4
conn = [(0, 1), (1, 2), (2, 3)]             # element connectivity (0-based)


# ---------------- Element matrices ----------------
def bar_stiffness(E, A, L):
    return (E * A / L) * np.array([[1, -1],
                                   [-1, 1]])


def bar_mass(rho, A, L):
    """Consistent mass matrix of a 2-node bar element."""
    return (rho * A * L / 6) * np.array([[2, 1],
                                         [1, 2]])


# ---------------- Global assembly ----------------
K = np.zeros((n_nodes, n_nodes))
M = np.zeros((n_nodes, n_nodes))
for e, (i, j) in enumerate(conn):
    ke = bar_stiffness(E, A[e], L[e])
    me = bar_mass(rho, A[e], L[e])
    idx = [i, j]
    for a in range(2):
        for b in range(2):
            K[idx[a], idx[b]] += ke[a, b]
            M[idx[a], idx[b]] += me[a, b]

print("Global stiffness matrix K (x1e8 N/m):")
print(np.round(K / 1e8, 4), "\n")

# ---------------- Boundary conditions ----------------
F = np.zeros(n_nodes)
F[2] = P                                    # P acts at node 3
fixed = [0, 3]                              # nodes 1 and 4 fixed (u = 0)
free = [n for n in range(n_nodes) if n not in fixed]

# ---------------- Static solution ----------------
U = np.zeros(n_nodes)
K_ff = K[np.ix_(free, free)]
U[free] = np.linalg.solve(K_ff, F[free])
R = K @ U - F

print("Nodal displacements (m):")
for i, u in enumerate(U, start=1):
    print(f"  u{i} = {u:.6e}")

print("\nReaction forces (N):")
for n in fixed:
    print(f"  R{n+1} = {R[n]:.3f}")

# Equilibrium check: reactions + applied load = 0
print(f"\nEquilibrium check: sum of reactions + P = {R[fixed].sum() + P:.6e} N (should be ~0)")

# ---------------- Element stresses ----------------
print("\nElement results:")
print(f"{'Element':<9}{'Strain':>14}{'Stress (MPa)':>16}{'Force (N)':>14}")
for e, (i, j) in enumerate(conn):
    strain = (U[j] - U[i]) / L[e]
    stress = E * strain
    force = stress * A[e]
    print(f"{e+1:<9}{strain:>14.4e}{stress/1e6:>16.4f}{force:>14.2f}")

# ---------------- Modal (free-vibration) analysis ----------------
# Generalized eigenvalue problem on free DOFs: (K - w^2 M) phi = 0
M_ff = M[np.ix_(free, free)]
vals, vecs = eigh(K_ff, M_ff)
omega = np.sqrt(vals)
freq = omega / (2 * np.pi)

print("\nNatural frequencies:")
for k, (w, f) in enumerate(zip(omega, freq), start=1):
    print(f"  Mode {k}: {w:.2f} rad/s = {f:.2f} Hz")

print("\nMode shapes (rows: node 2, node 3; columns: modes; mass-normalised):")
print(np.round(vecs, 4))
