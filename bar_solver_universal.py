"""
Universal 1D Finite Element Solver for Axial Bar Structures
-----------------------------------------------------------
Works for ANY number of bar elements connected in series (stepped bar),
with different length, area and Young's modulus for each element.

The program asks for all input values first, then prints:
  - global stiffness matrix
  - nodal displacements
  - support reactions + equilibrium check
  - element strain, stress and axial force
  - natural frequencies and mode shapes (optional)

Sign convention: +x direction to the right. Tension = positive stress.
"""

import numpy as np
from scipy.linalg import eigh


# ---------------------------------------------------------------
# Helper functions for safe input
# ---------------------------------------------------------------
def ask_float(prompt, positive=False):
    while True:
        try:
            value = float(input(prompt).strip())
            if positive and value <= 0:
                print("  Value must be greater than zero. Try again.")
                continue
            return value
        except ValueError:
            print("  Please enter a number.")


def ask_int(prompt, minimum=1):
    while True:
        try:
            value = int(input(prompt).strip())
            if value < minimum:
                print(f"  Enter a whole number >= {minimum}.")
                continue
            return value
        except ValueError:
            print("  Please enter a whole number.")


def ask_choice(prompt, options):
    options_lower = {o.lower(): o for o in options}
    while True:
        value = input(prompt).strip().lower()
        if value in options_lower:
            return options_lower[value]
        print(f"  Choose one of: {', '.join(options)}")


def ask_yes_no(prompt):
    return ask_choice(prompt + " (y/n): ", ["y", "n"]) == "y"


# ---------------------------------------------------------------
# Unit conversion tables (everything is converted to SI)
# ---------------------------------------------------------------
LENGTH_TO_M = {"m": 1.0, "cm": 1e-2, "mm": 1e-3}
AREA_TO_M2 = {"m2": 1.0, "cm2": 1e-4, "mm2": 1e-6}
E_TO_PA = {"Pa": 1.0, "MPa": 1e6, "GPa": 1e9}
FORCE_TO_N = {"N": 1.0, "kN": 1e3}


# ---------------------------------------------------------------
# Element matrices
# ---------------------------------------------------------------
def bar_stiffness(E, A, L):
    return (E * A / L) * np.array([[1, -1], [-1, 1]])


def bar_mass(rho, A, L):
    """Consistent mass matrix of a 2-node bar element."""
    return (rho * A * L / 6) * np.array([[2, 1], [1, 2]])


# ---------------------------------------------------------------
# Main program
# ---------------------------------------------------------------
def main():
    print("=" * 60)
    print("  UNIVERSAL 1D BAR FINITE ELEMENT SOLVER")
    print("=" * 60)

    # ---------- Step 1: number of elements ----------
    n_elem = ask_int("\nNumber of elements: ")
    n_nodes = n_elem + 1
    conn = [(e, e + 1) for e in range(n_elem)]
    print(f"  -> {n_elem} elements and {n_nodes} nodes "
          f"(element k joins node k and node k+1)")

    # ---------- Step 2: units ----------
    print("\nChoose units for your input values")
    len_unit = ask_choice("  Length unit (m / cm / mm): ", list(LENGTH_TO_M))
    area_unit = ask_choice("  Area unit (m2 / cm2 / mm2): ", list(AREA_TO_M2))
    e_unit = ask_choice("  Young's modulus unit (Pa / MPa / GPa): ", list(E_TO_PA))
    f_unit = ask_choice("  Force unit (N / kN): ", list(FORCE_TO_N))

    # ---------- Step 3: material and geometry ----------
    print("\nMaterial and geometry")
    same_E = ask_yes_no("  Same Young's modulus for all elements?")
    if same_E:
        E_common = ask_float(f"  Young's modulus E ({e_unit}): ", positive=True) * E_TO_PA[e_unit]

    L = np.zeros(n_elem)
    A = np.zeros(n_elem)
    E = np.zeros(n_elem)
    for e in range(n_elem):
        print(f"\n  Element {e+1} (node {e+1} to node {e+2})")
        L[e] = ask_float(f"    Length ({len_unit}): ", positive=True) * LENGTH_TO_M[len_unit]
        A[e] = ask_float(f"    Cross-section area ({area_unit}): ", positive=True) * AREA_TO_M2[area_unit]
        if same_E:
            E[e] = E_common
        else:
            E[e] = ask_float(f"    Young's modulus ({e_unit}): ", positive=True) * E_TO_PA[e_unit]

    # ---------- Step 4: loads ----------
    print(f"\nNodal loads (positive = pulling in +x direction, unit: {f_unit})")
    print("  Enter node number and force. Press Enter on an empty node to finish.")
    F = np.zeros(n_nodes)
    while True:
        text = input("  Load at node number (blank to finish): ").strip()
        if text == "":
            break
        try:
            node = int(text)
        except ValueError:
            print("  Please enter a whole number.")
            continue
        if not 1 <= node <= n_nodes:
            print(f"  Node must be between 1 and {n_nodes}.")
            continue
        F[node - 1] += ask_float(f"    Force at node {node} ({f_unit}): ") * FORCE_TO_N[f_unit]

    # ---------- Step 5: supports ----------
    print("\nSupports")
    while True:
        text = input(f"  Fixed node numbers, separated by commas (1 to {n_nodes}): ").strip()
        try:
            fixed = sorted({int(t) - 1 for t in text.split(",") if t.strip() != ""})
        except ValueError:
            print("  Enter numbers only, for example: 1,4")
            continue
        if not fixed:
            print("  At least one node must be fixed, otherwise the bar can move freely.")
            continue
        if any(n < 0 or n >= n_nodes for n in fixed):
            print(f"  Node numbers must be between 1 and {n_nodes}.")
            continue
        break
    free = [n for n in range(n_nodes) if n not in fixed]

    # ---------- Step 6: vibration option ----------
    do_modal = ask_yes_no("\nDo free-vibration (modal) analysis too?")
    rho = np.zeros(n_elem)
    if do_modal:
        same_rho = ask_yes_no("  Same density for all elements?")
        if same_rho:
            rho[:] = ask_float("  Density (kg/m3): ", positive=True)
        else:
            for e in range(n_elem):
                rho[e] = ask_float(f"  Density of element {e+1} (kg/m3): ", positive=True)

    # ===========================================================
    # All inputs collected. Now solve.
    # ===========================================================
    K = np.zeros((n_nodes, n_nodes))
    M = np.zeros((n_nodes, n_nodes))
    for e, (i, j) in enumerate(conn):
        ke = bar_stiffness(E[e], A[e], L[e])
        me = bar_mass(rho[e], A[e], L[e])
        idx = [i, j]
        for a in range(2):
            for b in range(2):
                K[idx[a], idx[b]] += ke[a, b]
                M[idx[a], idx[b]] += me[a, b]

    print("\n" + "=" * 60)
    print("  INPUT SUMMARY")
    print("=" * 60)
    print(f"{'Element':<9}{'L (m)':>12}{'A (m2)':>14}{'E (Pa)':>14}")
    for e in range(n_elem):
        print(f"{e+1:<9}{L[e]:>12.4f}{A[e]:>14.4e}{E[e]:>14.4e}")
    print("Loads (N):", {f"node {n+1}": float(F[n]) for n in range(n_nodes) if F[n] != 0})
    print("Fixed nodes:", [n + 1 for n in fixed])

    print("\nGlobal stiffness matrix K (N/m):")
    print(np.array2string(K, formatter={"float_kind": lambda x: f"{x:12.4e}"}))

    # ---------- Static solution ----------
    U = np.zeros(n_nodes)
    K_ff = K[np.ix_(free, free)] if free else np.zeros((0, 0))
    if free:
        try:
            U[free] = np.linalg.solve(K_ff, F[free])
        except np.linalg.LinAlgError:
            print("\nERROR: stiffness matrix is singular. Check supports and inputs.")
            return
    R = K @ U - F

    print("\n" + "=" * 60)
    print("  STATIC RESULTS")
    print("=" * 60)
    print("\nNodal displacements:")
    print(f"{'Node':<7}{'u (m)':>16}{'u (mm)':>14}")
    for n in range(n_nodes):
        print(f"{n+1:<7}{U[n]:>16.6e}{U[n]*1e3:>14.6f}")

    print("\nSupport reactions (N):")
    for n in fixed:
        print(f"  R{n+1} = {R[n]:.3f}")
    residual = R[fixed].sum() + F.sum()
    print(f"\nEquilibrium check: sum of reactions + applied loads = {residual:.3e} N (should be ~0)")

    print("\nElement results (tension +, compression -):")
    print(f"{'Element':<9}{'Strain':>14}{'Stress (MPa)':>16}{'Force (N)':>14}")
    for e, (i, j) in enumerate(conn):
        strain = (U[j] - U[i]) / L[e]
        stress = E[e] * strain
        force = stress * A[e]
        print(f"{e+1:<9}{strain:>14.4e}{stress/1e6:>16.4f}{force:>14.2f}")

    # ---------- Modal solution ----------
    if do_modal:
        print("\n" + "=" * 60)
        print("  MODAL (FREE-VIBRATION) RESULTS")
        print("=" * 60)
        if not free:
            print("No free nodes, so there are no vibration modes.")
        else:
            M_ff = M[np.ix_(free, free)]
            vals, vecs = eigh(K_ff, M_ff)
            vals = np.clip(vals, 0, None)
            omega = np.sqrt(vals)
            freq = omega / (2 * np.pi)
            print(f"\n{'Mode':<7}{'omega (rad/s)':>16}{'Frequency (Hz)':>18}")
            for k in range(len(omega)):
                print(f"{k+1:<7}{omega[k]:>16.2f}{freq[k]:>18.2f}")
            print("\nMode shapes (rows = free nodes "
                  f"{[n + 1 for n in free]}, columns = modes, mass-normalised):")
            print(np.round(vecs, 4))

    print("\nDone.")


if __name__ == "__main__":
    main()
