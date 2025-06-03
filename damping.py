from sympy import symbols, Function, diff, sin, cos, simplify, Matrix, Eq, solve, lambdify
import numpy as np
from numpy.linalg import eigvals

# === Symbolic setup ===
t = symbols('t', real=True)
g = symbols('g', real=True)
b1, b2, b3 = symbols('b1 b2 b3')  # Damping coefficients

# Joint angles
q1 = Function('q1')(t)
q2 = Function('q2')(t)
q3 = Function('q3')(t)

# Parameters
m1, m2, m3 = symbols('m1 m2 m3')
l1, l2, l3 = symbols('l1 l2 l3')
d1, d2, d3 = symbols('d1 d2 d3')

# Positions
x1 = d1*cos(q1); y1 = d1*sin(q1)
x2 = l1*cos(q1) + d2*cos(q1 + q2); y2 = l1*sin(q1) + d2*sin(q1 + q2)
x3 = l1*cos(q1) + l2*cos(q1 + q2) + d3*cos(q1 + q2 + q3)
y3 = l1*sin(q1) + l2*sin(q1 + q2) + d3*sin(q1 + q2 + q3)

# Velocities
v1_sq = simplify(diff(x1, t)**2 + diff(y1, t)**2)
v2_sq = simplify(diff(x2, t)**2 + diff(y2, t)**2)
v3_sq = simplify(diff(x3, t)**2 + diff(y3, t)**2)

# Energies
T = 0.5*m1*v1_sq + 0.5*m2*v2_sq + 0.5*m3*v3_sq
V = m1*g*y1 + m2*g*y2 + m3*g*y3
D = 0.5*b1*diff(q1, t)**2 + 0.5*b2*diff(q2, t)**2 + 0.5*b3*diff(q3, t)**2
L = T - V

# Inputs and derivatives
tau1, tau2, tau3 = symbols('tau1 tau2 tau3')
ddq1, ddq2, ddq3 = symbols('ddq1 ddq2 ddq3')
dq1, dq2, dq3 = symbols('dq1 dq2 dq3')
subs_derivs = {
    diff(q1, t): dq1, diff(q1, (t, 2)): ddq1,
    diff(q2, t): dq2, diff(q2, (t, 2)): ddq2,
    diff(q3, t): dq3, diff(q3, (t, 2)): ddq3,
}

# Substitutable parameters
param_values = {
    l1: 0.45, l2: 0.44, l3: 0.26,
    d1: 0.225, d2: 0.22, d3: 0.13,
    m1: 8.4, m2: 3.8, m3: 1.114,
    g: -9.81
}

# State + input definitions
x1_, x2_, x3_, x4_, x5_, x6_ = symbols('x1 x2 x3 x4 x5 x6')
u1_, u2_, u3_ = symbols('u1 u2 u3')

subs_state = {
    q1: x1_, q2: x2_, q3: x3_,
    dq1: x4_, dq2: x5_, dq3: x6_,
    tau1: u1_, tau2: u2_, tau3: u3_
}

# Initial conditions
x_eq = [1.75, 1.13, 2.18, 0, 0, 0]
u_eq = [19.2546, 12.0830, -0.4840]

# === Sweep damping values ===
for b_val in range(1, 500):
    damping_values = {b1: b_val, b2: b_val, b3: b_val}
    subs_all = {**param_values, **damping_values}

    # Recompute Lagrange + Dissipative forces
    EL1 = simplify(diff(diff(L, diff(q1, t)), t) - diff(L, q1) + diff(D, diff(q1, t))).subs(subs_derivs).subs(subs_all) - tau1
    EL2 = simplify(diff(diff(L, diff(q2, t)), t) - diff(L, q2) + diff(D, diff(q2, t))).subs(subs_derivs).subs(subs_all) - tau2
    EL3 = simplify(diff(diff(L, diff(q3, t)), t) - diff(L, q3) + diff(D, diff(q3, t))).subs(subs_derivs).subs(subs_all) - tau3

    # Solve for ddqi
    try:
        sol1 = solve(Eq(EL1, 0), ddq1)[0]
        sol2 = solve(Eq(EL2, 0), ddq2)[0]
        sol3 = solve(Eq(EL3, 0), ddq3)[0]
    except:
        continue

    # Substitute in state variables
    f1 = sol1.subs(subs_state)
    f2 = sol2.subs(subs_state)
    f3 = sol3.subs(subs_state)

    # Linearise
    stateFunction = Matrix([x4_, x5_, x6_, f1, f2, f3])
    A = stateFunction.jacobian(Matrix([x1_, x2_, x3_, x4_, x5_, x6_]))
    A_func = lambdify((x1_, x2_, x3_, x4_, x5_, x6_, u1_, u2_, u3_), A, modules='numpy')

    # Evaluate numerically
    try:
        A_eval = np.array(A_func(*x_eq, *u_eq)).astype(np.float64)
        eigs = eigvals(A_eval)
        if np.all(np.real(eigs) < 0):
            print(f"\n Stable A matrix at damping b = {b_val}")
            print("Eigenvalues:", eigs)
            break
    except Exception as e:
        print(f"b = {b_val} failed: {e}")

