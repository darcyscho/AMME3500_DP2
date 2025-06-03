from sympy import symbols, Function, diff, sin, cos, simplify, latex, pprint, Derivative as D, Matrix, N, solve, expand, collect, srepr, zeros 
from scipy.linalg import solve_continuous_are
import numpy as np 

t = symbols('t', real=True)
g = symbols('g', real=True)
b1, b2, b3 = symbols('b1 b2 b3')    # Damping coefficients

# === thigh link ===
q1 = Function('q1')(t)
m1 = symbols('m1', real=True)
d1 = symbols('d1', real=True)
x1 = d1*cos(q1)
y1 = d1*sin(q1)
x1_dt = diff(x1, t)
y1_dt = diff(y1, t)
v1_sqrd = simplify(x1_dt**2+y1_dt**2)

T1 = 0.5*m1*v1_sqrd
V1 = m1*g*y1
D1 = 0.5*b1*D(q1, t)**2

#print(V1)

# === shin link ===
q2 = Function('q2')(t)
m2 = symbols('m2', real=True)
l1 = symbols('l1', real=True)
d2 = symbols('d2', real=True)
x2 = l1*cos(q1)+d2*cos(q1 + q2) 
y2 = l1*sin(q1)+d2*sin(q1 + q2)
x2_dt = diff(x2, t)
y2_dt = diff(y2, t)
v2_sqrd = simplify(x2_dt**2+y2_dt**2)

T2 = 0.5*m2*v2_sqrd
V2 = m2*g*y2
D2 = 0.5*b2*D(q2, t)**2

#print(v2_sqrd)

# === foot link ===
q3 = Function('q3')(t)
m3 = symbols('m3', real=True)
l2 = symbols('l2', real=True)
d3 = symbols('d3', real=True)
x3 = l1*cos(q1)+l2*cos(q1+q2)+d3*cos(q1+q2+q3)
y3 = l1*sin(q1)+l2*sin(q1+q2)+d3*sin(q1+q2+q3) 
x3_dt = diff(x3, t)
y3_dt = diff(y3, t)
v3_sqrd = simplify(x3_dt**2+y3_dt**2)

T3 = 0.5*m3*v3_sqrd
V3 = m3*g*y3
D3 = 0.5*b3*D(q3, t)**2

#print(v3_sqrd)

l3 = symbols('l3', real=True)

# === total energies ===
T = T1 + T2 + T3
V = V1 + V2 + V3
Diss = D1 + D2 + D3
L = simplify(T-V)

# === Euler-Lagrange Equations ===
EL1 = diff(diff(L, D(q1, t)), t) - diff(L, q1) + diff(Diss, diff(q1, t))
EL2 = diff(diff(L, D(q2, t)), t) - diff(L, q2) + diff(Diss, diff(q2, t))
EL3 = diff(diff(L, D(q3, t)), t) - diff(L, q3) + diff(Diss, diff(q3, t))

# === Actuator Inputs ===
tau1 = symbols('tau1', real=True)
tau2 = symbols('tau2', real=True)
tau3 = symbols('tau3', real=True)

eq1 = simplify(EL1-tau1)
eq2 = simplify(EL2-tau2)
eq3 = simplify(EL3-tau3)

params = {
    l1: 0.45,
    l2: 0.44,
    l3: 0.26,
    d1: 0.225,
    d2: 0.22,
    d3: 0.13,
    m1: 8.4,
    m2: 3.8,
    m3: 1.114,
    g: -9.81
}

# Define equilibrium config 
equilibrium_values = {
    q1: 1.75,
    q2: 1.13,
    q3: 2.18,
    D(q1, t): 0,
    D(q2, t): 0,
    D(q3, t): 0,
    D(q1, (t, 2)): 0,
    D(q2, (t, 2)): 0,
    D(q3, (t, 2)): 0
}

# Combine with physical parameters
subs_all = {**params, **equilibrium_values}

# Solve E-L equations symbolically for tau_i at equilibrium
tau1_e = solve(eq1.subs(subs_all), tau1)[0]
tau2_e = solve(eq2.subs(subs_all), tau2)[0]
tau3_e = solve(eq3.subs(subs_all), tau3)[0]

# Display results
from sympy import N
print("\n--- Equilibrium Torques (3 DOF) ---")
print("tau1_e =", N(tau1_e))
print("tau2_e =", N(tau2_e))
print("tau3_e =", N(tau3_e))

