from sympy import symbols, Function, diff, sin, cos, simplify, Derivative as D

# this script makes sure the equilibrium inputs/joint angles are correct (makes equations = ~0)

# === Setup ===
t = symbols('t', real=True)
g = symbols('g', real=True)

# Generalized coordinates
q1 = Function('q1')(t)
q2 = Function('q2')(t)
q3 = Function('q3')(t)

# Physical parameters
m1, m2, m3 = symbols('m1 m2 m3')
d1, d2, d3 = symbols('d1 d2 d3')
l1, l2, l3 = symbols('l1 l2 l3')
b1, b2, b3 = symbols('b1 b2 b3')    # Damping coefficients

# === Energies ===
x1 = d1*cos(q1)
y1 = d1*sin(q1)
v1_sq = diff(x1, t)**2 + diff(y1, t)**2
T1 = 0.5 * m1 * v1_sq
V1 = m1 * g * y1
D1 = 0.5*b1*D(q1, t)**2

x2 = l1*cos(q1) + d2*cos(q1 + q2)
y2 = l1*sin(q1) + d2*sin(q1 + q2)
v2_sq = diff(x2, t)**2 + diff(y2, t)**2
T2 = 0.5 * m2 * v2_sq
V2 = m2 * g * y2
D2 = 0.5*b2*D(q2, t)**2

x3 = l1*cos(q1) + l2*cos(q1+q2) + d3*cos(q1+q2+q3)
y3 = l1*sin(q1) + l2*sin(q1+q2) + d3*sin(q1+q2+q3)
v3_sq = diff(x3, t)**2 + diff(y3, t)**2
T3 = 0.5 * m3 * v3_sq
V3 = m3 * g * y3
D3 = 0.5*b3*D(q3, t)**2

T = T1 + T2 + T3
V = V1 + V2 + V3
Diss = D1 + D2 + D3
L = T - V

# === E-L Equations ===
tau1, tau2, tau3 = symbols('tau1 tau2 tau3')
EL1 = diff(diff(L, D(q1, t)), t) - diff(L, q1) + diff(Diss, diff(q1, t)) 
EL2 = diff(diff(L, D(q2, t)), t) - diff(L, q2) + diff(Diss, diff(q2, t)) 
EL3 = diff(diff(L, D(q3, t)), t) - diff(L, q3) + diff(Diss, diff(q3, t)) 

# === Substitution Dictionaries ===

params = {
    l1: 0.45, l2: 0.44, l3: 0.26,
    d1: 0.225, d2: 0.22, d3: 0.13,
    m1: 8.4, m2: 3.8, m3: 1.114,
    g: -9.81
}

equilibrium = {
    q1: 2.62, q2: 1.75, q3: 1.40,
    D(q1, t): 0, D(q2, t): 0, D(q3, t): 0,
    D(D(q1, t), t): 0, D(D(q2, t), t): 0, D(D(q3, t), t): 0,
    tau1: 37.9203694366734, tau2: 3.03663048162315, tau3: -1.33120123436457
}

subs_all = {**params, **equilibrium}

# === Evaluate ===
res1 = (EL1-tau1).subs(subs_all).evalf()
res2 = (EL2-tau2).subs(subs_all).evalf()
res3 = (EL3-tau3).subs(subs_all).evalf()

print("Residual of EL1 at equilibrium:", res1)
print("Residual of EL2 at equilibrium:", res2)
print("Residual of EL3 at equilibrium:", res3)

