import numpy as np
from scipy.integrate import quad, odeint
from constant_23_5 import c, h, k, nu0, a_nu0
import matplotlib.pyplot as plt
from tabulate import tabulate
# solving the ionization equations for the hydrogen and proton densities for a quasar

def alphaB(T):
    if T == 5000:
        return 4.54e-19
    if T == 10000:
        return 2.59e-19
    if T == 20000:
        return 2.52e-19

def F(tau_nu0, T):
    constant = (10**24*a_nu0)/(8*np.pi*alphaB(T)*h)
    def integrand(x):
        return x**(-5.4)*np.exp(-(x/10)-(tau_nu0/(x**3)))
    I = quad(integrand, 1, 2)
    return constant*I[0]

def n_H0(tau_nu0, r, T, n_H):
    mess = (F(tau_nu0, T))/(r**2)
    answer = n_H + mess - np.sqrt(2*n_H*mess + mess**2)
    return answer

def neutral_fraction(tau_nu0_initial, stepsize, N, R, T, n_H):
    def function(tau_nu0, r):
        return n_H0(tau_nu0, r, T, n_H)*a_nu0
    r_vector = R + stepsize * np.array(list(range(N+1)))
    tau_nu0_vector = odeint(function, tau_nu0_initial, r_vector)
    # plt.plot(r_vector, tau_nu0_vector)
    # plt.show()
    n_H0_vector = [(n_H0(tau_nu0_vector[i], r_vector[i], T, n_H))/n_H for i in range(N + 1)]
    return n_H0_vector

# value of tau_nu0 at r=R
tau_nu0_initial = 0
#interstellar gas density
n_H = 1e6
#Interstellar gas temperatures
T = 1e4
R = 1e18

max_r = 3*1e20
N = 6000
stepsize = (max_r - R) / N
r_vector = R + stepsize * np.array(list(range(N + 1)))
n_H0_vector = neutral_fraction(tau_nu0_initial, stepsize, N, R, T, n_H)

# tabulate n_H0
sample_point = np.linspace(0, N, 16).astype(int)
table = np.zeros((len(sample_point),2))
for i in range(len(sample_point)):
    table[i][0] = r_vector[sample_point[i]]
    table[i][1] = n_H0_vector[sample_point[i]]

head = ["r /m", "Neutral fraction"]
print(tabulate(table, headers=head))