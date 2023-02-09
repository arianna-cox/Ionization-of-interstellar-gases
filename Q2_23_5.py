import numpy as np
from math import exp
from scipy.integrate import quad, odeint
from constant_23_5 import c, h, k, nu0, a_nu0

# solving the ionization equations for the hydrogen and proton densities
def alphaB(T):
    if T == 5000:
        return 4.54e-19
    if T == 10000:
        return 2.59e-19
    if T == 20000:
        return 2.52e-19

def integral(tau_nu0, R, Tstar, T):
    constant = (R**2*np.pi*nu0**3*a_nu0)/(c**2*alphaB(T))
    def integrand(x):
        return exp(-tau_nu0/(x**3))*exp(-(h*nu0*x)/(k*Tstar))/(x*(1-exp(-(h*nu0*x)/(k*Tstar))))
    I = quad(integrand, 1, 10)
    return constant*I[0]

# n_H0 as a function of tau_nu0 and r
def n_H0(tau_nu0, r, R, Tstar, T, n_H):
    mess = (integral(tau_nu0, R, Tstar, T))/(r**2)
    return n_H + mess - np.sqrt(2*n_H*mess + mess**2)

# neutral fraction
def n_H0_scaled(tau_nu0_initial, stepsize, N, R, Tstar, T, n_H):
    def function(tau_nu0, r):
        return n_H0(tau_nu0, r, R, Tstar, T, n_H)*a_nu0
    r_vector = R + stepsize * np.array(list(range(N+1)))
    tau_nu0_vector = odeint(function, tau_nu0_initial, r_vector)
    n_H0_vector = [(n_H0(tau_nu0_vector[i], r_vector[i], R, Tstar, T, n_H))/n_H for i in range(N + 1)]
    return n_H0_vector
