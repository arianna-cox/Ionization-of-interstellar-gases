import numpy as np
from scipy.integrate import quad
from Q1_23_5 import Radius, Tstar_sun
from constant_23_5 import c, h, k, nu0
from tabulate import tabulate
from Q2_23_5 import alphaB

#interstellar gas density
n_H = 1e6

#Interstellar gas temperatures
T = [5e3, 1e4, 2e4]

# Luminosity of the stars
L = np.array([4e29, 4e30, 3.9e26])
# Surface temperature of the stars
Tstar = np.array([2e4, 2.5e4, Tstar_sun])
# Radius of stars
R = [Radius(L[0], Tstar[0]), Radius(L[1], Tstar[1]), 6.96e8]

r_1_approx = []
Q = []
for i in range(3):
    constant = (8*np.pi**2*R[i]**2*nu0**3)/(c**2)
    def integrand(x):
        return x**2/(np.exp(h*nu0*x/(k*Tstar[i]))-1)
    QH = constant*quad(integrand, 1, 30)[0]
    Q.append(QH)
    # find r_1
    for j in range(3):
        r_1_approx.append((3*QH/(4*np.pi*alphaB(T[j])*n_H**2))**(1/3))


table = [['str',0.1,0.1, 0.1],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

r1_Q3 = np.load('Q3_r1.npy')
for i in range(3):
    table[i][0] = "7 solar mass"
    table[i+3][0] = "12 solar mass"
    table[i+6][0] = "sun"
    table[i][1] = T[i]
    table[i+3][1] = T[i]
    table[i+6][1] = T[i]
for i in range(9):
    table[i][2] = r_1_approx[i]
    table[i][3] = r1_Q3[i]


head = ["Star", "Q(H) /s^-1"]
print(tabulate( np.column_stack([["7 solar mass", "12 solar mass", "sun"], Q]), headers=head))

head = ["Star", "Gas temperature /K", "Question 4 r_1 /m", "Question 3 r_1 /m"]
print(tabulate(table, headers=head))