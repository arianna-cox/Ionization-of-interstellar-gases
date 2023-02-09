import numpy as np
import matplotlib.pyplot as plt
from Q1_23_5 import Radius, Tstar_sun
from Q2_23_5 import n_H0_scaled
from tabulate import tabulate

plt.rcParams.update({'font.size': 15})

# value of tau_nu0 at r=R
tau_nu0_initial = 0

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


def r_1(n_H0_vector, n_p_vector, r_vector):
    difference = np.abs(np.subtract(n_H0_vector, n_p_vector))
    return r_vector[np.argmin(difference)]

max_r = [4e17, 1e18, 7e13]
N = 10000
r_1_vector = []

for i in range(3):
    stepsize = (max_r[i] - R[i]) / N
    r_vector = R[i] + stepsize * np.array(list(range(N + 1)))
    for j in range(3):
        n_H0_vector= n_H0_scaled(tau_nu0_initial, stepsize, N, R[i], Tstar[i], T[j], n_H)
        n_p_vector = [1 - n_H0_vector[q] for q in range (N+1)]

        # plt.plot(r_vector, n_H0_vector)
        plt.plot(r_vector, n_p_vector)

        # find r_1
        r_1_vector.append(r_1(n_H0_vector, n_p_vector, r_vector))
    plt.xlabel('r')
    plt.ylabel(r'$n_{p}/n_H$')
    plt.xlim([0, r_vector[-1]])
    plt.ylim([0, 1])
    # plt.gca().legend([r'$n_{H^0}/n_H$', r'$n_{p}/n_H$'])
    plt.gca().legend(["T = 5000K", "T = 10000K", "T = 20000K"])
    plt.show()

np.save('Q3_r1.npy', r_1_vector)

table = [['str',0.1,0.1],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]

for i in range(3):
    table[i][0] = "7 solar mass"
    table[i+3][0] = "12 solar mass"
    table[i+6][0] = "sun"
    table[i][1] = T[i]
    table[i+3][1] = T[i]
    table[i+6][1] = T[i]
for i in range(9):
    table[i][2] = r_1_vector[i]

head = ["Star", "Gas temperature /K", "r_1 /m"]
print(tabulate(table, headers=head, tablefmt="grid"))