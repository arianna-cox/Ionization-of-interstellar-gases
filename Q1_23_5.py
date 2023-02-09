import numpy as np
from scipy.special import zeta
from constant_23_5 import c, h, k

def Radius(L,Tstar):
    return ((L*c**2*h**3)/(48*np.pi**2*k**4*Tstar**4*zeta(4)))**(0.5)

R_sun = 6.96e8
L_sun = 3.90e26
Tstar_sun = (h/k)*((L_sun*c**2)/(48*np.pi**2*h*R_sun**2*zeta(4)))**(1/4)

if __name__ == '__main__':
    print(f"T* = {Tstar_sun:.2e} K")

    for [L,Tstar] in [[4e29,20000], [4e30,25000]]:
        R = Radius(L,Tstar)
        print(f"R = {R:.2e} m")