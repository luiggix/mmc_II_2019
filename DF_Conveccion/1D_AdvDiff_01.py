#!/usr/bin/env pyhton
"""
    Construye y resuelve el sistema Ax = b, producto de aplicar \n
    el método de diferencias finitas a la ecuación de Poisson   \n
    en una dimensión, con condiciones de frontera tipo Dirchelet.

Ejecución:

    Solicita del usuario las dimensiones del dominio, \n
    el numero de nodos y el valor de las condiciones.
"""
import numpy as np

def Laplaciano1D(N, h, Gamma, rho, v):
    """
        Devuelve la matriz A, del sistema Ax = b. 
        Params:
            N: int
                Número de incognitas. \n
            diagonal: escalar
                Valor de la diagonal principal. \n
        Return:
            A: ndarray
                Matriz NxN, con tres bandas distintas a cero.
    """
    a =  Gamma / h**2
    b =  rho * v / (2*h)
    A = np.zeros((N,N))
    A[0,0] = 2 * a 
    A[0,1] = b - a
    for i in range(1,N-1):
        A[i,i] = 2 * a 
        A[i,i+1] = b - a
        A[i,i-1] = -b - a
    A[N-1,N-2] = -b - a
    A[N-1,N-1] = 2 * a
    return A

# Comienza la ejecución del programa
print()
print("+----------------------------------------------------+")
print("|      Solucion de la ecuacion de Laplace en 1D      |")
print("+----------------------------------------------------+")
print("| Autor: Luis M. de la Cruz S.                       |")
print("+----------------------------------------------------+")
print("|                 Datos de entrada                   |")
print("+----------------------------------------------------+")

def analyticSol(x):
    return (np.exp(rho * v * x / Gamma) - 1) / (np.exp(rho * v * L / Gamma) - 1) * (phiL - phi0) + phi0

L = 1.0 # m
rho = 1.0 # kg/m^3
v = 2.1 # m/s
Gamma = 0.1 # kg / m.s
phi0 = 1.0 #
phiL = 0.0 #
N = 20 # Número de incógnitas
h = L / (N+1)

print("| El tamanio de la malla es      : h = %g " % h)

# Definicion del sistema lineal de N x N
f = np.zeros(N)         # Lado derecho
A = Laplaciano1D(N, h, Gamma, rho, v) # Matriz del sistema

# Aplicacion de las cond. Dirichlet para el calculo correcto de la sol.
f[0]   =  phi0 * (rho * v / (2*h) + Gamma / h**2) 
f[N-1] = -phiL * (rho * v / (2*h) - Gamma / h**2)

# La solucion sera guardada en el arreglo u, que es de tamanio N+2, pues incluye las fronteras
u = np.zeros(N+2)

# Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
u[1:N+1] = np.linalg.solve(A,f)

# Los valores en los extremos son conocidos debido a las cond. Dirichlet
u[0] = phi0
u[N+1] = phiL

print("\n Lado derecho del sistema : \n", f)
print("\n Matriz del sistema : \n", A)
print("\n Solucion del sistema : \n", u)

input('\n Presione <enter> para generar la grafica ')

import matplotlib.pyplot as plt

x = np.linspace(0,L,N+2)
plt.plot(x,u,'bo--')

x1 = np.linspace(0,L,100)
ua = analyticSol(x1)
plt.plot(x1,ua,'k-')

plt.show()


