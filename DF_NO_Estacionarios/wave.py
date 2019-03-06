#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 17:03:18 2017

@author: luiggi
"""

import numpy as np
import matplotlib.pyplot as plt

def f(x):
    """
    Calcula la condición inicial del problema en todo el dominio.
    x : arreglo que contiene las coordenadas de la malla.
    """
    return np.sin(np.pi * x)

def g(x):
    """
    Calcula la velocidad inicial del problema.
    x : arreglo que contiene las coordenadas de la malla.
    """
    return 0

def solExacta(x, t):
    """
    Calcula la solución exacta del problema.
    x : arreglo que contiene las coordenadas de la malla.
    t : instante de tiempo en que se calcula la solución exacta.
    """
    return np.sin(np.pi * x) * np.cos(2 * np.pi * t)

def calcError(sol_n,sol_e):
    """
    Calcula el error de la solución aproximada con respecto a la exacta.
    sol_n : arreglo que contiene la solución numérica.
    sol_e : arreglo que contiene la solución exacta.
    """
    return np.abs(sol_n-sol_e)

def condicionesIniciales(l, ht, u, x):
    """
    Calcula la condición inicial en el instante n = 1
    l : representa a lambda = alpha * ht / h.
    ht : paso de tiempo.
    u : condición inicial en el instante correspondiente a n = 0.
    x : arreglo que contiene las coordenadas de la malla.
    """
    N = len(u)
    w = np.zeros(N)
    for i in range(1,N-1):
        w[i] = (1 - l**2) * u[i] + 0.5 * l**2 * (u[i+1] + u[i-1]) + ht * g(x[i])
#        w[i] = u[i] + ht * g(x[i])
    return w

def solver(u, w, N, x, Nt, l):
    """
    Calcula la solución del problema para n = 1, 0, ..., Nt
    u : condición inicial en el instante correspondiente a n = 0.
    w : condición inicial en el instante correspondiente a n = 1
    N : Número total de incógnitas internas.
    x : arreglo que contiene las coordenadas de la malla.
    Nt : Número total de pasos en el tiempo.
    l : representa a lambda = alpha * ht / h.
    """
    s = np.zeros(N+2)
    for n in range(1,Nt):
        for i in range(1,N+1):
            s[i] = 2 * (1 - l**2) * w[i] + l**2 * (w[i+1] + w[i-1]) - u[i]
        u = w.copy()
        w = s.copy()
        plt.plot(x,s,'--')
    return s



L = 1        # Longitud del dominio
N = 9        # Número de incógnitas internas 9
Tmax = 1.0   # Tiempo máximo de simulación
ht = 0.05    # Paso de tiempo .005
alpha = 2    # Dato físico

h = L / (N+1)        # Tamaño de la malla
Nt = int(Tmax / ht)  # Numero total de pasos
lamb = alpha * ht / h     # Parámetro lambda
Tmax = Nt * ht       # Tiempo total de simulación

print(" N = %d, h = %g, lambda = %g, Nt = %g, ht = %g, Tmax = %g" % (N,h,lamb,Nt,ht,Tmax))

x = np.linspace(0,L,N+2)  # Definición de la malla
u = f(x) # Condición inicial

#plt.plot(x, u,'^r-', label = "Cond. Inicial")  
#plt.plot(x, solExacta(x,Tmax),'^b--', label = "Sol. Exacta")  

w = condicionesIniciales(lamb, ht, u, x)

#plt.plot(x,u,'o-', label="Inicial 1")
#plt.plot(x,w,'s--', label="Inicial 2")

s = solver(u, w, N, x, Nt, lamb)

plt.plot(x, solExacta(x,Tmax),'^b-', label = "Exacta")  
plt.plot(x,calcError(s,solExacta(x,Tmax)),'xk-', label = "Error")

plt.grid()
plt.legend()
#plt.savefig("wave01.pdf")
plt.show()