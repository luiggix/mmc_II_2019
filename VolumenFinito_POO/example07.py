#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: luiggi

Problema 3 de la tarea

"""

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt

def solucionExacta(PA, PB, x, Lx, Gamma, t):
    N = 100
    p = np.zeros(len(x))
    PI = np.pi
    for i,xi in enumerate(x):
        xti = 0
        for n in range(1,N):
            xti += np.exp(-n*n*PI*PI/(Lx*Lx) * Gamma * t) * \
                   np.sin(n*PI*xi/Lx) / n
        p[i] = PA + (PB - PA) * (xi / Lx + 2 * xti / PI)
    return p

title_graf = '$\dfrac{\partial p}{\partial t} = \dfrac{\partial}{\partial x} (\Gamma \dfrac{\partial p}{\partial x})$'

# Propiedades petrofísicas
Lx = 100
poro = 0.2
mu = 1.0
k = 1.0
cT = 1e-4
PA = 2.0
PB = 1.0
Tmax = 0.8
delta_t = 0.0001 # Paso de tiempo
steps = int(Tmax/delta_t)
# Parámetros numéricos
m = 1
Nx = 10 * 2**m + 1
x = np.linspace(0,100,Nx)
dx = Lx / (Nx-1)
t = np.array([0.0005,0.0025,0.0075, 0.0150, 0.0375, 0.8000])
Gamma = k / (poro * mu * cT)

# L = 2.5 # m
# rho = 1.0 # kg/m^3
# u = 1.0 # m/s
# Gamma = 0.001 # kg / m.s
# phi0 = 1 #
# phiL = 0 #
# N = 200 # Número de nodos
# steps = 500
#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = Nx, length = Lx)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
x = malla.createMesh()    # Vector de coordenadas del dominio
xe = np.linspace(x[0], x[-1], 100)
#
# Se construye el arreglo donde se guardará la solución
#
phi = np.ones(nvx) # El arreglo contiene unos
phi[0]  = PA       # Condición de frontera izquierda
phi[-1] = PB       # Condición de frontera derecha
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = Lx,
              Porosidad = poro,
              Viscosidad = mu,
              Compresibilidad = cT,
              Coef_Diff = Gamma,
              Prop_0 = PA,
              Prop_L = PB,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta,
              Tmax = Tmax,
              Steps = steps,
              ht = delta_t)
#
# Se aloja memoria para los coeficientes
#
coef = fvm.Coefficients()
coef.alloc(nvx)
#
#  Definimos el tipo de coeficientes de la ecuación
#
dif = fvm.Diffusion1D(nvx, Gamma = Gamma, dx = delta)
tem = fvm.Temporal1D(nvx, rho = 1, dx = delta, dt = delta_t)

#for nt in t:
#    plt.plot(x,solucionExacta(PA,PB,x,Lx,Gamma,nt), label='t = %s'%nt)
#vfl.grafica(x,exac,label='Sol. Exac',kind='b-')

for i in range(1,steps+1):
#
# Iniciamos con coeficientes cero
#
    time_k = i * delta_t
#    print('i = {}, Time step = {}'.format(i, time_k), sep = '\t')
    coef.cleanCoefficients()
#
#  Calculamos los coeficientes de FVM de la Difusión
#
    dif.calcCoef()
#	dif.printCoefficients()
#	print('.'+'-'*70+'.')
#
#  Calculamos los coeficientes de FVM de la parte temporal
#
    tem.calcCoef(phi)
#	tem.printCoefficients()
#	print('dt = {}'.format(tem.deltaT()))
#	print('.'+'-'*70+'.') 
#
# Se aplican las condiciones de frontera
#
    coef.bcDirichlet('LEFT_WALL', PA)   # Se actualizan los coeficientes
    coef.bcDirichlet('RIGHT_WALL', PB) # de acuerdo a las cond. de frontera
#	coef.printCoefficients()
#	print('u = {}'.format(adv.u()))
#	print('.'+'-'*70+'.')
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
    Su = coef.Su()  # Vector del lado derecho
    A = fvm.Matrix(malla.volumes())  # Matriz del sistema
    A.build(coef) # Construcción de la matriz en la memoria
#print('A = ', A.mat(),
#      'b = {}'.format(Su[1:-1]), sep='\n')
#print('.'+'-'*70+'.')
#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
    phi[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
    
#	print('Solución = {}'.format(phi))
#	print('.'+'-'*70+'.')
#
# Usamos Viscoflow para graficar
# VISCOFLOW : Visualization of Complex Flows
#
    if (time_k in t or i == 75 or i == 150):
        print(i, time_k)
        plt.scatter(x,phi, label=' t = {:5.4f}'.format(time_k),
                    zorder=5)
        if time_k == 0.0005:
            plt.plot(xe,solucionExacta(PA,PB,xe,Lx,Gamma, time_k), 
                     'k-', lw=1, alpha=1.0, label='Exacta')
        else: 
            plt.plot(xe,solucionExacta(PA,PB,xe,Lx,Gamma, time_k), 
                     'k-', lw=1, alpha=1.0)
                 

plt.title(title_graf, fontsize=16, color='blue')
plt.xlabel('$x$ [cm]', fontsize=14)
plt.ylabel('$p$ [atm]', fontsize=14)
plt.legend()
plt.grid()
plt.savefig('example07.pdf')
plt.show()

