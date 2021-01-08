#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: luiggi

Example 4.1 from Malalasekera Book
----------------------------------
Consider the problem of source-free heat conduction in an insulated rod
whose ends are maintained at constant temperatures of 100°C and 500°C
respectively. The one-dimensional problem sketched in Figure 4.3 is governed
by
$
\dfrac{d}{dx} \left( k \dfrac{d T}{dx} \right) = 0
$

___________________________________________________________

 |<-------------- 0.5 m ------------->|
 A                                    B
 |====================================|
 |                                    |
 |====================================|
T_A = 100                            T_B = 500
___________________________________________________________
Figure 4.3

Calculate the steady state temperature distribution in the rod. Thermal conductivity
k equals 1000 W/m.K, cross-sectional area A is 10 × 10−3 m2.

Analytical solution:
    
T = 800 * x + 100
"""

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt

longitud = 0.5 # meters
TA = 100 # °C 
TB = 500 # °C 
k  = 1000 # W/m.K
N  = 6 # Número de nodos
#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = longitud)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = longitud,
              Temperatura_A = TA,
              Temperatura_B = TB,
              Conductividad = k,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
#  Creamos los coeficientes de FVM
#
df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta)
df1.alloc(nvx) # Se aloja memoria para los coeficientes
df1.calcCoef() # Se calculan los coeficientes
#
# Se construye el arreglo donde se guardará la solución
#
T = np.zeros(nvx) # El arreglo contiene ceros
T[0]  = TA        # Condición de frontera izquierda
T[-1] = TB        # Condición de frontera derecha
df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
df1.bcDirichlet('RIGHT_WALL', T[-1]) # de acuerdo a las cond. de frontera
print('aW = {}'.format(df1.aW()), 
      'aE = {}'.format(df1.aE()), 
      'Su = {}'.format(df1.Su()), 
      'aP = {}'.format(df1.aP()), sep='\n')
print('.'+'-'*70+'.')
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
Su = df1.Su()  # Vector del lado derecho
A = fvm.Matrix(malla.volumes())  # Matriz del sistema
A.build(df1) # Construcción de la matriz en la memoria
print('A = ', A.mat(),
      'b = {}'.format(Su[1:-1]), sep='\n')
print('.'+'-'*70+'.')
#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
print('Solución = {}'.format(T))
print('.'+'-'*70+'.')
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
#
# Calculamos la solución analítica
#
Ta = 800 * x + 100
#
#  Se grafica la solución
#
x *= 100 # Transformación a [cm]
plt.plot(x,Ta, '-', label = 'Sol. analítica') # Sol. analítica
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Solución de $k (\partial^2 T/\partial x^2) = 0$ con FVM')
plt.xlabel('$x$ [cm]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('example01.pdf')
plt.show()
