#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: luiggi

Example 4.2 from Malalasekera Book
----------------------------------
Now we discuss a problem that includes sources other than those arising
from boundary conditions. Figure 4.6 shows a large plate of thickness
L = 2 cm with constant thermal conductivity k = 0.5 W/m.K and uniform
heat generation q = 1000 kW/m3. The faces A and B are at temperatures
of 100°C and 200°C respectively. Assuming that the dimensions in the y- and
z-directions are so large that temperature gradients are significant in the 
x-direction only, calculate the steady state temperature distribution. Compare
the numerical result with the analytical solution. The governing equation is

$
\dfrac{d}{dx} \left( k \dfrac{d T}{dx} \right) + q = 0
$

___________________________________________________________

 |<-------------- 2.0 cm ------------>|
 A                                    B
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
 |                                    |
 |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|

T_A = 100                            T_B = 200
___________________________________________________________
Figure 4.6

Analytical solution:
    
$T = \left( \frac{T_B - T_A}{L} + \frac{q}{2k}(L - x) \right) x + T_A$

"""

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt

def analyticSol(x):
    return ((TB - TA) / longitud + q * (longitud - x) / (2 * k) ) * x + TA
    
longitud = 0.02 # meters
TA = 100  # °C 
TB = 200  # °C 
k  = 0.5  # W/m.K
q  = 1e+06 # 1e+06 W/m^3 = 1000 kW/m^3 Fuente uniforme
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
              Fuente = q,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
#  Creamos los coeficientes de FVM
#
df1 = fvm.Diffusion1D(nvx, Gamma = k, dx = delta)
df1.alloc(nvx) # Se aloja memoria para los coeficientes
df1.calcCoef() # Se calculan los coeficientes
df1.setSu(q)  # Se agrega la fuente uniforme
#
# Se construye el arreglo donde se guardará la solución
#
T = np.zeros(nvx) # El arreglo contiene ceros
T[0]  = TA        # Condición de frontera izquierda
T[-1] = TB        # Condición de frontera derecha
df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
df1.bcDirichlet('RIGHT_WALL', T[-1]) # de acuerdo a las cond. de frontera
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
Su = df1.Su()  # Vector del lado derecho
A = fvm.Matrix(malla.volumes())  # Matriz del sistema
A.build(df1) # Construcción de la matriz en la memoria
#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
print(A.mat())
print(Su[1:-1])
T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
#
# Calculamos la solución exacta y el error
#
Ta = analyticSol(x)
error = fvm.calcError(T, Ta)
datos = {'x(m)': x,
         'T(x)': T,
         'Analytic': Ta,
         'Error': error}
fvm.printFrame(datos)
print('||Error|| = ', np.linalg.norm(error))
print('.'+ '-'*70 + '.')
#
# Calculamos la solución exacta en una malla más fina para graficar
#
x1 = np.linspace(0,longitud,100)
Ta = analyticSol(x1)
#
#  Se grafica la solución
#
x *= 100 # Transformación a [cm]
x1 *= 100
plt.plot(x1,Ta, '-', label = 'Sol. analítica') 
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Solución de $k (\partial^2 T/\partial x^2)+q = 0$ con FVM')
plt.ylim(75,300)
plt.xlabel('$x$ [cm]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('example02.pdf')
plt.show()
