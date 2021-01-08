#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: luiggi

Example 4.3 from Malalasekera Book
----------------------------------
In the final worked example of this chapter we discuss the cooling of a
circular fin by means of convective heat transfer along its length. Convection
gives rise to a temperature-dependent heat loss or sink term in the governing
equation. Shown in Figure 4.9 is a cylindrical fin with uniform crosssectional
area A. The base is at a temperature of 100°C (T_A) and the end is
insulated. The fin is exposed to an ambient temperature of 20°C. Onedimensional
heat transfer in this situation is governed by

$
\dfrac{d}{dx} \left( k A \dfrac{d T}{dx} \right) - hP(T - T∞) = 0
$

where h is the convective heat transfer coefficient, P the perimeter, k the
thermal conductivity of the material and $T_\infty$ the ambient temperature.
Calculate the temperature distribution along the fin and compare the results
with the analytical solution given by

$\frac{T - T∞}{T_A - T∞} = \frac{\cosh[n(L-x)]}{\cosh(nL)}$

where $n^2 = hP/(kA)$, $L$ is the length of the fin and $x$ the distance along the
fin. Data: $L = 1$ m, $hP/(kA) = 25/m^2$ (note that $kA$ is constant).

The governing equation in the example contains a sink term, −hP(T − T∞),
the convective heat loss, which is a function of the local temperature T.

___________________________________________________________

 |<-------------- 1.0 m ------------->|
 A                                    
 |
 |
 |------------------------------------|
 |------------------------------------| Flujo igual a zero
 |        Temp. Ambiente (T∞)                           
 |                                    |
                                      |
T_A = 100                  \frac{\partial T}{\partial dx} = q = 0
___________________________________________________________
Figure 4.9

When $kA$ = constant, the governing equation can be written as

$
\dfrac{d}{dx} \left( \dfrac{d T}{dx} \right) - n^2(T - T∞) = 0
$

"""

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt

def analyticSol(x,n):
    return (TA - Tambiente) * np.cosh(n * (longitud - x)) / np.cosh(n * longitud) + Tambiente

longitud = 1.0 # metros
Tambiente = 20  # °C 
TA = 100  # °C 
n2 = 25 # /m^2
fluxB = 0 # Flujo igual a cero
N = 6 # Número de nodos
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
              Flujo_B = fluxB,
              n2 = n2,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
#  Creamos los coeficientes de FVM
#
df1 = fvm.Diffusion1D(nvx, Gamma = 1, dx = delta)
df1.alloc(nvx) # Se aloja memoria para los coeficientes
df1.calcCoef() # Se calculan los coeficientes
df1.setSu(n2 * Tambiente)  # Se agrega la fuente uniforme
df1.setSp(-n2)
#print(df1.aP(),df1.aW(), df1.aE(), df1.Su(), sep = '\n')
#
# Se construye el arreglo donde se guardará la solución
#
T = np.zeros(nvx) # El arreglo contiene ceros
T[0]  = TA        # Condición de frontera izquierda
df1.bcDirichlet('LEFT_WALL', T[0])   # Se actualizan los coeficientes
df1.bcNeumman('RIGHT_WALL', fluxB) # de acuerdo a las cond. de frontera
#
# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
Su = df1.Su()  # Vector del lado derecho
A = fvm.Matrix(malla.volumes())  # Matriz del sistema
A.build(df1) # Construcción de la matriz en la memoria
#
# Se resuelve el sistema usando un algoritmo del módulo linalg
#
T[1:-1] = np.linalg.solve(A.mat(),Su[1:-1])
T[-1] = T[-2] # Condición de frontera tipo Neumman
#
# Se construye un vector de coordenadas del dominio
#
x = malla.createMesh()
#
# Calculamos la solución exacta y el error
#
n = np.sqrt(n2)
Ta = analyticSol(x,n)
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
Ta = analyticSol(x1,n)
#
#  Se grafica la solución
#
plt.plot(x1,Ta, '-', label = 'Sol. analítica') 
plt.plot(x,T,'o', label = 'Sol. FVM')
plt.title('Solución de $\partial^2 T/\partial x^2 - hP(T-T_\infty) = 0$ con FVM')
plt.ylim(10,110)
plt.xlabel('$x$ [m]')
plt.ylabel('Temperatura [$^o$C]')
plt.grid()
plt.legend()
plt.savefig('example03.pdf')
plt.show()
