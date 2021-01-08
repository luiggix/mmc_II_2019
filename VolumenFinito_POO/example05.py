#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 15:32:38 2018

@author: luiggi

Problema 3 de la tarea

"""

import FiniteVolumeMethod as fvm
import viscoflow as vfl
import numpy as np
from scipy import special

def analyticSol(x, u, t, Gamma):
 	divisor = 2 * np.sqrt(Gamma * t)
 	sol = 0.5 * (special.erfc((x - u * t)/ divisor) + 
 		np.exp(u * x) * np.exp(-Gamma) * special.erfc((x + u * t)/divisor))
 	return sol

title_graf = '$\partial \phi / \partial t + \partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$'
L = 2.5 # m
rho = 1.0 # kg/m^3
u = 1.0 # m/s
Gamma = 0.001 # kg / m.s
phi0 = 1 #
phiL = 0 #
N = 200 # Número de nodos
delta_t = 0.002 # Paso de tiempo
steps = 500
#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = L)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
x = malla.createMesh()    # Vector de coordenadas del dominio
#
# Se construye el arreglo donde se guardará la solución
#
phi = np.zeros(nvx) # El arreglo contiene ceros
phi[0]  = phi0       # Condición de frontera izquierda
phi[-1] = phiL       # Condición de frontera derecha
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = L,
              Densidad = rho,
              Velocidad = u,
              Coef_Diff = Gamma,
              Prop_0 = phi0,
              Prop_L = phiL,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
# Se aloja memoria para los coeficientes
#
coef = fvm.Coefficients()
coef.alloc(nvx)
#
#  Definimos el tipo de coeficientes de la ecuación
#
dif = fvm.Diffusion1D(nvx, Gamma = Gamma, dx = delta)
adv = fvm.Advection1D(nvx, rho = rho, dx = delta)
tem = fvm.Temporal1D(nvx, rho = rho, dx = delta, dt = delta_t)

exac = analyticSol(x, u, delta_t * steps, Gamma)
vfl.grafica(x,exac,label='Sol. Exac',kind='b-')

for i in range(1,steps+1):
#
# Iniciamos con coeficientes cero
#
    print('Time step = {}'.format(i * delta_t), sep = '\t')
    coef.cleanCoefficients()
#
#  Calculamos los coeficientes de FVM de la Difusión
#
    dif.calcCoef()
#	dif.printCoefficients()
#	print('.'+'-'*70+'.')
#
#  Calculamos los coeficientes de FVM de la Advección
#
    adv.setU(u)
    adv.calcCoef()
#	adv.printCoefficients()
#	print('u = {}'.format(adv.u()))
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
    coef.bcDirichlet('LEFT_WALL', phi0)   # Se actualizan los coeficientes
    coef.bcDirichlet('RIGHT_WALL', phiL) # de acuerdo a las cond. de frontera
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
    if (i % 100 == 0):
        etiqueta = 'Step = {}'.format(i*delta_t)
        vfl.grafica(x,phi,title=title_graf, label=etiqueta)

vfl.show('example05.pdf')

