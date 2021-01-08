#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 13:21:50 2018

@author: luiggi
"""
import numpy as np
from pandas import DataFrame
from Mesh import Mesh
from Coefficients import Coefficients
from Diffusion import Diffusion1D
from Advection import Advection1D
from Temporal import Temporal1D
from Matrix import Matrix
import time

def crono(f):
 	"""
 	Regresa el tiempo que toma en ejecutarse la funcion.
 	"""
 	def eTime(A,b):
 		t1 = time.time()
 		f(A,b)
 		t2 = time.time()
 		return 'Elapsed time: ' + str((t2 - t1)) + "\n"
 	return eTime

def decorate(f):
    def nicePrint(**kargs):
        line = '-' * 70
        print('.'+ line + '.')
        print('|{:^70}|'.format('NoNacos : Numerical Objects for Natural Convection Systems'))
        print('.'+ line + '.')
        print('|{:^70}|'.format(' Ver. 0.1, Author LMCS, 2018, [GNU GPL License V3]'))
        print('.'+ line + '.')
        f(**kargs)
        print('.'+ line + '.')
    return nicePrint
 
@decorate
def printData(**kargs):
    for (key,value) in kargs.items():
        if (type(value) == str):
            print('|{:^70}|'.format('{0:>15s} = {1:11s}'.format(key, value)))
        elif (type(value) == int):
            print('|{:^70}|'.format('{0:>15s} = {1:<11d}'.format(key, value)))            
        else:
            print('|{:^70}|'.format('{0:>15s} = {1:10.5e}'.format(key, value)))

def printFrame(d):
    # Calculo el error porcentual y agrego al DataFrame
    # una columna con esos datos llamada 'Error %'
    d['Error %'] = d['Error'] / d['Analytic'] 
    print(DataFrame(d))
    print('.'+ '-'*70 + '.')

def calcError(phiA, phiN):
    return np.absolute(phiA - phiN)
        
if __name__ == '__main__':
 
    Coefficients.alloc(5)
    m = Mesh(nodes = 5)
    d = Diffusion1D(m.volumes())
    ma = Matrix(m.volumes())
    a = Advection1D(m.volumes())
    t = Temporal1D(m.volumes()) 

    print(m.delta(), d.aP(), a.aP(), t.aP(), ma.mat(), sep='\n')

    printData(Name='Laplace', nvx = 5, nx = 6, longitud = 1.3)
    
