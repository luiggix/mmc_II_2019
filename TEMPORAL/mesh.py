#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 17:00:46 2021

@author: luiggi
"""
import numpy as np

class Malla():
    
    def __init__(self, Lx, Nx, Ly = 1.0, Ny = 1):
        self.Nx = Nx
        self.Ny = Ny
        self.Lx = Lx
        self.Ly = Ly
        
#        print('Nx = {}, Ny = {}'.format(self.Nx, self.Ny))
#        print('Lx = {}, Ly = {}'.format(self.Lx, self.Ly))

    def calculaDelta(self):
        self.hx = self.Lx / (self.Nx+1)
        self.hy = self.Ly / (self.Ny+1)
        return (self.hx, self.hy)

if __name__ == '__main__':
        
    arreglo = np.array(10)

    malla = Malla(1.0, 5, 2.0, 8)

    hx, hy = malla.calculaDelta()

    print('hx = {}, hy = {}'.format(hx, hy))

    malla1D = Malla(2.0, 40)
    hx, _ = malla1D.calculaDelta()
    print('hx = {}'.format(hx))


