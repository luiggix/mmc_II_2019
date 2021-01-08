#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 17:45:32 2021

@author: luiggi
"""
import mesh

mesh1D = mesh.Malla(2.0, 21)

hx, _ = mesh1D.calculaDelta()

print('Conducción')
print('hx = {}'.format(hx))