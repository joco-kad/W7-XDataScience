# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:02:43 2016

@author: j.cosfeld
"""
"""

This is only a playground for the TestCase!

"""

import numpy as np # NumPy


NRadial = 101
ConnectionLength = 4000 
Area = 10 
R0 = 1e5
Z0 = 0
B0 = 1

NToroidal = 2
NPoloidal = 2

x_par = linspace( 0, 0.5 * ConnectionLength, NRadial )
phi_par = x_par * 180.0 / pi / R0

x_par_mid_index=len(x_par)/2

phi_perp=zeros(4)
phi_perp[0]=phi_par[x_par_mid_index]
phi_perp[1]=phi_par[x_par_mid_index+1]
phi_perp[2]=phi_par[x_par_mid_index+2]
phi_perp[3]=phi_par[x_par_mid_index+3]

phi=phi_perp

dR = sqrt(Area)
dZ = sqrt(Area)

R1 = linspace(-dR*(NRadial-1)/2, dR*(NRadial-1)/2, NRadial)

Z1= -1*ones(NRadial)
Z1=Z1*dZ*(NPoloidal-1)/2

Z=ones(2*NRadial)
R=ones(2*NRadial)

Z[:NRadial]=Z1
Z[NRadial:]=-Z1



#R, Z = mgrid[ -dR*(NRadial-1)/2 : dR*(NRadial-1)/2 : NRadial*1j,
#        -dZ*(NPoloidal-1)/2 : dZ*(NPoloidal-1)/2 : NPoloidal*1j ]
        
        
#for k in xrange( NToroidal ):
#    
#    print(phi[k])
