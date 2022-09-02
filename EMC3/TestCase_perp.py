# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 11:03:28 2016

@author: j.cosfeld
"""

from numpy import linspace, mgrid, pi, sqrt, savetxt, diff, zeros, ones 
import numpy as np # NumPy

from os import path, makedirs


"""Class that represents test cases for EMC3-EIRENE.

This class contains methods to create input files for EMC3-EIRENE that
can be used for testing the code.

"""

class TestCase2( object ):

  def __init__( self, path = '.' ):
    """Constructor

    Keyword arguments:
    path -- Path to the location where the test case will be created.
            (str, default '.')

    """
    self.path = path.rstrip('/')

  def Get1DGrid( self, NRadial = 101, ConnectionLength = 4000, Area = 10, 
                 R0 = 1e5, Z0 = 0, B0 = 1 ):
    """Creates the input files for an one dimensional grid.

    This methode creates the input files for an one dimensional grid 
    based on the given input parameters. The files are stored in a newly 
    created folder 'geometry'. The grid is stored in the file 
    'GRID_3D_DATA', the magnetic field data are stored in the file 
    'BFIELD_STRENGTH', and the EMC3 input file is named 'input.geo'.

    Keyword arguments:
    NToroidal -- Number of toroidal grid points. (int, default 101)
    ConnectionLength -- Target-target connection lenght of the grid in cm. 
                        (float, default 4000)
    Area -- Cross-section area of the grid in cm^2. (float, default 10)
    R0 -- Radial coordinate of the centre of the grid in cm. (float, 
          default 1e5)
    Z0 -- Horizontal coordinate of the centre of the grid in cm. (float, 
          default 0)
    B0 -- Magnetic field strength at the location of the grid in T.
          (float, default 1)

    """
    
    NToroidal = 4
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
    
    R[:NRadial]=R1
    R[NRadial:]=R1
     
    inputGeo = \
    "* input file for EMC3 code (1D parallel test)\n" + \
    "*** 1. fein grid mesh representing magnetic coordinates\n" + \
    "1                     !zone number\n" + \
    "%d %d %d     !core    !R-,P-,T_surfaces\n" % ( NRadial, NPoloidal, 
                                                    NToroidal ) + \
    "*** 2. non default non transparent, plate surface\n" + \
    "*** 2.1 non default surface\n" + \
    "* radial\n" + \
    "0                     !number of surfaces\n" + \
    "* poloidal\n" + \
    "* Type of surface=1:  periodic\n" + \
    "*                 2:  up/down symmetric\n" + \
    "*                 3:  Mapping\n" + \
    "2                     !number of surfaces\n" + \
    "0 0 1                 !IP  IZONE  Type of surface\n" + \
    "0 %d 0 %d             !IR1-2, IT1-2\n" % ( NRadial-2, NToroidal-2 ) + \
    "%d 0 1                !IP  IZONE  Type of surface\n" % ( NPoloidal-1 ) + \
    "0 %d 0 %d             !IR1-2, IT1-2\n" % ( NRadial-2, NToroidal-2 ) + \
    "* toroidal\n" + \
    "2                     !number of surfaces\n" + \
    "0 0 1                 !IT  IZONE  Type of surface\n" + \
    "0 %d 0 %d             !IR1-2, IP1-2\n" % ( NRadial-2, NPoloidal-2 ) + \
    "%d 0 1                 !IT  IZONE  Type of surface\n" % ( NToroidal-1 ) + \
    "0 %d 0 %d             !IR1-2, IP1-2\n" % ( NRadial-2, NPoloidal-2 ) + \
    "*** 2.2 non transparent (Boundary condition must be defined)\n" + \
    "* radial\n" + \
    "2                     !number of surfaces\n" + \
    "0 0 1                 !IR  IZONE  ISIDE\n" + \
    "0 %d 0 %d             !IP1-2, IT1-2\n" % ( NPoloidal-2, NToroidal-2 ) + \
    "%d 0 -1               !IR  IZONE  ISIDE\n" % ( NRadial-1 ) + \
    "0 %d 0 %d             !IP1-2, IT1-2\n" % ( NPoloidal-2, NToroidal-2 ) + \
    "* POLOIDAL\n" + \
    "0                     !number of surfaces\n" + \
    "* TOROIDAL\n" + \
    "0                     !number of surfaces\n" + \
    "*** 2.3 plate surface (Bohm Boundary condition)\n" + \
    "* radial\n" + \
    "0                     !number of surfaces\n" + \
    "* POLOIDAL\n" + \
    "0                     !number of surfaces\n" + \
    "* TOROIDAL\n" + \
    "1                     !number of surfaces\n" + \
    "%d 0 -1               !IT  IZONE  ISIDE\n" % ( NToroidal-1) + \
    "0 %d 0 %d             !IR1-2, IP1-2\n" % ( NRadial-2, NPoloidal-2 ) + \
    "*** 3. physical cell\n" + \
    "*<0: user defines cell\n" + \
    "* 0: default;\n" + \
    "*>0: read from file\n" + \
    "1 %d\n" % (NToroidal-1) + \
    "* check cell ?\n" + \
    "T"
    
    if not path.isdir( '%s/geometry' % self.path ):
      makedirs( '%s/geometry' % self.path )

    fileInput = open( '%s/geometry/input.geo' % self.path, 'w' )
    fileInput.write( inputGeo )
    fileInput.close()
    
    savetxt( '%s/geometry/x.dat' % self.path, x_par[:-1] + 0.5 * diff(x_par) )
    
    fileGrid = open( '%s/geometry/GRID_3D_DATA' % self.path, 'w' )
    fileField = open( '%s/geometry/BFIELD_STRENGTH' % self.path, 'w' )
    fileGrid.write( '%d %d %d %f %f\n' % ( NRadial, NPoloidal, NToroidal, R0, 
                                           Z0 ) )
    for k in xrange( NToroidal ):
      fileGrid.write( '%f\n' % phi[k] )
      fileGrid.write( NRadial * NPoloidal * ' %f' % tuple( R.flatten() ) )
      fileGrid.write( '\n' )
      fileGrid.write( NRadial * NPoloidal * ' %f' % tuple( Z.flatten() ) )
      fileGrid.write( '\n' )
      fileField.write( NRadial * NPoloidal * ' %f' % tuple( NRadial * 
                                                            NPoloidal * [B0] ) )
      fileField.write( '\n' )

    fileGrid.close()
    fileField.close()    
    
    fileTarget = open( '%s/geometry/TARGET_1D' % self.path, 'w' )
    fileTarget.write( ' Target for 1D test\n' )
    fileTarget.write( '2 2 1 %f %f\n' % ( R0, Z0 ) )
    fileTarget.write( '%f\n' % ( phi[-1] * (1-1e-7) ) ) # Reduce phi to ensure the plate is inside the grid
    fileTarget.write( '-10.    -10.\n' )
    fileTarget.write( ' 10.    -10.\n' )
    fileTarget.write( '%f\n' % ( phi[-1] * (1-1e-7) ) )
    fileTarget.write( '-10.     10.\n' )
    fileTarget.write( ' 10.     10.\n' )
    fileTarget.close()
    
    fileTarget = open( '%s/energy/EMC3_OUTPUT/E_SOURCE' % self.path, 'w' )

    #(NToroidal-1)*100 Sea entries are needed for (NToroidal-1) Gridboxes   
    
    linelength=10
    filelength=30
    Sea_entry=5.E-1
    Sea_file=zeros((10,linelength))
    Sea_file[:,:]=Sea_entry
    
    text_file = open('%s/energy/EMC3_OUTPUT/E_SOURCE' % self.path, "w")
    for i in range(1,filelength+1):
        text_file.write("5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\n")
        
    text_file.close()
    
    text_file = open('%s/energy/EMC3_OUTPUT/I_SOURCE' % self.path, "w")
    for i in range(1,filelength+1):
        text_file.write("5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\t5.E-1\n")
        
    text_file.close()


    fileTarget.close()
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    