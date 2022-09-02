# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 12:23:10 2016

@author: j.cosfeld
"""


class TestCase( object ):

  def __init__( self, path = '.' ):
    """Constructor

    Keyword arguments:
    path -- Path to the location where the test case will be created.
            (str, default '.')

    """
    self.path = path.rstrip('/')




from numpy import linspace, mgrid, pi, sqrt, savetxt, diff
from os import path, makedirs


NRadial = 101
ConnectionLength = 4000
Area = 10
R0 = 1e5
Z0 = 0
B0 = 1
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

NToroidal = 2
NPoloidal = 2

x = linspace( 0, 0.5 * ConnectionLength, NToroidal )
phi = x * 180.0 / pi / R0

dR = sqrt(Area)
dZ = sqrt(Area)

R, Z = mgrid[ -dR*(NToroidal-1)/2 : dR*(NToroidal-1)/2 : NToroidal*1j,
              -dZ*(NPoloidal-1)/2 : dZ*(NPoloidal-1)/2 : NPoloidal*1j ]
 
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
"1                     !number of surfaces\n" + \
"0 0 3                 !IT  IZONE  Type of surface\n" + \
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

savetxt( '%s/geometry/x.dat' % self.path, x[:-1] + 0.5 * diff(x) )

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

fileTarget = open( '/geometry/TARGET_1D' , 'w' )
fileTarget.write( ' Target for 1D test\n' )
fileTarget.write( '2 2 1 %f %f\n' % ( R0, Z0 ) )
fileTarget.write( '%f\n' % ( phi[-1] * (1-1e-7) ) ) # Reduce phi to ensure the plate is inside the grid
fileTarget.write( '-10.    -10.\n' )
fileTarget.write( ' 10.    -10.\n' )
fileTarget.write( '%f\n' % ( phi[-1] * (1-1e-7) ) )
fileTarget.write( '-10.     10.\n' )
fileTarget.write( ' 10.     10.\n' )
fileTarget.close()
