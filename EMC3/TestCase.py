from numpy import linspace, mgrid, pi, sqrt, savetxt, diff
from os import path, makedirs

class TestCase( object ):
  """Class that represents test cases for EMC3-EIRENE.
  
  This class contains methods to create input files for EMC3-EIRENE that
  can be used for testing the code.
  
  """

  def __init__( self, path = '.' ):
    """Constructor

    Keyword arguments:
    path -- Path to the location where the test case will be created.
            (str, default '.')

    """
    self.path = path.rstrip('/')

  def Get1DGrid_parallel( self, NToroidal = 101, ConnectionLength = 4000, Area = 10, 
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
    NRadial = 2
    NPoloidal = 2

    x = linspace( 0, 0.5 * ConnectionLength, NToroidal )
    phi = x * 180.0 / pi / R0

    dR = sqrt(Area)
    dZ = sqrt(Area)

    R, Z = mgrid[ -dR*(NRadial-1)/2 : dR*(NRadial-1)/2 : NRadial*1j,
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

  def Get1DGrid_perp( self, nomRADsurf = 101, nomPOLsurf = 2, nomTORsurf = 50, shift = 1e5, opt = 0):
    """
    create angular element by phi_begin phi_end and dphi 
    """
    phi_begin=0.572958 
    phi_end=1
    phi=np.linspace(phi_begin, phi_end, num=nomTORsurf)
    """
    create radial stringed boxes by rad_begin rad_end and dr
    """
    rad_begin=-1.113883
    rad_end=1.113883
    radial=np.linspace(rad_begin, rad_end, num=nomRADsurf)
    radial_dummy=ones(2*nomRADsurf)
    radial_dummy[0:len(radial_dummy)/2]=radial
    radial_dummy[len(radial_dummy)/2:len(radial_dummy)]=radial
    """
    create poloidal coordinates
    """
    pol_up=1.581139
    pol_bot=-1.581139
    pol=ones(nomRADsurf*nomPOLsurf)
    pol[0:len(pol)/2]=pol[0:len(pol)/2]*(-1)
    pol[:]=pol[:]*pol_up

    """
    write to file
    """        
    fh = open("./geometry/GRID_3D_DATA", "w")
    fh.write("%s " % nomRADsurf)
    fh.write("%s " % nomPOLsurf)
    fh.write("%s " % nomTORsurf)
    fh.write("%s " % shift)
    fh.write("%s \n" % opt)
    
    for i in range(0, nomTORsurf):
      fh.write("%s \n" % phi[i])
      for j in range(0, 2*nomRADsurf):
	  fh.write("%s " % radial_dummy[j])
      fh.write('\n')  
      for k in range(0, 2*nomRADsurf):
	  fh.write("%s " % pol[k])
      fh.write('\n')  
      
      
    bfield_dummy=(nomRADsurf*nomPOLsurf)    
    bfield_input=ones(bfield_dummy)
    fh2 = open("./geometry/BFIELD_STRENGTH", "w")
    for i in range(0, nomTORsurf):
	for j in range(0, 2*nomRADsurf):
	      fh2.write("%s " % bfield_input[j])
	fh2.write('\n')    


    NRadial=nomRADsurf
    NPoloidal=nomPOLsurf
    NToroidal=nomTORsurf
    
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
    

    fileInput = open( '%s/geometry/input.geo' % self.path, 'w' )
    fileInput.write( inputGeo )
    fileInput.close()


    savetxt( '%s/geometry/x.dat' % self.path, radial )
    

    
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
