from numpy import array, double 
from fortranformat import FortranRecordReader
from Grid import Grid
from Analysis import getAreaPolygon

class Depo_modded( object ):
  """Class that represents the deposition output of the EMC3-EIRENE code.
  
  This class can read the deposition output of streaming and energy.
  
  """
  
  def __init__( self, name, path = 'EMC3_OUTPUT', grid = None ):
    """Constructor

    Keyword arguments:
    name -- Name of the deposition to read. (str in {'PARTICLE', 'ENERGY'})
    path -- Path to the file. (str, default 'EMC3_OUTPUT')
    grid -- Grid. If None trying to initialize it. (Grid object, default 
            None)

    """
    self.name = name
    self.path = path.rstrip('/')

    if grid is None:
      self.grid = Grid()
    else:
      self.grid = grid

    self.__area = double( [] )
    self.__heatflux = double( [] )

    Formats = { 
      'PARTICLE': ( '(I6,3I3,3I5,E14.6)', '(I6,3I3,3I5,E12.4)' ),
      'ENERGY': ( '(I6,3I3,3I5,1pE12.4)', '(I6,3I3,3I5,1pE12.4)' ) }
 
    if self.name == 'PARTICLE':
      self.fluxUnit = 'A'
    elif self.name == 'ENERGY':
      self.fluxUnit = 'W'

    file = open( '%s/%s_DEPO' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content.reverse() # to allow a much more efficient pop from the end of the list

    Fformat = FortranRecordReader( Formats[name][0] )

    line = content.pop()
    if name == 'PARTICLE':
      Ns = line.split()
      self.NTargetSurfaces = int(Ns[0])
      self.NMappingSurfaces = int(Ns[1])

    isTarget = True
    self.isTarget = []
    self.index = []
    self.surfaceType = []
    self.index_toro = []
    self.index_polo = []
    self.index_radi = []
    self.phi = []
    self.R = []
    self.Z = []
    self.flux = []

    while len( content ) > 0:
      line = content.pop()
      if line.find( 'MAPPING' ) != -1:
        Fformat = FortranRecordReader( Formats[name][1] )
        isTarget = False
        continue
      self.isTarget.append( isTarget )
      values = Fformat.read( line )
      self.index.append( values[0] )
      self.surfaceType.append( values[2] )
      
      zone = self.grid.zones[ values[3] ]
      self.index_toro.append( [  values[6] ])
      self.index_polo.append( [  values[5] ])
      self.index_radi.append( [  values[4] ])
      self.R.append( [ zone.R[ values[6], values[5], values[4] ], \
                       zone.R[ values[6], values[5] + 1, values[4] ], \
                       zone.R[ values[6], values[5] + 1, values[4] + 1 ], \
                       zone.R[ values[6], values[5], values[4] + 1 ] ] )
      self.Z.append( [ zone.Z[ values[6], values[5], values[4] ], \
                       zone.Z[ values[6], values[5] + 1, values[4] ], \
                       zone.Z[ values[6], values[5] + 1, values[4] + 1 ], \
                       zone.Z[ values[6], values[5], values[4] + 1 ] ] )
      self.phi.append( zone.phi[ values[6] ] )
      self.flux.append( values[7] )

    self.isTarget = array( self.isTarget )
    self.index = array( self.index )
    self.surfaceType = array( self.surfaceType )
    self.index_toro=array(self.index_toro)
    self.index_polo=array(self.index_polo)
    self.index_radi=array(self.index_radi)
    self.phi = array( self.phi )
    self.flux = array( self.flux )
    self.R = array( self.R )
    self.Z = array( self.Z )
    
  def __calculateArea( self ):
    """ Calculates the area of each cell

    This method calculates the area of each cell using Gauss's area 
    formula.

    """
    area = []
    for R, Z in zip( self.R, self.Z ):
      area.append( getAreaPolygon( R, Z ) )

    self.__area = double( area )
    
  @property
  def heatflux( self ):
    """(Calculates and) returns the heat flux at each cell (W/cm**2)"""
    if self.__heatflux.size == 0:
      self.__heatflux = self.flux / self.area
    return self.__heatflux 
    
  @property
  def area( self ):
    """(Calculates and) returns the area of each cell (cm**2)"""
    if self.__area.size == 0:
      self.__calculateArea()
    return self.__area    
    
    
    
    
    
    
    
    
    
    
    
