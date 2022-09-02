from fortranformat import FortranRecordReader
import Surfaces

class GeometryParameters( object ):
  """Class that represents the geometry parameters of the magnetic grid.
  
  This class represents the input of EMC3 for geometry parameters. On 
  initialization it reads an existing input file.
  
  """

  def __init__( self, name = 'input.geo', path = '../../geometry' ):
    """Constructor

    Keyword arguments:
    name -- Name of the parameter input. (str, default 'input.geo')
    path -- Path to the file. (str, default '../../geometry')

    """
    self.name = name
    self.path = path.rstrip('/')
    
    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: '*' != s[0], content )

    self.NZones, \
    self.NDomains, \
    self.zoneIndices, \
    self.NSurfaces = self._ReadGridInformation( content )

    self.nonDefaultSurfaces = Surfaces.Read( content )
    
    self.nonTransparentSurfaces = Surfaces.Read( content )

    self.plateSurfaces = Surfaces.Read( content )

    # First int defines the case, all additional ints are used by the 
    # specific case. For a case number greater then 2 a file name is
    # read from the following line.
    physicalCellCase = map( int, content.pop(0).split('!')[0].split() )
    if physicalCellCase[0] <= 2:
      self.physicalCellCase = ( physicalCellCase[0], physicalCellCase[1:] )
    else:
      Fformat = FortranRecordReader( '(A72)' )
      self.physicalCellCase = ( physicalCellCase[0], 
                                Fformat.read( content.pop(0) )[0] )

    Fformat = FortranRecordReader( '(L1)' )
    self.check = Fformat.read( content.pop(0) )[0]

    if len(content) > 0:
      print 'There are still geometry input data to be read.'

  def _ReadGridInformation( self, content ):
    """Method that reads the grid information.

    This method reads the grid information and returns the number of 
    zones (int), number of domains (int), zone indices (list of int) and 
    number of radial, poloidal and toroidal surfaces for each zone (list 
    of list of int).
    
    Keyword arguments:
    content -- Content of the input file. (list of str)

    """
    Fformat = FortranRecordReader( '(i4,2x,i4)' )
    NZones, NDomains = Fformat.read( content.pop(0).split('!')[0] )
    zoneIndices = []
    for k in xrange( NDomains ):
      zoneIndices.append( int( content.pop(0).split('!')[0] ) )
    if NDomains <= 0:
      NDomains = NZones
    NSurfaces = []
    for k in xrange( NZones ):
      NSurfaces.append( tuple( map( int, content.pop(0).split('!')[0].split() ) ) )
    return NZones, NDomains, zoneIndices, NSurfaces

