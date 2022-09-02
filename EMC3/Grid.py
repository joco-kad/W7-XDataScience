from GeometryParameters import GeometryParameters
from Zone import Zone

class Grid( object ):
  """Class that represents the grid output from the grid generator.
  
  This class can read the grid coordinates (R, Z, Phi) produced by the 
  grid generator for the EMC3-EIRENE code. 
  
  """

  def __init__( self, name = 'GRID_3D_DATA', path = '../../geometry', 
                geometryParameters = None ):
    """Constructor

    Keyword arguments:
    name -- Name of the grid file. (str, default 'GRID_3D_DATA')
    path -- Path to the file. (str, default '../../geometry')
    geometryParameters -- Geometry parameters. If None trying to initialize 
                          it with data in 'path'. (GeometryParameters 
                          object, default None)

    """
    self.name = name
    self.path = path.rstrip('/')

    if geometryParameters is None:
      self.geometry = GeometryParameters()
    else:
      self.geometry = geometryParameters

    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content.reverse() # to allow a much more efficient pop from the end of the list

    self.zones = []
    for k in xrange( self.geometry.NZones ):
      zone = Zone( content )
      if self.geometry.NSurfaces[k] == (zone.NRadial, zone.NPoloidal, 
                                        zone.NToroidal):
        self.zones.append( zone )
      else:
        raise Exception( 'Information on number of surfaces does not ' + \
                         'match between geometry and grid input file.' )

    if len(content) > 0:
      print 'There are still grid data to be read.'

