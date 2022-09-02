from numpy import ndarray, nonzero
from fortranformat import FortranRecordReader
import FortranReader
from GeometryParameters import GeometryParameters
from Plot import InteractivePoloidalCrossSection

class PlateCells( object ):
  """Class that represents the plate cells in the magnetic grid.
  
  This class reads the cells that are considered to be within/behind plates 
  and sets up a list of ndarray stating which cells are plates.
  A plot method is implemented that enables a quick view on the data for 
  different toroidal angles.
  
  """

  def __init__( self, name = 'PLATES_MAG', path = '../../geometry', 
                geometryParameters = None ):
    """Constructor

    Keyword arguments:
    name -- Name of the plate cells file. (str, default 'PLATES_MAG')
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
    print name
    file.close()
    content.reverse()

    self.ids = self._InitiateEmptyIdsArray()
    self._FillIdsArray( content )

  def _InitiateEmptyIdsArray( self ):
    """Initiates an empty array for the cell ids. (list of ndarray)

    This method initiates an array for the cell ids in the correct 
    shape that fits to the grid. Default value in the array is False.

    """
    allIds = []
    for NRadial, NPoloidal, NToroidal in self.geometry.NSurfaces:
      ids = ndarray( (NToroidal-1, NPoloidal-1, NRadial-1), dtype = bool )
      ids[:] = False
      allIds.append( ids )
    return allIds

  def _FillIdsArray( self, content ):
    """Fills the array for the cell ids.

    This method fills the array for the cell ids with True if a plate is 
    present in that cell.

    Keyword arguments:
    content -- Content of the plates file in reversed order. (list of str)

    ATTENTION:
    Second case has not been tested.

    """
    if self.geometry.plateSurfaces[ 'radial' ] == -1:
      while len( content ) > 0:
        N = map( int, content[-1].split() )[3]
        indices = FortranReader.List( content, N+4, dtype = int ) 
        zoneIndex, radialIndex, poloidalIndex = indices[:3]
        for k in xrange(N/2):
          self.ids[zoneIndex][ indices[4+2*k]:indices[5+2*k]+1, 
                               poloidalIndex, 
                               radialIndex ] = True
    elif self.geometry.plateSurfaces[ 'radial' ] == -2:
      Fformat = FortranRecordReader( '(2(I3,1X),1x,80I1)' )
      while len( content ) > 0:
        indices = Fformat.read( content.pop() )
        radialIndex, poloidalIndex = indices[:2]
        self.ids[0][ nonzero(indices[2:]), poloidalIndex, radialIndex ] = True
    elif self.geometry.plateSurfaces[ 'radial' ] == -3:
      while len( content ) > 0:
        N = map( int, content[-1].split() )[3]
        indices = FortranReader.List( content, N+4, dtype = int )
        zoneIndex, radialIndex, poloidalIndex = indices[:3]
        for k in xrange(N):
          self.ids[zoneIndex][ indices[4+k], poloidalIndex, radialIndex ] = True

  def plot( self, grid, zoneIndex = 0, installations = [], xlim = [], 
            ylim = [], edgecolors = 'grey' ):
    """Plot method for a quick view on the plate cells.

    This is a plot method that provides a quick view on the plate cells. A 
    linear interpolation in toroidal direction is used for the grid. 
    
    Keyword arguments:
    grid -- Grid that correspond to the plates. (Grid object)
    zoneIndex -- Index of zone to plot. (int, default 0)
    installations -- Installations to be plotted. ((list of) Installation 
                     object(s))
    xlim -- Range of the x-axis. (list of int, default [])
    ylim -- Range of the y-axis. (list of int, default [])
    edgecolors -- Edge colour of mesh. (str, default 'grey')
    
    ATTENTION: 
    Interpolation of installations does not necessary fit the interpolation 
    of the grid.

    """
    zone = grid.zones[zoneIndex]
    R = 0.5 * ( zone.R[:-1,:,:] + zone.R[1:,:,:] )
    Z = 0.5 * ( zone.Z[:-1,:,:] + zone.Z[1:,:,:] )
    phi = 0.5 * ( zone.phi[:-1] + zone.phi[1:] )
    ids = self.ids[zoneIndex]

    IP = InteractivePoloidalCrossSection( phi, R, Z, ids, xlim, ylim, [0,1],
                                          edgecolors = edgecolors, 
                                          installations = installations )

