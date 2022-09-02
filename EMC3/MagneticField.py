from numpy import zeros
import FortranReader
from GeometryParameters import GeometryParameters
from Plot import InteractivePoloidalCrossSection

class MagneticField( list ):
  """Class that represents the magnetic field output from the grid generator.
  
  This class can read the magnetic field (B) produced by the grid generator 
  for the EMC3-EIRENE code and provides is as a list of ndarrays for the 
  different zones. 
  A plot method is implemented that enables a quick view on the data for 
  different toroidal angles.
  
  """

  def __init__( self, name = 'BFIELD_STRENGTH', path = '../../geometry', 
                geometryParameters = None ):
    """Constructor

    Keyword arguments:
    name -- Name of the magnetic field file. (str, default 
            'BFIELD_STRENGTH')
    path -- Path to the file. (str, default '../../geometry')
    geometryParameters -- Geometry parameters. If None trying to initialize 
                          it with data in 'path'. (GeometryParameters 
                          object, default None)

    """
    list.__init__( self )

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

    for NRadial, NPoloidal, NToroidal in self.geometry.NSurfaces:
      NCrossSection = NRadial * NPoloidal # Grid points in the poloidal cross-section
      B = zeros( (NToroidal, NCrossSection) ) # Magnetic field in Tesla
      for k in xrange( NToroidal ):
        B[k,:] = FortranReader.List( content, NCrossSection )
      B.shape = (NToroidal, NPoloidal, NRadial)
      self.append( B )

    if len(content) > 0:
      print 'There are still magnetic field data to be read.'

  def plot( self, grid, zoneIndex = 0, installations = [], xlim = [], 
            ylim = [], levels = None ):
    """Plot method for a quick view on the grid.

    This is a plot method that provides a quick view on the grid. 
    
    Keyword arguments:
    grid -- Grid that correspond to the magnetic field. (Grid object)
    zoneIndex -- Index of zone to plot. (int, default 0)
    installations -- Installations to be plotted. ((list of) Installation 
                     object(s))
    xlim -- Range of the x-axis. (list of int, default [])
    ylim -- Range of the y-axis. (list of int, default [])
    levels -- Level curves to draw. (list of float, default None)
    
    """
    zone = grid.zones[ zoneIndex ]
    B = self[ zoneIndex ]

    IP = InteractivePoloidalCrossSection( zone.phi, zone.R, zone.Z, B, 
                                          xlim, ylim, colourBarLabel = 'B [T]', 
                                          levels = levels, 
                                          installations = installations )

