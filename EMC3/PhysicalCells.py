from numpy import zeros, int32
import FortranReader
from GeometryParameters import GeometryParameters

class PhysicalCells( object ):
  """Class that represents the link between physical and magnetic cells.
  
  This class represents the cell geometry file that links the output 
  values of EMC3-EIRENE code to the cells of the magnetic grid used.
  
  """

  def __init__( self, name = 'CELL_GEO', path = '.', 
                geometryParameters = None ):
    """Constructor

    Keyword arguments:
    name -- Name of the physical cell file. (str, default 'CELL_GEO')
    path -- Path to the file. (str, default '.')
    geometryParameters -- Geometry parameters. If None trying to initialize 
                          it with data in 'path'. (GeometryParameters 
                          object, default None)

    """
    self.name = name
    self.path = path.rstrip('/')

    file = open( '%s/%s' % (self.path, self.name) )
    content = file.read().splitlines()
    file.close()
    content.reverse() # to allow a much more efficient pop from the end of the list

    if geometryParameters is None:
      self.geometry = GeometryParameters()
    else:
      self.geometry = geometryParameters
    
    Ns = content.pop().split()
    N = int(Ns[0]) # Number of data points in list
    self.NCellPlasma = int(Ns[1]) # Total cell number for plasma transport
    self.NCellNeutral = int(Ns[2]) # Total cell number for neutral transport

    idsAllZones = FortranReader.List( content, N, dtype = int ) # ids of the plasma values

    if len(content) > 0:
      print 'There are still data to be read.'

    self.ids = []
    Offset = 0
    for NRadial, NPoloidal, NToroidal in self.geometry.NSurfaces:
      NCells = (NRadial-1) * (NPoloidal-1) * (NToroidal-1)
      ids = idsAllZones[ Offset : Offset + NCells ]
      ids.shape = ( NToroidal-1, NPoloidal-1, NRadial-1 )
      self.ids.append( ids )
      Offset += NCells

