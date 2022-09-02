from numpy import zeros
import matplotlib.pyplot as plt
import FortranReader

class FluxConservation( object ):
  """Class that represents the output of the flux conservation check from the 
  grid generator.
  
  This class can read the output of the flux conservation check (R, Z, 
  diviation) from the grid generator for the EMC3-EIRENE code. 
  A plot method is implemented that enables a quick view on the data.
  
  """

  def __init__( self, name = 'FLUX_CONSERVATION', path = '../../geometry' ):
    """Constructor

    Keyword arguments:
    name -- Name of the flux conservation file. (str, default 
            'FLUX_CONSERVATION')
    path -- Path to the file. (str, default '../../geometry')

    """
    self.name = name
    self.path = path.rstrip('/')

    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content.reverse() # to allow a much more efficient pop from the end of the list
    
    Ns = content.pop().split()
    NR = int(Ns[0]) # Number of radial grid points 
    NZ = int(Ns[1]) # Number of vertical grid points
    NCrossSectionGrid = NR * NZ # Grid points in the poloidal cross-section
    NCrossSectionSurface = (NR-1) * (NZ-1) # Grid points in the poloidal cross-section
    
    self.R = FortranReader.List( content, NCrossSectionGrid , '(8f9.3)' ) # Radial positions in cm
    self.Z = FortranReader.List( content, NCrossSectionGrid , '(8f9.3)' ) # Horizontal positions in cm
    self.diviation = FortranReader.List( content, NCrossSectionSurface , '(1p6E12.4)' ) # Diviation in percent
    
    if len(content) > 0:
      print 'There are still data to be read.'
    
    self.R.shape = (NZ, NR)
    self.Z.shape = (NZ, NR)
    self.diviation.shape = (NZ-1, NR-1)
    
  def plot( self, xlim = [], ylim = [], edgecolors = 'grey' ):
    """Plot method for a quick view on the diviation.

    This is a plot method that provides a quick view on the diviation. 
    
    Keyword arguments:
    xlim -- Range of the x-axis. (list of int, default [])
    ylim -- Range of the y-axis. (list of int, default [])
    edgecolors -- Edge colour of mesh. (str, default 'grey')
    
    """
    plt.pcolor( self.R, self.Z, self.diviation, edgecolors = edgecolors )
    plt.colorbar()
    plt.axis( 'scaled' )
    plt.xlabel( 'x [cm]' )
    plt.ylabel( 'y [cm]' )
    plt.xlim( *xlim )
    plt.ylim( *ylim )
    plt.show()

