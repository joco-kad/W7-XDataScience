from numpy import zeros
from os.path import expanduser
from fortranformat import FortranRecordReader
from Plot import InteractivePoloidalCrossSection

class Installation( object ):
  """Class that represents the location of installations in the Kisslinger 
  format.
  
  This class reads the geometry of installations in the Kisslinger format 
  (R,Z-coorinates per toroidal angle), like limiters, divertors, or the 
  wall, and provides the coordinates of its nodes in the three ndarrays Phi 
  (degree), R (cm), and Z (cm).
  A plot method is implemented that enables a quick view on the data for 
  the different toroidal angles.
  
  """

  def __init__( self, name, path = './' ):
    """Constructor

    Keyword arguments:
    name -- Name of the installation file. (str)
    path -- Path to the file. (str, default './')

    """
    self.name = name
    self.path = expanduser( path.rstrip('/') )

    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content.reverse()

    Fformat = FortranRecordReader( '(A60)' )
    self.title = Fformat.read( content.pop() )
    values = content.pop().split()
    NToroidal, NPoloidal, self.NFoldPeriodicity = map( int, values[:3] )
    RReference, ZReference = map( float, values[3:5] )
    
    self.phi = zeros( NToroidal )
    self.R = zeros( ( NToroidal, NPoloidal ) )
    self.Z = zeros( ( NToroidal, NPoloidal ) )
    
    for kToroidal in xrange( NToroidal ):
      self.phi[kToroidal] = float( content.pop() )
      for kPoloidal in xrange( NPoloidal ):
        self.R[ kToroidal, kPoloidal ], \
        self.Z[ kToroidal, kPoloidal ] = map( float, content.pop().split()[:2] )

    self.R = self.R + RReference
    self.Z = self.Z + ZReference

  def plot( self, xlim = [], ylim = [] ):
    """Plot method for a quick view on the installation.

    This is a plot method that provides a quick view on the installation. 
    
    Keyword arguments:
    xlim -- Range of the x-axis. (list of int, default [])
    ylim -- Range of the y-axis. (list of int, default [])
    
    """
    IP = InteractivePoloidalCrossSection( self.phi, xlim = xlim, ylim = ylim, 
                                          installations = self )

