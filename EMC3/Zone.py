from numpy import zeros, double, abs, mean, diff, radians, pi, sqrt
import FortranReader
from Analysis import getAreaQuad, getAreaPolygon
from Plot import InteractivePoloidalCrossSection

class Zone( object ):
  """Class that represents the grid output from the grid generator for a zone.
  
  This class can read the grid coordinates (R, Z, Phi) produced by the grid 
  generator for the EMC3-EIRENE code for one zone. 
  A plot method is implemented that enables a quick view on the data for 
  different toroidal angles.
  
  """

  def __init__( self, content ):
    """Constructor

    Keyword arguments:
    content -- Content of the grid file in reversed order. (list of str)

    """
    Ns = content.pop().split()
    self.NRadial = int(Ns[0]) # Number of radial surfaces 
    self.NPoloidal = int(Ns[1]) # Number of poloidal surfaces 
    self.NToroidal = int(Ns[2]) # Number of toroidal surfaces
    self.RShift = 0 # Shift of zone in radial direction
    self.ZShift = 0 # Shift of zone in horizontal direction
    if len(Ns) > 3:
      self.RShift = int(Ns[3])
      if len(Ns) > 4:
        self.ZShift = int(Ns[4])
    self.shifted = abs(self.RShift) + abs(self.ZShift) > 0
    NCrossSection = self.NRadial * self.NPoloidal # Grid points in the poloidal cross-section
    
    self.phi = zeros( self.NToroidal ) # Toroidal angle in degree
    self.R = zeros( (self.NToroidal, NCrossSection) ) # Radial positions in cm
    self.Z = zeros( (self.NToroidal, NCrossSection) ) # Horizontal positions in cm
    
    for k in xrange( self.NToroidal ):
      self.phi[k] = float( content.pop() )
      self.R[k,:] = FortranReader.List( content, NCrossSection )
      self.Z[k,:] = FortranReader.List( content, NCrossSection )
    
    if self.shifted:
      self.R += self.RShift
      self.Z += self.ZShift

    self.R.shape = (self.NToroidal, self.NPoloidal, self.NRadial)
    self.Z.shape = (self.NToroidal, self.NPoloidal, self.NRadial)

    self.__volume = double( [] )
    self.__Reff = double( [] )
    
  def plot( self, installations = [], xlim = [], ylim = [], 
            edgecolors = 'grey' ):
    """Plot method for a quick view on the grid.

    This is a plot method that provides a quick view on the grid. 
    
    Keyword arguments:
    installations -- Installations to be plotted. ((list of) Installation 
                     object(s))
    xlim -- Range of the x-axis. (list of int, default [])
    ylim -- Range of the y-axis. (list of int, default [])
    edgecolors -- Edge colour of mesh. (str, default 'grey')
    
    """
    UntilABetterSolutionIsFound = zeros( ( self.Z.shape[0],
                                           self.Z.shape[1]-1, 
                                           self.Z.shape[2]-1 ) )
    IP = InteractivePoloidalCrossSection( self.phi, self.R, self.Z, 
                                          UntilABetterSolutionIsFound, 
                                          xlim, ylim, edgecolors = edgecolors,
                                          installations = installations )

  def __calculateVolume( self ):
    """Method to calculate the volume of each cell. (ndarray)

    This method calculates the volume of each cell assuming a linear
    change of the poloidal cell cross-section and radial location.

    """
    volume = zeros( ( self.NToroidal-1, self.NPoloidal-1, self.NRadial-1 ) )
    phi = radians( diff( self.phi ) )
    
    for I in xrange(self.NRadial-1):
      for J in xrange(self.NPoloidal-1):
        R = mean( self.R[:,J:J+2,I:I+2], (1,2) )
        A = getAreaQuad( self.R[:,J:J+2,I:I+2], self.Z[:,J:J+2,I:I+2] )
        volume[:,J,I] = 0.25 * ( R[:-1] + R[1:] ) * ( A[:-1] + A[1:] ) * phi

    self.__volume = volume
    
  def __calculateReff( self ):
    """Method to calculate the effective radius for each radial surface. (ndarray)

    This method calculates the corresponding effective radius for each radial 
    surface index.

    """
    area = zeros( self.NRadial )

    for k in xrange( self.NRadial ):
      for j in xrange( self.NToroidal ):
        area[k] += getAreaPolygon( self.R[j,:,k], self.Z[j,:,k] )
    
    area /= self.NToroidal 
            
    self.__Reff = sqrt( area / pi )

  @property
  def volume( self ):
    """(Calculates and) returns the volume of each cell (cm**3)"""
    if self.__volume.size == 0:
      self.__calculateVolume()
    return self.__volume

  @property
  def Reff( self ):
    """(Calculates and) returns the effective radius of each radial surface (cm)"""
    if self.__Reff.size == 0:
      self.__calculateReff()
    return self.__Reff

