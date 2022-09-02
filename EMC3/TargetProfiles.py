from numpy import linspace, meshgrid, pi, array, append
from fortranformat import FortranRecordReader
import FortranReader

class TargetProfiles( list ):
  """Class that represents the target profile ouputs of EMC3-EIRENE. (list)
  
  This class can read the target profiles output. 

  """

  def __init__( self, name = 'TARGET_PROFILES', path = 'EMC3_OUTPUT' ):
    """Constructor

    Keyword arguments:
    name -- Name of the target profiles file to read. (str, default 
            'TARGET_PROFILES') 
    path -- Path to the target profiles file. (str, default 'EMC3_OUTPUT')

    """
    self.name = name
    self.path = path.rstrip('/')

    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content.reverse() # to allow a much more efficient pop from the end of the list

    self.NTotal, self.NTargetElements = map( int, content.pop().split('!')[0].split() )
    
    for k in xrange( self.NTargetElements ):    
      self.append( TargetProfilesElement( content ) )

class TargetProfilesElement( object ):
  """Class that represents an element of the target profiles output of EMC3-EIRENE.

  This class reads the element parts of the target profiles and provides 
  features to perform the bilinear interpolation of the segments as 
  required for plotting the data.

  """

  def __init__( self, content ):
    """Constructor

    Keyword arguments:
    content -- Content to read. (list of str) 

    """
    FformatHeader = FortranRecordReader( 'I4,5x,A' )
    FformatTotal = FortranRecordReader( '1p2E12.4,A' )
    FformatDimensions = FortranRecordReader( '3I6' )

    elementId, fileName = FformatHeader.read( content.pop() )
    self.elementId = elementId
    self.fileName = fileName.strip()

    particle, power = FformatTotal.read( content.pop() )[:2] # ignore string
    self.totalParticleDeposition = ( 'A', particle )
    self.totalPowerDeposition = ( 'W', power )

    self.NSegments, \
    self.NPoloidal, \
    self.NToroidal = FformatDimensions.read( content.pop() )

    print self.NSegments
    print self.NPoloidal
    print self.NToroidal


    N = 4 * self.NSegments

    RCorner = FortranReader.List( content, N, '1p6E12.4' )
    ZCorner = FortranReader.List( content, N, '1p6E12.4' )
    phiCorner = FortranReader.List( content, N, '1p6E12.4' )

    RCorner.shape = ( self.NSegments, 4 )
    ZCorner.shape = ( self.NSegments, 4 )
    phiCorner.shape = ( self.NSegments, 4 )

    self.RCorner = ( 'cm', RCorner )
    self.ZCorner = ( 'cm', ZCorner )
    self.phiCorner = ( 'rad', phiCorner )

    N = self.NSegments * self.NPoloidal * self.NToroidal

    Gamma = FortranReader.List( content, N, '1p6E12.4' )
    q = FortranReader.List( content, N, '1p6E12.4' )
    n = FortranReader.List( content, N, '1p6E12.4' )
    Te = FortranReader.List( content, N, '1p6E12.4' )
    Ti = FortranReader.List( content, N, '1p6E12.4' )

    self.profiles = {}
    self.profiles['Gamma'] = ( 'A cm^{-2}', Gamma )
    self.profiles['q'] = ( 'W cm^{-2}', q )
    self.profiles['n'] = ( 'cm^{-3}', n )
    self.profiles['Te'] = ( 'eV', Te )
    self.profiles['Ti'] = ( 'eV', Ti )

    self.__interpolationParametersAlpha = []
    self.__interpolationParametersBeta = []

    self.__R = []
    self.__phi = []
    self.__Z = []

  def __calculateInterpolationParameters( self ):
    """Calculate parameters required for bilinear interpolation.

    This method calculates the parameters that are required for the bilinear
    interpolation from the target segments to the fine target segments.

    """
    self.__interpolationParametersAlpha = []
    self.__interpolationParametersBeta = []
    for R, Z in zip( self.RCorner[1], 
                     self.ZCorner[1] ):
      alpha = []
      alpha.append( 0.25 * (  R[0] + R[1] + R[2] + R[3] ) )
      alpha.append( 0.25 * ( -R[0] - R[1] + R[2] + R[3] ) )
      alpha.append( 0.25 * ( -R[0] + R[1] + R[2] - R[3] ) )
      alpha.append( 0.25 * (  R[0] - R[1] + R[2] - R[3] ) )

      beta = []
      beta.append( 0.25 * (  Z[0] + Z[1] + Z[2] + Z[3] ) )
      beta.append( 0.25 * ( -Z[0] - Z[1] + Z[2] + Z[3] ) )
      beta.append( 0.25 * ( -Z[0] + Z[1] + Z[2] - Z[3] ) )
      beta.append( 0.25 * (  Z[0] - Z[1] + Z[2] - Z[3] ) )

      self.__interpolationParametersAlpha.append( alpha )
      self.__interpolationParametersBeta.append( beta )

  @property
  def interpolationParametersAlpha( self ):
    """(Calculates and) return the interpolation parameters alpha"""
    if len( self.__interpolationParametersAlpha ) == 0:
      self.__calculateInterpolationParameters()
    return self.__interpolationParametersAlpha

  @property
  def interpolationParametersBeta( self ):
    """(Calculates and) return the interpolation parameters beta"""
    if len( self.__interpolationParametersBeta ) == 0:
      self.__calculateInterpolationParameters()
    return self.__interpolationParametersBeta

  def __interpolateGrid( self, iSegment, l, m ):
    """ Give interpolated point on target segment. (three floats)

    This method performes a bilinear interpolation on the target segment
    with the given index and returns the r-, phi-, z-coordinate of the point.

    Keyword arguements:
    iSegment -- Index of target segment. (int)
    l -- Interpolation in poloidal direction. (float in [-1,1])
    m -- Interpolation in toroidal direction. (float in [-1,1])

    """
    phiCorner = self.phiCorner[1][iSegment,:]
    phi = 0.5 * ( phiCorner[0] * ( 1 - m ) + phiCorner[1] * ( 1 + m ) )

    alpha = self.interpolationParametersAlpha[iSegment]
    beta = self.interpolationParametersBeta[iSegment]
    R = alpha[0] + l * alpha[1] + m * alpha[2] + l * m * alpha[3]
    Z = beta[0] + l * beta[1] + m * beta[2] + l * m * beta[3]

    return R, phi, Z

  def __calculateFineTargetSegments( self ):
    """Calculates and returns the fine target segments grid. (list of ndarray)

    This method calculates the fine target segments grid based on bilinear
    interpolation and returns the corner coordinates of the cells for each 
    target segment.
    The toroidal coordinated is translated to degree.

    """
    self.__R = []
    self.__phi = []
    self.__Z = []
    ls = linspace( -1, 1, self.NPoloidal + 1 )
    ms = linspace( -1, 1, self.NToroidal + 1 )
    lls, mms = meshgrid( ls, ms )
    for iSegment in xrange( self.NSegments ):
      RFine, phiFine, ZFine = self.__interpolateGrid( iSegment, lls, mms )
      for iToroidal in xrange( self.NToroidal ):
        for iPoloidal in xrange( self.NPoloidal ):
          self.__R.append( append( RFine[ iToroidal, 
                                          iPoloidal:iPoloidal + 2 ], 
                                   RFine[ iToroidal + 1, 
                                          iPoloidal:iPoloidal + 2 ][::-1] ) )
          self.__phi.append( append( phiFine[ iToroidal, 
                                              iPoloidal:iPoloidal + 2 ],
                                     phiFine[ iToroidal + 1,
                                              iPoloidal:iPoloidal + 2 ][::-1] ) )
          self.__Z.append( append( ZFine[ iToroidal,  
                                          iPoloidal:iPoloidal + 2 ],
                                   ZFine[ iToroidal + 1,
                                          iPoloidal:iPoloidal + 2 ][::-1] ) )
    self.__R = array( self.__R )
    self.__phi = array( self.__phi ) / pi * 180.0
    self.__Z = array( self.__Z )

  @property
  def R( self ):
    """(Calculates and) return the R-coordinates of the fine target segment (cm)"""
    if len( self.__R ) == 0:
      self.__calculateFineTargetSegments()
    return self.__R

  @property
  def phi( self ):
    """(Calculates and) return the phi-coordinates of the fine target segment (degree)"""
    if len( self.__phi ) == 0:
      self.__calculateFineTargetSegments()
    return self.__phi

  @property
  def Z( self ):
    """(Calculates and) return the Z-coordinates of the fine target segment (cm)"""
    if len( self.__Z ) == 0:
      self.__calculateFineTargetSegments()
    return self.__Z

