from numpy import zeros
from fortranformat import FortranRecordReader
import Surfaces

class NeutralParameters( object ):
  """Class that represents the neutral parameters.
  
  This class represents the input of EMC3 for neutral parameters. On 
  initialization it reads an existing input file. 
  After modifications of the parameter a new input file can be written.
  
  ATTENTSION:
  WriteAdditionalSurfaces and WriteLimiters is missing.
  
  """

  def __init__( self, name = 'input.n0g', path = '../run' ):
    """Constructor

    Keyword arguments:
    name -- Name of the parameter input. (str, default 'input.n0g')
    path -- Path to the file. (str, default '../run')

    """
    self.name = name
    self.path = path.rstrip('/')
    
    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: '*' != s[0], content )
    
    Fformat = FortranRecordReader( '(A72)' )

    self.nonTransparentSurfaces = Surfaces.Read( content )

    NSetsOfCells, \
    self.streamNumberPhysicalCells \
                        = map( int, content.pop(0).split('!')[0].split()[0:2] )

    self.cells = []
    for k in xrange( NSetsOfCells ):
      cell = {}
      cell['definition'], \
      cell['type'] = map( int, content.pop(0).split('!')[0].split()[0:2] )
      if cell['definition'] == 2:
        values = map( int, content.pop(0).split('!')[0].split()[0:10] )
        cell['zoneIndex'] = values[0]
        cell['indexRangeRadial'] = values[1:3]
        cell['stepSizeRadial'] = values[3]
        cell['indexRangePoloidal'] = values[4:6]
        cell['stepSizePoloidal'] = values[6]
        cell['indexRangeToroidal'] = values[7:9]
        cell['stepSizeToroidal'] = values[9]
      elif cell['definition'] == 3:
        cell['cellInformationFile'] = Fformat.read( content.pop(0) )[0]
      cell['n'], \
      cell['Te'], \
      cell['Ti'], \
      cell['M'] = map( float, content.pop(0).split('!')[0].split()[0:4] )
      self.cells.append( cell )

    self.sourceDefinition, \
    self.placeSourcePoints, \
    self.sideSourcePoints = map( int, 
                                 content.pop(0).split('!')[0].split()[0:3] )
    if self.sourceDefinition < 0:
      self.neutralSourceFile = Fformat.read( content.pop(0) )[0]
    self.additionalSurfacesFile = Fformat.read( content.pop(0) )[0].strip()

    if len(content) > 0:
      print 'There are still neutral parameter input data to be read.'

  def Write( self, name = None ):
    """Writes the input file.

    This methods creates the input file based on the data stored in this 
    object.

    Keyword arguments.
    name -- File name to write data to. If None, read-in file will be 
            overwritten. (str, default None)

    """
    if name is None:
      name = self.name
    f = open( name, 'w' )
    f.write( '***************** additional geometry and parameters for N0 ' + \
             '*********\n' )
    f.write( '*\n' )
    f.write( '* This file was created automaticly.\n' )
    f.write( '* Use the EMC3 python tool to modify it.\n' )
    f.write( '*\n' )
    f.write( '*** 1. non-transparent surfaces for neutral\n' )
    f.write( '*  non-transparent surfaces with informations about this ' + \
             'surface\n' )
    f.write( '*  being defined in Eirene. The surface number must be ' + \
             'indicated here.\n' )
    Surfaces.Write( f, self.nonTransparentSurfaces )
    f.write( '*** 2. DEFINE ADDITIONAL PHYSICAL CELLS FOR NEUTRALS\n' )
    f.write( '* the cells can be defined in following ways, depending on ' + \
             'IND_CELL\n' )
    f.write( '* IND_CELL=0: default case, geometric cell = physical cell\n' )
    f.write( '*          1: read cell file from fort.2 (defined by EMC3 ' + \
             'code)\n' )
    f.write( '*             (If it is the case, this kind of cells must ' + \
             'be first\n' )
    f.write( '*              defined since the cell number in EMC3 begins ' + \
             'always\n' )
    f.write( '*              with the number 1.)\n' )
    f.write( '*             ** not for the coupled case\n' )
    f.write( '*          2: restrict a range and give the steps\n' )
    f.write( '*          3: file provided by the user\n' )
    f.write( '* Define plasma: ne, Te, Ti, Mach_number\n' )
    f.write( '* number of different definitions\n' )
    f.write( '%d %d\n' % ( len( self.cells ), self.streamNumberPhysicalCells ) )
    f.write( '* One cell for plasma core (n0)\n' )
    f.write( '* 2.1. restrict a range and give the steps\n' )
    for cell in self.cells:
      f.write( '%d %d\n' % ( cell['definition'], cell['type'] ) )
      if cell['definition'] == 2:
        f.write( '* ZONE, R1, R2, DR, P1, P2, DP, T1, T2, DT\n' )
        f.write( '%d ' % cell['zoneIndex'] )
        f.write( '%d %d ' % tuple( cell['indexRangeRadial'] ) )
        f.write( '%d ' % cell['stepSizeRadial'] )
        f.write( '%d %d ' % tuple( cell['indexRangePoloidal'] ) )
        f.write( '%d ' % cell['stepSizePoloidal'] )
        f.write( '%d %d ' % tuple( cell['indexRangeToroidal'] ) )
        f.write( '%d\n' % cell['stepSizeToroidal'] )
      elif cell['definition'] == 3:
        f.write( '%s\n' % cell['cellInformationFile'] )
      f.write( '* plasma parameter\n' )
      f.write( '* ne, Te, Ti, M\n' )
      f.write( '%g %g %g %g\n' % ( cell['n'], cell['Te'], cell['Ti'], 
                                   cell['M'] ) )
    f.write( '*** 3 Neutral Source distribution\n' )
    f.write( '*   Neutral source distribution can either be readed from a ' + \
             'file\n' )
    f.write( '*   (N0S<0) or from EMC3 Code (N0S=0)\n' )
    f.write( '* NS_PLACE=0: the source points are on the additional ' + \
             'surfaces defined\n' )
    f.write( '*             in the geometric part\n' )
    f.write( '*          1: the source points are on the additional ' + \
             'surfaces defined\n' )
    f.write( '*             in the atomic part.\n' )
    f.write( '* NSSIDE=-1,1 (only for NS_PLACE=0)the source points are on ' + \
             'the\n' )
    f.write( '*             negative/positiside of a surface\n' )
    f.write( '* N0S NS_PLACE  NSSIDE\n' )
    f.write( '%d %d %d\n' % ( self.sourceDefinition, self.placeSourcePoints, 
                              self.sideSourcePoints ) )
    if self.sourceDefinition < 0:
      f.write( '* source.distr\n' )
      f.write( '%s\n' % self.neutralSourceFile )
    f.write( '*** 4 Additional surfaces\n' )
    f.write( '* The additional surfaces are represented by triangles.\n' )
    f.write( '* You can either give the triangles directly, or indirectly ' + \
             'give\n' )
    f.write( '* the file containing the plate with NTRIANG=0.\n' )
    f.write( '%s' % self.additionalSurfacesFile )
    f.close()

  def ReadAdditionalSurfaces( self, path = '.' ):
    """Reads the additional surfaces.

    This method reads and returns the additional surfaces from the file 
    specified in the neutral parameter input.

    Keyword arguments:
    path -- Path to the file. (str, default '.')

    """
    file = open( '%s/%s' % (path, self.additionalSurfacesFile) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: '*' != s[0], content )
    
    Fformat = FortranRecordReader( '(A72)' )

    NPlates = int( content.pop(0).split('!')[0].split()[0] )
    plates = []
    for k in xrange( NPlates ):
      plate = {}
      values = content.pop(0).split('!')[0].split()[0:13]
      NTriangles = int( values[0] )
      plate[ 'type' ] = int( values[1] )
      plate[ 'divertorBaffle' ] = int( values[2] )
      plate[ 'sputterCoefficients' ] = map( float, values[3:] )
      if NTriangles == 0:
        plate[ 'fileName' ] = Fformat.read( content.pop(0) )[0].strip()
      else:
        plate[ 'triangles' ] = []
        for k in xrange( NTriangles ):
          values = map( float, content.pop(0).split('!')[0].split()[0:6] )
          point1 = tuple( values[:3] )
          point2 = tuple( values[3:] )
          point3 = tuple( map( float, 
                               content.pop(0).split('!')[0].split()[0:3] ) )
          plate[ 'triangles' ].append( [ point1, point2, point3 ] )
      plates.append( plate )
    self.additionalSurfaces = plates

  def ReadLimiters( self, path = '.' ):
    """Reads the limiters.

    This method reads and returns the limiters from the file specified in the 
    additional surface file.

    Keyword arguments:
    path -- Path to the file. (str, default '.')

    REMARK:
    Name of method follows naming in EMC3, but is probably not a good 
    discription. Needs to be clarified at a later stage.

    """
    Fformat = FortranRecordReader( '(A60)' )

    for k in xrange( len( self.additionalSurfaces ) ):
      file = open( '%s/%s' % (path, self.additionalSurfaces[k]['fileName']) )
      content = file.read().splitlines()
      file.close()
      content = filter( lambda s: '*' != s[0], content )

      limiter = {}
      limiter[ 'name' ] = Fformat.read( content.pop(0) )[0].strip()
      values = content.pop(0).split('!')[0].split()[0:5]
      NPlates = int( values[0] )
      NPoints = int( values[1] )
      limiter[ 'periodicity' ] = int( values[2] )
      limiter[ 'RRef' ], limiter[ 'ZRef' ] = map( float, values[3:] )
      phi = zeros( NPlates )
      R = zeros( ( NPlates, NPoints ) )
      Z = zeros( ( NPlates, NPoints ) )
      for l in xrange( NPlates ):
        phi[l] = float( content.pop(0).split('!')[0].split()[0] )
        for m in xrange( NPoints ):
          R[l,m], Z[l,m] = map( float, 
                                content.pop(0).split('!')[0].split()[0:2] )
      limiter[ 'phi' ] = phi
      limiter[ 'R' ] = R
      limiter[ 'Z' ] = Z
      self.additionalSurfaces[k]['limiter'] = limiter

