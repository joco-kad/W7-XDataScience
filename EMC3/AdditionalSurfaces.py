from fortranformat import FortranRecordReader, FortranRecordWriter
from Installation import Installation

class AdditionalSurfaces( list ):
  """Class that represents the additional surfaces for neutrals.
  
  This class represents the input of EMC3 for additional surfaces for 
  neutrals. On initialization it reads an existing input file. 
  
  """

  def __init__( self, name, path = '../run' ):
    """Constructor

    Keyword arguments:
    name -- Name of the additional surface input. (str)
    path -- Path to the file. (str, default '../run')

    ATTENTION:
    Definition of plates via triangles is not tested!

    """
    self.name = name
    self.path = path.rstrip('/')
    
    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: len(s) != 0, content )
    content = filter( lambda s: s[0] != '*', content )

    NPlates = int( content.pop(0).split('!')[0] )
    for k in xrange( NPlates ):
      plate = {}
      line = content.pop(0)
      NTriangles, \
      plate['eireneId'], \
      plate['baffel'] = map( int, line.split()[:3] )
      plate['cSputter'] = map( float, line.split()[3:] )

      if NTriangles == 0:
       # plate given as Installation
       line = content.pop(0)
       installationFileName = line.split('/')[-1]
       installationPath = line.rstrip( installationFileName )
       plate['installation'] = Installation( installationFileName, 
                                             installationPath )
      else:
        plate['triangles'] = []
        for k in xrange( NTriangles ):
          line = content.pop(0)
          P1 = tuple( map( float, line.split()[:3] ) )
          P2 = tuple( map( float, line.split()[3:6] ) )
          P3 = tuple( map( float, content.pop(0).split()[:3] ) )
          plate['triangles'].append( [ P1, P2, P3 ] )

      plate['eireneId'] = -plate['eireneId']
      self.append( plate )

    if len(content) > 0:
      print 'There are still additional surface input data to be read.'

