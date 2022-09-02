from Installation import Installation

class TargetFile( list ):
  """Class that represents the target file for the deposition post-processing.
  
  This class represents the input of EMC3 for the target file for the
  deposition post-processing. On initialization it reads an existing input 
  file. 
  
  """

  def __init__( self, name, path = '../run' ):
    """Constructor

    Keyword arguments:
    name -- Name of the target file input. (str)
    path -- Path to the file. (str, default '../run')

    """
    self.name = name
    self.path = path.rstrip('/')
    
    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: len(s) != 0, content )
    content = filter( lambda s: s[0] != '*', content )

    NFiles = int( content.pop(0).split('!')[0] )
    for k in xrange( NFiles ):
      file = {}
      line = content.pop(0).split()
      line[0] = line[0].strip("'")
      installationFileName = line[0].split('/')[-1]
      installationPath = line[0].rstrip( installationFileName )
      file['installation'] = Installation( installationFileName, 
                                           installationPath )
      file['NPoloidal'], \
      file['NToroidal'] = map( int, line[1:3] )

      self.append( file )

    if len(content) > 0:
      print 'There are still target file input data to be read.'

