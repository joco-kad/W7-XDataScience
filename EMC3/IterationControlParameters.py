from fortranformat import FortranRecordReader, FortranRecordWriter

class IterationControlParameters( object ):
  """Class that represents the iteration control parameters.
  
  This class represents the input of EMC3 for iteration control parameters. 
  On initialization it reads an existing input file. 
  After modifications of the parameter a new input file can be written.
  
  """

  def __init__( self, name = 'input.ctr', path = '../run' ):
    """Constructor

    Keyword arguments:
    name -- Name of the parameter input. (str, default 'input.ctr')
    path -- Path to the file. (str, default '../run')

    """
    self.name = name
    self.path = path.rstrip('/')
    
    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: len(s) != 0, content )
    content = filter( lambda s: s[0] != '*', content )

    line = content.pop(0)
    self.globalIterations = []
    while line.find( 'FIN' ) == -1 and line.find( 'END' ) == -1:
      while line.find( 'FIN' ) == -1 and line.find( 'END' ) == -1:
        globalIteration = {}
        globalIteration[ 'N' ], \
        globalIteration[ 'MCSeed' ] = map( int, line.split('!')[0].split()[:2] )
        if globalIteration[ 'N' ] > 0:
          globalIteration[ 'NSubIterations' ], \
          NProcesses = map( int, content.pop(0).split('!')[0].split()[:2] )
          processes = []
          for k in xrange( NProcesses ):
            process = {}
            process[ 'name' ] = content.pop(0).split('!')[0]
            process[ 'integerParameters' ] \
                          = map( int, content.pop(0).split('!')[0].split()[:4] )
            process[ 'floatParameters' ] \
                        = map( float, content.pop(0).split('!')[0].split()[:4] )
            processes.append( process )
          globalIteration[ 'processes' ] = processes
        self.globalIterations.append( globalIteration )
        line = content.pop(0)
      line = content.pop(0)

    if len(content) > 0:
      print 'The iteration control parameter input contains unread data.'

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
    f.write( '* Program control\n' )
    f.write( '*\n' )
    f.write( '* This file was created automaticly.\n' )
    f.write( '* Use the EMC3 python tool to modify it.\n' )
    f.write( '*\n' )
    for globalIteration in self.globalIterations:
      f.write( '%d %d\n' % ( globalIteration[ 'N' ], 
                             globalIteration[ 'MCSeed' ] ) )
      if globalIteration[ 'N' ] > 0:
        f.write( '* SUB. ITER  TRANSPORT\n' )
        f.write( '%d %d\n' % ( globalIteration[ 'NSubIterations' ], 
                               len( globalIteration[ 'processes' ] ) ) )
        for process in globalIteration[ 'processes' ]:
          f.write( '%s\n' % process[ 'name' ] )
          f.write( '%d %d %d %d\n' % tuple( process[ 'integerParameters' ] ) )
          f.write( '%g %g %g %g\n' % tuple( process[ 'floatParameters' ] ) )
      f.write( 'END\n')
    f.write( 'FIN' )
    f.close()

