class ImpurityParameters( object ):
  """Class that represents the impurity parameters.
  
  This class represents the input of EMC3 for impurity parameters. On 
  initialization it reads an existing input file. 
  After modifications of the parameter a new input file can be written.
  
  ATTENTION:
  Not tested!
  
  """


  def __init__( self, name, path = '../run', transport = False ):
    """Constructor

    Keyword arguments:
    name -- Name of the parameter input. (str)
    path -- Path to the file. (str, default '../run')
    transport -- If True, additioanl transport parameter are read. (bool, 
                 default False)

    """
    self.name = name
    self.path = path.rstrip('/')
    self.transport = transport
    
    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: '*' != s[0], content )

    self.typeTransportModel = int( content.pop(0).split('!')[0].split()[0] )
    self.database = content.pop(0).split('!')[0].split()[0]
    if self.database.find( 'ADAS' ) != -1:
      NDatabaseFiles = int( content.pop(0).split('!')[0].split()[0] )
      self.databaseFiles = []
      for k in xrange( NDatabaseFiles ):
        self.databaseFiles.append( content.pop(0).split('!')[0].split()[0] )
    if transport:
      self.physicalSputtering, \
      self.chemicalSputtering, \
      self.selfSputtering = map( float, 
                                 content.pop(0).split('!')[0].split()[0:3] )
      values = content.pop(0).split('!')[0].split()[0:3]
      self.energy = float( values[0] )
      self.alpha = int( values[1] )
      self.beta = int( values[2] )
      values = content.pop(0).split('!')[0].split()[0:2]
      NSources = int( values[0] )
      self.betaSources = float( values[1] )
      self.sourceParameter = []
      for k in xrange( NSources ):
        self.sourceParameters.append( map( float, content.pop(0).split('!')[0].split()[0:8] ) )
      self.typeSolCoreInterface = int( content.pop(0).split('!')[0].split()[0] )
      if self.typeSolCoreInterface > 1:
        for k in xrange( self.typeSolCoreInterface ):
          self.solCoreInterfaceParameters = map( float, content.pop(0).split('!')[0].split()[0:5] )

    if len(content) > 0:
      print 'There are still impurity parameter input data to be read.'

  def Write( self, name = None, transport = None ):
    """Writes the input file.

    This methods creates the input file based on the data stored in this 
    object.

    Keyword arguments.
    name -- File name to write data to. If None, read-in file will be 
            overwritten. (str, default None)
    transport -- If True, additioanl transport parameter are written. If 
                 None, input value is used (bool, default None)

    """
    if name is None:
      name = self.name
    if transport is None:
      transport = self.transport
    f = open( name, 'w' )
    f.write( '*** 9 Additial inputs for impurity transport\n' )
    f.write( '*\n' )
    f.write( '* This file was created automaticly.\n' )
    f.write( '* Use the EMC3 python tool to modify it.\n' )
    f.write( '*\n' )
    f.write( '** 9.1 Transport model\n' )
    f.write( '%d\n' % self.typeTransportModel )
    f.write( '** 9.2 atomic-process Database\n' )
    f.write( '*STRAHL, ADAS, USER(not yet available)\n' )
    f.write( '*---------1.  STRAHL ------------------------------------------\n' )
    f.write( '* In case of using STRAHL database, you need simply to put the|\n' )
    f.write( '* keyword \'STRAHL\' bellow. At the same time, the file \'STRAHL\'|\n' )
    f.write( '* containing the atomic data must be provided.                |\n' )
    f.write( '*---------2. ADAS ---------------------------------------------\n' )
    f.write( '* Data block structure for using ADAS                         |\n' )
    f.write( '* ADAS      ! keyword                                         |\n' )
    f.write( '* NFILE     ! totoal number of adas data files                |\n' )
    f.write( '* PATH&files, where a class name like \'scd\',\'acd\' .... must   |\n' )
    f.write( '* appear at least once for classifying the preocesses.        |\n' )
    f.write( '* for example:                                                |\n' )
    f.write( '* ADAS      ! Call ADAS database                              |\n' )
    f.write( '* 3         ! 3 totoal processes                              |\n' )
    f.write( '*/afs/ipp/u/adas/adas/adf11/scd96/scd96_c.dat                 |\n' )
    f.write( '*/afs/ipp/u/adas/adas/adf11/acd96/acd96_c.dat                 |\n' )
    f.write( '*/afs/ipp/u/adas/adas/adf11/plt96/plt96_c.dat                 |\n' )
    f.write( '*---------3. USER supplied ------------------------------------\n' )
    f.write( '* USER (not yet available                                     |\n' )
    f.wirte( '%s\n' % self.database )
    if self.database.find( 'ADAS' ) != -1:
      f.write( '%d\n' % len( self.databaseFiles ) )
      for databaseFile in databaseFiles:
        f.write( '%s\n' % databaseFile )

    if transport:
      f.write( '** 9.3 impurity sources\n' )
      f.write( '%g %g %g\n' % ( self.physicalSputtering, 
                                self.chemicalSputtering, 
                                self.selfSputtering ) )
      f.write( '%g %d %d\n' % ( self.energy, self.alpha, self.beta ) )
      f.write( '%d %g\n' % ( len( self.sourceParameters ), self.betaSources ) )
      for sourceParameter in self.sourceParameters:
        f.write( '%g %g %g %g %g %g %g %g\n' % tuple( sourceParameter ) )
      f.write( '%d\n' % self.typeSolCoreInterface )
      if self.typeSolCoreInterface > 1:
        self.solCoreInterfaceParameters.sort()
        for solCoreInterfaceParameter in self.solCoreInterfaceParameters:
          f.write( '%g %g %g %g %g\n' % tuple( solCoreInterfaceParameter ) )

    f.close()

