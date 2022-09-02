from fortranformat import FortranRecordReader, FortranRecordWriter

"""Class that represents the plasma parameters.

This class represents the input of EMC3 for plasma parameters, boundary and 
initial conditions. On initialization it reads an existing input file. 
After modifications of the parameter a new input file can be written.

"""

class PlasmaParameters( object ):

  def __init__( self, name = 'input.par', path = '../run' ):
    """Constructor

    Keyword arguments:
    name -- Name of the parameter input. (str, default 'input.par')
    path -- Path to the file. (str, default '../run')

    """
    self.name = name
    self.path = path.rstrip('/')
    
    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content = filter( lambda s: '*' != s[0], content )

    NIonSpecies = int( content.pop(0).split('!')[0].split()[0] )
    self.ionSpecies = []
    for k in xrange( NIonSpecies ):
      line = content.pop(0).split('!')[0].split()
      self.ionSpecies.append( { 'name': line[0],
                                'excitationStages': int(line[1]),
                                'atomicMass': float(line[2]) } )
    
    self.crossFieldDiffusion = [] # in cm^2/s
    for k in xrange( NIonSpecies ):
      self.crossFieldDiffusion.append( float( content.pop(0).split('!')[0].split()[0] ) )

    line = content.pop(0).split('!')[0].split()
    self.crossFieldConductivityElectron = float( line[0] ) # in cm^2/s
    self.crossFieldConductivityIon = float( line[1] ) # in cm^2/s

    self.boundariesParticle = []
    for k in xrange( NIonSpecies ):
      line = content.pop(0).split('!')[0].split()
      boundaryType = int(line[0])
      coefficient = float(line[1])
      if k == 0 and boundaryType > 10:
        self.RSeparatrix = int( boundaryType / 10 )
        boundaryType %= 10
      boundaries = self._ReadBoundaries( content )
      self.boundariesParticle.append( { 'type': boundaryType, 
                                        'coefficient': coefficient, 
                                        'boundaries': boundaries } )

    self.boundariesEnergy = self._ReadBoundaries( content, 2 )

    self.boundariesMomentum = self._ReadBoundaries( content )

    self.inputFileUnits = {}

    line = self._ReadInputFileUnits( content, 'energy' )
    self.initialEnergyElectron = float( line[1] ) # in eV
    self.initialEnergyIon = float( line[2] ) # in eV

    self.initialDensity = float( self._ReadInputFileUnits( content, 'density' )[1] ) # in cm^-3

    self._ReadInputFileUnits( content, 'machNumber' )
    self._ReadInputFileUnits( content, 'particleVolumeSource' )
    self._ReadInputFileUnits( content, [ 'electronEnergySourceNeutrals',
                                         'electronEnergySourceImpurities' ] )
    self._ReadInputFileUnits( content, [ 'ionEnergySourceNeutrals',
                                         'ionEnergySourceImpurities' ] )
    self._ReadInputFileUnits( content, [ 'momentumSourceTransport',
                                         'momentumSourceChargeExchange' ] )

    if len(content) > 0:
      print 'There are still plasma parameter input data to be read.'

  def _ReadBoundaries( self, content, NSurfaceData = 1 ):
    """Reads the boundary information.

    This method reads the boundary information and returns in a dict for 
    the radial, poloidal and toroidal orientation a list of surfaces 
    containing a dict with the index of the corresponding non-transparent 
    surface, the type of the surface and the data for this surface.
    
    Type of surfaces, meaning of surface data: 
    -2: decay, length in cm, 
    -1: sink, negative flux through boundary surface in W (energy)/A (particle)/? (momentum), 
     1: source, positive flux through boundary surface in W (energy)/A (particle)/? (momentum) 

    Keyword arguments:
    content -- Content of the input file. (list of str)
    NSurfaceData -- Number of surface data to read. (int, default 1)

    """
    NSurfaces = int( content.pop(0).split('!')[0].split()[0] )
    orientations = [ 'radial', 'poloidal', 'toroidal' ]
    boundaries = {}
    for orientation in orientations:
      boundaries[ orientation ] = []
    for k in xrange( NSurfaces ):
      line = content.pop(0).split('!')[0].split()
      boundaries[ orientations[int(line[0])-1] ].append( { 'indexNonTransparentSurface': int(line[1]),
                                                         'type': int(line[2]), 
                                                         'data': map( float, line[3:3+NSurfaceData] ) } )
    return boundaries

  def _WriteBoundaries( self, fileHandle, boundaries ):
    """Writes the boundaries into the given file.

    This methods writes the boundaries into the given file.

    Keyword arguements:
    fileHandle -- File to write into. (file)
    boundaries -- Boundary conditions. (dict of list of dict)

    """
    orientations = [ 'radial', 'poloidal', 'toroidal' ]
    N = 0
    for orientation in orientations:
      N += len( boundaries[orientation] )
    fileHandle.write( '%d\n' % N )
    for k in xrange( len( orientations ) ):
      for b in boundaries[orientations[k]]:
        fileHandle.write( '%d %d %d' % ( k+1, 
                                         b['indexNonTransparentSurface'], 
                                         b['type'] ) )
        for data in b['data']:
          fileHandle.write( ' %g' % data )
        fileHandle.write( '\n' )
                
  def _ReadInputFileUnits( self, content, names ):
    """Reads input file units. (list of str)

    This methods reads one line from the content and extraces input file 
    units from it to add them to the dictionary 'InputFileUnits' based on 
    the given names.
    If the input file unit is less or equal zero it is ignored.
    The read line is return as a list of str.

    Keyword arguments:
    content -- Content of the input file. (list of str)
    names -- Names that describe the input file. ((list of) str)

    """
    if type(names) == str:
      names = [ names ]
    line = content.pop(0).split('!')[0].split()
    for k in xrange( len( names ) ):
      if int(line[k]) > 0:
        self.inputFileUnits[ names[k] ] = int(line[k])
    return line

  def _WriteInputFileUnits( self, fileHandle, names, additionalData = [] ):
    """Writes the input file units to the given file.

    This method writes the input file units to the given file as well as 
    additionally given floating points number(s).

    Keyword arguments:
    fileHandle -- File to write into. (file)
    names -- Names that describe the input file. ((list of) str)
    additionalData -- Additional float(s) to write. ((list of) float, 
                      default [])
    
    """
    if type(names) == str:
      names = [ names ]
    if type(additionalData) == float:
      additionalData = [ additionalData ]
    for name in names:
      fileHandle.write( '%d ' % self.inputFileUnits.get( name, 0 ) )
    for data in additionalData:
      fileHandle.write( '%g ' % data )
    fileHandle.write( '\n' )

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
    f.write( '* plasma parameters and boundary and initial conditions\n' )
    f.write( '*\n' )
    f.write( '* This file was created automaticly.\n' )
    f.write( '* Use the EMC3 python tool to modify it.\n' )
    f.write( '*\n' )
    f.write( '*** 4. ION SPECIES\n' )
    f.write( '** Nr.1 must be the main plasma ion\n' )
    f.write( '%d\n' % len( self.ionSpecies ) )
    for ionSpecies in self.ionSpecies:
      Fformat = FortranRecordWriter( '(A6,I6,F12.5)' )
      line = [ ionSpecies['name'], 
               ionSpecies['excitationStages'],
               ionSpecies['atomicMass'] ]
      f.write( '%s\n' % Fformat.write( line ) )

    f.write( '*** 5. transport parameters\n' )
    for crossFieldDiffusion in self.crossFieldDiffusion:
      f.write( '%g\n' % crossFieldDiffusion )
    f.write( '%g %g\n' % ( self.crossFieldConductivityElectron, 
                           self.crossFieldConductivityIon ) )

    f.write( '*** 6. boundary conditions\n' )
    f.write( '***6.1 particle transport for main ion + impurity\n' )
    for boundariesParticles in self.boundariesParticle:
      boundaryType = boundariesParticles['type']
      if hasattr( self, 'RSeparatrix' ):
        boundaryType += self.RSeparatrix * 10
      f.write( '%d %g\n' % ( boundaryType, boundariesParticles['coefficient'] ) )
      self._WriteBoundaries( f, boundariesParticles['boundaries'] )

    f.write( '***6.2 energy transport for e+i\n' )
    self._WriteBoundaries( f, self.boundariesEnergy )

    f.write( '***6.3 momentum transport\n' )
    self._WriteBoundaries( f, self.boundariesMomentum )

    f.write( '*** 7. intitial conditions\n' )
    self._WriteInputFileUnits( f, 'energy', 
                               [ getattr( self, 'initialEnergyElectron', 0.0 ),
                                 getattr( self, 'initialEnergyIon', 0.0 ) ] )
    self._WriteInputFileUnits( f, 'density', 
                               getattr( self, 'initialDensity', 0.0 ) )
    self._WriteInputFileUnits( f, 'machNumber' ) 

    f.write( '*** 8. Volume sources\n' )
    self._WriteInputFileUnits( f, 'particleVolumeSource' )
    self._WriteInputFileUnits( f, [ 'electronEnergySourceNeutrals',
                                    'electronEnergySourceImpurities' ] )
    self._WriteInputFileUnits( f, [ 'ionEnergySourceNeutrals',
                                    'ionEnergySourceImpurities' ] )
    self._WriteInputFileUnits( f, [ 'momentumSourceTransport',
                                    'momentumSourceChargeExchange' ] )
    f.close()

