"""Package containing functions to deal with surface definitions.

This package contains functions to deal with surface definitions in the 
parameter input files for EMC3. It contains functions to read and write 
the surface informations from and to input files.

"""

def Read( content ):
  """Function that reads surface information.

  This method reads the information on the e.g. non-default, non-
  transparent, or plate surfaces and returns the radial, poloidal and 
  toroidal surface informations as a dict of lists of dicts with the 
  required information that define the surfaces stored in the dictionary.

  The types of non-default surfaces are coded as follows: 
  1: periodic, 
  2: up/down symmetric, 
  3: Mapping.

  The types of non-transparent and plate surfaces are coded as follows:
  -1: ?
  1: ?
  
  Keyword arguments:
  content -- Content of the input file. (list of str)

  """
  directions = [ 'radial', 'poloidal', 'toroidal' ]
  allSurfaces = {}
  for direction in directions:
    label = filter( lambda s: s != direction, directions )
    N = int( content.pop(0).split('!')[0] )
    if N < 0:
      surfaces = N
    else:
      surfaces = []
      for k in xrange( N ):
        surface = {}
        info = map( int, content.pop(0).split('!')[0].split() )
        surface[ 'index' ] = info[0]
        surface[ 'zoneIndex' ] = info[1]
        surface[ 'typeOfSurface' ] = info[2]
        surface[ 'materialInfo' ] = info[3:8]
        indices = map( int, content.pop(0).split('!')[0].split() )
        surface[ '%sIndices' % label[0] ] = indices[:2]
        surface[ '%sIndices' % label[1] ] = indices[2:]
        surfaces.append( surface )
    allSurfaces[ direction ] = surfaces
  return allSurfaces

def Write( fileHandle, allSurfaces ):
  """Function that writes surface information.

  This method writes the information on the e.g. non-default, non-
  transparent, or plate surfaces in the parameter input file.

  Keyword arguments:
  fileHandle -- File handle to write data to. (file)
  allSurfaces -- Radial, poloidal and toroidal surface informations. 
                 (dict of lists of dicts)

  """
  directions = [ 'radial', 'poloidal', 'toroidal' ]
  for direction in directions:
    fileHandle.write( '* %s\n' % direction )
    label = filter( lambda s: s != direction, directions )
    if type( allSurfaces[ direction ] ) is int:
      if allSurfaces[ direction ] < 0:
        fileHandle.write( '%d\n' % allSurfaces[ direction ] )
      else:
        raise Exception( 'Unexpected integer value for surfaces in ' + \
                         '%s direction.' % direction )
    else:
      fileHandle.write( '%d\n' % len( allSurfaces[ direction ] ) )
      for surface in allSurfaces[ direction ]:
        fileHandle.write( '%d %d %d\n' % ( surface[ 'index' ], 
                                           surface[ 'zoneIndex' ], 
                                           surface[ 'typeOfSurface' ] ) )
        fileHandle.write( '%d %d ' % tuple( surface[ '%sIndices' % \
                                                     label[0] ] ) )
        fileHandle.write( '%d %d\n' % tuple( surface[ '%sIndices' % \
                                                      label[1] ] ) )

