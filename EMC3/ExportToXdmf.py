from xml.etree.cElementTree import Element, SubElement, ElementTree, parse
import h5py
from numpy import ravel_multi_index, ones, zeros, cos, sin, pi, int8, array
from Grid import Grid
from PlasmaField import PlasmaField

class ExportToXdmf( object ):
  """Class that describes the export into the XDMF format.

  This class describes the export of grid information and plasma field 
  values into the eXtensible Data Model and Format (XDMF). XDMF is for 
  instance readably by Paraview.

  """

  def __init__( self, filename, exportAll = False ):
    """Constructor

    Keyword arguments:
    filename -- Name of output files. (str)
    exportAll -- Whether or not all known plasma fields are exported at 
                 once. Only works if no file with given file name is 
                 present. (bool, default False)

    """
    self.filename = filename
    
    try:
      self.tree = parse( '%s.xmf' % self.filename )
      exportAll = False
    except:
      Xdmf = Element( 'Xdmf' )
      Xdmf.set( 'xmlns:xi', 'http://www.w3.org/2001/XInclude' )
      Xdmf.set( 'Version', '2.1' )
      Domain = SubElement( Xdmf, 'Domain' )
      self.tree = ElementTree( Xdmf )
    
    if exportAll:
      self.AddGrid()
      for name in ( 'TE_TI', 'DENSITY', 'MACH_NUMBER', 'SOURCE_RAD_E', 
                    'ALFA_SOURCE' ):
        try:
          PF = PlasmaField( name )
        except IOError, message:
          print message
        else:
          self.AddPlasmaField( PF )
      for name in ( 'SOURCE_P', 'SOURCE_E_E', 'SOURCE_E_I', 'SOURCE_M', 
                    'DENSITY_A', 'DENSITY_I', 'DENSITY_M', 'DENSITY_E_A', 
                    'DENSITY_E_I', 'DENSITY_E_M' ):
        try:
          PF = PlasmaField( name, path='EIRENE_OUTPUT' )
        except IOError, message:
          print message
        else:
          self.AddPlasmaField( PF )

  def AddGrid( self, grid = None, name = 'Grid', realCoordinates = True ):
    """Method to add a grid.

    This method adds a grid to the XDMF files.

    Keyword arguments:
    grid -- Grid. (Grid object, if None trying to initialize it with data 
            in current directory)
    name -- Name that describes the grid. (str, default 'Grid')
    realCoordinates -- Switch for real coordinates. If False, X = R and 
                       Y = Phi are chosen. (bool, default True)

    """
    if self.tree.find( 'Domain/Grid[@Name=\'%s\']' % name ) is None:
      if grid is None:
        grid = Grid()

      h = h5py.File( '%s.h5' % self.filename, 'a' )
      GridCollection = SubElement( self.tree.find( 'Domain' ), 'Grid', 
                                   Name = name, 
                                   GridType = 'Collection', 
                                   CollectionType = 'Spatial' )

      zoneIndex = 0
      for zone in grid.zones:
        x, y, z = self.__TransformCoordinates( zone, realCoordinates )
        topology = self.__GetHexahedralTopology( x.shape )

        GridElement = SubElement( GridCollection, 'Grid', 
                                  Name = 'Zone %d' % zoneIndex )
        Geometry = SubElement( GridElement, 'Geometry',
                               Type = 'X_Y_Z' )
        path = '%s/Zone%d' % ( name, zoneIndex )

        dset = h.create_dataset( '%s/Geometry/X' % path, data = x )
        DataItem = SubElement( Geometry, 'DataItem',
                               DataType = 'Float',
                               Dimensions = '%d %d %d' % x.shape, 
                               Format = 'HDF',
                               Precision = '8' )
        DataItem.text = '%s.h5:%s/Geometry/X' % ( self.filename, path )

        dset = h.create_dataset( '%s/Geometry/Y' % path, data = y )
        DataItem = SubElement( Geometry, 'DataItem',
                               DataType = 'Float',
                               Dimensions = '%d %d %d' % y.shape, 
                               Format = 'HDF',
                               Precision = '8' )
        DataItem.text = '%s.h5:%s/Geometry/Y' % ( self.filename, path )

        dset = h.create_dataset( '%s/Geometry/Z' % path, data = z )
        DataItem = SubElement( Geometry, 'DataItem',
                               DataType = 'Float',
                               Dimensions = '%d %d %d' % z.shape, 
                               Format = 'HDF',
                               Precision = '8' )
        DataItem.text = '%s.h5:%s/Geometry/Z' % ( self.filename, path )

        dset = h.create_dataset( '%s/Topology' % path, data = topology )
        Topology = SubElement( GridElement, 'Topology',
                               NumberOfElements = '%d' % topology[:,:,:,0].size,
                               Type = 'Hexahedron' )
        DataItem = SubElement( Topology, 'DataItem',
                               DataType = 'Int',
                               Dimensions = '%d %d %d %d' % topology.shape, 
                               Format = 'HDF' )
        DataItem.text = '%s.h5:%s/Topology' % ( self.filename, path )
        
        zoneIndex += 1

      h.close()
      self.__Write()
    else:
      print 'A grid named \'%s\' already exists.' % name

  def AddInstallation( self, installation, name = 'Installation', 
                       realCoordinates = True ):
    """Method to add an installation.

    This method adds an installation to the XDMF files.

    Keyword arguments:
    installation -- Installation. (Installation object)
    name -- Name that describes the installation. (str, default 'Installation')
    realCoordinates -- Switch for real coordinates. If False, X = R and 
                       Y = Phi are chosen. (bool, default True)

    """
    if self.tree.find( 'Domain/Grid[@Name=\'%s\']' % name ) is None:
      h = h5py.File( '%s.h5' % self.filename, 'a' )

      x, y, z = self.__TransformCoordinates( installation, realCoordinates )
      topology = self.__GetTriangleTopology( x.shape )

      GridElement = SubElement( self.tree.find( 'Domain' ), 'Grid', 
                                Name = name )
      Geometry = SubElement( GridElement, 'Geometry', Type = 'X_Y_Z' )

      dset = h.create_dataset( '%s/Geometry/X' % name, data = x )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % x.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/X' % ( self.filename, name )

      dset = h.create_dataset( '%s/Geometry/Y' % name, data = y )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % y.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/Y' % ( self.filename, name )

      dset = h.create_dataset( '%s/Geometry/Z' % name, data = z )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % z.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/Z' % ( self.filename, name )

      dset = h.create_dataset( '%s/Topology' % name, data = topology )
      Topology = SubElement( GridElement, 'Topology',
                             NumberOfElements = '%d' % topology[:,:,:,0].size,
                             Type = 'Triangle' )
      DataItem = SubElement( Topology, 'DataItem',
                             DataType = 'Int',
                             Dimensions = '%d %d %d %d' % topology.shape, 
                             Format = 'HDF' )
      DataItem.text = '%s.h5:%s/Topology' % ( self.filename, name )
        
      h.close()
      self.__Write()
    else:
      print 'An installation named \'%s\' already exists.' % name

  def AddTargetDeposition( self, depo, name = 'Target', 
                           realCoordinates = True ):
    """Method to add a target surface and the deposition flux onto this surface.

    This method adds a target surface and the deposition flux onto this 
    surface to the XDMF files.

    Keyword arguments:
    depo -- Deposition. (Depo object)
    name -- Name that describes the deposition target. The attribute name is
            chosen automatically based on the type of deposition. (str, 
            default 'Target')
    realCoordinates -- Switch for real coordinates. If False, X = R and 
                       Y = Phi are chosen. (bool, default True)

    """
    GridElement = self.tree.find( 'Domain/Grid[@Name=\'%s\']' % name )

    geometryMatch = True

    x, y, z = self.__TransformCoordinates( depo, realCoordinates )
    x = x[depo.isTarget == True]
    y = y[depo.isTarget == True]
    z = z[depo.isTarget == True]

    if GridElement is None:
      h = h5py.File( '%s.h5' % self.filename, 'a' )

      GridElement = SubElement( self.tree.find( 'Domain' ), 'Grid', 
                                Name = name )
      Geometry = SubElement( GridElement, 'Geometry', Type = 'X_Y_Z' )

      dset = h.create_dataset( '%s/Geometry/X' % name, data = x )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % x.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/X' % ( self.filename, name )

      dset = h.create_dataset( '%s/Geometry/Y' % name, data = y )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % y.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/Y' % ( self.filename, name )

      dset = h.create_dataset( '%s/Geometry/Z' % name, data = z )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % z.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/Z' % ( self.filename, name )

      Topology = SubElement( GridElement, 'Topology',
                             NumberOfElements = '%d' % x.shape[0],
                             Type = 'Quadrilateral' )

      h.close()
      self.__Write()
    else:
      Geometry = GridElement.find( 'Geometry' )
      xyzRead = []
      for DataItem in Geometry.getchildren():
        filename, path = DataItem.text.split(':')
        hRead = h5py.File( filename, 'r' )
        xyzRead.append( hRead.get( path ).value )
        hRead.close()
      if ( array(xyzRead) != array([x, y, z]) ).any():
        geometryMatch = False
        
    if geometryMatch:
      label = '%sFlux [%s]' % ( depo.name.capitalize(), depo.fluxUnit )

      if GridElement.find( 'Attribute[@Name=\'%s\']' % label ) is None:
        h = h5py.File( '%s.h5' % self.filename, 'a' )

        flux = depo.flux[depo.isTarget == True]

        path = '%s/Flux/%s' % ( name, depo.name.capitalize() )

        dset = h.create_dataset( path, data = flux )
        Attribute = SubElement( GridElement, 'Attribute',
                                Center = 'Cell',
                                Name = label,
                                Type = 'Scalar' )
        DataItem = SubElement( Attribute, 'DataItem',
                               DataType = 'Float',
                               Dimensions = '%d' % flux.size,
                               Format = 'HDF',
                               Precision = '8' )
        DataItem.text = '%s.h5:%s' % ( self.filename, path )

        h.close()
        self.__Write()
      else:
        print 'A \'%s\' deposition flux already exists for the target \'%s\'.' \
              % ( depo.name.capitalize(), name )
    else:
      print( 'A target named \'%s\' already exists, ' % name +
             'but does not match the given one.' )

  def AddMappingLoss( self, depo, name = 'MappingLoss', 
                      realCoordinates = True ):
    """Method to add the flux onto cells lost due to mapping.

    This method adds cell surfaces including the flux onto them which is lost
    due to not merging cells at mapping surfaces to the XDMF files.

    Keyword arguments:
    depo -- Deposition. (Depo object)
    name -- Tree name that describes the mapping loss deposition. The tree
            element will be named automatically based on type of deposition 
            flux ('Particle', 'Energy'). (str, default 'MappingLoss')
    realCoordinates -- Switch for real coordinates. If False, X = R and 
                       Y = Phi are chosen. (bool, default True)

    """
    GridTree = self.tree.find( 'Domain/Grid[@Name=\'%s\']' % name )
    if GridTree is None:
      GridTree = SubElement( self.tree.find( 'Domain' ), 'Grid', Name = name,
                             GridType = 'Tree' )
    if GridTree.find( 'Grid[@Name=\'%s\']' % depo.name.capitalize() ) is None:
      h = h5py.File( '%s.h5' % self.filename, 'a' )

      x, y, z = self.__TransformCoordinates( depo, realCoordinates )
      x = x[depo.isTarget == False]
      y = y[depo.isTarget == False]
      z = z[depo.isTarget == False]
      flux = depo.flux[depo.isTarget == False]

      GridElement = SubElement( GridTree, 'Grid', 
                                Name = depo.name.capitalize() )
      Geometry = SubElement( GridElement, 'Geometry', Type = 'X_Y_Z' )
      path = '%s/%s' % ( name, depo.name.capitalize() )

      dset = h.create_dataset( '%s/Geometry/X' % path, data = x )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % x.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/X' % ( self.filename, path )

      dset = h.create_dataset( '%s/Geometry/Y' % path, data = y )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % y.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/Y' % ( self.filename, path )

      dset = h.create_dataset( '%s/Geometry/Z' % path, data = z )
      DataItem = SubElement( Geometry, 'DataItem',
                             DataType = 'Float',
                             Dimensions = '%d %d' % z.shape, 
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Geometry/Z' % ( self.filename, path )

      Topology = SubElement( GridElement, 'Topology',
                             NumberOfElements = '%d' % x.shape[0],
                             Type = 'Quadrilateral' )
        
      dset = h.create_dataset( '%s/Flux' % path, data = flux )
      Attribute = SubElement( GridElement, 'Attribute',
                              Center = 'Cell',
                              Name = 'Lost%sFlux [%s]' 
                                     % ( depo.name.capitalize(), 
                                         depo.fluxUnit ), 
                              Type = 'Scalar' )
      DataItem = SubElement( Attribute, 'DataItem',
                             DataType = 'Float',

                             Dimensions = '%d' % flux.size,
                             Format = 'HDF',
                             Precision = '8' )
      DataItem.text = '%s.h5:%s/Flux' % ( self.filename, path )

      h.close()
      self.__Write()
    else:
      print( 'A mapping loss deposition named ' + 
             '\'%s\' already exists in the tree \'%s\'.' 
             % ( depo.name.capitalize(), name ) )

  def AddTargetProfiles( self, targetProfiles, name = 'TargetProfiles', 
                         realCoordinates = True ):
    """Method to add a target surface and the deposition profiles onto this 
    surface.

    This method adds a target surface and the deposition profiles onto this 
    surface to the XDMF files.

    Keyword arguments:
    targetProfiles -- Target profiles. (TargetProfiles object)
    name -- Name that describes the target profiles. The attribute name is
            chosen automatically based on the type of profile. (str, 
            default 'Target Profiles')
    realCoordinates -- Switch for real coordinates. If False, X = R and 
                       Y = Phi are chosen. (bool, default True)

    ATTENTION:
    Option realCoordinates is not yet supported

    """
    if self.tree.find( 'Domain/Grid[@Name=\'%s\']' % name ) == None:
      h = h5py.File( '%s.h5' % self.filename, 'a' )
      GridCollection = SubElement( self.tree.find( 'Domain' ), 'Grid', 
                                   Name = name, 
                                   GridType = 'Collection', 
                                   CollectionType = 'Spatial' )

      for targetProfilesElement in targetProfiles:
        x, y, z = self.__TransformCoordinates( targetProfilesElement, realCoordinates )

        GridElement = SubElement( GridCollection, 'Grid', 
                                  Name = 'Element %d' % targetProfilesElement.elementId )
        Geometry = SubElement( GridElement, 'Geometry',
                               Type = 'X_Y_Z' )
        path = '%s/Element%d' % ( name, targetProfilesElement.elementId )

        dset = h.create_dataset( '%s/Geometry/X' % path, data = x )
        DataItem = SubElement( Geometry, 'DataItem',
                               DataType = 'Float',
                               Dimensions = '%d %d' % x.shape, 
                               Format = 'HDF',
                               Precision = '8' )
        DataItem.text = '%s.h5:%s/Geometry/X' % ( self.filename, path )

        dset = h.create_dataset( '%s/Geometry/Y' % path, data = y )
        DataItem = SubElement( Geometry, 'DataItem',
                               DataType = 'Float',
                               Dimensions = '%d %d' % y.shape, 
                               Format = 'HDF',
                               Precision = '8' )
        DataItem.text = '%s.h5:%s/Geometry/Y' % ( self.filename, path )

        dset = h.create_dataset( '%s/Geometry/Z' % path, data = z )
        DataItem = SubElement( Geometry, 'DataItem',
                               DataType = 'Float',
                               Dimensions = '%d %d' % z.shape, 
                               Format = 'HDF',
                               Precision = '8' )
        DataItem.text = '%s.h5:%s/Geometry/Z' % ( self.filename, path )

        Topology = SubElement( GridElement, 'Topology',
                               NumberOfElements = '%d' % x.shape[0],
                               Type = 'Quadrilateral' )
        
        for quantity in [ 'Gamma', 'q', 'n', 'Te', 'Ti' ]:
          label = '%s [%s]' % ( quantity, targetProfilesElement.profiles[quantity][0] )

          if GridElement.find( 'Attribute[@Name=\'%s\']' % label ) is None:
            dset = h.create_dataset( '%s/Profiles/%s' % ( path, quantity ), 
                                     data = targetProfilesElement.profiles[quantity][1] )
            Attribute = SubElement( GridElement, 'Attribute',
                                    Center = 'Cell',
                                    Name = label,
                                    Type = 'Scalar' )
            DataItem = SubElement( Attribute, 'DataItem',
                                   DataType = 'Float',
                                   Dimensions = '%d' % targetProfilesElement.profiles[quantity][1].size,
                                   Format = 'HDF',
                                   Precision = '8' )
            DataItem.text = '%s.h5:%s/Profiles/%s' % ( self.filename, path, 
                                                       quantity )
          else:
            print( 'A \'%s\' target profile already exists ' % quantity 
                   + 'for the target \'%s\' (element %d).' 
                   % ( name, targetProfilesElement.elementId ) )
      h.close()
      self.__Write()
    else:
      print 'Target profiles named \'%s\' already exists.' % name

  def __TransformCoordinates( self, gridlike, realCoordinates ):
    """Transforms coordinates into the required form.

    This method transforms the coordinates stored in the grid like object 
    into the required format. If real coordinates are not requested X = R 
    and Y = Phi is chosen.

    Keyword arguments:
    gridlike -- Grid like object containing R, Z, phi coordinates to 
                transform. (Grid/Installation/Depo/TargetProfilesElement)
    realCoordinates -- Switch for real coordinates. If False, X = R and 
                       Y = Phi are chosen. (bool, default False)

    """
    phi = gridlike.phi
    for k in xrange( gridlike.R.ndim - phi.ndim ):
      phi = phi[:,None]
    if not realCoordinates:
      x = gridlike.R
      y = ones( x.shape ) * phi
    else:
      x = gridlike.R * cos( phi / 180. * pi )
      y = gridlike.R * sin( phi / 180. * pi )
    z = gridlike.Z
    return x, y, z

  def __Write( self ):
    self.tree.write( '%s.xmf' % self.filename, encoding = 'utf-8', 
                     xml_declaration = True )

  def __CellTopology( self, nToroidal, nPoloidal, nRadial, shape ):
    """Method that returns the topology/connectivity of the cell. (list)

    This method constructs the topology/connectivity of the cell and 
    returns it as a list of raveled indixes.

    Keyword arguments:
    nToroidal -- Index in toroidal direction. (int)
    nPoloidal -- Index in poloidal direction. (int)
    nRadial -- Index in radial direction. (int)
    shape -- Shape of the grid array. (tuple of ints)

    """
    surface1 = [ (nToroidal, nPoloidal, nRadial), 
                 (nToroidal, nPoloidal+1, nRadial), 
                 (nToroidal, nPoloidal+1, nRadial+1), 
                 (nToroidal, nPoloidal, nRadial+1) ]
    surface2 = [ (nToroidal+1, nPoloidal, nRadial), 
                 (nToroidal+1, nPoloidal+1, nRadial), 
                 (nToroidal+1, nPoloidal+1, nRadial+1), 
                 (nToroidal+1, nPoloidal, nRadial+1) ]
    cell = surface1 + surface2
    cellRaveled = []
    for node in cell:
      cellRaveled.append( ravel_multi_index( node, shape ) )
    return cellRaveled

  def __TriangleTopology( self, nToroidal, nPoloidal, shape ):
    """Method that returns the topology/connectivity of the two triangles of 
    on limiter/divertor surface in the Kisslinger format. (list)

    This method constructs the topology/connectivity of the two triangle that 
    result from the limiter/divertor surface given in the Kisslinger format 
    and returns it as a list of raveled indixes.

    Keyword arguments:
    nToroidal -- Index in toroidal direction. (int)
    nPoloidal -- Index in poloidal direction. (int)
    shape -- Shape of the grid array. (tuple of ints)

    """
    triangle1 = [ (nToroidal, nPoloidal), 
                  (nToroidal+1, nPoloidal+1), 
                  (nToroidal+1, nPoloidal) ]
    triangle2 = [ (nToroidal, nPoloidal), 
                  (nToroidal, nPoloidal+1), 
                  (nToroidal+1, nPoloidal+1) ]
    triangle1Raveled = []
    for node in triangle1:
      triangle1Raveled.append( ravel_multi_index( node, shape ) )
    triangle2Raveled = []
    for node in triangle2:
      triangle2Raveled.append( ravel_multi_index( node, shape ) )
    return triangle1Raveled, triangle2Raveled

  def __GetHexahedralTopology( self, shape ):
    """Method that returns the topology of a hexahedral grid.

    This method constructs and returns the topology of a hexahedral grid.

    Keyword arguments:
    shape -- Shape of the grid array. (tuple of ints)

    ATTENTION:
    Slow! Needs to be speeded-up by using vector operations.

    """
    topology = zeros( ( shape[0]-1, shape[1]-1, shape[2]-1, 8 ), dtype=int )
    for nToroidal in xrange( topology.shape[0] ):
      for nPoloidal in xrange( topology.shape[1] ):
        for nRadial in xrange( topology.shape[2] ):
          topology[nToroidal, 
                   nPoloidal, 
                   nRadial, :] = self.__CellTopology( nToroidal, 
                                                      nPoloidal, 
                                                      nRadial, shape )
    return topology

  def __GetTriangleTopology( self, shape ):
    """Method that returns the topology of a triangle grid.

    This method constructs and returns the topology of a triangle grid.

    Keyword arguments:
    shape -- Shape of the grid array. (tuple of ints)

    ATTENTION:
    Slow! Needs to be speeded-up by using vector operations.

    """
    topology = zeros( ( shape[0]-1, shape[1]-1, 2, 3 ), dtype=int )
    for nToroidal in xrange( topology.shape[0] ):
      for nPoloidal in xrange( topology.shape[1] ):
        trianglesRaveled = self.__TriangleTopology( nToroidal, nPoloidal, 
                                                    shape )
        topology[nToroidal, nPoloidal, 0, :] = trianglesRaveled[0]
        topology[nToroidal, nPoloidal, 1, :] = trianglesRaveled[1]
    return topology

  def AddPlasmaField( self, plasmaField, quantity = None, gridName = 'Grid' ):
    """Method to add plasma field quantities.

    This method adds plasma field quantities to the XDMF files.

    Keyword arguments:
    plasmaField -- Plasma field object. (PlasmaField)
    quantity -- Quantity of plasma field object to be added. If None given
                all available are added. (str, default None)
    gridName -- Name of the grid to add the plasma field to. (str, default 
                'Grid')

    """
    if quantity is None:
      quantity = plasmaField.keys()

    GridCollection = self.tree.find( 'Domain/Grid[@Name=\'%s\']' % gridName )

    for q in quantity:
      label = q
      if ( plasmaField[q][0] != '' ):
        label += ' [%s]' % plasmaField[q][0]

      zoneIndex = 0
      for zone in plasmaField(q):
        GridElement = GridCollection.find( 'Grid[@Name=\'Zone %d\']' % zoneIndex )
        if GridElement.find( 'Attribute[@Name=\'%s\']' % label ) is None:
          h = h5py.File( '%s.h5' % self.filename, 'a' )

          path = '%s/Zone%d' % ( gridName, zoneIndex )
          dset = h.create_dataset( '%s/PlasmaField/%s' % ( path, q ), 
                                   data = zone )
          Attribute = SubElement( GridElement, 'Attribute',
                                  Center = 'Cell',
                                  Name = label, 
                                  Type = 'Scalar' )
          DataItem = SubElement( Attribute, 'DataItem',
                                 DataType = 'Float',
                                 Dimensions = '%d %d %d' % ( zone.shape ),
                                 Format = 'HDF',
                                 Precision = '8' )
          DataItem.text = '%s.h5:%s/PlasmaField/%s' % ( self.filename, path, q )

          h.close()
          self.__Write()
        else:
          print( 'An attribute named ' +
                 '\'%s\' already exists in zone %d.' % ( label, zoneIndex ) )
        zoneIndex += 1

  def AddCellProperty( self, cellProperty, name = 'CellProperty', gridName = 'Grid' ):
    """Method to add cell properties.

    This method adds cell properties to the XDMF files.

    Keyword arguments:
    cellProperty -- Cell property object. (PlateCells or PhysicalCells)
    name -- Name that describes the cell property. (str, default 
            'CellProperty')
    gridName -- Name of the grid to add the cell property to. (str, default 
                'Grid')

    """
    GridCollection = self.tree.find( 'Domain/Grid[@Name=\'%s\']' % gridName )

    zoneIndex = 0
    for zone in cellProperty.ids:
      GridElement = GridCollection.find( 'Grid[@Name=\'Zone %d\']' % zoneIndex )
      if GridElement.find( 'Attribute[@Name=\'%s\']' % name ) is None:
        h = h5py.File( '%s.h5' % self.filename, 'a' )

        path = '%s/Zone%d' % ( gridName, zoneIndex )
        if zone.dtype == bool:
          zone = int8( zone )
        dset = h.create_dataset( '%s/%s' % ( path, name ), data = zone )
        Attribute = SubElement( GridElement, 'Attribute',
                                Center = 'Cell',
                                Name = name, 
                                Type = 'Scalar' )
        DataItem = SubElement( Attribute, 'DataItem',
                               DataType = 'Int',
                               Dimensions = '%d %d %d' % ( zone.shape ),
                               Format = 'HDF',
                               Precision = '%d' % zone.itemsize )
        DataItem.text = '%s.h5:%s/%s' % ( self.filename, path, name )

        h.close()
        self.__Write()
      else:
        print( 'An attribute for the cell property ' + 
               '\'%s\' already exists in zone %d.' % ( name, zoneIndex ) )
      zoneIndex += 1

