from numpy import zeros, where, int32, abs
import FortranReader
from PhysicalCells import PhysicalCells
from Plot import InteractivePoloidalCrossSection

class PlasmaField( dict ):
  """Class that represents the plasma field ouputs of EMC3-EIRENE. (dict)
  
  This class can read the field output of energy, streaming and impurity 
  and links it to the grid used. The plasma field is given for each cell.
  A plot method is implemented that enables a quick view on the data for 
  different toroidal angles.
  
  """

  def __init__( self, name, path = 'EMC3_OUTPUT', physicalCells = None ):
    """Constructor

    Keyword arguments:
    name -- Name of the plasma field file to read. (str in {'TE_TI', 
            'DENSITY', 'MACH_NUMBER', 'ALFA_SOURCE', 'SOURCE_P', 
            'SOURCE_E_E', 'SOURCE_E_I', 'SOURCE_RAD_E, 'SOURCE_M', 
            'DENSITY_A', 'DENSITY_I', 'DENSITY_M', 'DENSITY_E_A', 
            'DENSITY_E_I', 'DENSITY_E_M', 'IMPURITY_NEUTRAL', 
            'IMP_RADIATION'})
    path -- Path to the plasma field file. (str, default 'EMC3_OUTPUT')
    physicalCells -- Physical cells reference. (PhysicalCells object, if None
                    trying to initialize it with data in 'path')

    """
    self.name = name
    self.path = path.rstrip('/')

    outputs = {
      'TE_TI': ( [ 'Te', 'Ti' ], 'eV' ),
      'DENSITY': ( [ 'n' ], 'cm^-3' ),
      'MACH_NUMBER': ( [ 'M' ], '' ),
      'ALFA_SOURCE': ( [ 'ploss' ], 'A g cm/(s cm^3) ?' ),
      'SOURCE_P': ( [ 'Sion' ], 'A/cm^3' ),
      'SOURCE_E_E': ( [ 'See' ], 'W/cm^3' ),
      'SOURCE_E_I': ( [ 'Sei' ], 'W/cm^3' ),
      'SOURCE_RAD_E': ( [ 'Simp' ], 'W/cm^3 ?' ),
      'SOURCE_M': ( [ 'Sm', 'Smlin' ], 'A g cm/(s cm^3)' ), 
      'DENSITY_A': ( [ 'natm' ], 'cm^-3' ), 
      'DENSITY_I': ( [ 'nion' ], 'cm^-3' ), 
      'DENSITY_M': ( [ 'nmol' ], 'cm^-3' ), 
      'DENSITY_E_A': ( [ 'uatm' ], 'eV/cm^3' ), 
      'DENSITY_E_I': ( [ 'uion' ], 'eV/cm^3' ), 
      'DENSITY_E_M': ( [ 'umol' ], 'eV/cm^3' ), 
      'TEMPERATURE_A': ( [ 'Tatm' ], 'eV' ), 
      'TEMPERATURE_I': ( [ 'Tion' ], 'eV' ), 
      'TEMPERATURE_M': ( [ 'Tmol' ], 'eV' ),
      'IMPURITY_NEUTRAL': ( [ 'nimp' ], 'cm^-3' ), 
      'IMP_RADIATION': ( [ 'Prad' ], 'W/cm^3' ),
      'CONNECTION_LENGTH': ( [ 'Lc' ], 'cm' )
    }

    neutrals = ( 'DENSITY_A', 'DENSITY_I', 'DENSITY_M', 
                 'DENSITY_E_A', 'DENSITY_E_I', 'DENSITY_E_M',
                 'TEMPERATURE_A', 'TEMPERATURE_I', 'TEMPERATURE_M' )
    multispecies = ( 'DENSITY', 'SOURCE_P' ) + neutrals

    file = open( '%s/%s' % ( self.path, self.name ) )
    content = file.read().splitlines()
    file.close()
    content.reverse() # to allow a much more efficient pop from the end of the list

    if physicalCells is None:
      self.physicalCells = PhysicalCells()
    else:
      self.physicalCells = physicalCells

    if name not in neutrals:
      self.NCell = self.physicalCells.NCellPlasma
    else:
      self.NCell = self.physicalCells.NCellNeutral

    if name in multispecies:
      outputTemp = outputs[name][0][0]
      output = []
    else:
      output = outputs[name][0]

    valuesList = [] # list of plasma field values
    k = 0
    while len(content) > 0:
      for quantity in outputs[name][0]:
        if name in multispecies:
          if name in 'DENSITY':
            k = int( content.pop() )
          else:
            k += 1
          output.append( '%s%d' % ( outputTemp, k ) )
        values = FortranReader.List( content, self.NCell, '(1p6E12.4)' ) # plasma field values

        valuesList.append( ( outputs[name][1], values ) )

    self.update( dict( zip( output, valuesList ) ) )

  def __call__( self, quantity ):
    """Returns the data of the quantity in each zone. (list of ndarrays)

    The caller links the plasma field values on the physical (irregular) 
    grid to the cells of the geometrical (regular) gird and returns the 
    data of the requested quantity for each zone.

    Keyword arguments:
    quantity -- Physical quantity to return. (str)
    """
    fieldAllZones = []
    for ids in self.physicalCells.ids:
      field = zeros( ids.shape )
      indices = where( ids-1 < self.NCell )
      values = self[quantity][1]
      field[indices] = values[ (ids-1)[indices] ]
      fieldAllZones.append( field )

    return fieldAllZones

  def plot( self, quantity, grid, zoneIndex = 0, installations = [], xlim = [], 
            ylim = [], clim = [] ):
    """Plot method for a quick view on the plasma field values.

    This is a plot method that provides a quick view on the plasma field 
    values. A linear interpolation in toroidal direction is used for the 
    grid. 
    
    Keyword arguments:
    quantity -- Physical quantity to plot. (str)
    grid -- Grid that correspond to the quantity. (Grid object)
    zoneIndex -- Index of zone to plot. (int, default 0)
    installations -- Installations to be plotted. ((list of) Installation 
                     object(s))
    xlim -- Range of the x-axis. (list of int, default [])
    ylim -- Range of the y-axis. (list of int, default [])
    clim -- Range of the colourbar. No effect for quantity 'M' and 'ploss'.
            (list of int, default [])
    edgecolors -- Edge colour of mesh. (str, default 'grey')
    
    """
    colourBarLabel = quantity
    if ( self[quantity][0] != '' ):
      colourBarLabel += ' [%s]' % self[quantity][0]

    zone = grid.zones[zoneIndex]
    R = 0.5 * ( zone.R[:-1,:,:] + zone.R[1:,:,:] )
    Z = 0.5 * ( zone.Z[:-1,:,:] + zone.Z[1:,:,:] )
    phi = 0.5 * ( zone.phi[:-1] + zone.phi[1:] )

    data = self(quantity)[zoneIndex]

    if ( quantity in ['M', 'ploss'] ):
      cmap = 'RdGy'
      AbsMax = abs(data).max()
      clim = [-AbsMax,AbsMax] # to ensure that 0 gets the colour white
    else:
      cmap = 'hot_r'

    IP = InteractivePoloidalCrossSection( phi, R, Z, data, xlim, ylim, clim, cmap, colourBarLabel,
                                          installations = installations )

