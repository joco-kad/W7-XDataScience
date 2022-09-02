from fortranformat import FortranRecordReader
import matplotlib.pyplot as plt

class Info( dict ):
  """Class that represents the info output of the EMC3-EIRENE code. (dict)
  
  This class can read the info output of streaming, energy, impurity and 
  neutral. The class inherits from a dictonary and stores the different 
  info parameter as key, value pairs.
  Additionally, a plot method is implemented.
  
  """
  
  def __init__( self, name, path = 'EMC3_OUTPUT' ):
    """Constructor

    Keyword arguments:
    name -- Name of the info to read. (str in {'STREAMING', 'ENERGY', 
            'IMPURITY', 'NEUTRAL'})
    path -- Path to the file. (str, default 'EMC3_OUTPUT')

    """
    self.name = name
    self.path = path.rstrip('/')

    Formats = { 
      'STREAMING': ( '(F6.2,F6.3,1p5E11.3)', 
                     [ r'$\delta n \; (\%)$', r'$\delta M$', 
                       r'$\Gamma_\mathrm{total} \; (\mathrm{A})$', 
                       r'$n_\mathrm{e, up} \; (\mathrm{cm}^{-3})$', 
                       r'$n_\mathrm{e, down, left} \; (\mathrm{cm}^{-3})$', 
                       r'$<n_\mathrm{e, down}> \; (\mathrm{cm}^{-3})$',
                       r'$n_\mathrm{e, down, right} \; (\mathrm{cm}^{-3})$' ] ),
      'ENERGY': ( '(F6.1,1p4E12.4,/,F6.1,1p4E12.4,/,18x,1p3E12.4)', 
                  [ r'$\delta T_\mathrm{e} \; (\%)$', 
                    r'$T_\mathrm{e, up} \; (\mathrm{eV})$', 
                    r'$T_\mathrm{e, down, left} \; (\mathrm{eV})$', 
                    r'$<T_\mathrm{e, down}> \; (\mathrm{eV})$', 
                    r'$T_\mathrm{e, down, right} \; (\mathrm{eV})$', 
                    r'$\delta T_\mathrm{i} \; (\%)$', 
                    r'$T_\mathrm{i, up} \; (\mathrm{eV})$', 
                    r'$T_\mathrm{i, down, left} \; (\mathrm{eV})$', 
                    r'$<T_\mathrm{i, down}> \; (\mathrm{eV})$', 
                    r'$T_\mathrm{i, down, right} \; (\mathrm{eV})$', 
                    r'$P_\mathrm{loss, neutral} \; (\mathrm{W})$',
                    r'$P_\mathrm{loss, radiation} \; (\mathrm{W})$',
                    r'$P_\mathrm{loss, target} \; (\mathrm{W})$' ] ),
      'IMPURITY': ( '(1p2E12.4)', 
                    [ 'CRAD', 'PRAD' ] ), 
      'NEUTRAL': ( '(1p6E12.4)', 
                   [ 'SP_MAIN_NEW', 'SP', 'SE_SUM_N0_E', 'SE_SUM_N0_I', 
                     'SM_SUMP_ALL', 'SM_SUMN_ALL' ] ) }

    file = open( '%s/%s_INFO' % (self.path, self.name) )
    Fformat = FortranRecordReader( Formats[name][0] )
    content = file.read().splitlines()
    file.close()

    dlines = Fformat.format.count('/') + 1

    values = []
    for k in xrange( 0, len(content), dlines ):
      line = content[k]
      for l in xrange( 1, dlines ):
        line += '\n' + content[k+l]
      values.append( Fformat.read( line ) )

    self.update( dict( zip( Formats[name][1], zip(*values) ) ) )

  def plot( self, parameter = '', xlim = [] ):
    """Plot method for a quick view on the data.

    This is a plot method that provides a quick view on the data. If a 
    parameter is given only this parameter will be plotted, otherwise 
    all parameters are plotted in subplots. The x-range can be 
    restricted.

    Keyword arguments:
    parameter -- Parameter to plot. (str, default '')
    xlim -- Range of the x-axis. (list of int, default [])
    
    """
    if len(xlim) == 0:
      xlim.append( 0 )
    if len(xlim) == 1:
      xlim.append( len(self.values()[0]) )
    if parameter == '':
      fig, ax = plt.subplots( ncols = 2, 
                              nrows = int( 0.5 * len(self.keys()) + 0.5 ), 
                              sharex = True )
      ax = ax.flatten()
      k = 0
      for parameter in sorted( self.keys() ):
        ax[k].plot( self.get(parameter)[xlim[0]:xlim[1]] )
        ax[k].set_ylabel( parameter )
        k += 1
      fig.subplots_adjust( hspace = 0 )
      ax[-2].set_xlabel( '# iterations' )
      ax[-1].set_xlabel( '# iterations' )
      ax[-1].set_xlim( xlim[0], xlim[1] )
      fig.suptitle( self.name.lower() )
    else:
      plt.plot( self.get(parameter)[xlim[0]:xlim[1]] )
      plt.xlabel( '# iterations' )
      plt.xlim( xlim[0], xlim[1] )
      plt.ylabel( parameter )
      plt.title( self.name.lower() )
    plt.show()

