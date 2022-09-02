from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Button, Slider
from matplotlib.pyplot import subplot, show, draw, colorbar
from numpy import where

"""Package for more advanced plot routines.

This package contains more advanced plot routines for providing a first 
view on the simulation results.

"""

class InteractivePoloidalCrossSection():
  """Class that describes an interactive plot for cell data in the poloidal 
  cross-section.

  This class describes an interactive plot for cell data in the poloidal 
  cross-section. 
  The class is designed with two buttons (previous, next) and a slider to 
  scroll through the poloidal cross-sections at the differenet toroidal 
  locations.

  """
    
  def __init__( self, phi, R = None, Z = None, data = None, 
                xlim = [], ylim = [], clim = [], cmap = 'hot_r', 
                colourBarLabel = '', levels = None,
                installations = [] ):
    """Constructor.
      
    Keyword arguments:
    phi -- Toroidal angle of poloidal cross-sections. (list of floats)
    R -- Radial position of grid nodes. (ndarray of floats, default None)
    Z -- Horizontal position of grid nodes. (ndarray floats, default None)
    data -- Cell data to plot. (ndarray of floats, default None)
    xlim -- Range of the x-axis. (list of ints, default [])
    ylim -- Range of the y-axis. (list of ints, default [])
    clim -- Range of the colourbar. (list of ints, default [])
    edgecolors -- Edge colour of mesh. (str, default 'grey')
    colourBarLabel -- Label of the colourbar. (str, default '')
    levels -- Level curves of contour plot to draw. (list of float, default 
              None)
    installations -- Installations to be plotted. ((list of) Installation 
                     object(s), default [])

    ATTENTION: 
    Interpolation of installations does not necessary fit the interpolation 
    of the grid

    """
    self.k = 0
    self.hColorbar = None
    self.R = R
    self.Z = Z
    self.phi = phi
    self.data = data
    self.xlim = xlim
    self.ylim = ylim
    self.clim = clim     
    self.cmap = cmap
    self.colourBarLabel = colourBarLabel
    self.levels = levels
    if type( installations ) is not list:
      installations = [ installations ]
    self.installations = installations

    gs = GridSpec( 2, 3, height_ratios = [14,1], width_ratios = [10,1,1],
                   hspace = 0.3, wspace = 0.1 )
    self.ax = subplot( gs[0,:] )
    axSlider = subplot( gs[1,0] )
    axPrevious = subplot( gs[1,1] )
    axNext = subplot( gs[1,2] )

    self.sAngle = Slider( axSlider, '', phi[0], phi[-1], phi[0], valfmt='' )
    bPrev = Button( axPrevious, '<' )
    bNext = Button( axNext, '>' )

    self.sAngle.on_changed( self.update )
    bPrev.on_clicked( self.previous )
    bNext.on_clicked( self.next )

    self.draw()
    show()

  def next( self, event ):
    """Event catcher of the next button.

    This method catches the event of the next button, increases the 
    toroidal position and redraw the figure.

    Keyword argument:
    event -- Event to catch. (event)

    """
    if self.k < len(self.phi)-1:
      self.k += 1
      self.draw()
      draw()

  def previous( self, event ):
    """Event catcher of the previous button.

    This method catches the event of the previous button, decreases the 
    toroidal position and redraw the figure.

    Keyword argument:
    event -- Event to catch. (event)

    """
    if self.k > 0:
      self.k -= 1
      self.draw()
      draw()

  def update( self, val ):
    """Event catcher of the slider.

    This method catches the event of the slider, calculates the index of the
    nearest toroidal position that correspond to the slider value and redraw 
    the figure.

    Keyword argument:
    val -- Value of the slider. (float)

    """
    k = abs( self.phi - val ).argmin()
    if k != self.k:
      self.k = k
      self.draw()
      draw()

  def draw( self ):
    """Draws the poloidal cross-section.

    This method draws the poloidal cross-section.

    """
    self.ax.cla()

    if self.R is not None and self.Z is not None and self.data is not None:
      if self.R.shape == self.data.shape:
        hPlot = self.ax.contour( self.R[self.k,:,:], self.Z[self.k,:,:],
                                 self.data[self.k,:,:], levels = self.levels )
      else:
        hPlot = self.ax.pcolor( self.R[self.k,:,:],
                                self.Z[self.k,:,:],
                                self.data[self.k,:,:])
      if self.data.min() != self.data.max():
        if self.hColorbar:
          self.hColorbar.ax.cla()
          self.hColorbar = colorbar( hPlot, cax=self.hColorbar.ax )
        else:
          self.hColorbar = colorbar( hPlot, ax=self.ax )
        self.hColorbar.set_label( self.colourBarLabel )
      hPlot.set_cmap( self.cmap )

    for installation in self.installations:
      ms = where( installation.phi == self.phi[self.k] )[0]
      if ms.size > 0:
        # toroidal angle of installation correspond to plot angle
        for m in ms:
          self.ax.plot( installation.R[m,:], installation.Z[m,:], 'k' )
      else:
        # find toroidal interval of installation that contains plot angle
        ms = where( ( ( installation.phi[:-1] < self.phi[self.k] ) &
                      ( installation.phi[1:] > self.phi[self.k] ) ) | 
                    ( ( installation.phi[:-1] > self.phi[self.k] ) &
                      ( installation.phi[1:] < self.phi[self.k] ) ) )[0]
        for m in ms:
          f = ( self.phi[self.k] - installation.phi[m] ) / \
              ( installation.phi[m+1] - installation.phi[m] )
          self.ax.plot( installation.R[m,:] + \
                        f * ( installation.R[m+1,:] - installation.R[m,:] ),
                        installation.Z[m,:] + \
                        f * ( installation.Z[m+1,:] - installation.Z[m,:] ),
                        'k' )

    self.ax.set_aspect( 'equal' )
    self.ax.set_title( '%6.3f degree' % self.phi[self.k] )
    self.ax.set_xlabel( 'R [cm]' )
    self.ax.set_ylabel( 'Z [cm]' )
    self.ax.set_xlim( *self.xlim )
    self.ax.set_ylim( *self.ylim )
    if self.clim != []:
      hPlot.set_clim( *self.clim )

    self.sAngle.set_val( self.phi[self.k] )

