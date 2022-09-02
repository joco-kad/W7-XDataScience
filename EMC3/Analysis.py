from numpy import diff, append, sqrt, arctan2, radians, mod, pi, linspace, \
                  where, abs, dot, roll, array, isnan, ones, int
from matplotlib.path import Path
from PlateCells import PlateCells
import sys
from Installation import Installation

"""Contains functions for analysis of EMC3 results.

Contains functions that support the user during the analysis of EMC3 
simulation results.

"""

def getFluxTubeIndices( i, j, k, Grid, Plates, direction=1 ,phiSegment = 36.0 ):
  """ Returns the indices that belong to one flux tube. [list of three lists, 
  two ndarrys]

  This function returns the indices that belong to on flux tube starting at 
  (i,j,k). Additionally, the toroidal angle and radius along the flux tube is 
  returned.

  Keyword arguments:
  i -- Start index toroidal. (int)
  j -- Start index poloidal. (int)
  k -- Start index radial. (int)
  Grid -- Grid to use. (Grid object)
  Plates -- Plates to consider. (PlateCells object)
  phiSegment -- Toroidal angle of one segment. (float, default 36.0)
  
  TODO:
  Need to be extendet to search in both direction and detect closed loops 
  automatically. Need to check for plates along the flight.
  
  A value for the direction of iss has to be implemented, for the evaluation of 
  Machnumber profiles. 
  """
  print direction
  counter=1  
  phis = Grid.zones[0].phi[:-1] + 0.5 * diff( Grid.zones[0].phi )
  phiss = phis[i:]
  Rs = Grid.zones[0].R[:,j:j+2,k:k+2].mean((1,2))
  Rss = Rs[:-1] + 0.5 * diff( Rs )
  Rss = Rss[i:]
  phiOffset = 0
  found = True
  if i < 0:
    iss = range( i, 0 )
    isd = []
    isd.extend(-1*ones(len(iss),dtype=int))
    
    minus
  else:
    if direction==1:  
        iss = range( i, len(phis) )
        i = -1
    else:
        iss = range( 0, i )
        iss = iss[::-1]
        i = 1
    isd = []
    isd.extend(ones(len(phis),dtype=int))
  jss = [j] * len(iss)
  kss = [k] * len(iss)  
  while found:
    #Get center of cell at toroidal surface i
    R = Grid.zones[0].R[i,j:j+2,k:k+2].mean()
    Z = Grid.zones[0].Z[i,j:j+2,k:k+2].mean()
    # Get mapping cell
    for j in xrange( Grid.zones[0].Z.shape[1]-1 ):
      for k in xrange( Grid.zones[0].Z.shape[2]-1 ):
        p = [ [Grid.zones[0].R[i,j  ,k  ], Grid.zones[0].Z[i,j  ,k  ]],
              [Grid.zones[0].R[i,j+1,k  ], Grid.zones[0].Z[i,j+1,k  ]],
              [Grid.zones[0].R[i,j+1,k+1], Grid.zones[0].Z[i,j+1,k+1]],
              [Grid.zones[0].R[i,j  ,k+1], Grid.zones[0].Z[i,j  ,k+1]],
              [Grid.zones[0].R[i,j  ,k  ], Grid.zones[0].Z[i,j  ,k  ]] ]
        poly = Path( p )
        found = poly.contains_point( (R, -Z) )
        if found:
          break
      if found:
        break
    if not found:
      print 'No mapping possible'
    else:
      phiOffset -= phiSegment
      phiss = append( phiss, phis + phiOffset )
      if i == -1:
        i = 0
        iss.extend( range( len(phis)-1, -1, -1) )
        isd.extend(-1*ones(len(phis),dtype=int))
        jss.extend( [j] * len(phis) )
        kss.extend( [k] * len(phis) )
        Rs = Grid.zones[0].R[::-1,j:j+2,k:k+2].mean((1,2))
        Rss = append( Rss, Rs[:-1] + 0.5 * diff( Rs ) )
        if Plates.ids[0][i,j,k]:
          print 'Downstream'
          iss = iss[:-1]
          jss = jss[:-1]
          kss = kss[:-1]
          phiss = phiss[:-1]
          Rss = Rss[:-1]
          break
      else:
        print counter
        if counter == 1000:
            break
            print 'running on and on.'
        counter=counter+1
        i = -1
        iss.extend( range( 0, len(phis) ) )
        isd.extend(ones(len(phis),dtype=int))
        jss.extend( [j] * len(phis) )
        kss.extend( [k] * len(phis) )
        Rs = Grid.zones[0].R[:,j:j+2,k:k+2].mean((1,2))
        Rss = append( Rss, Rs[:-1] + 0.5 * diff( Rs ) )
  return iss, jss, kss, phiss, Rss, isd[:-1]

def Cartesian2Grid( xs, ys, zs, PhiSimulated = 36.0, Periodicity = 5 ):
  """Returns the corresponding grid coordinates. (list of three ndarrays)

  This function transfers Cartesian coordinates into cylindrical 
  coordinates and applies the mirror symmetry at the grid boundary if
  necessary.

  Keyword arguments:
  xs -- x-coordinates. (ndarray or float)
  ys -- y-coordinates. (ndarray or float)
  zs -- z-coordinates. (ndarray or float)
  PhiSimulated -- Simulated angle. (float, default 36.0)
  Periodicity -- Periodicity of device. (int)

  """
  pToro = [ sqrt( xs**2 + ys**2), arctan2( ys, xs ), zs ]

  pToro[1] = mod( pToro[1], radians( 360.0 / Periodicity ) )

  for k in xrange( len( pToro[1] ) ):
    if radians( PhiSimulated ) < pToro[1][k]:
      # Assuming mirror symmetry at grid boundary
      pToro[1][k] = radians( 360.0 / Periodicity ) - pToro[1][k]
      pToro[2][k] *= -1.0

  pToro[1] = pToro[1] / pi * 180.0
  
  return pToro

def getCellIndex( points, Grid ):
  """Returns cell indices corresponding to the given points. (three lists of int)

  This function finds and returns all cell indices of the given points 
  specified in cylinder coordinates. 
  The indices are returned in the order toroidal, poloidal, radial.
  "NaN" is return, if no cell could be found.

  Keyword arguments:
  points -- Points given in R-, phi- and Z-coordinates ((list) of tuple(s) of 
            three float)
  Grid -- Grid to use. (Grid object)

  TODO:
  Needs speed-up by e.g. using VTK routines for finding cell index.

  ATTENTION:
  Only valid for first zone at the moment.
  To be generalised.

  """
  if isinstance( points, tuple ):
    points = [ points ]

  iTor = []
  iPol = []
  iRad = []
  phiOld = float( 'NaN' )
  for point in points:
    R, phi, Z = point

    if phi != phiOld:
      # Find toroidal angles of interest
      i = where( ( Grid.zones[0].phi[:-1] <= phi ) & \
                 ( phi < Grid.zones[0].phi[1:] ) )[0][0]
      Rinterp =  Grid.zones[0].R[i,:,:] + \
                 ( Grid.zones[0].R[i+1,:,:] - Grid.zones[0].R[i,:,:] ) * \
                 ( phi - Grid.zones[0].phi[i] ) / ( Grid.zones[0].phi[i+1] - Grid.zones[0].phi[i] )
      Zinterp =  Grid.zones[0].Z[i,:,:] + \
                 ( Grid.zones[0].Z[i+1,:,:] - Grid.zones[0].Z[i,:,:] ) * \
                 ( phi - Grid.zones[0].phi[i] ) / ( Grid.zones[0].phi[i+1] - Grid.zones[0].phi[i] )
      phiOld = phi

    # Get poloidal cross-section and find cell in cross-section
    for j in xrange( Grid.zones[0].Z.shape[1]-1 ):
      for k in xrange( Grid.zones[0].Z.shape[2]-1 ):
        p = [ [Rinterp[j  ,k  ], Zinterp[j  ,k  ]],
              [Rinterp[j+1,k  ], Zinterp[j+1,k  ]],
              [Rinterp[j+1,k+1], Zinterp[j+1,k+1]],
              [Rinterp[j  ,k+1], Zinterp[j  ,k+1]],
              [Rinterp[j  ,k  ], Zinterp[j  ,k  ]] ]
        poly = Path( p )
        found = poly.contains_point( (R, Z) )
        if found:
          iTor.append( i )
          iPol.append( j )
          iRad.append( k )
          break
      if found:
        break
    else:
      iTor.append( float('NaN') )
      iPol.append( float('NaN') )
      iRad.append( float('NaN') )
 
  return iTor, iPol, iRad

def getAreaQuad( R, Z ):
  """ Calculates the area of quads. (ndarray)

  This method calculates the area of quads via the cross-product
  by summation of the area of the two parallelograms.
  
  Keyword arguments:
  R -- Radial coordinates of the quad ((Nx)2x2-ndarray)
  Z -- Horizontal coordinates of the quad ((Nx)2x2-ndarray)

  """
  if R.ndim < 3:
    R = R[None,:,:]
  if Z.ndim < 3:
    Z = Z[None,:,:]
    
  area = 0.5 * abs( ( R[:,1,1] - R[:,1,0] ) * ( Z[:,1,0] - Z[:,0,0] )
                   -( R[:,1,0] - R[:,0,0] ) * ( Z[:,1,1] - Z[:,1,0] ) ) \
        +0.5 * abs( ( R[:,0,1] - R[:,0,0] ) * ( Z[:,1,1] - Z[:,0,1] ) 
                   -( R[:,1,1] - R[:,0,1] ) * ( Z[:,0,1] - Z[:,0,0] ) )
  return area

def getAreaPolygon( R, Z ):
  """ Calculates the area of a polygon. (float)

  This method calculates the area of a polygon using Gauss's area 
  formula.

  Keyword arguments:
  R -- Radial coordinates of the polygon (ndarray)
  Z -- Horizontal coordinates of the polygon (ndarray)

  """
  return 0.5 * abs( dot( R, roll( Z, 1 ) ) - dot( Z, roll( R, 1 ) ) )

def getLCFS( Grid, Limiter ):
  """Calculates and returns the effective radius of the last closed 
  flux surface (LCFS). (float)

  This function calculates and returns the effective radius of the LCFS 
  in a limiter configuration.
  It can concider a PlateCells or Installation object as a representative
  of the limiter.

  Keyword arguments:
  Grid -- Grid. (Grid object)
  Limiter -- Limiter defined via the plate cells or an installation. 
             (PlateCells or Installation object)

  ATTENTION:
  Only valid for first zone at the moment.
  To be generalised.

  """
  zone = Grid.zones[0]

  if isinstance( Limiter, PlateCells ):
    indicesLimiter = where( Limiter.ids[0] )[2]
  elif isinstance( Limiter, Installation ):
    indicesLimiter = []
    for R, phi, Z in zip( Limiter.R, Limiter.phi, Limiter.Z ):                         
      indicesLimiter.extend( getCellIndex( zip( R, [phi]*len(R), Z ), Grid )[2] )
    indicesLimiter = array( indicesLimiter )[~isnan( indicesLimiter )]
    indicesLimiter = array( indicesLimiter, int )
  else:
    raise TypeError( "Limiter type %s unknown." % type( Limiter ) )

  return zone.Reff[ indicesLimiter ].min()

