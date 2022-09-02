import numpy as np # NumPy
import scipy as sp # SciPy
import time
from scipy.optimize import curve_fit
import os
import math
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import EMC3
import json
from matplotlib.colors import LinearSegmentedColormap
import scipy.io
from itertools import product
import matplotlib.lines as lines


plt.ioff()





class limiter_hf_rad(object):
    """ Returns the radial profile onto the limiter starting from the coordinate givin by the start index

  Keyword arguments:
  timeframe -- time frame of the exp. data
  shift -- shift of the profiles below the original profile
  
  Keyword arguments for return:
  hf_line -- indicees of poloidal hf line onto the limiter 
    """  
  
    def __init__(self,timeframe,spot):
        #plt.ion()
        g=EMC3.Grid() 
        
        pathH = '.'
        Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
        PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )

        z=g.zones[0]
        rad=z.R
        tor=z.phi
        pol=z.Z
        
        tp=EMC3.TargetProfiles()
        tp0=tp[1]
        
        mat = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/marcin/160308022_3Dframes_emic.mat')
        r_line_exp=np.sqrt(mat['X']**2+mat['Y']**2)*100 

        r_line_idx=np.where(np.logical_and( mat['Z'][0]*100>=spot-1 , mat['Z'][0]*100 <= spot+1 ))
        r_line_coords=r_line_exp[0][r_line_idx]
        
        
        nearest=[]
        nearest_idx=[]
        for l in range( 0,len(r_line_coords)  ):
            nearest.append(find_nearest(tp0.R[:,0],r_line_coords[l])[0])
            nearest_idx.append(find_nearest(tp0.R[:,0],r_line_coords[l])[1])   

        hf_sim=tp0.profiles['q'][1][nearest_idx]
        hf_exp=mat['heatflux'][timeframe]
        hf_exp=hf_exp[r_line_idx]

        hf_sim=hf_sim
        hf_exp=hf_exp*1e2
        
        hf_sim[hf_sim==0] = np.nan
        hf_exp[hf_exp==0] = np.nan
             
        f, ax = plt.subplots()
        plt.plot(r_line_coords,hf_sim,'ro')
        plt.plot(r_line_coords,hf_exp,'bo')
        ax.set_ylabel( r'$q[W/cm^2]$')
        ax.set_xlabel( r'$R\mathrm{[cm]}$')
        ax.grid() 
        ax.set_ylim(0,250)
        #ax.set_ylim(0,200)
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()
        
        
        

        self.r_line_coords=r_line_coords
        self.nearest=nearest
        self.nearest_idx=nearest_idx
        self.hf_sim=hf_sim
        self.hf_exp=hf_exp
        self.r_line_exp=r_line_exp





class fluxtube_area(object):
    """ Returns the crosssection of a fluxtube for each toroidal 

  Keyword arguments:
  i -- Start index toroidal. (int)
  j -- Start index poloidal. (int)
  k -- Start index radial. (int)tive charge state as toroidal plasma property
  
  Keyword arguments for return:
  area -- indicees of the fluxtube  
    """  
  
    def __init__(self,i,j,k):
        #plt.ion()
        g=EMC3.Grid() 
        
        pathH = '.'
        Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
        PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )

        z=g.zones[0]
        rad=z.R
        tor=z.phi
        pol=z.Z

        fluxtube_tl=EMC3.Analysis.getFluxTubeIndices(i, j, k, g, PlatesH)
        
        meanR=[]
        for i in range(1,len(fluxtube_tl[4])):
            meanR.append((fluxtube_tl[4][i-1]+fluxtube_tl[4][i])/2)
        Lc=np.array(meanR)*0.5*1e-2*np.pi/180    
        Lc=np.sum(Lc)

        dummy1=rad[:-1,:-1,:-1]
        dummy2=tor[:-1]
        dummy3=pol[:-1,:-1,:-1]
        
        R_flux_tl=dummy1[fluxtube_tl[0],fluxtube_tl[1],fluxtube_tl[2]]
        Phi_flux_tl=dummy2[fluxtube_tl[0]]
        Z_flux_tl=dummy3[fluxtube_tl[0],fluxtube_tl[1],fluxtube_tl[2]]
        
        R_flux_tr=dummy1[fluxtube_tl[0],fluxtube_tl[1],(np.array(fluxtube_tl[2])-1)]
        Phi_flux_tr=dummy2[fluxtube_tl[0]]
        Z_flux_tr=dummy3[fluxtube_tl[0],fluxtube_tl[1],np.array(fluxtube_tl[2])-1]
        
        R_flux_bl=dummy1[fluxtube_tl[0],np.array(fluxtube_tl[1])+1,fluxtube_tl[2]]
        Phi_flux_bl=dummy2[fluxtube_tl[0]]
        Z_flux_bl=dummy3[fluxtube_tl[0],np.array(fluxtube_tl[1])+1,fluxtube_tl[2]]
        
        R_flux_br=dummy1[fluxtube_tl[0],np.array(fluxtube_tl[1])+1,np.array(fluxtube_tl[2])-1]
        Phi_flux_br=dummy2[fluxtube_tl[0]]
        Z_flux_br=dummy3[fluxtube_tl[0],np.array(fluxtube_tl[1])+1,np.array(fluxtube_tl[2])-1]
        
        
        dummyR=(2,2)
        dummyR=np.zeros(dummyR)
        dummyZ=(2,2)
        dummyZ=np.zeros(dummyZ)
        #self.R_flux_tl=R_flux_tl
        #self.R_flux_tr=R_flux_tr
        #self.R_flux_bl=R_flux_bl
        #self.R_flux_br=R_flux_br
        #self.Z_flux_tl=Z_flux_tl
        #self.Z_flux_tr=Z_flux_tr
        #self.Z_flux_bl=Z_flux_bl
        #self.Z_flux_br=Z_flux_br
        dummy=[]
        for i in range(0,len(R_flux_tl)):    
            dummyR[0,0]=R_flux_tl[i]
            dummyR[0,1]=R_flux_tr[i]
            dummyR[1,0]=R_flux_bl[i]
            dummyR[1,1]=R_flux_br[i]
            dummyZ[0,0]=Z_flux_tl[i]
            dummyZ[0,1]=Z_flux_tr[i]
            dummyZ[1,0]=Z_flux_bl[i]
            dummyZ[1,1]=Z_flux_br[i]
            dummy.append(EMC3.Analysis.getAreaQuad(dummyR,dummyZ))
        
        
        f, ax = plt.subplots()
        plt.plot(np.linspace(0,Lc,len(dummy)),dummy)
        ax.set_ylabel( r'$A_\mathrm{fluxtube}\mathrm{[cm^2]}$')
        ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
        ax.grid() 
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
#        ax.set_xlim(0,120)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()
        
        f, ax = plt.subplots()
        for i in range(340,380):    
            plt.plot([R_flux_tl[i], R_flux_tr[i]], [Z_flux_tl[i], Z_flux_tr[i]],'r')
            ax.text(R_flux_tl[i]+2, Z_flux_tl[i]+2, str(i), fontsize=2)
            plt.plot([R_flux_tr[i], R_flux_br[i]], [Z_flux_tr[i], Z_flux_br[i]],'b')
            plt.plot([R_flux_br[i], R_flux_bl[i]], [Z_flux_br[i], Z_flux_bl[i]],'k')
            plt.plot([R_flux_bl[i], R_flux_tl[i]], [Z_flux_bl[i], Z_flux_tl[i]],'g')


        ax.set_ylabel( r'$R_\mathrm{corners}[\mathrm{cm}]$')
        ax.set_xlabel( r'$Z[\mathrm{cm}]$')
        ax.grid() 
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()
    

        f, ax = plt.subplots()
        i=358
        plt.plot([R_flux_tl[i], R_flux_tr[i]], [Z_flux_tl[i], Z_flux_tr[i]],'r')
        ax.text(R_flux_tl[i]+2, Z_flux_tl[i]+2, str(i), fontsize=2)
        plt.plot([R_flux_tr[i], R_flux_br[i]], [Z_flux_tr[i], Z_flux_br[i]],'b')
        plt.plot([R_flux_br[i], R_flux_bl[i]], [Z_flux_br[i], Z_flux_bl[i]],'k')
        plt.plot([R_flux_bl[i], R_flux_tl[i]], [Z_flux_bl[i], Z_flux_tl[i]],'g')


        ax.set_ylabel( r'$R_\mathrm{corners}[\mathrm{cm}]$')
        ax.set_xlabel( r'$Z[\mathrm{cm}]$')
        ax.grid() 
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()


        f, ax = plt.subplots()
        i=359
        plt.plot([R_flux_tl[i], R_flux_tr[i]], [Z_flux_tl[i], Z_flux_tr[i]],'r')
        ax.text(R_flux_tl[i]+2, Z_flux_tl[i]+2, str(i), fontsize=2)
        plt.plot([R_flux_tr[i], R_flux_br[i]], [Z_flux_tr[i], Z_flux_br[i]],'b')
        plt.plot([R_flux_br[i], R_flux_bl[i]], [Z_flux_br[i], Z_flux_bl[i]],'k')
        plt.plot([R_flux_bl[i], R_flux_tl[i]], [Z_flux_bl[i], Z_flux_tl[i]],'g')


        ax.set_ylabel( r'$R_\mathrm{corners}[\mathrm{cm}]$')
        ax.set_xlabel( r'$Z[\mathrm{cm}]$')
        ax.grid() 
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()



        self.dummy1=dummy1
        self.dummy2=dummy2
        self.dummy3=dummy3
        self.R_flux_tl=R_flux_tl
        self.R_flux_tr=R_flux_tr
        self.R_flux_bl=R_flux_bl
        self.R_flux_br=R_flux_br
        self.Z_flux_tl=Z_flux_tl
        self.Z_flux_tr=Z_flux_tr
        self.Z_flux_bl=Z_flux_bl
        self.Z_flux_br=Z_flux_br
        self.dummy=dummy



class toroidal_profile(object):
    """ Returns the toroidal Profile as plot and as np.arrays

  Keyword arguments:
  i -- Start index toroidal. (int)
  j -- Start index poloidal. (int)
  k -- Start index radial. (int)
  measurand -- which measurand should be plotted (string)
  
  Keyword arguments for return:
  fluxtube -- indicees of the fluxtube [iss,jss,kss,phiss,Rss]
  dens -- density profile along the fluxtube
  temp -- temperature profile along the fluxtube
  
  TODO:
  
  the Zeff profile along the fluxtube must be implemented
  
  ATTENTION:
  function d_sum is not usable since the function only sums up the density over the radial axis
  
    """  
  
    def __init__(self,i,j,k,measurand,direction=1):
        #plt.ion()
        g=EMC3.Grid() 

        pathH = '.'
        Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
        PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )

        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')
        Mach=EMC3.PlasmaField('MACH_NUMBER')

        fluxtube=EMC3.Analysis.getFluxTubeIndices(i, j, k, g, PlatesH, direction)

        dummy1=d('n1')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy1_2=d('n2')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy1_3=d('n3')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy1_4=d('n4')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy1_5=d('n5')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy1_6=d('n6')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy1_7=d('n7')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy2=T('Te')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy3=T('Ti')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        dummy4=Mach('M')[0][fluxtube[0],fluxtube[1],fluxtube[2]]
        
        mi=1.66*1e-27*1.00794
        cs=np.sqrt((dummy2*1.602*1e-19+dummy3*1.602*1e-19)/mi)

        counter=0
        meanR=[]
        for i in range(1,len(fluxtube[4])):
            meanR.append((fluxtube[4][i-1]+fluxtube[4][i])/2)
        Lc=np.array(meanR)*0.5*1e-2*np.pi/180    
        Lc=np.sum(Lc)

        x_plot=np.linspace(0,Lc,len(dummy1))
        
        if measurand=='ne':
            f, ax = plt.subplots()
            n_e=dummy1+1*dummy1_2+2*dummy1_3+3*dummy1_4+4*dummy1_5+5*dummy1_6+6*dummy1_7
            plt.plot(x_plot,n_e/1e12,color='#0070bd')
            ax.set_ylabel( r'$n_\mathrm{e}[10^{12}\mathrm{m^3}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.grid() 
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
            
        if measurand=='ni':
            f, ax = plt.subplots()
            n_e=dummy1
            plt.plot(x_plot,n_e/1e12,color='#0070bd')
            ax.set_ylabel( r'$n_\mathrm{i}[10^{18}\mathrm{m^{-3}}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.grid() 
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
            
        if measurand=='Te':
            f, ax = plt.subplots()
            plt.plot(x_plot,dummy2,color='#0070bd')
            ax.set_ylabel( r'$T_\mathrm{e}[\mathrm{eV}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.grid()
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)

        if measurand=='Ti':
            f, ax = plt.subplots()
            plt.plot(x_plot,dummy3,color='#0070bd')
            ax.set_ylabel( r'$T_\mathrm{i}[\mathrm{eV}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.grid()
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)

        if measurand=='pe':
            f, ax = plt.subplots()
            n_e=dummy1+1*dummy1_2+2*dummy1_3+3*dummy1_4+4*dummy1_5+5*dummy1_6+6*dummy1_7
            plt.plot(x_plot,n_e*dummy2*1.602*1e-19*1e6,color='#0070bd')
            ax.set_ylabel( r'$p_\mathrm{e}[\mathrm{10^{-6}Pa}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.set_ylim(0.2,1.8)
            ax.grid()
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
            
        if measurand=='pi':
            f, ax = plt.subplots()
            n_i=dummy1+dummy1_2+dummy1_3+dummy1_4+dummy1_5+dummy1_6+dummy1_7
            plt.plot(x_plot,n_i*dummy3*1.602*1e-19*1e6,color='#0070bd')
            ax.set_ylabel( r'$p_\mathrm{i}[\mathrm{10^{-6}Pa}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.set_ylim(0.2,1.8)
            ax.grid()
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)                    
            
        if measurand=='ptot':
            f, ax = plt.subplots()
            n_i=dummy1+dummy1_2+dummy1_3+dummy1_4+dummy1_5+dummy1_6+dummy1_7
            n_e=dummy1+1*dummy1_2+2*dummy1_3+3*dummy1_4+4*dummy1_5+5*dummy1_6+6*dummy1_7
            pi=n_i*dummy3*1.602*1e-19
            pe=n_e*dummy2*1.602*1e-19
            pkin=n_i*dummy4*dummy4*mi*cs*cs
            plt.plot(x_plot, (pe+pi+pkin)*1e6 ,color='#0070bd')
            ax.set_ylabel( r'$p_\mathrm{total}[\mathrm{10^{-6}Pa}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.grid()
            ax.set_ylim(0.2,1.8)
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)        

        if measurand=='pkin':
            f, ax = plt.subplots()
            n_i=dummy1+dummy1_2+dummy1_3+dummy1_4+dummy1_5+dummy1_6+dummy1_7
            n_e=dummy1+1*dummy1_2+2*dummy1_3+3*dummy1_4+4*dummy1_5+5*dummy1_6+6*dummy1_7
            pi=n_i*dummy3*1.602*1e-19
            pe=n_e*dummy2*1.602*1e-19
            pkin=n_i*dummy4*dummy4*mi*cs*cs
            plt.plot(x_plot, (pkin)*1e6 ,color='#0070bd')
            ax.set_ylabel( r'$p_\mathrm{kin}[\mathrm{10^{-6}Pa}]$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.grid()
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)        
        
        if measurand=='M':
            f, ax = plt.subplots()
            plt.plot(x_plot,dummy4*fluxtube[-1][1:],color='#ffa500')
            ax.set_ylabel( r'$M$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.set_xlim(0,120)
            ax.grid()
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            #plt.setp(ax.get_yticklabels()[::2], visible=False)        


        if measurand=='Msquared':
            f, ax = plt.subplots()
            plt.plot(x_plot,dummy4*fluxtube[-1][1:]*dummy4*fluxtube[-1][1:],color='#0070bd')
            ax.set_ylabel( r'$M^2$')
            ax.set_xlabel( r'$L_\mathrm{c}\mathrm{[m]}$')
            ax.grid()
            plt.setp(ax.get_xticklabels()[::2], visible=False)

        print 'end point R'
        print fluxtube[4][-1]
        print 'end point Phi'
        print fluxtube[3][-1]

        plt.rcParams.update({'font.size': 18})
        f.tight_layout()

        self.fluxtube=fluxtube
        self.d=dummy1
        self.Te=dummy2
        self.Ti=dummy3
        self.M=dummy4
        self.Lc=x_plot
        self.cs=cs
        self.mi=mi



class limdata_reff(object):
    """ Returns the radial Profile as plot and as np.arrays

  Keyword arguments:
  measurand -- which measurand should be plotted (string)
  
  Stings are:
  'ne' -- input string to get density of the electrons plotted
  'Te' -- input string to get temperature of the electrons plotted
  'ne_up_down' -- input string to get density of the electrons plotted with seperation of up and down probeheads 
  'Te_up_down' -- input string to get temperature of the electrons plotted with seperation of up and down probeheads 
  
  TODO:
  Legend for the up and down probehead seperation must be implemented 

  ATTENTION:
  Red exp_results are corresponding to the probehead below z==0.21[m]
    """


    def __init__(self,measurerand):
        #plt.ion()
        g=EMC3.Grid() 


        a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

#Boyd Data Load
        data_json=json.load(open('/home/j.cosfeld/runs/W7X_limiter/exp_results/boyd/LP20160308_23_L53_2k2.json'))
        if measurerand=='ne':
            data=np.array(data_json['ne18'])
        if measurerand=='Te':
            data=np.array(data_json['Te'])
        if measurerand=='ne_sep':
            data=np.array(data_json['ne18'])
        if measurerand=='Te_sep':
            data=np.array(data_json['Te'])


            
        max_data_boyd=np.zeros(data_json.shape[1])        
        errorbar_boyd=np.zeros(data_json.shape[1])   
        
        
        for i in range (0,data.shape[1]):
            max_data_boyd[i]=np.nanmax(data[250,i])
            errorbar_boyd[i]=np.average(data[240:275,7]) 

        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')

        ne_sum=d_sum(g,d,0,0)

               
        f, ax = plt.subplots()
        if measurerand=='Te':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(T('Te')[0][1,0,:]),'x',linewidth=2,color="blue")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd,yerr=errorbar_boyd-max_data_boyd, color='orange',fmt='o')
            ax.set_ylabel( r'$T_e[eV]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,180)
            ax.grid()
        if measurerand=='ne':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6),'x',linewidth=2,color="r")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*1e18,yerr=errorbar_boyd*1e18-max_data_boyd*1e18, color='darkviolet',fmt='o')
            ax.set_ylabel( r'$n_e[1/m^3]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,2*1e19)
            ax.grid()
        if measurerand=='Te_sep':
            coords=np.array(data_json['info']['coords'])
            coordsR=np.sqrt(coords[:,0]**2+coords[:,1]**2)
            coordsZup=np.where(coordsR*100-coordsR[10]*100>0)[0]
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(T('Te')[0][1,0,:]),'x',linewidth=2,color="blue")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd,yerr=errorbar_boyd-max_data_boyd, color='orange',fmt='o')
            for i in range (0,coordsZup.shape[0]):
                ax.errorbar(r_eff_boyd[i]/a_eff_exp,max_data_boyd[i],yerr=errorbar_boyd[i]-max_data_boyd[i], color='red',fmt='o')
            ax.set_ylabel( r'$T_e[eV]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,180)
            ax.grid()
        if measurerand=='ne_sep':
            coords=np.array(data_json['info']['coords'])
            coordsR=np.sqrt(coords[:,0]**2+coords[:,1]**2)
            coordsZup=np.where(coordsR*100-coordsR[10]*100>0)[0]
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6),'x',linewidth=2,color="r")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*1e18,yerr=errorbar_boyd*1e18-max_data_boyd*1e18, color='darkviolet',fmt='o')
            for i in range (0,coordsZup.shape[0]):
                ax.errorbar(r_eff_boyd[i]/a_eff_exp,max_data_boyd[i]*1e18,yerr=errorbar_boyd[i]*1e18-max_data_boyd[i]*1e18, color='blue',fmt='o')
            ax.set_ylabel( r'$n_e[1/m^3]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,2*1e19)
            ax.grid()

        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()

class divertor_hf(object):
    """ Returns the radial Profile as plot and as np.arrays

  Keyword arguments:
  measurand -- which measurand should be plotted (string)
  
  Stings are:
  'ne' -- input string to get density of the electrons plotted
  'Te' -- input string to get temperature of the electrons plotted
  'ne_up_down' -- input string to get density of the electrons plotted with seperation of up and down probeheads 
  'Te_up_down' -- input string to get temperature of the electrons plotted with seperation of up and down probeheads 
  
  TODO:
  Legend for the up and down probehead seperation must be implemented 

  ATTENTION:
  Red exp_results are corresponding to the probehead below z==0.21[m]
    """

    def __init__(self):
        #plt.ion()
        g=EMC3.Grid() 

        radius,r_eff_ken=get_a_r_eff_div(g)
        a_eff=radius[32]
        
        hf = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_divertor/exp_results/heatflux/20171026.038_AEF20_heatflux_V4/module20/20171026.038_AEF20_heatflux_V4.mat')
        
        hf00=hf['heat_20_0']
        X00=hf['plot_X_20_0']
        Y00=hf['plot_Y_20_0']
        Z00=hf['plot_Z_20_0']
        
        pCart2=[X00[0,:]*100,Y00[0,:]*100,Z00[0,:]*100]
        Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart2 )
        iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), g )
        
        self.iTor=iTor
        self.iPol=iPol
        self.iRad=iRad
        
        tp=EMC3.TargetProfiles()
        
        #finding alg. for R
        dummyR=np.zeros(shape=(4,2))
        dummyPhi=np.zeros(shape=(4,2))
        dummyZ=np.zeros(shape=(4,2))

        m=7
        idxR=np.zeros(len(Rs))
        idxP=np.zeros(len(Rs))
        idxZ=np.zeros(len(Rs))
        for i in range(0,len(Rs)):
            idxR[i]=int(find_nearest(tp[m].R[:,0],Rs[i])[1])
            idxP[i]=int(find_nearest(tp[m].phi[:,0],phis[i])[1])
            idxZ[i]=int(find_nearest(tp[m].Z[:,0],Zs[i])[1])
            
        self.idxR=idxR
        self.idxP=idxP
        self.idxZ=idxZ
        
        hfR=tp[7].profiles['q'][1][idxR.astype(int)]
        self.hfR=hfR
        
        f, ax = plt.subplots()
        plt.plot(Rs,hfR*1e4/1e6)
        plt.plot(Rs,hf00[50,:]/1e6)
        ax.set_ylabel( r'$q[MW/m^2]$')
        ax.set_xlabel( r'$R\mathrm{[cm]}$')
        ax.grid() 
        ax.set_xlim(520,540)
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()
        
        
        
        #for i in range(0,4):
            #arr1=tp[m].R[:,i]
            #arr2=Rs
            #dummyR[i,:]=sorted(product(arr1, arr2), key=lambda t: abs(t[0]-t[1]))[0]
        #for i in range(0,4):
            #arr1=tp[m].phi[:,i]
            #arr2=phis
            #dummyPhi[i,:]=sorted(product(arr1, arr2), key=lambda t: abs(t[0]-t[1]))[0]
        #for i in range(0,4):
            #arr1=tp[m].Z[:,i]
            #arr2=Zs
            #dummyZ[i,:]=sorted(product(arr1, arr2), key=lambda t: abs(t[0]-t[1]))[0]            

        #self.dummyR=dummyR
        #self.dummyPhi=dummyPhi
        #self.dummyZ=dummyZ
        
        #iTorEMC3, iPolEMC3, iRadEMC3 = EMC3.Analysis.getCellIndex( zip(dummyR[:,0], dummyPhi[:,0], dummyZ[:,0]), g )
        #iTorEXP, iPolEXP, iRadEXP    = EMC3.Analysis.getCellIndex( zip(dummyR[:,1], dummyPhi[:,1], dummyZ[:,1]), g )        
        
        #self.iTorEMC3=iTorEMC3
        #self.iTorEXP=iTorEXP
        #self.iPolEMC3=iPolEMC3      
        #self.iPolEXP=iPolEXP
        #self.iRadEMC3=iRadEMC3      
        #self.iRadEXP=iRadEXP


class divertor_lp(object):
    """ Returns the radial Profile as plot and as np.arrays

  Keyword arguments:
  measurand -- which measurand should be plotted (string)
  
  Stings are:
  'ne' -- input string to get density of the electrons plotted
  'Te' -- input string to get temperature of the electrons plotted
  'ne_up_down' -- input string to get density of the electrons plotted with seperation of up and down probeheads 
  'Te_up_down' -- input string to get temperature of the electrons plotted with seperation of up and down probeheads 
  
  TODO:
  Legend for the up and down probehead seperation must be implemented 

  ATTENTION:
  Red exp_results are corresponding to the probehead below z==0.21[m]
    """


    def __init__(self,measurerand,z_eff_downstream = 'none', z_eff_upstream = 'none'):
        #plt.ion()
        g=EMC3.Grid() 


        radius,r_eff_ken=get_a_r_eff_div(g)
        a_eff=radius[32]

#Boyd Data Load
        ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/coordinatesLD.txt')
        ken=ken/10.0
        if measurerand=='ne':
            ne_ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_ne038.txt')
            ne_ken_err=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_err_ne038.txt')
            ne_ken=ne_ken[0:10,422]
            ne_ken_err=ne_ken_err[0:10,422]
        if measurerand=='ne2':
            ne_ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_ne038.txt')
            ne_ken_err=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_err_ne038.txt')
            ne_ken=ne_ken[0:10,422]
            ne_ken_err=ne_ken_err[0:10,422]
        if measurerand=='ne_sorted':
            ne_ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_ne038.txt')
            ne_ken_err=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_err_ne038.txt')
            ne_ken=ne_ken[0:20,422]
            ne_ken_err=ne_ken_err[0:20,422]            
        if measurerand=='Te':
            Te_ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_Te038.txt')
            Te_ken_err=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_err_Te038.txt')
            Te_ken=Te_ken[0:10,422]
            Te_ken_err=Te_ken_err[0:10,422]
        if measurerand=='Te2':
            Te_ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_Te038.txt')
            Te_ken_err=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_err_Te038.txt')
            Te_ken=Te_ken[0:10,422]
            Te_ken_err=Te_ken_err[0:10,422]            
        if measurerand=='Te_sorted':
            Te_ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_Te038.txt')
            Te_ken_err=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/LD_err_Te038.txt')
            Te_ken=Te_ken[0:20,422]
            Te_ken_err=Te_ken_err[0:20,422]

            
        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')

        ne_sum=d_sum(g,d,0,0)

        x_emc3=radius
        x_ken=r_eff_ken
        
        self.x_emc3=x_emc3
        self.x_ken=x_ken
        
        overlap=[]
        for i in range(0,x_ken.shape[0]):
            overlap.append(find_nearest(x_emc3,x_ken[i])[1])
        
        self.overlap=overlap
        
        f, ax = plt.subplots()         

        if measurerand=='ne':
            ax.plot((radius)[0:-1]/a_eff
            ,zero_to_nan(np.mean(d('n1')[0][17,354:393,:]*1e6,axis=0))
            ,'--',linewidth=4,color='#d91219')
            ax.errorbar(radius[48]/a_eff,ne_ken[0],yerr=ne_ken_err[0], color='darkviolet',fmt='o')
            ax.errorbar(radius[44]/a_eff,ne_ken[1],yerr=ne_ken_err[1], color='darkviolet',fmt='o')
            ax.errorbar(radius[39]/a_eff,ne_ken[2],yerr=ne_ken_err[2], color='darkviolet',fmt='o')
            ax.errorbar(radius[37]/a_eff,ne_ken[3],yerr=ne_ken_err[3], color='darkviolet',fmt='o')
            ax.errorbar(radius[36]/a_eff,ne_ken[4],yerr=ne_ken_err[4], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,ne_ken[5],yerr=ne_ken_err[5], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,ne_ken[6],yerr=ne_ken_err[6], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,ne_ken[7],yerr=ne_ken_err[7], color='darkviolet',fmt='o')
            ax.errorbar(radius[34]/a_eff,ne_ken[8],yerr=ne_ken_err[8], color='darkviolet',fmt='o')
            ax.errorbar(radius[35]/a_eff,ne_ken[9],yerr=ne_ken_err[9], color='darkviolet',fmt='o')
            ax.set_ylabel( r'$n_\mathrm{e}[\mathrm{1/m^3}]$')
            ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
            ax.set_ylim(0,1.5*1e18)
            ax.set_xlim(0.9,1.2)
            ax.grid()

        if measurerand=='ne2':
            pCart2=[ken[:,0],ken[:,1],ken[:,2]]
            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart2 )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), g )
            sorting=np.zeros([np.shape(iTor)[0]/2.0,4])
            sorting[:,0]=ne_ken
            sorting[:,1]=iTor[0:10]
            sorting[:,2]=iPol[0:10]
            sorting[:,3]=iRad[0:10]
            ne_sorted=sorting[np.lexsort((sorting[:, 0], ))]
            a=0

            densities=np.zeros(10)
            for i in range(0,10):
                print (d('n1')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6+a)
                densities[i]=d('n1')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6+\
                            d('n2')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6*1+\
                            d('n3')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6*2+\
                            d('n4')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6*3+\
                            d('n5')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6*4+\
                            d('n6')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6*5+\
                            d('n7')[0][int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3])]*1e6*6
                            
            print densities
            
            
            ax.plot(radius[48]/a_eff,densities[0],'o',color='#d91219')
            ax.plot(radius[44]/a_eff,densities[1],'o',color='#d91219')
            ax.plot(radius[39]/a_eff,densities[2],'o',color='#d91219')
            ax.plot(radius[37]/a_eff,densities[3],'o',color='#d91219')
            ax.plot(radius[36]/a_eff,densities[4],'o',color='#d91219')
            ax.plot(radius[33]/a_eff,densities[5],'o',color='#d91219')
            ax.plot(radius[33]/a_eff,densities[6],'o',color='#d91219')
            ax.plot(radius[33]/a_eff,densities[7],'o',color='#d91219')
            ax.plot(radius[34]/a_eff,densities[8],'o',color='#d91219')
            ax.plot(radius[35]/a_eff,densities[9],'o',color='#d91219')
            z_eff_downstream=0
            counter=0
            for i in range(0,10):
                ne_z_sum_downstream=d_z_sum_impu_point(g,d,int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3]),'C')
                ne_sum_downstream=d_sum_impu_point(g,d,int(ne_sorted[i,1]),int(ne_sorted[i,2]),int(ne_sorted[i,3]),'C')
                if ne_z_sum_downstream[0]==0:
                    counter+=1
                #print np.type(np.min(ne_z_sum_downstream/ne_sum_downstream))
                z_eff_downstream+=np.nan_to_num(np.min(ne_z_sum_downstream/ne_sum_downstream))            
            z_eff_downstream=z_eff_downstream/(10-counter)
            ax.errorbar(radius[48]/a_eff,ne_ken[0]*np.min(z_eff_downstream),yerr=ne_ken_err[0], color='darkviolet',fmt='o')         
            ax.errorbar(radius[44]/a_eff,ne_ken[1]*np.min(z_eff_downstream),yerr=ne_ken_err[1], color='darkviolet',fmt='o')
            ax.errorbar(radius[39]/a_eff,ne_ken[2]*np.min(z_eff_downstream),yerr=ne_ken_err[2], color='darkviolet',fmt='o')
            ax.errorbar(radius[37]/a_eff,ne_ken[3]*np.min(z_eff_downstream),yerr=ne_ken_err[3], color='darkviolet',fmt='o')
            ax.errorbar(radius[36]/a_eff,ne_ken[4]*np.min(z_eff_downstream),yerr=ne_ken_err[4], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,ne_ken[5]*np.min(z_eff_downstream),yerr=ne_ken_err[5], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,ne_ken[6]*np.min(z_eff_downstream),yerr=ne_ken_err[6], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,ne_ken[7]*np.min(z_eff_downstream),yerr=ne_ken_err[7], color='darkviolet',fmt='o')
            ax.errorbar(radius[34]/a_eff,ne_ken[8]*np.min(z_eff_downstream),yerr=ne_ken_err[8], color='darkviolet',fmt='o')
            ax.errorbar(radius[35]/a_eff,ne_ken[9]*np.min(z_eff_downstream),yerr=ne_ken_err[9], color='darkviolet',fmt='o')
            ax.set_ylabel( r'$n_\mathrm{e}[\mathrm{1/m^3}]$')
            ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
            #ax.set_ylim(0,1.5*1e18)
            #ax.set_xlim(0.9,1.2)
            ax.grid()

            
        if measurerand=='ne_sorted':
            a=0
            pCart2=[ken[:,0],ken[:,1],ken[:,2]]
            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart2 )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), g )
            sorting=np.zeros([np.shape(iTor)[0],4])
            sorting[:,0]=ne_ken
            sorting[:,1]=iTor[0:20]
            sorting[:,2]=iPol[0:20]
            sorting[:,3]=iRad[0:20]
            ne_sorted=sorting[np.lexsort((sorting[:, 0], ))]
            
            nes1=ne_sorted[:,1]
            nes2=ne_sorted[:,2]
            nes3=ne_sorted[:,3]            
            dummy_dens=d('n1')[0][zip(nes1),zip(nes2),zip(nes3)]*1e6
            wasdalos=zero_to_nan(dummy_dens)
            self.wasdalos=wasdalos
            self.ne_sorted=ne_sorted
            for i in range(0, 20):
                ax.plot(ne_sorted[i,0]/1e18,
                (wasdalos[i])/1e18,'o'
                ,color='darkviolet')
                
            self.nes1=ne_sorted[:,1]
            self.nes2=ne_sorted[:,2]
            self.nes3=ne_sorted[:,3]            
            self.d=d
                
            f2, ax2 = plt.subplots()    
            ax2.plot(range(1,21),ne_sorted[:,0]/1e18,'o',color='darkviolet')
            dummy_dens=d('n1')[0][zip(nes1),zip(nes2),zip(nes3)]*1e6
            a=zero_to_nan(dummy_dens)
            self.a=a
            ax2.plot(np.array(a)/1e18,'s',color='darkviolet')
            dots_x=np.array([0,2])
            dots_y=np.array([0,2])
            ax.plot(dots_x,dots_y, '-',color='darkviolet')
            ax.set_ylim(0,2)
            ax.set_xlim(0,2)
            ax.set_ylabel( r'$n_\mathrm{e,EMC3}[\mathrm{10^{18}/m^3}]$')
            ax.set_xlabel( r'$n_\mathrm{e,EXP}[\mathrm{10^{18}/m^3}]$')
            ax.grid()            
            ax2.set_ylabel( r'$n_\mathrm{e,EMC3}[1/m^3}]$')
            ax2.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
            ax2.grid()            
            
            
        if measurerand=='Te':
            ax.plot((radius)[0:-1]/a_eff
            ,zero_to_nan(np.mean(T('Te')[0][17,354:393,:],axis=0))
            ,'--',linewidth=4,color='#0070bd')
            ax.plot((radius)[0:-1]/a_eff
            ,zero_to_nan(np.mean(T('Ti')[0][17,354:393,:],axis=0))
            ,'--',linewidth=4,color='#d91219')
            ax.errorbar(radius[48]/a_eff,Te_ken[0],yerr=Te_ken_err[0], color='orange',fmt='o')
            ax.errorbar(radius[44]/a_eff,Te_ken[1],yerr=Te_ken_err[1], color='orange',fmt='o')
            ax.errorbar(radius[39]/a_eff,Te_ken[2],yerr=Te_ken_err[2], color='orange',fmt='o')
            ax.errorbar(radius[37]/a_eff,Te_ken[3],yerr=Te_ken_err[3], color='orange',fmt='o')
            ax.errorbar(radius[36]/a_eff,Te_ken[4],yerr=Te_ken_err[4], color='orange',fmt='o')
            ax.errorbar(radius[33]/a_eff,Te_ken[5],yerr=Te_ken_err[5], color='orange',fmt='o')
            ax.errorbar(radius[33]/a_eff,Te_ken[6],yerr=Te_ken_err[6], color='orange',fmt='o')
            ax.errorbar(radius[33]/a_eff,Te_ken[7],yerr=Te_ken_err[7], color='orange',fmt='o')
            ax.errorbar(radius[34]/a_eff,Te_ken[8],yerr=Te_ken_err[8], color='orange',fmt='o')
            ax.errorbar(radius[35]/a_eff,Te_ken[9],yerr=Te_ken_err[9], color='orange',fmt='o')
            ax.set_ylabel( r'$T_\mathrm{e}[\mathrm{eV}]$')
            ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
            ax.set_ylim(0,100)
            ax.set_xlim(0.9,1.2)
            ax.grid()
     
        if measurerand=='Te2':
            pCart2=[ken[:,0],ken[:,1],ken[:,2]]
            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart2 )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), g )
            sorting=np.zeros([np.shape(iTor)[0]/2.0,4])
            sorting[:,0]=Te_ken
            sorting[:,1]=iTor[0:10]
            sorting[:,2]=iPol[0:10]
            sorting[:,3]=iRad[0:10]
            ne_sorted=sorting[np.lexsort((sorting[:, 0], ))]
            a=0
            ax.plot(radius[48]/a_eff,T('Te')[0]   [int(ne_sorted[0,1]),int(ne_sorted[0,2]),int(ne_sorted[0,3])]+a,'o',color='#d91219')                        
            ax.plot(radius[44]/a_eff,T('Te')[0][int(ne_sorted[1,1]),int(ne_sorted[1,2]),int(ne_sorted[1,3])]+a,'o',color='#d91219')
            ax.plot(radius[39]/a_eff,T('Te')[0][int(ne_sorted[2,1]),int(ne_sorted[2,2]),int(ne_sorted[2,3])]+a,'o',color='#d91219')
            ax.plot(radius[37]/a_eff,T('Te')[0][int(ne_sorted[3,1]),int(ne_sorted[3,2]),int(ne_sorted[3,3])]+a,'o',color='#d91219')
            ax.plot(radius[36]/a_eff,T('Te')[0][int(ne_sorted[4,1]),int(ne_sorted[4,2]),int(ne_sorted[4,3])]+a,'o',color='#d91219')
            ax.plot(radius[33]/a_eff,T('Te')[0][int(ne_sorted[5,1]),int(ne_sorted[5,2]),int(ne_sorted[5,3])]+a,'o',color='#d91219')
            ax.plot(radius[33]/a_eff,T('Te')[0][int(ne_sorted[6,1]),int(ne_sorted[6,2]),int(ne_sorted[6,3])]+a,'o',color='#d91219')
            ax.plot(radius[33]/a_eff,T('Te')[0][int(ne_sorted[6,1]),int(ne_sorted[6,2]),int(ne_sorted[6,3])]+a,'o',color='#d91219')
            ax.plot(radius[33]/a_eff,T('Te')[0][int(ne_sorted[7,1]),int(ne_sorted[7,2]),int(ne_sorted[7,3])]+a,'o',color='#d91219')
            ax.plot(radius[34]/a_eff,T('Te')[0][int(ne_sorted[8,1]),int(ne_sorted[8,2]),int(ne_sorted[8,3])]+a,'o',color='#d91219')
            ax.plot(radius[35]/a_eff,T('Te')[0][int(ne_sorted[9,1]),int(ne_sorted[9,2]),int(ne_sorted[9,3])]+a,'o',color='#d91219')
            ax.errorbar(radius[48]/a_eff,Te_ken[0],yerr=Te_ken_err[0], color='darkviolet',fmt='o')
            ax.errorbar(radius[44]/a_eff,Te_ken[1],yerr=Te_ken_err[1], color='darkviolet',fmt='o')
            ax.errorbar(radius[39]/a_eff,Te_ken[2],yerr=Te_ken_err[2], color='darkviolet',fmt='o')
            ax.errorbar(radius[37]/a_eff,Te_ken[3],yerr=Te_ken_err[3], color='darkviolet',fmt='o')
            ax.errorbar(radius[36]/a_eff,Te_ken[4],yerr=Te_ken_err[4], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,Te_ken[5],yerr=Te_ken_err[5], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,Te_ken[6],yerr=Te_ken_err[6], color='darkviolet',fmt='o')
            ax.errorbar(radius[33]/a_eff,Te_ken[7],yerr=Te_ken_err[7], color='darkviolet',fmt='o')
            ax.errorbar(radius[34]/a_eff,Te_ken[8],yerr=Te_ken_err[8], color='darkviolet',fmt='o')
            ax.errorbar(radius[35]/a_eff,Te_ken[9],yerr=Te_ken_err[9], color='darkviolet',fmt='o')
            ax.set_ylabel( r'$T_\mathrm{e}[\mathrm{eV}]$')
            ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
            ax.set_ylim(0,100)
            ax.set_xlim(0.9,1.2)
            ax.grid()    
     
        if measurerand=='Te_sorted':
            pCart2=[ken[:,0],ken[:,1],ken[:,2]]
            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart2 )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), g )
            sorting=np.zeros([np.shape(iTor)[0],4])
            sorting[:,0]=Te_ken
            sorting[:,1]=iTor[0:20]
            sorting[:,2]=iPol[0:20]
            sorting[:,3]=iRad[0:20]
            Te_sorted=sorting[np.lexsort((sorting[:, 0], ))]
            
            Tes1=Te_sorted[:,1]
            Tes2=Te_sorted[:,2]
            Tes3=Te_sorted[:,3]            
            dummy_Tes=T('Te')[0][zip(Tes1),zip(Tes2),zip(Tes3)]
            wasdalos=zero_to_nan(dummy_Tes)
            
            for i in range(0, 20):
                ax.plot(Te_sorted[i,0],
                (wasdalos[i]),'o',color='#ffa500')
            dots_x=np.array([0,100])
            dots_y=np.array([0,100])
            ax.plot(dots_x,dots_y, '-',color='#ffa500')
            ax.set_ylim(0,100)
            ax.set_xlim(0,100)
            ax.set_ylabel( r'$T_\mathrm{e,EMC3}[\mathrm{eV}]$')
            ax.set_xlabel( r'$T_\mathrm{e,EXP}[\mathrm{eV}]$')
            ax.grid()              
            
            
            f2, ax2 = plt.subplots()    
            ax2.plot(range(1,21),Te_sorted[:,0],'o',color='#ffa500')
            dummy_dens=T('Te')[0][zip(Tes1),zip(Tes2),zip(Tes3)]
            a=zero_to_nan(dummy_dens)
            self.a=a
            ax2.plot(np.array(a),'s',color='#ffa500')
            ax2.set_ylabel( r'$n_\mathrm{e,EMC3}[\mathrm{10^{18}/m^3}]$')
            ax2.set_xlabel( r'LP probes')       
            ax2.grid()  
            
            
            

        self.radius=radius    







class limdata_reff_converge(object):
    """ Returns the radial Profile as plot and as np.arrays

  Keyword arguments:
  measurand -- which measurand should be plotted (string)
  
  Stings are:
  'ne' -- input string to get density of the electrons plotted
  'Te' -- input string to get temperature of the electrons plotted
  'ne_up_down' -- input string to get density of the electrons plotted with seperation of up and down probeheads 
  'Te_up_down' -- input string to get temperature of the electrons plotted with seperation of up and down probeheads 
  
  TODO:
  Legend for the up and down probehead seperation must be implemented 

  ATTENTION:
  Red exp_results are corresponding to the probehead below z==0.21[m]
    """


    def __init__(self,measurerand,z_eff_downstream = 'none', z_eff_upstream = 'none'):
        #plt.ion()
        g=EMC3.Grid() 


        a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

#Boyd Data Load
        data_json=json.load(open('/home/j.cosfeld/runs/W7X_limiter/exp_results/boyd/LP20160308_23_L53_2k2.json'))
        if measurerand=='ne':
            data=np.array(data_json['ne18'])
        if measurerand=='Te':
            data=np.array(data_json['Te'])
        if measurerand=='ne_sep':
            data=np.array(data_json['ne18'])
        if measurerand=='Te_sep':
            data=np.array(data_json['Te'])
        if measurerand=='ne_modded':
            data=np.array(data_json['ne18'])
        if measurerand=='Te_modded':
            data=np.array(data_json['Te'])
        if measurerand=='ne_comb':
            data=np.array(data_json['ne18'])
        if measurerand=='delta':
            data=np.array(data_json['ne18'])
            
        max_data_boyd=np.zeros(data.shape[1])        
        errorbar_boyd=np.zeros(data.shape[1])   
        
        
        for i in range (0,data.shape[1]):
            max_data_boyd[i]=np.nanmax(data[250,i])
            errorbar_boyd[i]=np.average(data[240:275,7]) 
            

        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')

        ne_sum=d_sum(g,d,0,0)
        
        ne_z_sum_downstream = d_z_sum(g,d,1,0)
        ne_z_sum_upstream   = d_z_sum(g,d,3,247)
        ne_sum_downstream=d_sum(g,d,1,0) 
        ne_sum_upstream=d_sum(g,d,3,247)
        
        ne_z_sum_downstream = d_z_sum_impu(g,d,1,0,'C')
        ne_sum_downstream=d_sum_impu(g,d,1,0,'C')
        m_eff_downstream=sum_eff_mass(g,d,1,0,'C')
        
        self.m_eff_downstream=m_eff_downstream
        
        z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
        
        x_emc3=radius/a_eff_emc3
        x_boyd=r_eff_boyd/a_eff_exp
        
        self.x_emc3=x_emc3
        self.x_boyd=x_boyd
        overlap=[]
        for i in range(0,x_boyd.shape[0]):
            overlap.append(find_nearest(x_emc3,x_boyd[i])[1])
        
        self.overlap=overlap
        self.max_data_boyd=max_data_boyd
        
        f, ax = plt.subplots()

        if measurerand=='Te_modded':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(T('Te')[0][1,0,:]),'--',linewidth=2,color="blue")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*z_eff_downstream[overlap],yerr=errorbar_boyd-max_data_boyd, color='orange',fmt='o')
            ax.set_ylabel( r'$T_\mathrm{e}[\mathrm{eV}]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,120)
            ax.set_xlim(0.8,1.3)
            ax.grid()
        
        if measurerand=='Te':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(T('Te')[0][1,0,:]),'--',linewidth=2,color="blue")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd,yerr=errorbar_boyd-max_data_boyd, color='orange',fmt='o')
            ax.set_ylabel( r'$T_\mathrm{e}[\mathrm{eV}]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,75)
            ax.set_xlim(0.8,1.3)
            ax.grid()

        if measurerand=='ne_modded':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6),'--',linewidth=2,color="r")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*1e18*z_eff_downstream[overlap],yerr=errorbar_boyd*1e18-max_data_boyd*1e18, color='darkviolet',fmt='s')
            ax.set_ylabel( r'$n_\mathrm{e}[\mathrm{1/m^3}]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,1*1e19)
            #ax.set_xlim(0.97,1.25)
            ax.grid()
            ne_boyd=max_data_boyd*1e18              

        if measurerand=='ne':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6),'--',linewidth=2,color="r")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*1e18,yerr=errorbar_boyd*1e18-max_data_boyd*1e18, color='darkviolet',fmt='o')
            ax.set_ylabel( r'$n_\mathrm{e}[\mathrm{1/m^3}]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,1.5*1e19)
            ax.set_xlim(0.8,1.3)
            ax.grid()
            ne_boyd=max_data_boyd*1e18              
        if measurerand=='Te_sep':
            coords=np.array(data_json['info']['coords'])
            coordsR=np.sqrt(coords[:,0]**2+coords[:,1]**2)
            coordsZup=np.where(coordsR*100-coordsR[10]*100>0)[0]
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(T('Te')[0][1,0,:]),'x',linewidth=2,color="blue")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd,yerr=errorbar_boyd-max_data_boyd, color='orange',fmt='o')
            for i in range (0,coordsZup.shape[0]):
                ax.errorbar(r_eff_boyd[i]/a_eff_exp,max_data_boyd[i],yerr=errorbar_boyd[i]-max_data_boyd[i], color='red',fmt='o')
            ax.set_ylabel( r'$T_e[eV]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,180)
            ax.grid()
        if measurerand=='ne_sep':
            coords=np.array(data_json['info']['coords'])
            coordsR=np.sqrt(coords[:,0]**2+coords[:,1]**2)
            coordsZup=np.where(coordsR*100-coordsR[10]*100>0)[0]
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6),'x',linewidth=2,color="r")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*1e18,yerr=errorbar_boyd*1e18-max_data_boyd*1e18, color='darkviolet',fmt='o')
            for i in range (0,coordsZup.shape[0]):
                ax.errorbar(r_eff_boyd[i]/a_eff_exp,max_data_boyd[i]*1e18,yerr=errorbar_boyd[i]*1e18-max_data_boyd[i]*1e18, color='blue',fmt='o')
            ax.set_ylabel( r'$n_e[1/m^3]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,2*1e19)
            ax.grid()
        if measurerand=='ne_comb':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan((d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6)),'x',linewidth=2,color="r")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*1e18*z_eff_downstream[overlap],yerr=errorbar_boyd*1e18*z_eff_downstream[overlap]-max_data_boyd*1e18*z_eff_downstream[overlap], color='darkviolet',fmt='o')
            ax.set_ylabel( r'$n_e[1/m^3]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,2*1e19)
            ax.grid()
            ne_boyd=max_data_boyd*1e18*z_eff_downstream[overlap]
        if measurerand=='delta':
            dummy_error=errorbar_boyd*1e18*z_eff_downstream[overlap]-max_data_boyd*1e18*z_eff_downstream[overlap]
            dummy=zero_to_nan((d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6))
            mean_emc3=np.nanmean(dummy[13:25])
            mean_boyd=np.nanmean(max_data_boyd*1e18*z_eff_downstream[overlap])
            mean_error=np.nanmean(dummy_error)/mean_boyd
            delta=abs(mean_emc3-mean_boyd)/mean_boyd
            print mean_emc3
            print mean_boyd
            print delta
            print mean_error

        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()
        
        ne=zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6)
        self.z_eff_downstream=z_eff_downstream
        self.z_eff_upstream=z_eff_upstream
        self.ne_boyd=max_data_boyd*1e18
        print np.nanmax(max_data_boyd*1e18)
        #self.ne_boyd=ne_boyd
        self.radius=radius
        self.r_eff_boyd=r_eff_boyd
        self.a_eff_emc3=a_eff_emc3
        self.a_eff_exp=a_eff_exp
        self.ne=ne
        self.ne_emc3=(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6)

        self.z_eff_downstream=z_eff_downstream
        self.z_eff_upstream=z_eff_upstream
        self.TeEMC3=zero_to_nan(T('Te')[0][1,0,:])
        self.TiEMC3=zero_to_nan(T('Ti')[0][1,0,:])
        
        
class limdata_reff_converge2(object):
    """ Returns the radial Profile as plot and as np.arrays

  Keyword arguments:
  measurand -- which measurand should be plotted (string)
  
  Stings are:
  'ne' -- input string to get density of the electrons plotted
  'Te' -- input string to get temperature of the electrons plotted
  'ne_up_down' -- input string to get density of the electrons plotted with seperation of up and down probeheads 
  'Te_up_down' -- input string to get temperature of the electrons plotted with seperation of up and down probeheads 
  
  TODO:
  Legend for the up and down probehead seperation must be implemented 

  ATTENTION:
  Red exp_results are corresponding to the probehead below z==0.21[m]
    """


    def __init__(self,measurerand,z_eff_downstream_old = 'none', TeEMC3old='none', TiEMC3old='none',m_eff_downstream_old = 'none'):
        #plt.ion()
        g=EMC3.Grid() 


        a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

#Boyd Data Load
        data_json=json.load(open('/home/j.cosfeld/runs/W7X_limiter/exp_results/boyd/LP20160308_23_L53_2k2.json'))

        if measurerand=='ne_modded':
            data=np.array(data_json['ne18'])

            
        max_data_boyd=np.zeros(data.shape[1])        
        errorbar_boyd=np.zeros(data.shape[1])   
        
        
        for i in range (0,data.shape[1]):
            max_data_boyd[i]=np.nanmax(data[250,i])
            errorbar_boyd[i]=np.average(data[240:275,7]) 
            

        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')

        ne_sum=d_sum(g,d,0,0)
        
        ne_z_sum_downstream = d_z_sum(g,d,1,0)
        ne_z_sum_upstream   = d_z_sum(g,d,3,247)
        ne_sum_downstream=d_sum(g,d,1,0) 
        ne_sum_upstream=d_sum(g,d,3,247)
        
        ne_z_sum_downstream = d_z_sum_impu(g,d,1,0,'C')
        ne_sum_downstream=d_sum_impu(g,d,1,0,'C')
        m_eff_downstream=sum_eff_mass(g,d,1,0,'C')
        
        z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream

        TeEMC3old=np.array(TeEMC3old)
        TiEMC3old=np.array(TiEMC3old)
        
        
        
        x_emc3=radius/a_eff_emc3
        x_boyd=r_eff_boyd/a_eff_exp
        
        self.x_emc3=x_emc3
        self.x_boyd=x_boyd
        overlap=[]
        for i in range(0,x_boyd.shape[0]):
            overlap.append(find_nearest(x_emc3,x_boyd[i])[1])
        
        if m_eff_downstream_old=='none':
            m_eff_downstream_old=np.zeros(len(m_eff_downstream))
            m_eff_downstream_old[:]=1.66*1e-27      
        m_eff_downstream_old=np.array(m_eff_downstream_old)
        self.m_eff_downstream=m_eff_downstream
        eps_oben2=(min(TeEMC3old[overlap])+min(z_eff_downstream_old[overlap])*min(TiEMC3old[overlap]))   *min(m_eff_downstream[overlap])
        self.eps_oben2=eps_oben2
        self.z_eff_downstream=z_eff_downstream
        self.TeEMC3=zero_to_nan(T('Te')[0][1,0,:])
        self.TiEMC3=zero_to_nan(T('Ti')[0][1,0,:])
        Te=T('Te')[0][1,0,:]
        Ti=T('Ti')[0][1,0,:]
        eps_unten2=(Te[overlap]+max(z_eff_downstream[overlap])*(Ti[overlap])) *max(m_eff_downstream_old[overlap])
        eps2=eps_oben2/eps_unten2
        self.eps_unten2=eps_unten2
        self.eps2=eps2

        
        z_eff_downstream[overlap]
        
        
        
        self.overlap=overlap
        self.max_data_boyd=max_data_boyd
        
        f, ax = plt.subplots()

        if measurerand=='ne_modded':
            ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6),'--',linewidth=2,color="r")
            ax.errorbar(r_eff_boyd/a_eff_exp,max_data_boyd*1e18*eps2,yerr=errorbar_boyd*1e18-max_data_boyd*1e18, color='darkviolet',fmt='s')
            ax.set_ylabel( r'$n_\mathrm{e}[\mathrm{1/m^3}]$')
            ax.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            ax.set_ylim(0,0.5*1e19)
            ax.set_xlim(0.97,1.25)
            ax.grid()
            ne_boyd=max_data_boyd*1e18              
       
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()
        
        ne=zero_to_nan(d('n1')[0][1,0,:]*1e6+d('n2')[0][1,0,:]*1e6)
        self.z_eff_downstream=z_eff_downstream
        self.ne_boyd=max_data_boyd*1e18*eps2
        print np.nanmax(max_data_boyd*1e18*eps2)
        self.radius=radius
        self.r_eff_boyd=r_eff_boyd
        self.a_eff_emc3=a_eff_emc3
        self.a_eff_exp=a_eff_exp
        self.ne=ne


        
        
        
class z_eff(object):
    def __init__(self):
        #plt.ion()
        g=EMC3.Grid() 

        a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')

        if len(d)>1:
            ne_sum=d_sum(g,d,0,0)

            print 'downstream'
            ne_z_sum_downstream = d_z_sum_impu(g,d,1,0,'C')
            ne_sum_downstream=d_sum_impu(g,d,1,0,'C')
            m_eff_downstream=sum_eff_mass(g,d,1,0,'C')
            print 'upstream'
            ne_z_sum_upstream   = d_z_sum_impu(g,d,3,247,'C')
            ne_sum_upstream=d_sum_impu(g,d,3,247,'C')                
            m_eff_upstream=sum_eff_mass(g,d,3,247,'C')
            if len(d)>7:
                print 'downstream'
                ne_z_sum_downstream_C = ne_z_sum_downstream
                m_eff_downstream_C=m_eff_downstream
                ne_z_sum_downstream_O = d_z_sum_impu(g,d,1,0,'O')
                m_eff_downstream_O=sum_eff_mass(g,d,1,0,'O')
                ne_sum_downstream_C = ne_sum_downstream
                ne_sum_downstream_O = d_sum_impu(g,d,1,0,'O')
                z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                print 'upstream'
                ne_z_sum_upstream_C   = ne_z_sum_upstream
                m_eff_upstream_C=m_eff_upstream
                ne_z_sum_upstream_O   = d_z_sum_impu(g,d,3,247,'O')
                m_eff_upstream_O=sum_eff_mass(g,d,3,247,'O')
                ne_sum_upstream_C   = ne_sum_upstream
                ne_sum_upstream_O   = d_sum_impu(g,d,3,247,'O')
                z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                
                ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                self.dummy1=ne_z_sum_downstream
                ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                self.dummy2=ne_sum_downstream
                ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                
            z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
            z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream   


        pathH = '.'

        Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
        GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
        PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
        PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
        niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
        TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


        # Multi-purpose manipulator tip plunge minimum and maximum location
        p1Cart = [ -556.276339, -210.053078, -17.736684 ]
        p2Cart = [ -587.407899, -222.425466, -17.736684 ]
        pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                    np.linspace( p1Cart[1], p2Cart[1], 40 ),
                    np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

        Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
        iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )
        RLimTip = 0.5 * ( GH.zones[0].R[iTor,iPol,iRad][PlatesH.ids[0][0,0,iRad]].min() +
                            GH.zones[0].R[iTor,iPol,iRad][PlatesH.ids[0][0,0,iRad]==False].max() )

            
            
        f, ax = plt.subplots()
        ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(z_eff_downstream),'-',linewidth=2,color="r",label='down-stream')
        ax.plot((radius/a_eff_emc3)[0:-1],zero_to_nan(z_eff_upstream),'--',linewidth=2,color="b",label='up-stream')
        ax.set_ylabel( r'$Z_\mathrm{eff}$')
        ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        ax.set_ylim(bottom=0)
        ax.set_xlim(1,1.25)
        ax.grid()
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()

        f, ax = plt.subplots()
        ax.plot((radius/a_eff_emc3)[0:-1],np.array((m_eff_downstream))/1.66e-27,'-',linewidth=2,color="r",label='down-stream')
        ax.plot((radius/a_eff_emc3)[0:-1],np.array(zero_to_nan(m_eff_upstream))/1.66e-27,'--',linewidth=2,color="b",label='up-stream')
        ax.set_ylabel( r'$m_\mathrm{eff}/m_\mathrm{H}$')
        ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        ax.set_ylim(bottom=0)
        ax.set_xlim(1,1.25)
        ax.set_ylim(1,3)
        ax.grid()
        plt.setp(ax.get_xticklabels()[1::2], visible=False)
        plt.setp(ax.get_yticklabels()[1::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()



        self.z_eff_downstream=z_eff_downstream
        self.z_eff_upstream=z_eff_upstream
        self.m_eff_downstream=m_eff_downstream
        self.m_eff_upstream=m_eff_upstream
        self.radius=radius
        self.a_eff_emc3=a_eff_emc3
    

        markerarry=['.','1','2','3','4','x','1','2','3','4','x'];
        colorarray=['b','k','k','k','k','k','k','k','k','k','k'];
        labelarry=[r'$n_\mathrm{H^{+}}$',r'$n_\mathrm{C^{1+}}$',r'$n_\mathrm{C^{2+}}$',r'$n_\mathrm{C^{3+}}$',r'$n_\mathrm{C^{4+}}$',r'$n_\mathrm{C^{5+}}$',r'$n_\mathrm{C^{6+}}$'
                   ,r'$n_\mathrm{C^{7+}}$',r'$n_\mathrm{C^{8+}}$',r'$n_\mathrm{C^{9+}}$',r'$n_\mathrm{C^{10+}}$',r'$n_\mathrm{C^{11+}}$',r'$n_\mathrm{C^{12+}}$'
                   ]

        max_d=np.zeros(len(d))
        f, ax = plt.subplots()
        for m in range(2,len(d)):
            rank=str(m)
            max_d[m]=max((d('n'+rank)[0][1,0,:])/(d('n1')[0][1,0,:]))
            max_d_pos=np.argmax(max_d)
        rank=str(max_d_pos)
        ax.plot((radius/a_eff_emc3)[0:-1],(d('n' + rank)[0][1,0,:])/(d('n1')[0][1,0,:]),linewidth=2, label=labelarry[max_d_pos],color='blue')            
            
        ax.set_ylabel( r'$\mathrm{max}(n_\mathrm{x})/n_\mathrm{H}$')
        ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        #ax.set_title('maximum contribution by impurity state x')
        ax.set_ylim(bottom=0)
        ax.set_xlim(1,1.25)
        ax.grid()
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        ax.legend(fontsize=12,loc=2)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()


        f, ax = plt.subplots()
        ax.plot((radius/a_eff_emc3)[0:-1],(d_sum_avg(g,d))/(np.average(d('n1')[0],(0,1),g.zones[0].volume)),linewidth=2, label=labelarry[max_d_pos],color='blue')            
        ax.set_ylabel( r'$\overline{n_\mathrm{e,imp}}/\overline{n_\mathrm{e}}$')
        ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        ax.set_ylim(bottom=0)
        ax.set_xlim(1,1.25)
        ax.grid()
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()



      
class z_eff_div(object):
    def __init__(self):
        #plt.ion()
        g=EMC3.Grid() 

        a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)
        a_eff_emc3=radius[11]
        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')

        if len(d)>1:
            ne_sum=d_sum(g,d,0,0)

            print 'downstream'
            ne_z_sum_downstream = d_z_sum_impu(g,d,17,102,'C')
            ne_sum_downstream=d_sum_impu(g,d,17,102,'C')
            m_eff_downstream=sum_eff_mass(g,d,17,102,'C')
            print 'upstream'
            ne_z_sum_upstream   = d_z_sum_impu(g,d,30,237,'C')
            ne_sum_upstream=d_sum_impu(g,d,30,237,'C')                
            m_eff_upstream=sum_eff_mass(g,d,30,237,'C')
            if len(d)>7:
                print 'downstream'
                ne_z_sum_downstream_C = ne_z_sum_downstream
                m_eff_downstream_C=m_eff_downstream
                ne_z_sum_downstream_O = d_z_sum_impu(g,d,17,102,'O')
                m_eff_downstream_O=sum_eff_mass(g,d,17,102,'O')
                ne_sum_downstream_C = ne_sum_downstream
                ne_sum_downstream_O = d_sum_impu(g,d,17,102,'O')
                z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                print 'upstream'
                ne_z_sum_upstream_C   = ne_z_sum_upstream
                m_eff_upstream_C=m_eff_upstream
                ne_z_sum_upstream_O   = d_z_sum_impu(g,d,30,237,'O')
                m_eff_upstream_O=sum_eff_mass(g,d,30,237,'O')
                ne_sum_upstream_C   = ne_sum_upstream
                ne_sum_upstream_O   = d_sum_impu(g,d,30,237,'O')
                z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                
                ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                self.dummy1=ne_z_sum_downstream
                ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                self.dummy2=ne_sum_downstream
                ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                
            z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
            z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream   


        pathH = '.'

        Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
        GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
        PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
        PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
        niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
        TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


            
            
        f, ax = plt.subplots()
        #ax.plot((radius/a_eff_emc3)[::-1][:-1],zero_to_nan(z_eff_downstream),'-',linewidth=2,color="r",label='down-stream')
        ax.plot((radius/a_eff_emc3)[::-1][:-1],zero_to_nan(z_eff_upstream),'--',linewidth=2,color="b",label='up-stream')
        ax.set_ylabel( r'$Z_\mathrm{eff}$')
        ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        ax.set_ylim(bottom=0)
        #ax.set_xlim(1,1.17)
        ax.grid()
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()

        self.z_eff_downstream=z_eff_downstream
        self.z_eff_upstream=z_eff_upstream
        self.m_eff_downstream=m_eff_downstream
        self.m_eff_upstream=m_eff_upstream
        self.radius=radius
        self.a_eff_emc3=a_eff_emc3
    

        markerarry=['.','1','2','3','4','x','1','2','3','4','x'];
        colorarray=['b','k','k','k','k','k','k','k','k','k','k'];
        labelarry=[r'$n_\mathrm{H^{+}}$',r'$n_\mathrm{C^{1+}}$',r'$n_\mathrm{C^{2+}}$',r'$n_\mathrm{C^{3+}}$',r'$n_\mathrm{C^{4+}}$',r'$n_\mathrm{C^{5+}}$',r'$n_\mathrm{C^{6+}}$'
                   ,r'$n_\mathrm{C^{7+}}$',r'$n_\mathrm{C^{8+}}$',r'$n_\mathrm{C^{9+}}$',r'$n_\mathrm{C^{10+}}$',r'$n_\mathrm{C^{11+}}$',r'$n_\mathrm{C^{12+}}$'
                   ]

        max_d=np.zeros(len(d))
        f, ax = plt.subplots()
        for m in range(2,len(d)):
            rank=str(m)
            max_d[m]=max((d('n'+rank)[0][17,102,:])/(d('n1')[0][17,102,:]))
            max_d_pos=np.argmax(max_d)
        rank=str(max_d_pos)
        ax.plot((radius/a_eff_emc3)[0:-1][::-1],(d('n' + rank)[0][17,102,:])/(d('n1')[0][17,102,:]),linewidth=2, label=labelarry[max_d_pos],color='blue')            
            
        ax.set_ylabel( r'$\mathrm{max}(n_\mathrm{x})/n_\mathrm{H}$')
        ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        #ax.set_title('maximum contribution by impurity state x')
        ax.set_ylim(bottom=0)
        #ax.set_xlim(1,1.17)
        ax.grid()
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        ax.legend(fontsize=12,loc=2)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()


        #f, ax = plt.subplots()
        #ax.plot((radius/a_eff_emc3)[0:-1],(d_sum_avg(g,d))/(np.average(d('n1')[0],(0,1),g.zones[0].volume)),linewidth=2, label=labelarry[max_d_pos],color='blue')            
        #ax.set_ylabel( r'$\overline{n_\mathrm{e,imp}}/\overline{n_\mathrm{e}}$')
        #ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        #ax.set_ylim(bottom=0)
        ##ax.set_xlim(1,1.17)
        #ax.grid()
        #plt.setp(ax.get_xticklabels()[::2], visible=False)
        #plt.setp(ax.get_yticklabels()[::2], visible=False)
        #plt.rcParams.update({'font.size': 18})
        #f.tight_layout()

        f, ax = plt.subplots()
        #ax.plot((radius/a_eff_emc3)[11:32],np.array((m_eff_downstream[11:32]))/1.66e-27,'-',linewidth=2,color="r",label='down-stream')
        ax.plot((radius/a_eff_emc3)[::-1][:-1],np.array(zero_to_nan(m_eff_upstream))/1.66e-27,'--',linewidth=2,color="b",label='up-stream')
        ax.set_ylabel( r'$m_\mathrm{eff}/m_\mathrm{H}$')
        ax.set_xlabel( r'$r_\mathrm{eff}/a_\mathrm{eff}$')
        ax.set_ylim(bottom=0)
        #ax.set_xlim(1,1.17)
        #ax.set_ylim(0,8)
        ax.grid()
        plt.setp(ax.get_xticklabels()[1::2], visible=False)
        plt.setp(ax.get_yticklabels()[1::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()

class z_eff_pavone(object):
    def __init__(self):
        #plt.ion()
        g=EMC3.Grid() 

        a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)
        T=EMC3.PlasmaField('TE_TI')
        d=EMC3.PlasmaField('DENSITY')
        z=g.zones[0]
        reff=z.Reff[0:-1]

        if len(d)>1:
            ne_sum=d_sum(g,d,0,0)

            print 'upstream'
            ne_z_sum_upstream   = d_z_sum_impu(g,d,55,300,'C')
            ne_sum_upstream=d_sum_impu(g,d,55,300,'C')                
            m_eff_upstream=sum_eff_mass(g,d,55,300,'C')
            if len(d)>7:
                print 'upstream'
                ne_z_sum_upstream_C   = ne_z_sum_upstream
                m_eff_upstream_C=m_eff_upstream
                ne_z_sum_upstream_O   = d_z_sum_impu(g,d,55,300,'O')
                m_eff_upstream_O=sum_eff_mass(g,d,55,300,'O')
                ne_sum_upstream_C   = ne_sum_upstream
                ne_sum_upstream_O   = d_sum_impu(g,d,55,300,'O')
                z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                
                ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                
            z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream   


        pathH = '.'

        Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
        GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
        PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
        PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
        niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
        TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


            
            
        f, ax = plt.subplots()
        ax.plot((reff),zero_to_nan(z_eff_upstream),'--',linewidth=2,color="b",label='up-stream')
        ax.set_ylabel( r'$Z_\mathrm{eff}$')
        ax.set_xlabel( r'$r_\mathrm{eff}$')
        ax.set_ylim(bottom=0)
        #ax.set_xlim(1,1.17)
        ax.grid()
        plt.setp(ax.get_xticklabels()[::2], visible=False)
        plt.setp(ax.get_yticklabels()[::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()

        self.z_eff_upstream=z_eff_upstream
        self.radius=radius
        self.a_eff_emc3=a_eff_emc3
    

        markerarry=['.','1','2','3','4','x','1','2','3','4','x'];
        colorarray=['b','k','k','k','k','k','k','k','k','k','k'];
        labelarry=[r'$n_\mathrm{H^{+}}$',r'$n_\mathrm{C^{1+}}$',r'$n_\mathrm{C^{2+}}$',r'$n_\mathrm{C^{3+}}$',r'$n_\mathrm{C^{4+}}$',r'$n_\mathrm{C^{5+}}$',r'$n_\mathrm{C^{6+}}$'
                   ,r'$n_\mathrm{C^{7+}}$',r'$n_\mathrm{C^{8+}}$',r'$n_\mathrm{C^{9+}}$',r'$n_\mathrm{C^{10+}}$',r'$n_\mathrm{C^{11+}}$',r'$n_\mathrm{C^{12+}}$'
                   ]

        f, ax = plt.subplots()
        ax.plot((reff),np.array(zero_to_nan(m_eff_upstream))/1.66e-27,'--',linewidth=2,color="b",label='up-stream')
        ax.set_ylabel( r'$m_\mathrm{eff}/m_\mathrm{H}$')
        ax.set_xlabel( r'$r_\mathrm{eff}$')
        ax.set_ylim(bottom=0)
        #ax.set_xlim(1,1.17)
        #ax.set_ylim(0,8)
        ax.grid()
        plt.setp(ax.get_xticklabels()[1::2], visible=False)
        plt.setp(ax.get_yticklabels()[1::2], visible=False)
        plt.rcParams.update({'font.size': 18})
        f.tight_layout()


class Cut_Pol_Avg_Logscale_Dens(object):
    def __init__(self, measurerand):
    #
    #Gives the normed Density of Carbon impurities (normed at every radial position)
    #
        #plt.ion()

        plt.matplotlib.rcParams.update({'font.size': 18})

        g=EMC3.Grid()
        
        a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

    #
    # Plotting
    #
        markerarry=['.','1','2','3','4','x','I'];
        colorarray=['b','k','k','k','k','k','k'];
        labelarry=[r'$n_\mathrm{H^{+}}$',r'$n_\mathrm{C^{1+}}$',r'$n_\mathrm{C^{2+}}$',r'$n_\mathrm{C^{3+}}$',r'$n_\mathrm{C^{4+}}$',r'$n_\mathrm{C^{5+}}$',r'$n_\mathrm{C^{6+}}$']

    #
    # Case check if array index in plot is wanted or not
    # Case check which measureand is wanted by the user 
#
        if measurerand=='n':
            data=EMC3.NeutralParameters("input_modded.n0g")
            d_eirene=np.zeros(7)
            for i in range(0, 7):
                d_eirene[i]=data.cells[i]['n']
            d=EMC3.PlasmaField('DENSITY')
            matrix_length=g.geometry.NSurfaces[0][0]
            matrix_high=len(d)
            avg=np.zeros(shape=(matrix_high,matrix_length-1))
#
# Calculate the average density for each ionization state along the radial axis 
#
            for m in range(1,len(d)):
                rank=str(m)
                avg[m,:]=np.average(d('n'+ rank)[0],(0,1),g.zones[0].volume)
#
# Plot is generated in figure f 
#
            f = plt.figure()
            ax1 = f.add_axes([0.2, 0.2, 0.55, 0.75])


            for m in range(1,len(d)):
                    ax1.plot((radius/a_eff_emc3)[8:-2],avg[m,:][8:-1]/avg[1,:][8:-1],'x',linewidth=2,marker=markerarry[m-1], label=labelarry[m-1],color=colorarray[m-1])

            plt.yscale('log', nonposy='clip')
            ax1.set_xlim( 0.7, 1.3 )
            ax1.set_ylabel( r'$n_i/\mathrm{max}(n_\mathrm{H^+})$')
            ax1.set_xlabel( r'$r_\mathrm{eff.}/a_\mathrm{eff.}$')
            lines1, labels1 = ax1.get_legend_handles_labels()
            ax1.grid()
            ax1.legend( lines1, labels1, bbox_to_anchor=(1.05, 1), borderaxespad=0.,
                        loc = 2, labelspacing = 0.0, handletextpad = 1.0, handlelength = 1.5)
            topbot=ax1.get_ylim()
            #plt.tight_layout()


        self.avg=avg
        self.radius=radius
        self.a_eff_emc3=a_eff_emc3
        self.a_eff_exp=a_eff_exp

class Plot_all_Profiles_Experiment(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  

  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self,impu='C',case='modded',z_eff_downstream = 'none', z_eff_upstream = 'none'):
            plt.rcParams.update({'font.size': 24})
            g=EMC3.Grid() 
            
            a_eff_emc3,a_eff_exp,radius,r_eff_drews=get_a_r_eff2(g)
            
            
            self.a_eff_emc3=a_eff_emc3
            self.a_eff_exp=a_eff_exp
            self.radius=radius
            self.r_eff_drews=r_eff_drews

            pathH = '.'

            Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
            GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
            PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
            PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
            niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
            TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


            # Multi-purpose manipulator tip plunge minimum and maximum location
            p1Cart = [ -556.276339, -210.053078, -17.736684 ]
            p2Cart = [ -587.407899, -222.425466, -17.736684 ]
            pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                      np.linspace( p1Cart[1], p2Cart[1], 40 ),
                      np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )

            RLimTip = 0.5 * ( GH.zones[0].R[iTor,iPol,iRad][PlatesH.ids[0][0,0,iRad]].min() +
                              GH.zones[0].R[iTor,iPol,iRad][PlatesH.ids[0][0,0,iRad]==False].max() )
            n_e_imp_C=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            n_e_imp_O=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            if impu=='C':
                for m in range(1,len(niH.keys())):
                    dummy_false='n'+ str(m) in niH.keys()
                    if dummy_false==False:
                        n_e_imp_C+=n_imp*0                        
                    if dummy_false==True:    
                        n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                        n_e_imp_C+=n_imp*m
                n_e_imp=n_e_imp_C       
                print 'Carbon'
            counter=1    
            if impu!='C':
                for m in range(1,8):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                for m in range(8,len(niH.keys())+1):
                    counter+=1
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_O+=n_imp*counter
                n_e_imp=n_e_imp_C+n_e_imp_O
                print 'Carbon plus Oxygen'
                self.n_e_imp_C=n_e_imp_C
                self.n_e_imp_O=n_e_imp_O
            for Te, Ti, n, TiLabel in zip( [ TeTiH('Te')[0] ],
                                                  [ TeTiH('Ti')[0] ],
                                                  [ niH('n1')[0] ],
                                                  [ r'$T_\mathrm{H^{+}}$']):
              fig, ax1 = plt.subplots(  )
#a radial shift is implemented here, because of the radial offset of the
#langmuir probe head              
              shift=0
              
              #self.n1=niH('n1')[0][iTor,iPol,iRad]
              #self.n2=niH('n2')[0][iTor,iPol,iRad]
              #self.n3=niH('n3')[0][iTor,iPol,iRad]
              #self.n4=niH('n4')[0][iTor,iPol,iRad]
              #self.n5=niH('n5')[0][iTor,iPol,iRad]
              #self.n6=niH('n6')[0][iTor,iPol,iRad]
              #self.n7=niH('n7')[0][iTor,iPol,iRad]
              
              data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/base/profiles.mat')


              Te_drews=data_drews['te']
              Te22=Te_drews[:,0]
              Te23=Te_drews[:,1]
              Te24=Te_drews[:,2]
              ne_drews=data_drews['ne']
              ne22=ne_drews[:,0]*1e18
              ne23=ne_drews[:,1]*1e18
              ne24=ne_drews[:,2]*1e18
              x=data_drews['x']
              x22=x[:,0]/10
              x23=x[:,1]/10
              x24=x[:,2]/10
              y=data_drews['y']
              y22=y[:,0]/10
              y23=y[:,1]/10
              y24=y[:,2]/10
              radius22=np.sqrt(x22**2+y22**2)
              radius23=np.sqrt(x23**2+y23**2)
              radius24=np.sqrt(x24**2+y24**2)
              ne_up_err=data_drews['ne_error_1']
              ne_up_err22=ne_up_err[:,0]
              ne_up_err23=ne_up_err[:,1]
              ne_up_err24=ne_up_err[:,2]
              ne_down_err=data_drews['ne_error_2']
              ne_down_err22=ne_down_err[:,0]
              ne_down_err23=ne_down_err[:,1]              
              ne_down_err24=ne_down_err[:,2]              
              error_ne_22=np.array([ne_down_err22,ne_up_err22])
              error_ne_23=np.array([ne_down_err23,ne_up_err23])
              error_ne_24=np.array([ne_down_err24,ne_up_err24])
              te_err=data_drews['te_error']
              te_err22=te_err[:,0]
              te_err23=te_err[:,1]              
              te_err24=te_err[:,2]             
              
              self.radius22=radius22
              self.radius23=radius23
              self.radius24=radius24
              
              U22 = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/Upin/160308022voltages.mat')
              U23 = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/Upin/160308023voltages.mat')
              U24 = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/Upin/160308024voltages.mat')
              
              mi=1.66*1e-27*1.00794
              e=1.602e-19
              Aeff=2e-3*3e-3*2

              self.U22=U22
              self.U23=U23
              self.U24=U24
              
              self.ne_22_unmodded=ne22*1e-18
              self.ne_23_unmodded=ne23*1e-18
              self.ne_24_unmodded=ne24*1e-18

              T=EMC3.PlasmaField('TE_TI')
              d=EMC3.PlasmaField('DENSITY')

              if len(d)>1:
                ne_sum=d_sum(g,d,0,0)

                print 'downstream'
                ne_z_sum_downstream = d_z_sum_impu(g,d,1,0,'C')
                ne_sum_downstream=d_sum_impu(g,d,1,0,'C')
                m_eff_downstream=sum_eff_mass(g,d,1,0,'C')
                print 'upstream'
                ne_z_sum_upstream   = d_z_sum_impu(g,d,28,247,'C')
                ne_sum_upstream=d_sum_impu(g,d,28,247,'C')                
                m_eff_upstream=sum_eff_mass(g,d,28,247,'C')
                if len(d)>7:
                    print 'downstream'
                    ne_z_sum_downstream_C = ne_z_sum_downstream
                    m_eff_downstream_C=m_eff_downstream
                    ne_z_sum_downstream_O = d_z_sum_impu(g,d,1,0,'O')
                    m_eff_downstream_O=sum_eff_mass(g,d,1,0,'O')
                    ne_sum_downstream_C = ne_sum_downstream
                    ne_sum_downstream_O = d_sum_impu(g,d,1,0,'O')
                    z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                    z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                    print 'upstream'
                    ne_z_sum_upstream_C   = ne_z_sum_upstream
                    m_eff_upstream_C=m_eff_upstream
                    ne_z_sum_upstream_O   = d_z_sum_impu(g,d,3,247,'O')
                    m_eff_upstream_O=sum_eff_mass(g,d,3,247,'O')
                    ne_sum_upstream_C   = ne_sum_upstream
                    ne_sum_upstream_O   = d_sum_impu(g,d,3,247,'O')
                    z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                    z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                    
                    ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                    self.dummy1=ne_z_sum_downstream
                    ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                    self.dummy2=ne_sum_downstream
                    ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                    ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                    
                z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
                z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream        

                self.z_eff_downstream=z_eff_downstream
                self.z_eff_upstream=z_eff_upstream
                self.m_eff_upstream=m_eff_upstream
                self.m_eff_downstream=m_eff_downstream
                self.Rs=Rs
                self.iTor=iTor
                self.iPol=iPol
                self.iRad=iRad

              self.Rs=Rs
              self.RLimTip=RLimTip

              overlap=[]
              for i in range(0,radius22.shape[0]):
                  overlap.append(find_nearest(Rs,radius22[i])[1])
            
              self.overlap=overlap

              if case=='unmodded':
                ax1.plot( radius22-602, Te22,
                            '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                ax1.errorbar( radius22[::8]-602, Te22[::8],yerr=te_err22[::8], color='orange')
                ax1.plot( radius23-602, Te23,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius23[::8]-602, Te23[::8],yerr=te_err23[::8], color='orange')
                ax1.plot( radius24-602, Te24,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius24[::8]-602, Te24[::8],yerr=te_err24[::8], color='orange')
                
                ax1.plot( Rs - RLimTip, zero_to_nan(Te[iTor,iPol,iRad]),
                            'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )              
                
                ax2 = ax1.twinx()

                ax2.plot( radius22-602, ne22*1e-18,
                            '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                ax2.errorbar( radius22[::8]-602, ne22[::8]*1e-18,yerr=error_ne_22[:,::8], color='darkviolet')
                ax2.plot( radius23-602, ne23*1e-18,
                            '-', linewidth = 3,color='darkviolet' )
                ax2.errorbar( radius23[::8]-602, ne23[::8]*1e-18,yerr=error_ne_23[:,::8], color='darkviolet')              
                ax2.plot( radius24-602, ne24*1e-18,
                            '-', linewidth = 3,color='darkviolet' )
                ax2.errorbar( radius24[::8]-602, ne24[::8]*1e-18,yerr=error_ne_24[:,::8], color='darkviolet')              
                ax2.plot( Rs - RLimTip, (zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12,
                            'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ (cm)' )
                ax1.set_ylabel( r'$T$ (eV)' )
                ax1.set_xlim( 3, 10 )
                ax1.set_ylim( 0, 30 )
                ax2.set_ylabel( r'$n_\mathrm{e}$ ($10^{18}$ m$^{-3}$)' )
                ax2.set_ylim( 0, 3.5 )

                ax1.grid()
                lines1, labels1 = ax1.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                ax2.legend( lines1 + lines2, labels1 + labels2,
                            loc = 0, labelspacing = 0.0, handletextpad = 0.0, handlelength = 1.5 )

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )
                
                self.ne_22_unmodded=ne22*1e-18
                self.ne_23_unmodded=ne23*1e-18
                self.ne_24_unmodded=ne24*1e-18


              if case=='unmodded_single':
                f, ax1 = plt.subplots()  
                ax1.plot( radius22-602, Te22,
                            '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                ax1.errorbar( radius22[::8]-602, Te22[::8],yerr=te_err22[::8], color='orange')
                ax1.plot( radius23-602, Te23,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius23[::8]-602, Te23[::8],yerr=te_err23[::8], color='orange')
                ax1.plot( radius24-602, Te24,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius24[::8]-602, Te24[::8],yerr=te_err24[::8], color='orange')
                ax1.plot( Rs - RLimTip, zero_to_nan( Te [iTor,iPol,iRad]),
                            'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )           
                ax1.plot( Rs - RLimTip, zero_to_nan( Ti [iTor,iPol,iRad]),
                'r--', linewidth = 3, label = r'$T_\mathrm{e}$' )   
                
                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ (cm)' )
                ax1.set_ylabel( r'$T$ (eV)' )
                ax1.set_xlim( 3, 7 )
                ax1.set_ylim( 0, 40 )

                ax1.grid()

                plt.tight_layout( pad = 0.25 )
                
                f, ax2 = plt.subplots()

                ax2.plot( radius22-602, ne22*1e-18,
                            '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                ax2.errorbar( radius22[::8]-602, ne22[::8]*1e-18,yerr=error_ne_22[:,::8], color='darkviolet')
                ax2.plot( radius23-602, ne23*1e-18,
                            '-', linewidth = 3,color='darkviolet' )
                ax2.errorbar( radius23[::8]-602, ne23[::8]*1e-18,yerr=error_ne_23[:,::8], color='darkviolet')              
                ax2.plot( radius24-602, ne24*1e-18,
                            '-', linewidth = 3,color='darkviolet' )
                ax2.errorbar( radius24[::8]-602, ne24[::8]*1e-18,yerr=error_ne_24[:,::8], color='darkviolet')              
                ax2.plot( Rs - RLimTip, (zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12,
                            'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )

                ax2.set_xlim( 3, 7 )
                ax2.set_ylabel( r'$n_\mathrm{e}$ ($10^{18}$ m$^{-3}$)' )
                ax2.set_ylim( 0, 1.5 )
                ax2.grid()
                self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                self.TiEMC3=zero_to_nan(Ti[iTor,iPol,iRad])



                plt.tight_layout( pad = 0.25 )

                self.ne_22_unmodded=ne22*1e-18
                self.ne_23_unmodded=ne23*1e-18
                self.ne_24_unmodded=ne24*1e-18
				
              if case=='modded_single':
                ne_22_modded=ne22*1e-18*z_eff_upstream[overlap] 
                ne_23_modded=ne23*1e-18*z_eff_upstream[overlap] 
                ne_24_modded=ne24*1e-18*z_eff_upstream[overlap] 
                
                
                xaxis_plot=np.linspace((r_eff_drews/a_eff_exp)[-1],(r_eff_drews/a_eff_exp)[0],len(Te22))
                ax1.plot( xaxis_plot, Te22,
                            '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                #ax1.errorbar( xaxis_plot[::-8], Te22[::8],yerr=te_err22[::8], color='orange')
                #ax1.plot( xaxis_plot[::-1], Te23,
                            #'-', linewidth = 3,color='orange' )
                #ax1.errorbar( xaxis_plot[::-8], Te23[::8],yerr=te_err23[::8], color='orange')
                #ax1.plot( xaxis_plot[::-1], Te24,
                            #'-', linewidth = 3,color='orange' )
                #ax1.errorbar( xaxis_plot[::-8], Te24[::8],yerr=te_err24[::8], color='orange')
                
                #ax1.plot( (radius/a_eff_emc3)[1:-1], zero_to_nan(Te[iTor,iPol,iRad]),
                            #'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )              
                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ [cm]' )
                ax1.set_ylabel( r'$T$ [eV]' )
                ax1.set_xlim( 1, 1.1 )
                ax1.set_ylim( 0, 20 )
                ax1.grid()
                self.Te22=Te22
                plt.tight_layout( pad = 0.25 )

				
                f, ax2 = plt.subplots()
                
                ax2.plot( xaxis_plot[::-1], ne_22_modded,
                            '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                ax2.errorbar( xaxis_plot[::-8], ne_22_modded[::8],yerr=error_ne_22[:,::8], color='darkviolet')
                ax2.plot( xaxis_plot[::-1], ne_23_modded,
                            '-', linewidth = 3,color='darkviolet' )
                ax2.errorbar( xaxis_plot[::-8], ne_23_modded[::8],yerr=error_ne_23[:,::8], color='darkviolet')              
                ax2.plot( xaxis_plot[::-1], ne_24_modded,
                            '-', linewidth = 3,color='darkviolet' )
                ax2.errorbar( xaxis_plot[::-8], ne_24_modded[::8],yerr=error_ne_24[:,::8], color='darkviolet')              
                ax2.plot( (radius/a_eff_emc3)[1:-1], (zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12,
                            'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                ax2.set_xlabel( r'$R-R_\mathtm{LCFS}$ [cm] ' )
                ax2.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )		
                ax2.set_xlim( 1, 1.1 )
                ax2.set_ylim( 0, 2.5 )
                ax2.grid()

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )
                self.ne_22_modded=ne_22_modded
                self.ne_23_modded=ne_23_modded
                self.ne_24_modded=ne_24_modded
                self.neEMC3=(zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12

              if case=='mod_sep':
                ne_22_modded=ne22*1e-18*z_eff_upstream[overlap] 
                ne_23_modded=ne23*1e-18*z_eff_upstream[overlap] 
                ne_24_modded=ne24*1e-18*z_eff_upstream[overlap] 

                f, ax1 = plt.subplots()
                ax1.plot( radius22-602, Te22,
                            '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                ax1.errorbar( radius22[::8]-602, Te22[::8],yerr=te_err22[::8], color='orange')
                ax1.plot( radius23-602, Te23,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius23[::8]-602, Te23[::8],yerr=te_err23[::8], color='orange')
                ax1.plot( radius24-602, Te24,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius24[::8]-602, Te24[::8],yerr=te_err24[::8], color='orange')
                
                ax1.plot( Rs - RLimTip, zero_to_nan(Te[iTor,iPol,iRad]),
                            'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )              
                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ [cm]' )
                ax1.set_ylabel( r'$T$ [eV]' )
                ax1.set_xlim( 3, 8 )
                ax1.set_ylim( 0, 20 )

                ax1.grid()
                lines1, labels1 = ax1.get_legend_handles_labels()

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )

        
        
                f, ax1 = plt.subplots()
                ax1.plot( radius22-602, ne_22_modded,
                            '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                ax1.errorbar( radius22[::8]-602, ne_22_modded[::8],yerr=error_ne_22[:,::8], color='darkviolet')
                ax1.plot( radius23-602, ne_23_modded,
                            '-', linewidth = 3,color='darkviolet' )
                ax1.errorbar( radius23[::8]-602, ne_23_modded[::8],yerr=error_ne_23[:,::8], color='darkviolet')              
                ax1.plot( radius24-602, ne_24_modded,
                            '-', linewidth = 3,color='darkviolet' )
                ax1.errorbar( radius24[::8]-602, ne_24_modded[::8],yerr=error_ne_24[:,::8], color='darkviolet')              
                ax1.plot( Rs - RLimTip, (zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12,
                            'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )

                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ [cm]' )
                ax1.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}/\mathrm{m}^3$]' )
                #ax1.set_xlim( 3, 8 )
                #ax1.set_ylim( 0, 2.5 )

                ax1.grid()
                lines1, labels1 = ax1.get_legend_handles_labels()

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )

                self.ne_22_modded=ne_22_modded
                self.ne_23_modded=ne_23_modded
                self.ne_24_modded=ne_24_modded
                self.neemc3=(zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12
                
                self.radius22=radius22
                self.Rs=Rs
                self.RLimTip=RLimTip
                
                self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                self.TiEMC3=zero_to_nan(Ti[iTor,iPol,iRad])


              if case=='unmod_sep':
                ne_22_modded=ne22*1e-18*z_eff_upstream[overlap] 
                ne_23_modded=ne23*1e-18*z_eff_upstream[overlap] 
                ne_24_modded=ne24*1e-18*z_eff_upstream[overlap] 

                f, ax1 = plt.subplots()
                ax1.plot( radius22-602, Te22,
                            '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                ax1.errorbar( radius22[::8]-602, Te22[::8],yerr=te_err22[::8], color='orange')
                ax1.plot( radius23-602, Te23,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius23[::8]-602, Te23[::8],yerr=te_err23[::8], color='orange')
                ax1.plot( radius24-602, Te24,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius24[::8]-602, Te24[::8],yerr=te_err24[::8], color='orange')
                
                ax1.plot( Rs - RLimTip, zero_to_nan(Te[iTor,iPol,iRad]),
                            'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )              
                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ [cm]' )
                ax1.set_ylabel( r'$T$ [eV]' )
                ax1.set_xlim( 3, 8 )
                ax1.set_ylim( 0, 20 )

                ax1.grid()
                lines1, labels1 = ax1.get_legend_handles_labels()

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )

                f, ax1 = plt.subplots()
                ax1.plot( radius22-602, ne22,
                            '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                ax1.errorbar( radius22[::8]-602, ne22[::8]*1e-18,yerr=error_ne_22[:,::8], color='darkviolet')
                ax1.plot( radius23-602, ne23,
                            '-', linewidth = 3,color='darkviolet' )
                ax1.errorbar( radius23[::8]-602, ne23[::8]*1e-18,yerr=error_ne_23[:,::8], color='darkviolet')              
                ax1.plot( radius24-602, ne24,
                            '-', linewidth = 3,color='darkviolet' )
                ax1.errorbar( radius24[::8]-602, ne24[::8]*1e-18,yerr=error_ne_24[:,::8], color='darkviolet')              
                ax1.plot( Rs - RLimTip, (zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12,
                            'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ [cm]' )
                ax1.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}/\mathrm{m}^3$]' )
                ax1.set_xlim( 3, 8 )
                ax1.set_ylim( 0, 2.5 )

                ax1.grid()
                lines1, labels1 = ax1.get_legend_handles_labels()

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )

                self.ne22=ne22
                self.ne23=ne23
                self.ne24=ne24


            self.radius22=radius22
            self.radius23=radius23
            self.radius22=radius24
            self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
            self.Te_emc3=Te[iTor,iPol,iRad]
            self.Ti_emc3=Ti[iTor,iPol,iRad]
            self.TeEMC3complete=Te
            self.TiEMC3complete=Ti
            self.TEMC3complete=T
            self.dEMC3complete=d
            


class Plot_all_Profiles_Experiment2(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  
  
  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self,impu='C',case='modded',z_eff_downstream_old = 'none', z_eff_upstream_old = 'none', TeEMC3old='none', TiEMC3old='none',m_eff_downstream_old = 'none', m_eff_upstream_old = 'none', Told='none', dold='none',state='none'):
            plt.rcParams.update({'font.size': 24})
            g=EMC3.Grid() 
            z=g.zones[0]
            radial=z.R
            poloidal=z.Z
            toroidal=z.phi
            a_eff_emc3,a_eff_exp,radius,r_eff_drews=get_a_r_eff2(g)
            
            
            self.a_eff_emc3=a_eff_emc3
            self.a_eff_exp=a_eff_exp
            self.radius=radius
            self.r_eff_drews=r_eff_drews
            g=EMC3.Grid() 

            pathH = '.'

            Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
            GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
            PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
            PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
            niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
            TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


            # Multi-purpose manipulator tip plunge minimum and maximum location
            p1Cart = [ -556.276339, -210.053078, -17.736684 ]
            p2Cart = [ -587.407899, -222.425466, -17.736684 ]
            pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                      np.linspace( p1Cart[1], p2Cart[1], 40 ),
                      np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )

            RLimTip = 0.5 * ( GH.zones[0].R[iTor,iPol,iRad][PlatesH.ids[0][0,0,iRad]].min() +
                              GH.zones[0].R[iTor,iPol,iRad][PlatesH.ids[0][0,0,iRad]==False].max() )
            n_e_imp_C=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            n_e_imp_O=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            if impu=='C':
                for m in range(1,len(niH.keys())):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                n_e_imp=n_e_imp_C       
                print 'Carbon'
            counter=1    
            if impu!='C':
                for m in range(1,8):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                for m in range(8,len(niH.keys())+1):
                    counter+=1
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_O+=n_imp*counter
                n_e_imp=n_e_imp_C+n_e_imp_O
                print 'Carbon plus Oxygen'
                self.n_e_imp_C=n_e_imp_C
                self.n_e_imp_O=n_e_imp_O
            for Te, Ti, n, TiLabel in zip( [ TeTiH('Te')[0] ],
                                                  [ TeTiH('Ti')[0] ],
                                                  [ niH('n1')[0] ],
                                                  [ r'$T_\mathrm{H^{+}}$']):
              fig, ax1 = plt.subplots(  )
#a radial shift is implemented here, because of the radial offset of the
#langmuir probe head              
              shift=0
              
              data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/base/profiles.mat')

              Te_drews=data_drews['te']
              Te22=Te_drews[:,0]
              Te23=Te_drews[:,1]
              Te24=Te_drews[:,2]
              ne_drews=data_drews['ne']
              ne22=ne_drews[:,0]*1e18
              ne23=ne_drews[:,1]*1e18
              ne24=ne_drews[:,2]*1e18
              x=data_drews['x']
              x22=x[:,0]/10
              x23=x[:,1]/10
              x24=x[:,2]/10
              y=data_drews['y']
              y22=y[:,0]/10
              y23=y[:,1]/10
              y24=y[:,2]/10
              radius22=np.sqrt(x22**2+y22**2)
              radius23=np.sqrt(x23**2+y23**2)
              radius24=np.sqrt(x24**2+y24**2)
              ne_up_err=data_drews['ne_error_1']
              ne_up_err22=ne_up_err[:,0]
              ne_up_err23=ne_up_err[:,1]
              ne_up_err24=ne_up_err[:,2]
              ne_down_err=data_drews['ne_error_2']
              ne_down_err22=ne_down_err[:,0]
              ne_down_err23=ne_down_err[:,1]              
              ne_down_err24=ne_down_err[:,2]              
              error_ne_22=np.array([ne_down_err22,ne_up_err22])
              error_ne_23=np.array([ne_down_err23,ne_up_err23])
              error_ne_24=np.array([ne_down_err24,ne_up_err24])
              te_err=data_drews['te_error']
              te_err22=te_err[:,0]
              te_err23=te_err[:,1]              
              te_err24=te_err[:,2]             
              
              self.radius22=radius22
              self.radius23=radius23
              self.radius24=radius24
              
              U22 = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/Upin/160308022voltages.mat')
              U23 = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/Upin/160308023voltages.mat')
              U24 = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/Upin/160308024voltages.mat')
              
              mi=1.66*1e-27*1.00794
              e=1.602e-19
              Aeff=2e-3*3e-3*2

              self.U22=U22
              self.U23=U23
              self.U24=U24
              
              self.ne_22_unmodded=ne22*1e-18
              self.ne_23_unmodded=ne23*1e-18
              self.ne_24_unmodded=ne24*1e-18

              T=EMC3.PlasmaField('TE_TI')
              d=EMC3.PlasmaField('DENSITY')

              if len(d)>1:
                ne_sum=d_sum(g,d,0,0)

                print 'downstream'
                ne_z_sum_downstream = d_z_sum_impu(g,d,1,0,'C')
                ne_sum_downstream=d_sum_impu(g,d,1,0,'C')
                m_eff_downstream=sum_eff_mass(g,d,1,0,'C')
                print 'upstream'
                ne_z_sum_upstream   = d_z_sum_impu(g,d,28,247,'C')
                ne_sum_upstream=d_sum_impu(g,d,28,247,'C')                
                m_eff_upstream=sum_eff_mass(g,d,28,247,'C')
                cs_oben=cs3oben(g,d,T('Te'),T('Ti'),28,247,'C')
                cs_unten=cs3unten(g,d,T('Te'),28,247,'C')
                
                if Told=='none':
                    cs_oben_old=cs3oben(g,d,T('Te'),T('Ti'),30,237,'C')
                    cs_unten_old=cs3unten(g,d,T('Te'),30,237,'C') 
                else:
                    cs_oben_old=cs3oben(g,dold,Told('Te'),Told('Ti'),30,237,'C')
                    cs_unten_old=cs3unten(g,dold,Told('Te'),30,237,'C')
                self.cs_oben=cs_oben
                self.cs_oben_old=cs_oben_old
                self.cs_unten=cs_unten
                self.cs_unten_old=cs_unten_old      
                
                
                if len(d)>7:
                    print 'downstream'
                    ne_z_sum_downstream_C = ne_z_sum_downstream
                    m_eff_downstream_C=meff_downstream
                    ne_z_sum_downstream_O = d_z_sum_impu(g,d,1,0,'O')
                    m_eff_downstream_O=sum_eff_mass(g,d,1,0,'O')
                    ne_sum_downstream_C = ne_sum_downstream
                    ne_sum_downstream_O = d_sum_impu(g,d,1,0,'O')
                    z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                    z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                    print 'upstream'
                    ne_z_sum_upstream_C   = ne_z_sum_upstream
                    m_eff_upstream_C=meff_upstream
                    ne_z_sum_upstream_O   = d_z_sum_impu(g,d,3,247,'O')
                    m_eff_upstream_O=sum_eff_mass(g,d,3,247,'O')
                    ne_sum_upstream_C   = ne_sum_upstream
                    ne_sum_upstream_O   = d_sum_impu(g,d,3,247,'O')
                    z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                    z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                    
                    ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                    self.dummy1=ne_z_sum_downstream
                    ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                    self.dummy2=ne_sum_downstream
                    ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                    ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                    
                z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
                z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream        

                self.z_eff_downstream=z_eff_downstream
                self.z_eff_upstream=z_eff_upstream
                self.m_eff_upstream=m_eff_upstream
                self.m_eff_downstream=m_eff_downstream
                self.Rs=Rs
                self.iTor=iTor
                self.iPol=iPol
                self.iRad=iRad


                overlap=[]
                for i in range(0,radius22.shape[0]):
                    overlap.append(find_nearest(Rs,radius22[i])[1])
            
                self.Te=Te
                if m_eff_upstream_old=='none':
                    m_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_upstream_old[:]=1.66*1e-27
                    m_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_downstream_old[:]=1.66*1e-27        
                if z_eff_upstream_old=='none':
                    z_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    z_eff_upstream_old[:]=1.0
                    z_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    z_eff_downstream_old[:]=1.0
                if TeEMC3old=='none':
                    TeEMC3old=np.zeros(len(Te[iTor,iPol,iRad]))
                    TeEMC3old[:]=Te[iTor,iPol,iRad]
                if TiEMC3old=='none':
                    TiEMC3old=np.zeros(len(Ti[iTor,iPol,iRad]))
                    TiEMC3old[:]=Ti[iTor,iPol,iRad]      
                TeEMC3old=np.array(TeEMC3old)
                TiEMC3old=np.array(TiEMC3old)
                self.overlap=overlap
                self.TeEMC3old=TeEMC3old
                self.TiEMC3old=TiEMC3old
                self.z_eff_upstream_old=z_eff_upstream_old
                self.Te=Te[iTor,iPol,iRad]
                self.Ti=Ti[iTor,iPol,iRad]
                self.z_eff_upstream=z_eff_upstream
                eps_oben=(TiEMC3old[overlap]+z_eff_upstream_old[overlap]*TeEMC3old[overlap])   
                eps_unten=(Ti[iTor,iPol,iRad][overlap]+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap])   
                eps=eps_oben/eps_unten

                
                if m_eff_upstream_old=='none':
                    m_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_upstream_old[:]=1.66*1e-27
                    m_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_downstream_old[:]=1.66*1e-27      
                m_eff_upstream_old=np.array(m_eff_upstream_old)
                eps_oben2=(TiEMC3old[overlap]+z_eff_upstream_old[overlap]*TeEMC3old[overlap])   *m_eff_upstream[overlap]   
                eps_unten2=(Ti[iTor,iPol,iRad][overlap]+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap])  *m_eff_upstream_old[overlap]
                eps2=np.sqrt(eps_oben2/eps_unten2)

                self.TiEMC3old=TiEMC3old[overlap]
                self.TeEMC3old=TeEMC3old[overlap]
                self.z_effEMC3old=z_eff_upstream_old[overlap]
                self.TeEMC3new=Te[iTor,iPol,iRad][overlap]
                self.TiEMC3new=Ti[iTor,iPol,iRad][overlap]
                self.z_effEMC3new=z_eff_upstream[overlap]

                cs2_old=(z_eff_upstream_old[overlap]*TeEMC3old[overlap]*1.602e-19+TiEMC3old[overlap]*1.602e-19)/m_eff_upstream_old[overlap]
                cs2_old=np.sqrt(cs2_old)
                cs2_modded=(z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap]*1.602e-19+Ti[iTor,iPol,iRad][overlap]*1.602e-19)/m_eff_upstream[overlap]
                cs2_modded=np.sqrt(cs2_modded)
                eps2_with_cs=cs2_old/cs2_modded

                TeEMC3old_long=Told('Te')[0][28,247,:]
                TiEMC3old_long=Told('Ti')[0][28,247,:]


                cs2_old_long=(z_eff_upstream_old*TeEMC3old_long*1.602e-19+TiEMC3old_long*1.602e-19)/m_eff_upstream_old
                cs2_old_long=np.sqrt(cs2_old_long)
                cs2_modded_long=(z_eff_upstream_old*TeEMC3old_long*1.602e-19+TiEMC3old_long*1.602e-19)/m_eff_upstream_old
                cs2_modded_long=np.sqrt(cs2_modded_long)

                cs2_old_long_down=(z_eff_downstream_old*TeEMC3old_long*1.602e-19*2+TiEMC3old_long*2*1.602e-19)/m_eff_downstream_old
                cs2_old_long_down=np.sqrt(cs2_old_long_down)
                cs2_modded_long_down=(z_eff_downstream_old*TeEMC3old_long*2*1.602e-19+TiEMC3old_long*2*1.602e-19)/m_eff_downstream_old
                cs2_modded_long_down=np.sqrt(cs2_modded_long_down)                


                cs2_old_down=(z_eff_downstream_old[overlap]*TeEMC3old[overlap]*1.602e-19+TiEMC3old[overlap]*1.602e-19)/m_eff_downstream_old[overlap]
                cs2_old_down=np.sqrt(cs2_old_down)
                cs2_modded_down=(z_eff_downstream[overlap]*Te[iTor,iPol,iRad][overlap]*1.602e-19+Ti[iTor,iPol,iRad][overlap]*1.602e-19)/m_eff_downstream[overlap]
                cs2_modded_down=np.sqrt(cs2_modded_down)
                eps2_with_cs_down=cs2_old_down/cs2_modded_down



                
                if state=='clean':
                    cs_oben_old=cs3oben_clean(g,d,T('Te'),T('Ti'),28,247,'C')
                    cs_unten_old=cs3unten_clean(g,d,T('Te'),28,247,'C')
                    print 'piep'
                
                
                eps3_oben=cs_oben_old*cs_unten
                eps3_unten=cs_oben*cs_unten_old
                eps3=np.sqrt(eps3_oben/eps3_unten)
                self.cs_oben_old=cs_oben_old
                self.cs_unten=cs_unten
                self.cs_oben=cs_oben
                self.cs_unten_old=cs_unten_old
                
                
                self.eps3=eps3
                self.eps2=eps2
                self.eps=eps
                cs3_modded=np.sqrt(cs_oben/cs_unten)
                cs3_old=np.sqrt(cs_oben_old/cs_unten_old)
                self.cs3_modded=cs3_modded
                self.cs3_old=cs3_old
                eps3_with_cs=cs3_old/cs3_modded
                
                print 'this is eps ** min then max **'
                print min(eps)
                print max(eps)
                print 'this is eps2 ** min then max **'
                print min(eps2)
                print max(eps2)
                print 'this is eps2.1 ** min then max **'
                print min(eps2_with_cs)                
                print max(eps2_with_cs)                
                print 'this is eps3 ** min then max **'
                print np.nanmin(eps3)
                print np.nanmax(eps3)
                print 'this is eps2.1 ** min then max **'
                print min(eps2_with_cs)                
                print max(eps2_with_cs)                                
                self.overlap=overlap
                self.m_eff_upstream=m_eff_upstream[overlap]
                self.z_eff_upstream=z_eff_upstream[overlap]
                self.Rs=Rs

              
              
              if case=='mod_sep':
                ne_22_modded=ne22*1e-18*min(eps2)
                ne_23_modded=ne23*1e-18*min(eps2)
                ne_24_modded=ne24*1e-18*min(eps2)

                f, ax1 = plt.subplots()
                ax1.plot( radius22-602, Te22,
                            '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                ax1.errorbar( radius22[::8]-602, Te22[::8],yerr=te_err22[::8], color='orange')
                ax1.plot( radius23-602, Te23,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius23[::8]-602, Te23[::8],yerr=te_err23[::8], color='orange')
                ax1.plot( radius24-602, Te24,
                            '-', linewidth = 3,color='orange' )
                ax1.errorbar( radius24[::8]-602, Te24[::8],yerr=te_err24[::8], color='orange')
                
                ax1.plot( Rs - RLimTip, zero_to_nan(Te[iTor,iPol,iRad]),
                            'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )              
                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ [cm]' )
                ax1.set_ylabel( r'$T$ [eV]' )
                ax1.set_xlim( 3, 8 )
                ax1.set_ylim( 0, 20 )

                ax1.grid()
                lines1, labels1 = ax1.get_legend_handles_labels()

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )

        
        
                f, ax1 = plt.subplots()

                ax1.plot( radius22-602, ne_22_modded,
                            '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                ax1.errorbar( radius22[::8]-602, ne_22_modded[::8],yerr=error_ne_22[:,::8], color='darkviolet')
                ax1.plot( radius23-602, ne_23_modded,
                            '-', linewidth = 3,color='darkviolet' )
                ax1.errorbar( radius23[::8]-602, ne_23_modded[::8],yerr=error_ne_23[:,::8], color='darkviolet')              
                ax1.plot( radius24-602, ne_24_modded,
                            '-', linewidth = 3,color='darkviolet' )
                ax1.errorbar( radius24[::8]-602, ne_24_modded[::8],yerr=error_ne_24[:,::8], color='darkviolet')              
                ax1.plot( Rs - RLimTip, (zero_to_nan(n[iTor,iPol,iRad])+n_e_imp) * 1e-12,
                            'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )

                ax1.set_xlabel( r'$R - R_\mathrm{LCFS}$ [cm]' )
                ax1.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}/\mathrm{m}^3$]' )
                ax1.set_xlim( 3, 8 )
                ax1.set_ylim( 0, 1.5 )

                ax1.grid()
                lines1, labels1 = ax1.get_legend_handles_labels()

                for label in ax1.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

                plt.tight_layout( pad = 0.25 )

                self.ne_22_modded=ne_22_modded
                self.ne_23_modded=ne_23_modded
                self.ne_24_modded=ne_24_modded
                self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                self.TiEMC3=zero_to_nan(Ti[iTor,iPol,iRad])

                f, ax = plt.subplots()
                ax.plot( radial[28,247,:-1], cs2_modded_long_down/1e6,
                            '-', linewidth = 3, label = r'$c_\mathrm{s}^2$',color='#7e2f8e' )
                ax.plot( radial[28,247,:-1], cs2_modded_long/1e6,
                            '-', linewidth = 3, label = r'$c_\mathrm{s}^2$',color='#7e2f8e' )
                #ax.plot( radial[28,247,:-1], cs3_modded/1e6,
                            #'--',color='#d91219', linewidth = 3, label = r'$c_\mathrm{s}^3$' )
                ax.set_xlabel( r'$R$ [cm]' )
                ax.set_ylabel( r'$c_\mathrm{s}$ [$\mathrm{10^6m/s}$]' )
                ax.set_xlim(600,615)
                #ax.set_ylim(0,6)
                ax.grid()

                plt.tight_layout()
                for label in ax.xaxis.get_ticklabels()[::2]:
                    label.set_visible(False)
                for label in ax.yaxis.get_ticklabels()[::2]:
                    label.set_visible(False)

            self.radius22=radius22
            self.radius23=radius23
            self.radius22=radius24
            self.ne_EMC3=(n[iTor,iPol,iRad]+n_e_imp)* 1e-12





            
            
class Plot_all_Profiles_Divertor(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  
  
  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self,impu='C',case='modded',z_eff_downstream = 'none', z_eff_upstream = 'none', minus_shift=0):
            plt.rcParams.update({'font.size': 18})
            g=EMC3.Grid() 
            
            a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

            pathH = '.'

            Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
            GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
            PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
            PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
            niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
            TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


            # Multi-purpose manipulator tip plunge minimum and maximum location
            p1Cart = [ -560.482757, -212.503934, -17.065822 ]
            p2Cart = [ -593.011862, -225.421244, -17.035162 ]
            pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                      np.linspace( p1Cart[1], p2Cart[1], 40 ),
                      np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
            Zs=Zs-minus_shift
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )

            n_e_imp_C=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            n_e_imp_O=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            if impu=='C':
                for m in range(1,len(niH.keys())):
                    dummy_false='n'+ str(m) in niH.keys()
                    if dummy_false==False:
                        n_e_imp_C+=n_imp*0                        
                    if dummy_false==True:    
                        n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                        n_e_imp_C+=n_imp*m
                n_e_imp=n_e_imp_C       
                print 'Carbon'
            counter=1    
            if impu!='C':
                for m in range(1,8):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                for m in range(8,len(niH.keys())+1):
                    counter+=1
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_O+=n_imp*counter
                n_e_imp=n_e_imp_C+n_e_imp_O
                print 'Carbon plus Oxygen'
                self.n_e_imp_C=n_e_imp_C
                self.n_e_imp_O=n_e_imp_O
            for Te, Ti, n, TiLabel in zip( [ TeTiH('Te')[0] ],
                                                  [ TeTiH('Ti')[0] ],
                                                  [ niH('n1')[0] ],
                                                  [ r'$T_\mathrm{H^{+}}$']):
#a radial shift is implemented here, because of the radial offset of the
#langmuir probe head              
              shift=0
              
              data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_divertor/exp_results/drews/21_26/profiles_analysis_joerg171026025.mat')

              Te_drews=data_drews['Te']
              Te_drews=Te_drews[:,0]
              ne_drews=data_drews['Ne']
              ne=ne_drews[:,0]*1e18
              x=data_drews['x']
              x=x[:,0]/10
              y=data_drews['y']
              y=y[:,0]/10
              radius_drews=np.sqrt(x**2+y**2)*1000
         
              
              self.radius_drews=radius_drews
              
              mi=1.66*1e-27*1.00794
              e=1.602e-19
              Aeff=2e-3*3e-3*2

              T=EMC3.PlasmaField('TE_TI')
              d=EMC3.PlasmaField('DENSITY')
              
              if len(d)==1:
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                'r--', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    #ax.set_xlim(608-6+1.3,616)
                    ax.set_ylim(0,1000)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    #ax.set_xlim(608-6+1.3,613-6+1.3)
                    #ax.set_ylim(0,1)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)            

              if len(d)>1:
                ne_sum=d_sum(g,d,0,0)

                print 'downstream'
                ne_z_sum_downstream = d_z_sum_impu(g,d,17,102,'C')
                ne_sum_downstream=d_sum_impu(g,d,17,102,'C')
                m_eff_downstream=sum_eff_mass(g,d,17,102,'C')
                print 'upstream'
                ne_z_sum_upstream   = d_z_sum_impu(g,d,30,237,'C')
                ne_sum_upstream=d_sum_impu(g,d,30,237,'C')                
                m_eff_upstream=sum_eff_mass(g,d,30,237,'C')
                if len(d)>7:
                    print 'downstream'
                    ne_z_sum_downstream_C = ne_z_sum_downstream
                    m_eff_downstream_C=m_eff_downstream
                    ne_z_sum_downstream_O = d_z_sum_impu(g,d,17,102,'O')
                    m_eff_downstream_O=sum_eff_mass(g,d,17,102,'O')
                    ne_sum_downstream_C = ne_sum_downstream
                    ne_sum_downstream_O = d_sum_impu(g,d,17,102,'O')
                    z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                    z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                    print 'upstream'
                    ne_z_sum_upstream_C   = ne_z_sum_upstream
                    m_eff_upstream_C=m_eff_upstream
                    ne_z_sum_upstream_O   = d_z_sum_impu(g,d,30,237,'O')
                    m_eff_upstream_O=sum_eff_mass(g,d,30,237,'O')
                    ne_sum_upstream_C   = ne_sum_upstream
                    ne_sum_upstream_O   = d_sum_impu(g,d,30,237,'O')
                    z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                    z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                    
                    ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                    self.dummy1=ne_z_sum_downstream
                    ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                    self.dummy2=ne_sum_downstream
                    ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                    ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                    
                z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
                z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream        

                self.z_eff_downstream=z_eff_downstream
                self.z_eff_upstream=z_eff_upstream
                self.m_eff_upstream=m_eff_upstream
                self.m_eff_downstream=m_eff_downstream
                self.Rs=Rs
                self.iTor=iTor
                self.iPol=iPol
                self.iRad=iRad

                overlap=[]
                for i in range(0,radius_drews.shape[0]):
                    overlap.append(find_nearest(Rs,radius_drews[i])[1])
                
                self.overlap=overlap
                self.Te=Te

                if case=='unmodded':
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, zero_to_nan((n[iTor,iPol,iRad]+n_e_imp) * 1e-12),
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    #ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.radius_drews=radius_drews
                    self.Te_drews=Te_drews
                    self.ne_drews=ne_drews
                    self.Rs=Rs
                    self.neEMC3=(n[iTor,iPol,iRad]+n_e_imp)*1e-12
                    self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                    
                    
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                    self.TiEMC3=zero_to_nan(Ti[iTor,iPol,iRad])
                    self.TeEMC3complete=Te
                    self.TiEMC3complete=Ti
                    self.TEMC3complete=T
                    self.dEMC3complete=d

                if case=='modded':
                    ne_drews_modded=ne_drews[:,0]*z_eff_upstream[overlap] 
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews_modded,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                    self.TiEMC3=zero_to_nan(Ti[iTor,iPol,iRad])
                    self.ne_drews_modded=ne_drews_modded
                    self.TeEMC3complete=Te
                    self.TiEMC3complete=Ti
                    self.TEMC3complete=T
                    self.dEMC3complete=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)






class Plot_all_Profiles_Divertor_H(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  
  
  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self,impu='C',case='unmodded',z_eff_downstream = 'none', z_eff_upstream = 'none', minus_shift=0):
            plt.rcParams.update({'font.size': 18})
            g=EMC3.Grid() 
            
            a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

            pathH = '.'

            Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
            GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
            PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
            PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
            niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
            TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


            # Multi-purpose manipulator tip plunge minimum and maximum location
            p1Cart = [ -560.482757, -212.503934, -17.065822 ]
            p2Cart = [ -593.011862, -225.421244, -17.035162 ]
            pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                      np.linspace( p1Cart[1], p2Cart[1], 40 ),
                      np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
            Zs=Zs-minus_shift
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )

            n_e_imp_C=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            n_e_imp_O=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            if impu=='C':
                for m in range(1,len(niH.keys())):
                    dummy_false='n'+ str(m) in niH.keys()
                    if dummy_false==False:
                        n_e_imp_C+=n_imp*0                        
                    if dummy_false==True:    
                        n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                        n_e_imp_C+=n_imp*m
                n_e_imp=n_e_imp_C       
                print 'Carbon'
            counter=1    
            if impu!='C':
                for m in range(1,8):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                for m in range(8,len(niH.keys())+1):
                    counter+=1
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_O+=n_imp*counter
                n_e_imp=n_e_imp_C+n_e_imp_O
                print 'Carbon plus Oxygen'
                self.n_e_imp_C=n_e_imp_C
                self.n_e_imp_O=n_e_imp_O
            for Te, Ti, n, TiLabel in zip( [ TeTiH('Te')[0] ],
                                                  [ TeTiH('Ti')[0] ],
                                                  [ niH('n1')[0] ],
                                                  [ r'$T_\mathrm{H^{+}}$']):
#a radial shift is implemented here, because of the radial offset of the
#langmuir probe head              
              shift=0
              
              data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_divertor/exp_results/drews/2018/mpm_profiles_joerg180814036.mat')

              Te_drews=data_drews['Te_mpm2']
              Te_drews=Te_drews[:,0]
              ne_drews=data_drews['ne_mpm2']
              ne=ne_drews[:,0]*1e18
              x=data_drews['x_mpm2']
              x=x[:,0]/10
              y=data_drews['y_mpm2']
              y=y[:,0]/10
              radius_drews=np.sqrt(x**2+y**2)*1000
         
              
              self.radius_drews=radius_drews
              
              mi=1.66*1e-27*1.00794
              e=1.602e-19
              Aeff=2e-3*3e-3*2

              T=EMC3.PlasmaField('TE_TI')
              d=EMC3.PlasmaField('DENSITY')
              
              if len(d)==1:
                    f, ax = plt.subplots()
                    ax.plot( radius_drews[0:-len(radius_drews)/2], Te_drews[0:-len(Te_drews)/2],
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                'r--', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    #ax.set_xlim(608-6+1.3,616)
                    ax.set_ylim(0,1000)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews[0:-len(radius_drews)/2], ne_drews[0:-len(ne_drews)/2],
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    #ax.set_xlim(608-6+1.3,613-6+1.3)
                    #ax.set_ylim(0,1)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)            

              if len(d)>1:
                ne_sum=d_sum(g,d,0,0)

                print 'downstream'
                ne_z_sum_downstream = d_z_sum_impu(g,d,17,102,'C')
                ne_sum_downstream=d_sum_impu(g,d,17,102,'C')
                m_eff_downstream=sum_eff_mass(g,d,17,102,'C')
                print 'upstream'
                ne_z_sum_upstream   = d_z_sum_impu(g,d,30,237,'C')
                ne_sum_upstream=d_sum_impu(g,d,30,237,'C')                
                m_eff_upstream=sum_eff_mass(g,d,30,237,'C')
                if len(d)>7:
                    print 'downstream'
                    ne_z_sum_downstream_C = ne_z_sum_downstream
                    m_eff_downstream_C=m_eff_downstream
                    ne_z_sum_downstream_O = d_z_sum_impu(g,d,17,102,'O')
                    m_eff_downstream_O=sum_eff_mass(g,d,17,102,'O')
                    ne_sum_downstream_C = ne_sum_downstream
                    ne_sum_downstream_O = d_sum_impu(g,d,17,102,'O')
                    z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                    z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                    print 'upstream'
                    ne_z_sum_upstream_C   = ne_z_sum_upstream
                    m_eff_upstream_C=m_eff_upstream
                    ne_z_sum_upstream_O   = d_z_sum_impu(g,d,30,237,'O')
                    m_eff_upstream_O=sum_eff_mass(g,d,30,237,'O')
                    ne_sum_upstream_C   = ne_sum_upstream
                    ne_sum_upstream_O   = d_sum_impu(g,d,30,237,'O')
                    z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                    z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                    
                    ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                    self.dummy1=ne_z_sum_downstream
                    ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                    self.dummy2=ne_sum_downstream
                    ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                    ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                    
                z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
                z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream        

                self.z_eff_downstream=z_eff_downstream
                self.z_eff_upstream=z_eff_upstream
                self.m_eff_upstream=m_eff_upstream
                self.m_eff_downstream=m_eff_downstream
                self.Rs=Rs
                self.iTor=iTor
                self.iPol=iPol
                self.iRad=iRad


                overlap=[]
                for i in range(0,radius_drews.shape[0]):
                    overlap.append(find_nearest(Rs,radius_drews[i])[1])
                
                self.overlap=overlap
                self.Te=Te

                if case=='unmodded':
                    f, ax = plt.subplots()
                    ax.plot( radius_drews[0:-len(radius_drews)/2-30], Te_drews[0:-len(Te_drews)/2-30],
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(604,612)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews[0:-len(radius_drews)/2-30], ne_drews[0:-len(Te_drews)/2-30],
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, zero_to_nan((n[iTor,iPol,iRad]+n_e_imp) * 1e-12),
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(604,612)
                    ax.set_ylim(0,8)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    self.Rs=Rs
                    self.ne=zero_to_nan((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.Te=zero_to_nan(Te[iTor,iPol,iRad])
                    self.ne_drews_unmodded=ne_drews
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                    self.TiEMC3=zero_to_nan(Ti[iTor,iPol,iRad])
                    #self.ne_drews_modded=ne_drews_modded
                    self.TeEMC3complete=Te
                    self.TiEMC3complete=Ti
                    self.TEMC3complete=T
                    self.dEMC3complete=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.dEMC3complete=d





class Plot_all_Profiles_Divertor2(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  
  
  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self,impu='C',case='modded',z_eff_downstream_old = 'none', z_eff_upstream_old = 'none', TeEMC3old='none', TiEMC3old='none',m_eff_downstream_old = 'none', m_eff_upstream_old = 'none', Told='none', dold='none',state='none'):
            plt.rcParams.update({'font.size': 18})
            g=EMC3.Grid() 
            z=g.zones[0]
            radial=z.R
            poloidal=z.Z
            toroidal=z.phi
            
            a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

            pathH = '.'

            Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
            GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
            PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
            PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
            niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
            TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )

            # Multi-purpose manipulator tip plunge minimum and maximum location
            p1Cart = [ -560.482757, -212.503934, -17.065822 ]
            p2Cart = [ -593.011862, -225.421244, -17.035162 ]
            pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                      np.linspace( p1Cart[1], p2Cart[1], 40 ),
                      np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )

            n_e_imp_C=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            n_e_imp_O=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            if impu=='C':
                for m in range(1,len(niH.keys())):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                n_e_imp=n_e_imp_C       
                print 'Carbon'
            if impu!='C':
                for m in range(1,7):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                for m in range(7,len(niH.keys())+1):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_O+=n_imp*m
                    print m
                n_e_imp=n_e_imp_C+n_e_imp_O
                print 'Carbon plus Oxygen'
                self.n_e_imp_C=n_e_imp_C
                self.n_e_imp_O=n_e_imp_O
            for Te, Ti, n, TiLabel in zip( [ TeTiH('Te')[0] ],
                                                  [ TeTiH('Ti')[0] ],
                                                  [ niH('n1')[0] ],
                                                  [ r'$T_\mathrm{H^{+}}$']):
#a radial shift is implemented here, because of the radial offset of the
#langmuir probe head              
              shift=0
              
              data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_divertor/exp_results/drews/profiles_analysis_mod171026038.mat')

              Te_drews=data_drews['Te']
              Te_drews=Te_drews[:,0]
              self.Te_drews=Te_drews
              ne_drews=data_drews['Ne']
              ne=ne_drews[:,0]*1e18
              x=data_drews['x']
              x=x[:,0]/10
              y=data_drews['y']
              y=y[:,0]/10
              radius_drews=np.sqrt(x**2+y**2)*1000
              self.radius_drews=radius_drews
              
              
              mi=1.66*1e-27*1.00794
              e=1.602e-19
              Aeff=2e-3*3e-3*2

              T=EMC3.PlasmaField('TE_TI')
              d=EMC3.PlasmaField('DENSITY')
              

              if len(d)==1:
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                'r--', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    #ax.set_xlim(608-6+1.3,616)
                    ax.set_ylim(0,1000)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    #ax.set_xlim(608-6+1.3,613-6+1.3)
                    #ax.set_ylim(0,1)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)            

              if len(d)>1:
                ne_sum=d_sum(g,d,0,0)

                print 'downstream'
                ne_z_sum_downstream = d_z_sum_impu(g,d,17,102,'C')
                ne_sum_downstream=d_sum_impu(g,d,17,102,'C')
                m_eff_downstream=sum_eff_mass(g,d,17,102,'C')
                print 'upstream'
                ne_z_sum_upstream=  d_z_sum_impu(g,d,30,237,'C')
                ne_sum_upstream=    d_sum_impu(g,d,30,237,'C')                
                m_eff_upstream=sum_eff_mass(g,d,30,237,'C')
                cs_oben=cs3oben(g,d,T('Te'),T('Ti'),30,237,'C')
                cs_unten=cs3unten(g,d,T('Te'),30,237,'C')
                print np.shape(T)
                print np.shape(Told)
                if Told=='none':
                    cs_oben_old=cs3oben(g,d,T('Te'),T('Ti'),30,237,'C')
                    cs_unten_old=cs3unten(g,d,T('Te'),30,237,'C') 
                else:
                    cs_oben_old=cs3oben(g,dold,Told('Te'),Told('Ti'),30,237,'C')
                    cs_unten_old=cs3unten(g,dold,Told('Te'),30,237,'C')
                self.cs_oben=cs_oben
                self.cs_oben_old=cs_oben_old
                self.cs_unten=cs_unten
                self.cs_unten_old=cs_unten_old                
                if len(d)>7:
                    print 'downstream'
                    ne_z_sum_downstream_C = ne_z_sum_downstream
                    m_eff_downstream_C=meff_downstream
                    ne_z_sum_downstream_O = d_z_sum_impu(g,d,17,102,'O')
                    m_eff_downstream_O=sum_eff_mass(g,d,17,102,'O')
                    ne_sum_downstream_C = ne_sum_downstream
                    ne_sum_downstream_O = d_sum_impu(g,d,17,102,'O')
                    z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                    z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                    print 'upstream'
                    ne_z_sum_upstream_C   = ne_z_sum_upstream
                    m_eff_upstream_C=meff_upstream
                    ne_z_sum_upstream_O   = d_z_sum_impu(g,d,30,237,'O')
                    m_eff_upstream_O=sum_eff_mass(g,d,30,237,'O')
                    ne_sum_upstream_C   = ne_sum_upstream
                    ne_sum_upstream_O   = d_sum_impu(g,d,30,237,'O')
                    z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                    z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                    
                    ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                    self.dummy1=ne_z_sum_downstream
                    ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                    self.dummy2=ne_sum_downstream
                    ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                    ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                        
                self.ne_z_sum_upstream=ne_z_sum_upstream
                self.ne_sum_upstream=ne_sum_upstream
                z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
                z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream        
                overlap=[]
                for i in range(0,radius_drews.shape[0]):
                    overlap.append(find_nearest(Rs,radius_drews[i])[1])

                if m_eff_upstream_old=='none':
                    m_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_upstream_old[:]=1.66*1e-27
                    m_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_downstream_old[:]=1.66*1e-27        
                if z_eff_upstream_old=='none':
                    z_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    z_eff_upstream_old[:]=1.0
                    z_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    z_eff_downstream_old[:]=1.0
                if TeEMC3old=='none':
                    TeEMC3old=Te[iTor,iPol,iRad]
                if TiEMC3old=='none':
                    TiEMC3old=Ti[iTor,iPol,iRad]                 
                    
                TiEMC3old=np.array(TiEMC3old)
                TeEMC3old=np.array(TeEMC3old)

                self.z_eff_upstream=z_eff_upstream

                eps_oben=(TiEMC3old[overlap]+z_eff_upstream_old[overlap]*TeEMC3old[overlap])   
                eps_unten=(Ti[iTor,iPol,iRad][overlap]+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap])   
                eps=eps_oben/eps_unten
                
                m_eff_upstream_old=np.array(m_eff_upstream_old)
                
                if m_eff_upstream_old=='none':
                    m_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_upstream_old[:]=1.66*1e-27
                    m_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_downstream_old[:]=1.66*1e-27      
                m_eff_upstream_old=np.array(m_eff_upstream_old)
                eps_oben2=(TiEMC3old[overlap]*1.602e-19+z_eff_upstream_old[overlap]*TeEMC3old[overlap]*1.602e-19)*m_eff_upstream[overlap]   
                eps_unten2=(Ti[iTor,iPol,iRad][overlap]*1.602e-19+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap]*1.602e-19)  *m_eff_upstream_old[overlap]
                self.eps_oben2=eps_oben2
                self.eps_unten2=eps_unten2
                self.TeEMC3old=TeEMC3old[overlap]
                self.z_eff_upstream_old=z_eff_upstream_old[overlap]
                self.TiEMC3old=TiEMC3old[overlap]
                self.m_eff_upstream=m_eff_upstream[overlap]   
                self.Te=Te[iTor,iPol,iRad][overlap]
                self.z_eff_upstream=z_eff_upstream[overlap]
                self.Ti=Ti[iTor,iPol,iRad][overlap]
                self.m_eff_upstream_old=m_eff_upstream_old[overlap]
                eps2=np.sqrt(eps_oben2/eps_unten2)
                
                
                
                cs2_old=(TiEMC3old[overlap]*1.602e-19+z_eff_upstream_old[overlap]*TeEMC3old[overlap]*1.602e-19)/m_eff_upstream_old[overlap]
                cs2_old=np.sqrt(cs2_old)
                cs2_modded=(Ti[iTor,iPol,iRad][overlap]*1.602e-19+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap]*1.602e-19)/m_eff_upstream[overlap]
                cs2_modded=np.sqrt(cs2_modded)
                eps2_with_cs=cs2_old/cs2_modded

                TeEMC3old_long=Told('Te')[0][30,237,:]
                TiEMC3old_long=Told('Ti')[0][30,237,:]

                cs2_old_long=(TiEMC3old_long*1.602e-19+z_eff_upstream_old*TeEMC3old_long*1.602e-19)/m_eff_upstream_old
                cs2_old_long=np.sqrt(cs2_old_long)
                cs2_modded_long=(TiEMC3old_long*1.602e-19+z_eff_upstream_old*TeEMC3old_long*1.602e-19)/m_eff_upstream_old
                cs2_modded_long=np.sqrt(cs2_modded_long)
                
                
                self.cs2_modded_long=cs2_modded_long
                self.cs2_old_long=cs2_old_long
                
                if state=='clean':
                    cs_oben_old=cs3oben_clean(g,d,T('Te'),T('Ti'),30,237,'C')
                    cs_unten_old=cs3unten_clean(g,d,T('Te'),30,237,'C')
                    print 'piep'
                
                eps3_oben=cs_oben_old*cs_unten
                eps3_unten=cs_oben*cs_unten_old
                eps3=np.sqrt(eps3_oben/eps3_unten)
                self.cs_oben_old=cs_oben_old
                self.cs_unten=cs_unten
                self.cs_oben=cs_oben
                self.cs_unten_old=cs_unten_old
                
                
                self.eps3=eps3
                self.eps2=eps2
                self.eps=eps
                cs3_modded=np.sqrt(cs_oben/cs_unten)
                cs3_old=np.sqrt(cs_oben_old/cs_unten_old)
                self.cs3_modded=cs3_modded
                self.cs3_old=cs3_old
                eps3_with_cs=cs3_old/cs3_modded
                
                print 'this is eps ** min then max **'
                print min(eps)
                print max(eps)
                print 'this is eps2 ** min then max **'
                print min(eps2)
                print max(eps2)
                print 'this is eps2.1 ** min then max **'
                print min(eps2_with_cs)                
                print max(eps2_with_cs)                
                print 'this is eps3 ** min then mean **'
                print np.nanmin(eps3)
                print np.nanmean(eps3)
                
                self.overlap=overlap
                self.m_eff_upstream=m_eff_upstream[overlap]
                self.z_eff_upstream=z_eff_upstream[overlap]
                self.Rs=Rs
                if case=='unmodded':
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.Te_emc3=Te[iTor,iPol,iRad]
                    self.Ti_emc3=Ti[iTor,iPol,iRad]
                
                if case=='modded':
                    ne_drews_modded=ne_drews[:,0]*min(eps2)
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radial[30,237,:-1], cs2_modded_long/1e6,
                                '-', linewidth = 3, label = r'$c_\mathrm{s}^2$',color='#7e2f8e' )
                    ax.plot( radial[30,237,:-1], cs3_modded/1e6,
                                '--',color='#d91219', linewidth = 3, label = r'$c_\mathrm{s}^3$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$c_\mathrm{s}$ [$\mathrm{10^6m/s}$]' )
                    ax.set_xlim(597,617)
                    #ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_drews_modded=ne_drews_modded
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.Te_emc3=Te[iTor,iPol,iRad]
                    self.Ti_emc3=Ti[iTor,iPol,iRad]


                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews_modded,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_drews_modded=ne_drews_modded
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.Te_emc3=Te[iTor,iPol,iRad]
                    self.Ti_emc3=Ti[iTor,iPol,iRad]




class Plot_all_Profiles_Divertor2_H(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  
  
  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self,impu='C',case='modded',z_eff_downstream_old = 'none', z_eff_upstream_old = 'none', TeEMC3old='none', TiEMC3old='none',m_eff_downstream_old = 'none', m_eff_upstream_old = 'none', Told='none', dold='none',state='none',cs='eff'):
            plt.rcParams.update({'font.size': 18})
            g=EMC3.Grid() 
            z=g.zones[0]
            radial=z.R
            poloidal=z.Z
            toroidal=z.phi
            
            a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

            pathH = '.'

            Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
            GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
            PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
            PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
            niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
            TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )

            # Multi-purpose manipulator tip plunge minimum and maximum location
            p1Cart = [ -560.482757, -212.503934, -17.065822 ]
            p2Cart = [ -593.011862, -225.421244, -17.035162 ]
            pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                      np.linspace( p1Cart[1], p2Cart[1], 40 ),
                      np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )

            n_e_imp_C=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            n_e_imp_O=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            if impu=='C':
                for m in range(1,len(niH.keys())):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                n_e_imp=n_e_imp_C       
                print 'Carbon'
            if impu!='C':
                for m in range(1,7):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                for m in range(7,len(niH.keys())+1):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_O+=n_imp*m
                    print m
                n_e_imp=n_e_imp_C+n_e_imp_O
                print 'Carbon plus Oxygen'
                self.n_e_imp_C=n_e_imp_C
                self.n_e_imp_O=n_e_imp_O
            for Te, Ti, n, TiLabel in zip( [ TeTiH('Te')[0] ],
                                                  [ TeTiH('Ti')[0] ],
                                                  [ niH('n1')[0] ],
                                                  [ r'$T_\mathrm{H^{+}}$']):
#a radial shift is implemented here, because of the radial offset of the
#langmuir probe head              
              shift=0
              
              data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_divertor/exp_results/drews/2018/mpm_profiles_joerg180814036.mat')

              Te_drews=data_drews['Te_mpm2']
              Te_drews=Te_drews[:,0]
              ne_drews=data_drews['ne_mpm2']
              ne=ne_drews[:,0]*1e18
              x=data_drews['x_mpm2']
              x=x[:,0]/10
              y=data_drews['y_mpm2']
              y=y[:,0]/10
              radius_drews=np.sqrt(x**2+y**2)*1000
              self.radius_drews=radius_drews
              
              
              mi=1.66*1e-27*1.00794
              e=1.602e-19
              Aeff=2e-3*3e-3*2

              T=EMC3.PlasmaField('TE_TI')
              d=EMC3.PlasmaField('DENSITY')
              

              if len(d)==1:
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                'r--', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    #ax.set_xlim(608-6+1.3,616)
                    ax.set_ylim(0,1000)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    #ax.set_xlim(608-6+1.3,613-6+1.3)
                    #ax.set_ylim(0,1)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)            

              if len(d)>1:
                ne_sum=d_sum(g,d,0,0)

                print 'downstream'
                ne_z_sum_downstream = d_z_sum_impu(g,d,17,102,'C')
                ne_sum_downstream=d_sum_impu(g,d,17,102,'C')
                m_eff_downstream=sum_eff_mass(g,d,17,102,'C')
                print 'upstream'
                ne_z_sum_upstream=  d_z_sum_impu(g,d,30,237,'C')
                ne_sum_upstream=    d_sum_impu(g,d,30,237,'C')                
                m_eff_upstream=sum_eff_mass(g,d,30,237,'C')
                cs_oben=cs3oben(g,d,T('Te'),T('Ti'),30,237,'C')
                cs_unten=cs3unten(g,d,T('Te'),30,237,'C')
                print np.shape(T)
                print np.shape(Told)
                if Told=='none':
                    cs_oben_old=cs3oben(g,d,T('Te'),T('Ti'),30,237,'C')
                    cs_unten_old=cs3unten(g,d,T('Te'),30,237,'C') 
                else:
                    cs_oben_old=cs3oben(g,dold,Told('Te'),Told('Ti'),30,237,'C')
                    cs_unten_old=cs3unten(g,dold,Told('Te'),30,237,'C')
                self.cs_oben=cs_oben
                self.cs_oben_old=cs_oben_old
                self.cs_unten=cs_unten
                self.cs_unten_old=cs_unten_old                
                if len(d)>7:
                    print 'downstream'
                    ne_z_sum_downstream_C = ne_z_sum_downstream
                    m_eff_downstream_C= m_eff_downstream
                    ne_z_sum_downstream_O = d_z_sum_impu(g,d,17,102,'O')
                    m_eff_downstream_O=sum_eff_mass(g,d,17,102,'O')
                    ne_sum_downstream_C = ne_sum_downstream
                    ne_sum_downstream_O = d_sum_impu(g,d,17,102,'O')
                    z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                    z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                    print 'upstream'
                    ne_z_sum_upstream_C   = ne_z_sum_upstream
                    m_eff_upstream_C=m_eff_upstream
                    ne_z_sum_upstream_O   = d_z_sum_impu(g,d,30,237,'O')
                    m_eff_upstream_O=sum_eff_mass(g,d,30,237,'O')
                    ne_sum_upstream_C   = ne_sum_upstream
                    ne_sum_upstream_O   = d_sum_impu(g,d,30,237,'O')
                    z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                    z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                    
                    ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                    self.dummy1=ne_z_sum_downstream
                    ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                    self.dummy2=ne_sum_downstream
                    ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                    ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                        
                self.ne_z_sum_upstream=ne_z_sum_upstream
                self.ne_sum_upstream=ne_sum_upstream
                z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
                z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream        
                overlap=[]
                for i in range(0,radius_drews.shape[0]):
                    overlap.append(find_nearest(Rs,radius_drews[i])[1])

                if m_eff_upstream_old=='none':
                    m_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_upstream_old[:]=1.66*1e-27
                    m_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_downstream_old[:]=1.66*1e-27        
                if z_eff_upstream_old=='none':
                    z_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    z_eff_upstream_old[:]=1.0
                    z_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    z_eff_downstream_old[:]=1.0
                if TeEMC3old=='none':
                    TeEMC3old=Te[iTor,iPol,iRad]
                if TiEMC3old=='none':
                    TiEMC3old=Ti[iTor,iPol,iRad]                 
                    
                TiEMC3old=np.array(TiEMC3old)
                TeEMC3old=np.array(TeEMC3old)

                self.z_eff_upstream=z_eff_upstream

                eps_oben=(TiEMC3old[overlap]+z_eff_upstream_old[overlap]*TeEMC3old[overlap])   
                eps_unten=(Ti[iTor,iPol,iRad][overlap]+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap])   
                eps=eps_oben/eps_unten
                
                m_eff_upstream_old=np.array(m_eff_upstream_old)
                
                if m_eff_upstream_old=='none':
                    m_eff_upstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_upstream_old[:]=1.66*1e-27
                    m_eff_downstream_old=np.zeros(len(m_eff_upstream))
                    m_eff_downstream_old[:]=1.66*1e-27      
                m_eff_upstream_old=np.array(m_eff_upstream_old)
                eps_oben2=(TiEMC3old[overlap]*1.602e-19+z_eff_upstream_old[overlap]*TeEMC3old[overlap]*1.602e-19)*m_eff_upstream[overlap]   
                eps_unten2=(Ti[iTor,iPol,iRad][overlap]*1.602e-19+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap]*1.602e-19)  *m_eff_upstream_old[overlap]
                self.eps_oben2=eps_oben2
                self.eps_unten2=eps_unten2
                self.TeEMC3old=TeEMC3old[overlap]
                self.z_eff_upstream_old=z_eff_upstream_old[overlap]
                self.TiEMC3old=TiEMC3old[overlap]
                self.m_eff_upstream=m_eff_upstream[overlap]   
                self.Te=Te[iTor,iPol,iRad][overlap]
                self.z_eff_upstream=z_eff_upstream[overlap]
                self.Ti=Ti[iTor,iPol,iRad][overlap]
                self.m_eff_upstream_old=m_eff_upstream_old[overlap]
                eps2=np.sqrt(eps_oben2/eps_unten2)
                
                
                
                cs2_old=(TiEMC3old[overlap]*1.602e-19+z_eff_upstream_old[overlap]*TeEMC3old[overlap]*1.602e-19)/m_eff_upstream_old[overlap]
                cs2_old=np.sqrt(cs2_old)
                cs2_modded=(Ti[iTor,iPol,iRad][overlap]*1.602e-19+z_eff_upstream[overlap]*Te[iTor,iPol,iRad][overlap]*1.602e-19)/m_eff_upstream[overlap]
                cs2_modded=np.sqrt(cs2_modded)
                eps2_with_cs=cs2_old/cs2_modded

                TeEMC3old_long=Told('Te')[0][30,237,:]
                TiEMC3old_long=Told('Ti')[0][30,237,:]

                cs2_old_long=(TiEMC3old_long*1.602e-19+z_eff_upstream_old*TeEMC3old_long*1.602e-19)/m_eff_upstream_old
                cs2_old_long=np.sqrt(cs2_old_long)
                cs2_modded_long=(TiEMC3old_long*1.602e-19+z_eff_upstream_old*TeEMC3old_long*1.602e-19)/m_eff_upstream_old
                cs2_modded_long=np.sqrt(cs2_modded_long)
                
                
                self.cs2_modded_long=cs2_modded_long
                self.cs2_old_long=cs2_old_long
                
                if state=='clean':
                    cs_oben_old=cs3oben_clean(g,d,T('Te'),T('Ti'),30,237,'C')
                    cs_unten_old=cs3unten_clean(g,d,T('Te'),30,237,'C')
                    print 'piep'
                
                eps3_oben=cs_oben_old*cs_unten
                eps3_unten=cs_oben*cs_unten_old
                eps3=np.sqrt(eps3_oben/eps3_unten)
                self.cs_oben_old=cs_oben_old
                self.cs_unten=cs_unten
                self.cs_oben=cs_oben
                self.cs_unten_old=cs_unten_old
                
                
                self.eps3=eps3
                self.eps2=eps2
                self.eps=eps
                cs3_modded=np.sqrt(cs_oben/cs_unten)
                cs3_old=np.sqrt(cs_oben_old/cs_unten_old)
                self.cs3_modded=cs3_modded
                self.cs3_old=cs3_old
                eps3_with_cs=cs3_old/cs3_modded
                
                print 'this is eps ** min then max **'
                print min(eps)
                print max(eps)
                print 'this is eps2 ** min then max **'
                print min(eps2)
                print max(eps2)
                print 'this is eps2.1 ** min then max **'
                print min(eps2_with_cs)                
                print max(eps2_with_cs)                
                print 'this is eps3 ** min then mean **'
                print np.nanmin(eps3)
                print np.nanmean(eps3)
                
                self.overlap=overlap
                self.m_eff_upstream=m_eff_upstream[overlap]
                self.z_eff_upstream=z_eff_upstream[overlap]
                self.Rs=Rs
                if case=='unmodded':
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.Te_emc3=Te[iTor,iPol,iRad]
                    self.Ti_emc3=Ti[iTor,iPol,iRad]
                
                if case=='modded':
                    if cs=='tok':
                        ne_drews_modded=ne_drews[:,0]*np.nanmin(eps3)
                    if cs=='eff':
                        ne_drews_modded=ne_drews[:,0]*np.nanmin(eps2)
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radial[30,237,:-1], cs2_modded_long/1e6,
                                '-', linewidth = 3, label = r'$c_\mathrm{s}^2$',color='#7e2f8e' )
                    ax.plot( radial[30,237,:-1], cs3_modded/1e6,
                                '--',color='#d91219', linewidth = 3, label = r'$c_\mathrm{s}^3$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$c_\mathrm{s}$ [$\mathrm{10^6m/s}$]' )
                    ax.set_xlim(597,617)
                    #ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_drews_modded=ne_drews_modded
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.Te_emc3=Te[iTor,iPol,iRad]
                    self.Ti_emc3=Ti[iTor,iPol,iRad]


                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews_modded,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.ne_drews_unmodded=ne_drews
                    self.ne_drews_modded=ne_drews_modded
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.Te_emc3=Te[iTor,iPol,iRad]
                    self.Ti_emc3=Ti[iTor,iPol,iRad]




class Plot_all_Profiles_Divertor(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  
  
  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self,impu='C',case='modded',z_eff_downstream = 'none', z_eff_upstream = 'none', minus_shift=0):
            plt.rcParams.update({'font.size': 18})
            g=EMC3.Grid() 
            
            a_eff_emc3,a_eff_exp,radius,r_eff_boyd=get_a_r_eff(g)

            pathH = '.'

            Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
            GH = EMC3.Grid( path = pathH + '/../../geometry', geometryParameters = Geo )
            PC = EMC3.PhysicalCells( path = pathH, geometryParameters = Geo )
            PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
            niH = EMC3.PlasmaField( 'DENSITY', pathH + '/EMC3_OUTPUT', physicalCells = PC )
            TeTiH = EMC3.PlasmaField( 'TE_TI', pathH + '/EMC3_OUTPUT', physicalCells = PC )


            # Multi-purpose manipulator tip plunge minimum and maximum location
            p1Cart = [ -560.482757, -212.503934, -17.065822 ]
            p2Cart = [ -593.011862, -225.421244, -17.035162 ]
            pCart = [ np.linspace( p1Cart[0], p2Cart[0], 40 ),
                      np.linspace( p1Cart[1], p2Cart[1], 40 ),
                      np.linspace( p1Cart[2], p2Cart[2], 40 ) ]

            Rs, phis, Zs = EMC3.Analysis.Cartesian2Grid( *pCart )
            Zs=Zs-minus_shift
            iTor, iPol, iRad = EMC3.Analysis.getCellIndex( zip(Rs, phis, Zs), GH )

            n_e_imp_C=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            n_e_imp_O=np.zeros( len( niH('n1')[0][iTor,iPol,iRad] ) )
            if impu=='C':
                for m in range(1,len(niH.keys())):
                    dummy_false='n'+ str(m) in niH.keys()
                    if dummy_false==False:
                        n_e_imp_C+=n_imp*0                        
                    if dummy_false==True:    
                        n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                        n_e_imp_C+=n_imp*m
                n_e_imp=n_e_imp_C       
                print 'Carbon'
            counter=1    
            if impu!='C':
                for m in range(1,8):
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_C+=n_imp*m
                for m in range(8,len(niH.keys())+1):
                    counter+=1
                    n_imp=niH('n'+ str(m))[0][iTor,iPol,iRad]
                    n_e_imp_O+=n_imp*counter
                n_e_imp=n_e_imp_C+n_e_imp_O
                print 'Carbon plus Oxygen'
                self.n_e_imp_C=n_e_imp_C
                self.n_e_imp_O=n_e_imp_O
            for Te, Ti, n, TiLabel in zip( [ TeTiH('Te')[0] ],
                                                  [ TeTiH('Ti')[0] ],
                                                  [ niH('n1')[0] ],
                                                  [ r'$T_\mathrm{H^{+}}$']):
#a radial shift is implemented here, because of the radial offset of the
#langmuir probe head              
              shift=0
              
              data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_divertor/exp_results/drews/2018/mpm_profiles_joerg180814036.mat')

              Te_drews=data_drews['Te_mpm1']
              Te_drews=Te_drews[:,0]
              ne_drews=data_drews['ne_mpm1']
              ne=ne_drews[:,0]*1e18
              x=data_drews['x_mpm1']
              x=x[:,0]/10
              y=data_drews['y_mpm1']
              y=y[:,0]/10
              radius_drews=np.sqrt(x**2+y**2)*1000
         
              
              self.radius_drews=radius_drews
              
              mi=1.66*1e-27*1.00794
              e=1.602e-19
              Aeff=2e-3*3e-3*2

              T=EMC3.PlasmaField('TE_TI')
              d=EMC3.PlasmaField('DENSITY')
              
              if len(d)==1:
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='orange' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                'b--', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                'r--', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    #ax.set_xlim(608-6+1.3,616)
                    ax.set_ylim(0,1000)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='darkviolet' )
                    ax.plot( Rs, (n[iTor,iPol,iRad]+n_e_imp) * 1e-12,
                                'r--', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    #ax.set_xlim(608-6+1.3,613-6+1.3)
                    #ax.set_ylim(0,1)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)            

              if len(d)>1:
                ne_sum=d_sum(g,d,0,0)

                print 'downstream'
                ne_z_sum_downstream = d_z_sum_impu(g,d,17,102,'C')
                ne_sum_downstream=d_sum_impu(g,d,17,102,'C')
                m_eff_downstream=sum_eff_mass(g,d,17,102,'C')
                print 'upstream'
                ne_z_sum_upstream   = d_z_sum_impu(g,d,30,237,'C')
                ne_sum_upstream=d_sum_impu(g,d,30,237,'C')                
                m_eff_upstream=sum_eff_mass(g,d,30,237,'C')
                if len(d)>7:
                    print 'downstream'
                    ne_z_sum_downstream_C = ne_z_sum_downstream
                    m_eff_downstream_C=m_eff_downstream
                    ne_z_sum_downstream_O = d_z_sum_impu(g,d,17,102,'O')
                    m_eff_downstream_O=sum_eff_mass(g,d,17,102,'O')
                    ne_sum_downstream_C = ne_sum_downstream
                    ne_sum_downstream_O = d_sum_impu(g,d,17,102,'O')
                    z_eff_downstream_C=ne_z_sum_downstream_C/ne_sum_downstream_C
                    z_eff_downstream_O=ne_z_sum_downstream_O/ne_sum_downstream_O
                    print 'upstream'
                    ne_z_sum_upstream_C   = ne_z_sum_upstream
                    m_eff_upstream_C=m_eff_upstream
                    ne_z_sum_upstream_O   = d_z_sum_impu(g,d,30,237,'O')
                    m_eff_upstream_O=sum_eff_mass(g,d,30,237,'O')
                    ne_sum_upstream_C   = ne_sum_upstream
                    ne_sum_upstream_O   = d_sum_impu(g,d,30,237,'O')
                    z_eff_upstream_C=ne_z_sum_upstream_C/ne_sum_upstream_C
                    z_eff_upstream_O=ne_z_sum_upstream_O/ne_sum_upstream_O
                    
                    ne_z_sum_downstream=ne_z_sum_downstream_C+ne_z_sum_downstream_O
                    self.dummy1=ne_z_sum_downstream
                    ne_sum_downstream=ne_sum_downstream_C+ne_sum_downstream_O
                    self.dummy2=ne_sum_downstream
                    ne_z_sum_upstream=ne_z_sum_upstream_C+ne_z_sum_upstream_O
                    ne_sum_upstream=ne_sum_upstream_C+ne_sum_upstream_O
                    
                z_eff_downstream=ne_z_sum_downstream/ne_sum_downstream
                z_eff_upstream=ne_z_sum_upstream/ne_sum_upstream        

                self.z_eff_downstream=z_eff_downstream
                self.z_eff_upstream=z_eff_upstream
                self.m_eff_upstream=m_eff_upstream
                self.m_eff_downstream=m_eff_downstream
                self.Rs=Rs
                self.iTor=iTor
                self.iPol=iPol
                self.iRad=iRad

                overlap=[]
                for i in range(0,radius_drews.shape[0]):
                    overlap.append(find_nearest(Rs,radius_drews[i])[1])
                
                self.overlap=overlap
                self.Te=Te

                if case=='unmodded':
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, Te_drews,
                                '-', linewidth = 3, label = r'$T_\mathrm{e,drews}$',color='#ffa500' )
                    ax.plot( Rs, zero_to_nan(Te[iTor,iPol,iRad]),
                                '--',color='#0070bd', linewidth = 3, label = r'$T_\mathrm{e}$' )
                    ax.plot( Rs, zero_to_nan(Ti[iTor,iPol,iRad]),
                                '--',color='#d91219', linewidth = 3, label = r'$T_\mathrm{i}$' )              
                    ax.set_xlabel( r'$R~[\mathrm{cm}]$' )
                    ax.set_ylabel( r'$T~[\mathrm{eV}]$' ) 
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    ax.set_ylim(0,100)
                    ax.grid()
                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    f, ax = plt.subplots()
                    ax.plot( radius_drews, ne_drews,
                                '-', linewidth = 3, label = r'$n_\mathrm{e,drews}$',color='#7e2f8e' )
                    ax.plot( Rs, zero_to_nan((n[iTor,iPol,iRad]+n_e_imp) * 1e-12),
                                '--',color='#d91219', linewidth = 3, label = r'$n_\mathrm{e}$' )
                    ax.set_xlabel( r'$R$ [cm]' )
                    ax.set_ylabel( r'$n_\mathrm{e}$ [$10^{18}$ m$^{-3}$]' )
                    ax.set_xlim(608-6+1.3,613-6+1.3)
                    #ax.set_ylim(0,6)
                    ax.grid()

                    plt.tight_layout()
                    for label in ax.xaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    for label in ax.yaxis.get_ticklabels()[::2]:
                        label.set_visible(False)
                    
                    self.radius_drews=radius_drews
                    self.Te_drews=Te_drews
                    self.ne_drews=ne_drews
                    self.Rs=Rs
                    self.neEMC3=(n[iTor,iPol,iRad]+n_e_imp)*1e-12
                    self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                    self.ne_drews_unmodded=ne_drews
                    self.ne_emc3=((n[iTor,iPol,iRad]+n_e_imp) * 1e-12)
                    self.TeEMC3=zero_to_nan(Te[iTor,iPol,iRad])
                    self.TiEMC3=zero_to_nan(Ti[iTor,iPol,iRad])
                    self.TeEMC3complete=Te
                    self.TiEMC3complete=Ti
                    self.TEMC3complete=T
                    self.dEMC3complete=d









class plot_thomson_damm(object):
    def __init__(self,time):
            plt.rcParams.update({'font.size': 18})

            g=EMC3.Grid()
            z=g.zones[0]
            d=EMC3.PlasmaField('DENSITY')
            T=EMC3.PlasmaField('TE_TI')
            n1=d('n1')[0][55,300,:]
            TeEMC3=T('Te')[0][55,300,:]

            filename=('/home/j.cosfeld/runs/W7X_divertor/exp_results/damm/thomson_data_20171109.045_new_spectral_cal_to_Int.csv')    
            
            # usecols select columns
            data = np.loadtxt(filename, delimiter=';', skiprows=1)
            t=data[1:-1,1]
            r_eff_damm=data[1:-1,2]
            ne=data[1:-1,3]
            Te=data[1:-1,6]
            shift=time*16
            
            
            f, ax = plt.subplots()

            plt.plot(r_eff_damm[0+shift:16+shift],ne[0+shift:16+shift],'o',linewidth=3.0)
            plt.plot(r_eff_damm[10+shift:10+16+shift],ne[10+shift:10+16+shift],'o',linewidth=3.0)
            plt.plot(r_eff_damm[20+shift:20+16+shift],ne[20+shift:20+16+shift],'o',linewidth=3.0)
            plt.plot(z.Reff[0:-1]/100,n1/1e12,'--',linewidth=3.0)
            
            ax.set_ylabel( r'$n_\mathrm{e}[10^{18}/\mathrm{m^3}]$')
            ax.set_xlabel( r'$r_\mathrm{eff}$')
            ax.grid() 
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
        ##        ax.set_xlim(0,120)
            #plt.rcParams.update({'font.size': 18})
            f.tight_layout()

            f, ax = plt.subplots()

            plt.plot(r_eff_damm[0+shift:16+shift],Te[0+shift:16+shift],'o',linewidth=3.0)
            plt.plot(r_eff_damm[10+shift:10+16+shift],Te[10+shift:10+16+shift],'o',linewidth=3.0)
            plt.plot(r_eff_damm[20+shift:20+16+shift],Te[20+shift:20+16+shift],'o',linewidth=3.0)
            plt.plot(z.Reff[0:-1]/100,TeEMC3/1e3,'--',linewidth=3.0)
            
            ax.set_ylabel( r'$T_\mathrm{e}[\mathrm{eV}]$')
            ax.set_xlabel( r'$r_\mathrm{eff}$')
            ax.grid() 
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
            ax.set_ylim(0,5)
            #plt.rcParams.update({'font.size': 18})
            f.tight_layout()

            self.data=data
            self.t=t
            self.r_eff_damm=r_eff_damm
            self.ne=ne
            self.Te=Te
            



class Plot_thomson_compare(object):
    """ Returns the a comparison plot of the MLP profiles 

  Keyword arguments:
  impu -- if impu string is not C a simulation with more ion specied than carbon impurities is considered
  P_sol -- string for labeling the SOL power 
  P_rad -- string for labeling the radiated power
  
  
  TODO:
  
  
  
  ATTENTION:
  Te=Ti is assumed for the claculation of csi_drews!
  
    """  


    def __init__(self):
            plt.rcParams.update({'font.size': 24})
            g=EMC3.Grid() 
            
            z=g.zones[0]
            rad=z.R
            tor=z.phi
            pol=z.Z
            radius=z.Reff

            data_flecken3 = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_divertor/exp_results/flecken/171026038neTS.mat')    
            T=EMC3.PlasmaField('TE_TI')
            d=EMC3.PlasmaField('DENSITY')

            f, ax = plt.subplots()
            for m in range(0, np.shape(data_flecken3['ne']['y'][0][0][:])[0] ):
                plt.plot(data_flecken3['ne']['x1'][0][0][0],data_flecken3['ne']['y'][0][0][m],'b')
            plt.plot(radius[:-1]/100,zero_to_nan(d('n1')[0][0,255,:]/1e12),'r--',linewidth=3.0)
            ax.set_ylabel( r'$n_\mathrm{e}[10^{18}/\mathrm{m^3}]$')
            ax.set_xlabel( r'$r_\mathrm{eff}$')
            ax.grid() 
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
        #        ax.set_xlim(0,120)
            plt.rcParams.update({'font.size': 18})
            f.tight_layout()





class Cut_Pol_Avg_Sum(object):
    def __init__(self):
    #
    #Gives the normed Density of Carbon impurities (normed at every radial position)
    #
        #plt.ion()

        plt.matplotlib.rcParams.update({'font.size': 18})

        g=EMC3.Grid()
        G=EMC3.Grid('GRID_3D_DATA')
        z=G.zones[0]
        radial=z.R
        poloidal=z.Z

        pathH = '.'


        Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
        PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )
    #
    #Calculate the area of a Grid Cell and thereby the effective raidus
    #
        area=np.zeros(poloidal.shape[0])
        radius=np.zeros(poloidal.shape[2])

        for k in range(0,poloidal.shape[2]):
            for j in range(0,poloidal.shape[0]):
                area[j]=EMC3.Analysis.getAreaPolygon(radial[j,:,k],poloidal[j,:,k])

            area_mean=np.mean(area)
            radius[k]=np.sqrt(area_mean/math.pi)

        dummy=0
        lim_index=np.zeros(poloidal.shape[0])
        lcfs_index=0


        radius_mid=np.zeros(poloidal.shape[2])
        radius_plot=np.zeros(poloidal.shape[2])
        for i in range(1, poloidal.shape[2]):
            radius_mid[i]=(radius[i]-radius[i-1])/2
        radius_plot=radius-radius_mid

    #        
    #lcfs_index gives the index in an radial array, at which the limiter last close flux surface is placed
    #
        for j in range(0,poloidal.shape[0]):
            dummy = np.append(PlatesH.ids[0][0,j,:],[True])
            lim_index[j]=dummy.tolist().index(True)
        lcfs_index=int(min(lim_index[:]))

        self.g = g
        self.G = G
        self.radial = radial
        self.poloidal = poloidal
        self.area = area
        self.radius = radius
        self.PlatesH=PlatesH
        self.lcfs_index=lcfs_index
    #
    # Case check if array index in plot is wanted or not
    # Case check which measureand is wanted by the user 
    #
        data=EMC3.NeutralParameters('./../run/input_modded.n0g')
        d_eirene=np.zeros(7)
        for i in range(0, 7):
            d_eirene[i]=data.cells[i]['n']
        d=EMC3.PlasmaField('DENSITY')
        matrix_length=g.geometry.NSurfaces[0][0]
        matrix_high=len(d)
        avg=np.zeros(shape=(matrix_high,matrix_length-1))
#
# Calculate the average density for each ionization state along the radial axis 
#

        dens_sum=np.zeros(matrix_length-1)

        avg[1,:]=np.average(d('n1')[0],(0,1),g.zones[0].volume)
        for m in range(2,len(d)):
            rank=str(m)
            avg[m,:]=np.average(d('n'+ rank)[0],(0,1),g.zones[0].volume)
            dens_sum[:]+=avg[m,:]*m

        plot_arry=np.zeros(np.size(avg[1,:])-8)
        counter=0
        for j in range(8, np.size(avg[1,:])):
            plot_arry[counter]=j
            counter=counter+1

        f, ax = plt.subplots()
#        ax.plot(radius_plot[7:-3]/radius_plot[lcfs_index],dens_sum[8:-1],linewidth=2,color="blue")
#        ax.plot(radius_plot[7:-3]/radius_plot[lcfs_index],avg[1,8:-1],linewidth=2,color="red")
        ax.plot(radius_plot[7:-3]/radius_plot[lcfs_index],dens_sum[8:-1]/avg[1,8:-1],linewidth=2,color="blue")
        ax.set_ylabel( r'$n_e/n_{e,\mathrm{imp}}$')
        ax.set_xlabel( r'$r/r_{LCFS}$')
        ax.grid()
        f.tight_layout()
        plt.setp(ax.get_xticklabels()[::2], visible=False)


        f, ax = plt.subplots()
        ax.plot(radius_plot[7:-3]/radius_plot[lcfs_index],dens_sum[8:-1],linewidth=2,color="blue")
        ax.plot(radius_plot[7:-3]/radius_plot[lcfs_index],avg[1,8:-1],linewidth=2,color="red")
        ax.set_ylabel( r'$n[1/m^3]$')
        ax.set_xlabel( r'$r/r_{LCFS}$')
        ax.grid()
        f.tight_layout()
        plt.setp(ax.get_xticklabels()[::2], visible=False)

        self.avg=avg
        self.dens_sum=dens_sum


def fit_func(x, a, b):
    return a*arctan(-b*x)

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]



def get_a_r_eff_div(g):
    z=g.zones[0]
    radial=z.R
    poloidal=z.Z
    toroidal=z.phi
#getting effective radius of the machine
    radius=z.Reff

#Boyd Data Load
    ken=np.loadtxt('/home/j.cosfeld/runs/W7X_divertor/exp_results/Kenneth/coordinatesLD.txt')


    x_ken=ken[:,0]/10
    y_ken=ken[:,1]/10
    z_ken=ken[:,2]/10
    r_ken=np.sqrt(x_ken**2+y_ken**2)

    index_ken=np.zeros(len(r_ken))
    r_eff_ken=np.zeros(len(r_ken))
    
    index_ken=np.asarray(EMC3.Analysis.getCellIndex(zip(r_ken,[8.39999]*len(r_ken),z_ken),g)[2])
    index_ken=np.nan_to_num(index_ken)
    r_eff_ken=radius[index_ken.astype(int)]

    return radius,r_eff_ken

def get_a_r_eff(g):
    z=g.zones[0]
    radial=z.R
    poloidal=z.Z
    toroidal=z.phi
#getting effective radius of the machine
    radius=z.Reff

#getting the index of the limiter RLimTip
    pathH = '.'

    Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
    PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )

    lim_index=np.zeros(poloidal.shape[0])
    lcfs_index_min=np.zeros(toroidal.shape[0]-1)

    for i in range(0,toroidal.shape[0]-1):
        for j in range(0,poloidal.shape[0]):
            dummy = np.append(PlatesH.ids[0][i,j,:],[True])
            lim_index[j]=dummy.tolist().index(True)
        lcfs_index=min(lim_index)
        lcfs_index_min[i]=lcfs_index
        
    a_eff=radius[min(lcfs_index_min)]
    a_eff_emc3=a_eff
#Limiter coords of the W7X device


    i=EMC3.Installation('LIM_R570.56_Z43','~/devices/W7X/Limiter')
    lim_r=(i.R[0]+0.1).tolist()
    lim_z=i.Z[0].tolist()
    lim_phi=i.phi.tolist()

    limindex_device=np.asarray(EMC3.Analysis.getCellIndex(zip(lim_r,[0]*len(lim_r),[0]*len(lim_r)),g)[2])
    a_eff_exp=radius[np.nanmin(limindex_device[:])]

#Boyd Data Load
    data_json=json.load(open('/home/j.cosfeld/runs/W7X_limiter/exp_results/boyd/LP20160308_24_L53_2k2.json'))

    x_boyd=np.zeros(np.array(data_json['ne18']).shape[1])
    y_boyd=np.zeros(np.array(data_json['ne18']).shape[1])
    z_boyd=np.zeros(np.array(data_json['ne18']).shape[1])
    r_boyd=np.zeros(np.array(data_json['ne18']).shape[1])
    max_data_boyd=np.zeros(np.array(data_json['ne18']).shape[1])
    
    for i in range (0,np.array(data_json['ne18']).shape[1]):

        x_boyd[i]=data_json['info']['coords'][i][1]*(-100)
        y_boyd[i]=data_json['info']['coords'][i][0]*100
        z_boyd[i]=data_json['info']['coords'][i][2]*100
        r_boyd[i]=np.sqrt(x_boyd[i]*x_boyd[i]+y_boyd[i]*y_boyd[i])

    index_boyd=np.zeros(len(r_boyd))
    r_eff_boyd=np.zeros(len(r_boyd))
    
    index_boyd=np.asarray(EMC3.Analysis.getCellIndex(zip(r_boyd,[0]*len(r_boyd),z_boyd),g)[2])        
    r_eff_boyd=radius[index_boyd]

    


    return a_eff_emc3,a_eff_exp,radius,r_eff_boyd

def get_a_r_eff2(g):
    z=g.zones[0]
    radial=z.R
    poloidal=z.Z
    toroidal=z.phi
#getting effective radius of the machine
    radius=z.Reff

#getting the index of the limiter RLimTip
    pathH = '.'

    Geo = EMC3.GeometryParameters( path = pathH + '/../../geometry' )
    PlatesH = EMC3.PlateCells( path = pathH + '/../../geometry', geometryParameters = Geo )

    lim_index=np.zeros(poloidal.shape[0])
    lcfs_index_min=np.zeros(toroidal.shape[0]-1)

    for i in range(0,toroidal.shape[0]-1):
        for j in range(0,poloidal.shape[0]):
            dummy = np.append(PlatesH.ids[0][i,j,:],[True])
            lim_index[j]=dummy.tolist().index(True)
        lcfs_index=min(lim_index)
        lcfs_index_min[i]=lcfs_index
        
    a_eff=radius[min(lcfs_index_min)]
    a_eff_emc3=a_eff
#Limiter coords of the W7X device


    i=EMC3.Installation('LIM_R570.56_Z43','~/devices/W7X/Limiter')
    lim_r=(i.R[0]+0.1).tolist()
    lim_z=i.Z[0].tolist()
    lim_phi=i.phi.tolist()

    limindex_device=np.asarray(EMC3.Analysis.getCellIndex(zip(lim_r,[0]*len(lim_r),[0]*len(lim_r)),g)[2])
    a_eff_exp=radius[np.nanmin(limindex_device[:])]


    data_drews = scipy.io.loadmat('/home/j.cosfeld/runs/W7X_limiter/exp_results/drews/base/profiles.mat')
    x=data_drews['x']
    x22=x[:,0]/10
    y=data_drews['y']
    y22=y[:,0]/10
    z=data_drews['z']
    z22=z[:,0]/10
    radius22=np.sqrt(x22**2+y22**2)

    index_drews=np.zeros(len(radius22))
    r_eff_boyd=np.zeros(len(radius22))
    
    index_drews=np.asarray(EMC3.Analysis.getCellIndex(zip(radius22,[16]*len(radius22),z22),g)[2])        
    r_eff_drews=radius[index_drews]

    


    return a_eff_emc3,a_eff_exp,radius,r_eff_drews


def d_sum_avg(g,d):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    avg=np.zeros(shape=(matrix_high,matrix_length-1))
#
# Calculate the average density for each ionization state along the radial axis 
#

    dens_sum=np.zeros(matrix_length-1)

    for m in range(2,len(d)):
        rank=str(m)
        avg[m,:]=np.average(d('n'+ rank)[0],(0,1),g.zones[0].volume)
        dens_sum[:]+=avg[m,:]*m
               

    return dens_sum

def d_z_sum(g,d,phi_slice,z):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high,matrix_length-1))
#
# Calculate the average density for each ionization state along the radial axis 
#

    dens_sum=np.zeros(matrix_length-1)
    dummy[1,:]=d('n1')[0][phi_slice,z,:]
    for m in range(2,len(d)):
        rank=str(m)
        dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
        dens_sum[:]+=dummy[m,:]*m*m
    dens_sum[:]+=dummy[1,:]

    return dens_sum



def d_sum_paraview(d):
#
# Calculate the average density for each ionization state along the radial axis 
#
    dummy_save=0
    for m in range(1,len(d)):
        rank=str(m)
        dummy=np.array(d('n'+ rank))
        dummy_save+=dummy
    dens_sum=np.array(d('n1'))+dummy_save

    return dens_sum


def d_sum_impu(g,d,phi_slice,z,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        print 'belonging to C, Nenner'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,:]
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy_false='n'+ rank in d.keys()
            if dummy_false==False:
                dummy[m,:]=0
            if dummy_false==True:
                dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
                print np.max(d('n'+ rank)[0][phi_slice,z,:])/1e13*(m-1)
            dens_sum[:]+=dummy[m,:]*(m-1)
        dens_sum[:]+=dummy[1,:]
        print max(dens_sum)/1e13
    
    if impu=='O':
        print 'belonging to O, Nenner'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,:]=0
        for m in range(2+Z_C-1,len(d)):
            rank=str(m)
            dummy_false='n'+ rank in d.keys()
            if dummy_false==False:
                dummy[m,:]=0
            if dummy_false==True:    
                dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
                print np.max(d('n'+ rank)[0][phi_slice,z,:])/1e13*(m-Z_C)
            dens_sum[:]+=dummy[m,:]*(m-Z_C)
        dens_sum[:]+=dummy[1,:]
        print max(dens_sum)/1e13


    return dens_sum


def cs3oben(g,d,Te,Ti,phi_slice,z,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        cs=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,:]*1e6*Ti[0][phi_slice,z,:]*1.602e-19+d('n1')[0][phi_slice,z,:]*1e6*Te[0][phi_slice,z,:]*1.602e-19
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]*1e6*Ti[0][phi_slice,z,:]*1.602e-19+d('n'+ rank)[0][phi_slice,z,:]*1e6*Te[0][phi_slice,z,:]*(m-1)*1.602e-19
            cs[:]+=dummy[m,:]
        cs[:]+=dummy[1,:]

    return cs


def cs3oben_clean(g,d,Te,Ti,phi_slice,z,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        cs=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,:]*1e6*Ti[0][phi_slice,z,:]*1.602e-19+d('n1')[0][phi_slice,z,:]*1e6*Te[0][phi_slice,z,:]*1.602e-19
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]*1e6*Ti[0][phi_slice,z,:]*1.602e-19+d('n'+ rank)[0][phi_slice,z,:]*1e6*Te[0][phi_slice,z,:]*1.602e-19
            cs[:]+=dummy[m,:]
        cs[:]+=dummy[1,:]

    return cs


def cs3unten(g,d,T,phi_slice,z,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
    u=1.6605e-27
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        cs=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,:]*1e6*1.6605e-27
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]*1e6*2*Z_C*1.6605e-27
            cs[:]+=dummy[m,:]
        cs[:]+=dummy[1,:]

    return cs

def cs3unten_clean(g,d,T,phi_slice,z,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
    u=1.6605e-27
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        cs=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,:]*1e6*1.6605e-27
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]*1e6*1.6605e-27
            print max(dummy[m,:])
            cs[:]+=dummy[m,:]
        cs[:]+=dummy[1,:]

    return cs


def sum_eff_mass(g,d,phi_slice,z,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    dummy2=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_H=1
    Z_C=6
    Z_O=16
    u=1.66*1e-27
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        dens_sum=np.zeros(matrix_length-1)
        dens_sum2=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,:]*Z_H*u
        dummy2[1,:]=d('n1')[0][phi_slice,z,:]
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy_false='n'+ rank in d.keys()
            if dummy_false==False:
                dummy[m,:]=0
            if dummy_false==True:
                dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
                dummy2[m,:]=d('n'+ rank)[0][phi_slice,z,:]
            dens_sum[:]+=dummy[m,:]*Z_C*2*u
        dens_sum[:]+=dummy[1,:]
        dens_sum2[:]+=dummy2[1,:]
    
    if impu=='O':
        dens_sum=np.zeros(matrix_length-1)
        dens_sum2=np.zeros(matrix_length-1)
        dummy[1,:]=0
        for m in range(2+Z_C-1,len(d)):
            rank=str(m)
            dummy_false='n'+ rank in d.keys()
            if dummy_false==False:
                dummy[m,:]=0
            if dummy_false==True:    
                dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
                dummy2[m,:]=d('n'+ rank)[0][phi_slice,z,:]
            dens_sum[:]+=dummy[m,:]*Z_O*2*u
        dens_sum[:]+=dummy[1,:]
        dens_sum2[:]+=dummy2[1,:]

    dens_sum=dens_sum/dens_sum2
    return dens_sum



def d_sum_impu_point(g,d,phi_slice,z,ir,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        print 'belonging to C, Nenner'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,ir]=d('n1')[0][phi_slice,z,ir]
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy[m,ir]=d('n'+ rank)[0][phi_slice,z,ir]
            print np.max(d('n'+ rank)[0][phi_slice,z,ir])/1e13*(m-1)
            dens_sum[ir]+=dummy[m,ir]*(m-1)
        dens_sum[ir]+=dummy[1,ir]
        print max(dens_sum)/1e13
    
    if impu=='O':
        print 'belonging to O, Nenner'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,:]=0
        for m in range(2+Z_C-1,len(d)):
            rank=str(m)
            dummy[m,ir]=d('n'+ rank)[0][phi_slice,z,ir]
            print np.max(d('n'+ rank)[0][phi_slice,z,ir])/1e13*(m-Z_C)
            dens_sum[ir]+=dummy[m,ir]*(m-Z_C)
        dens_sum[ir]+=dummy[1,ir]
        print max(dens_sum)/1e13


    return dens_sum



def d_z_sum_impu(g,d,phi_slice,z,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        print 'belonging to C, Zaehler'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,:]
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy_false='n'+ rank in d.keys()
            if dummy_false==False:
                dummy[m,:]=0            
            if dummy_false==True:
                dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
                print np.max(d('n'+ rank)[0][phi_slice,z,:])/1e13*(m-1)*(m-1)
            dens_sum[:]+=dummy[m,:]*(m-1)*(m-1)
        dens_sum[:]+=dummy[1,:]
        print max(dens_sum)/1e13
    
    if impu=='O':
        print 'belonging to O, Zaehler'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,:]=0
        for m in range(2+Z_C-1,len(d)):
            rank=str(m)
            dummy_false='n'+ rank in d.keys()
            if dummy_false==False:
                dummy[m,:]=0   
            if dummy_false==True:                
                dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
                print np.max(d('n'+ rank)[0][phi_slice,z,:])/1e13*(m-Z_C)*(m-Z_C)
            dens_sum[:]+=dummy[m,:]*(m-Z_C)*(m-Z_C)
        dens_sum[:]+=dummy[1,:]
        print max(dens_sum)/1e13


    return dens_sum


def d_z_sum_impu_point(g,d,phi_slice,z,ir,impu):

    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high+1,matrix_length-1))
    
    Z_C=6
    Z_O=16
#
# Calculate the average density for each ionization state along the radial axis 
#
    if impu=='C':
        print 'belonging to C, Zaehler'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,:]=d('n1')[0][phi_slice,z,ir]
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy[m,:]=d('n'+ rank)[0][phi_slice,z,ir]
            print np.max(d('n'+ rank)[0][phi_slice,z,ir])/1e13*(m-1)*(m-1)
            dens_sum[:]+=dummy[m,:]*(m-1)*(m-1)
        dens_sum[:]+=dummy[1,:]
        print max(dens_sum)/1e13
    
    if impu=='O':
        print 'belonging to O, Zaehler'
        dens_sum=np.zeros(matrix_length-1)
        dummy[1,:]=0
        for m in range(2+Z_C-1,len(d)):
            rank=str(m)
            dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
            print np.max(d('n'+ rank)[0][phi_slice,z,:])/1e13*(m-Z_C)*(m-Z_C)
            dens_sum[:]+=dummy[m,:]*(m-Z_C)*(m-Z_C)
        dens_sum[:]+=dummy[1,:]
        print max(dens_sum)/1e13


    return dens_sum




def z_eff_paraview(d,impu):

    number_of_Z=len(d)-1
    len_of_Z=len(d['n1'][1])
    Z_C=6
    Z_O=16
#    
# Calculate Z_eff for each GridCell
#
    if impu=='C':
        dummy_save=0
        for m in range(2,2+Z_C):
            rank=str(m)
            dummy=np.array(d('n'+ rank))
            dummy_save+=dummy*m*m
            print np.max(dummy_save[0][0,0])
        dens_sum=np.array(d('n1'))+dummy_save
    
    if impu=='O':
        dummy_save=0
        for m in range(2+Z_C,len(d)):
            rank=str(m)
            dummy=np.array(d('n'+ rank))
            dummy_save+=dummy*m*m
            print np.max(dummy_save[0][0,0])
        dens_sum=np.array(d('n1'))+dummy_save

    return dens_sum

def d_sum(g,d,phi_slice,z):
    
    matrix_length=g.geometry.NSurfaces[0][0]
    matrix_high=len(d)
    dummy=np.zeros(shape=(matrix_high,matrix_length-1))
#
# Calculate the average density for each ionization state along the radial axis 
#
    dens_sum=np.zeros(matrix_length-1)
    dummy[1,:]=d('n1')[0][phi_slice,z,:]
    for m in range(2,len(d)):
        rank=str(m)
        dummy_false='n'+ rank in d.keys()
        if dummy_false==False:
            dummy[m,:]=0
        if dummy_false==True:
            dummy[m,:]=d('n'+ rank)[0][phi_slice,z,:]
        dens_sum[:]+=dummy[m,:]
    dens_sum[:]+=dummy[1,:]

    return dens_sum

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx
