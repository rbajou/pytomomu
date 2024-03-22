#!/usr/bin/python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.integrate import quad
from scipy import interpolate

class FluxModel:
    
    
    def __init__(self, altitude:float=0., corsika_flux:np.ndarray=None, corsika_std:np.ndarray=None, energy_bins:np.ndarray=None, theta_bins:np.ndarray=None):
        self.altitude = altitude
        ####Common parameters
        self.par = {
            'mu': 105.658,#MeV
            "B" : 0.054,
            "Epi":115/1.1,
            "Ek" : 850/1.1, #or 810??? (Lesparre et al. 2010)
            "gamma" : 2.70 
        }
        self.available = ["gaisser_std","gaisser", "tang",  "guan", "shukla"]#, ""]
        self.is_corsika = False
        if all(item is not None for item in [corsika_flux,corsika_std, energy_bins, theta_bins]): 
            self.available.append("corsika_soufriere")
            self.is_corsika = True
            theta_bw = abs(theta_bins[:-1]-theta_bins[1:])
            x, y = energy_bins, theta_bins #x=column coordinate, y=row
            ###cut on energy 
            emax = 1e5#1.5e3 #GeV
            arg_emax = np.min(np.where(energy_bins > emax)[0])
            xslice = slice(1, arg_emax)
            ###interpolate empty flux values (nan surrounded by non-nan values, does not work for nan at the edges of the grid)
            mean= corsika_flux[:, xslice]
            std= corsika_std[:, xslice]
            X, Y = np.meshgrid(x,y)
            X, Y  = X[:, xslice], Y[:, xslice]
            xnew  = np.logspace(np.log10(np.nanmin(x[xslice])), np.log10(np.nanmax(x[xslice])), 100)
            ynew = np.linspace(np.nanmin(y), np.nanmax(y), 100)
            Xnew, Ynew =  np.meshgrid(xnew,ynew)
            for i, v in enumerate([mean, std]):
                nnan = ((~np.isnan(v)) & (v!=0.))
                points = np.zeros(shape=(len(v[nnan].flatten()),2))
                points[:,0] = X[nnan].flatten()
                points[:,1]= Y[nnan].flatten()
                values = v[nnan].flatten()
                grid_v = interpolate.griddata(points, values, (Xnew, Ynew), method='linear')#, *args, **kwargs)
                # print("grid_v=", grid_v)
                # fig, ax = plt.subplots()
                # ax.imshow(grid_v)
                # plt.show()
                
                ###fill nan value on top right of the flux matrix (thetamin, emax) 
                #grid_v[np.where(np.isnan(grid_v))] = np.mean([grid_v[np.argwhere(np.isnan(grid_v))[0][0]+1, np.argwhere(np.isnan(grid_v))[0][1]], grid_v[np.argwhere(np.isnan(grid_v))[0][0], np.argwhere(np.isnan(grid_v))[0][1]-1]])
                ##cut on theta
                #yslice = slice(0, len(theta_bins)) #all values
                if i==0:corsika_flux = grid_v
                if i==1:corsika_std = grid_v

            #x, y, z, zerr = energy_bins[xslice], theta_bins[yslice],  corsika_flux, corsika_std
            x, y, z, zerr = xnew, ynew,  corsika_flux, corsika_std
            #print(f"x, y, z = {x}, {y}, {z}")
            self.DF_corsika = interpolate.interp2d(x, y, z, fill_value="extrapolate", kind='linear')
            self.DF_corsika_std = interpolate.interp2d(x, y,  zerr, fill_value="extrapolate", kind='linear')
            #print("test_interp = ", self.DF_corsika(100, 0.5)) 
            #exit()
    
        ###For Guan (2015) and Tang et al. (2006) models
    def cosThetaStar(self, theta):
        """
        'theta' in rad
        """
        p_1 = 0.102573 
        p_2 = -0.068287
        p_3 = 0.958633 
        p_4 = 0.0407253 
        p_5 = 0.817285
        cosThetaStar = np.sqrt( (np.cos(theta)**2 + p_1**2 + p_2 * np.cos(theta)**p_3 + p_4 * np.cos(theta)**p_5) / (1 + p_1**2 + p_2 + p_4) )
        return cosThetaStar
    
    def gaisser_std(self, energy, theta):
        """
        Gaisser 1990
        muon 'energy' at the surface in GeV
        'theta' in rad
        """
        B,Epi,Ek,gamma=self.par["B"],self.par["Epi"],self.par["Ek"],self.par["gamma"]
        P = 1/(1+energy*np.cos(theta)/Epi)+B/(1+energy*np.cos(theta)/Ek)
        diffFluxAtSeaLevel = 0.14 * energy**(-gamma) * P
        return diffFluxAtSeaLevel 
      
    def gaisser(self, energy, theta):
        """
        Modified Gaisser formula with energy loss in atm D_E (50 GeV -> 99% ionisation (Bichsel et al. 2010)) and  W term = probability for muons to reach sea level (Dar, 1983)
        muon 'energy' at the surface in GeV
        'theta' in rad
        """
        B,Epi,Ek,gamma=self.par["B"],self.par["Epi"],self.par["Ek"],self.par["gamma"]
        T = np.arccos(np.sqrt(1-(1-(np.cos(theta))**2)/(1+32/6371)**2)) # %H_atm = 32km
        D_E = (2.06*1e-3)*(1030./np.cos(T)-120)  # Muon energy loss through the atmosphere [GeV] 
        E_top = energy + D_E  # Muon energy at production level
        ##Production
        r_c = 1e-4 # ratio of the prompt muons to pions
        P = 1/(1+E_top*np.cos(T)/Epi)+B/(1+E_top*np.cos(T)/Ek)+r_c
        ##Survival
        E_M = energy + D_E/2  # Muon mean energy through the atmosphere
        W = (120 * np.cos(T) / 1030)**(1.04/(np.cos(T)*E_M)) 
        ## Differential Flux at Sea Level
        diffFluxAtSeaLevel = 0.14 * E_top**(-gamma) * P * W 
        return diffFluxAtSeaLevel
    
    def guan(self, energy, theta):
        """ 
        Guan et al. (2015)
        Modified Gaisser parametrization
        """     
        B,Epi,Ek,gamma=self.par["B"],self.par["Epi"],self.par["Ek"],self.par["gamma"]
        cosThetaStar = self.cosThetaStar(theta)
        E_star = energy * (1 + 3.64/(energy * (cosThetaStar**1.29)))
        #####Production
        P = 1./(1+energy*cosThetaStar/Epi)+B/(1+energy*cosThetaStar/Ek)
        ####Differential Flux at Sea Level
        diffFluxAtSeaLevel = 0.14 * E_star**(-gamma)* P 
        return diffFluxAtSeaLevel

    def tang(self, energy, theta):
        """
        Tang et al. (2006)
        Modified Gaisser parametrization
        """
        B,Epi,Ek,gamma=self.par["B"],self.par["Epi"],self.par["Ek"],self.par["gamma"]
        cosThetaStar = self.cosThetaStar(theta)      
        a_atm =   2.06*1e-3 ### GeV/mwe
        #a_atm = a_atm * 1e2 #GeV/(g/cm^-2)
        h_f = 950 #atmosphère opacity (g/cm^-2)
        lambda_N = 90#opacité moyenne entre point d'entrée primary et point de production des muons (g/cm^-2)
        D_E =  a_atm*(h_f/cosThetaStar - lambda_N) ### Muon's energy loss through the atmosphere [GeV] 
        E_top = lambda e : e + D_E ### Muon's energy at production level
        A = lambda e : 1.1 * (lambda_N * np.sqrt(np.cos(theta)+0.001)/h_f)**(4.5/(E_top(e)*cosThetaStar))  #### Emu0 (ground) -> E_top(e) at production ERROR in TANG article ????
        r_c = 1e-4### ratio of the prompt muons to pions
        P = lambda e : 1/(1+E_top(e)*cosThetaStar/Epi) + B/(1+E_top(e)*cosThetaStar/Ek) + r_c
        

        #####
        # 3 regimes : 
        #1) Eμ0 > (100/cosθ⋆) GeV (the standard Gaisser formula),
        if 100/cosThetaStar <= energy:
            diffFluxAtSeaLevel = self.gaisser(energy, theta) 
        #2) (1/cosθ⋆) < Eμ0 ≤ (100/cosθ⋆) GeV (Eqs. (4)–(7)) 
            return diffFluxAtSeaLevel
        elif (1/cosThetaStar < energy ) and ( energy < 100/cosThetaStar ):
            pass
        #3) and Eμ0 ≤ (1/cosθ⋆) GeV (Eq. (9))
        else :
            #print("Eμ0 ≤ (1/cosθ⋆)  ",energy, "GeV") 
            energy = (3*energy+7/cosThetaStar)/10  #### sign difference in https://www.frontiersin.org/articles/10.3389/fenrg.2021.750159/full
            #print("new energy  ",energy,"GeV") 
            #theta = np.arccos(cosThetaStar)
            #diffFluxAtSeaLevel = self.gaisser_std(energy, theta) 
            #print(f"1/cosThetaStar,  100/cosThetaStar={1/cosThetaStar}, {100/cosThetaStar}GeV")

        diffFluxAtSeaLevel = A(energy) * 0.14 * energy**(-gamma) * P(energy)
        #diffFluxAtSeaLevel = self.gaisser(energy, theta) 
        return diffFluxAtSeaLevel


    def shukla(self, energy, theta, E0, Ec, I0, n, eps):
        '''
        E0 : energy loss due to hadronic and electromagnetic interactions in the atm [GeV]
        Ec : cut-off value of the data [GeV]
        I0 : vertical muon flux [cm^2.s.sr]^-1
        n : power index
        eps : parameter that modifies the high energy part and should account for the finite lifetime of pi and K [GeV]
        return shukla differential flux [GeV.cm^2.s.sr]^-1 in  
        '''
        R = 6371*1e3 #m earth radius
        d = 1.5e4 #m  typical muon production altitude (?)
        D = lambda theta : np.sqrt(R**2/d**2 * np.cos(theta)**2 + 2*R/d + 1) - R/d * np.cos(theta)   #column density for inclined trajectory in curved Earth
        N = lambda n : (n-1)*(E0+Ec)**(n-1)
        I = I0*N(n)*(E0+energy)**-n * ( 1+ energy/eps)**-1 * D(theta)**(-(n-1)) 
        return I
    
    
    def ComputeDiffFlux(self, energy:float, theta:float, model:str="guan", sys_unc_df=None, altitude_correction:bool=False ):
        """
        energy in GeV
        theta in rad
        """
        #print("model=", model)
        self.diffFluxAtSeaLevel,  self.diffFluxAtSeaLevel_std = None, None
        
        if model == "gaisser_std": self.diffFluxAtSeaLevel = self.gaisser_std(energy, theta)
        elif model == "gaisser": self.diffFluxAtSeaLevel = self.gaisser(energy, theta)
        elif model == "guan": self.diffFluxAtSeaLevel = self.guan(energy, theta)
        elif model == "tang": self.diffFluxAtSeaLevel = self.tang(energy, theta)
        elif model == "shukla": 
            ###TEST param at sea level (theta=0°)
            if theta == 0 :
                I0 = 70.7 *1e-4 #(cm^2.s.sr)^-1
                n = 3.0
                E0 = 4.29 #GeV
                Ec = 0.5 #GeV
                eps = 854 #GeV
            elif theta == 75*np.pi/180 : 
                #print(("here shukla 75"))
                I0 = 65.2 *1e-4 #(cm^2.s.sr)^-1
                n = 3.0
                E0 = 23.78 #GeV
                Ec = 1.0 #GeV
                eps = 2e3 #GeV
            else: raise Exception("Shukla parameters are only known for theta=0 or 75°")
            self.diffFluxAtSeaLevel = self.shukla(energy=energy,
                                                    theta=theta,
                                                    E0=E0, Ec=Ec, I0=I0,
                                                    n=n, eps=eps)
        elif model == "corsika_soufriere": 
            if self.is_corsika:  
                self.diffFluxAtSeaLevel = self.DF_corsika(energy, theta)
                self.diffFluxAtSeaLevel_std = self.DF_corsika_std(energy, theta)
                
        else : raise ValueError(f"Unknown flux model, please choose among: {self.available}")
        c = 1
        mu = self.par["mu"]*1e-3 #in GeV
        #lorentz_gamma = energy/(mu*c**2)
        #beta  = np.sqrt(1 - 1/lorentz_gamma**2)
        p = np.sqrt(energy**2 + mu**2)#lorentz_gamma*beta*mu
        self.diffFluxAtTelescopeLevel = self.diffFluxAtSeaLevel 
        if  self.diffFluxAtSeaLevel_std is not None:  self.diffFluxAtTelescopeLevel_std = self.diffFluxAtSeaLevel_std
        #if (self.altitude < 1000) & (p > 10)  & (model != "tang"): 
        if altitude_correction & (model != "corsika_soufriere"):
            #if p > 10 : 
            ##Lesparre GJI 2010
            #altitude correction by Hebbeker & Timmermans (2002)
            h0 = (4900 + 750 * p)   
            print(p, energy, np.exp(self.altitude/h0))
            self.diffFluxAtTelescopeLevel *= np.exp(self.altitude/h0)
        if sys_unc_df is not None: self.diffFluxAtTelescopeLevel *= sys_unc_df
        return self.diffFluxAtTelescopeLevel
    
    def ComputeOpenSkyFlux(self, emin:float, emax:float, theta:float, model:str="gaisser",  sys_unc_df=None, **kwargs):

        flux = np.zeros(shape=(theta.shape[0])) 
        for i in range(theta.shape[0]):
            res, err = quad(self.ComputeDiffFlux, emin, emax, args=(theta[i], model, sys_unc_df), **kwargs)
            flux[i] = res
        return flux
    
    def get_transmitted_flux(self, emin:float, theta:float, model:str="gaisser", sys_unc_df=None, emax=np.inf, **kwargs):
        """
        Integrate flux upon energy range [emin, emax] 
        emin (float) in GeV : minimal crossing energy at given opacity
        theta (float) in rad : incident vertical angle (vertical axis: theta=0, horizon: theta=pi/2)
        out: tuple (result, error) integration
        """
        res, err = quad(self.ComputeDiffFlux, emin, emax, args=(theta, model, sys_unc_df), **kwargs)
        return res, err
    
    def ComputeTransmittedFlux(self, emin:np.ndarray, theta:np.ndarray, model:str="gaisser", sys_unc_df=None, emax=np.inf, **kwargs):
        """
        Integrate flux upon energy range [emin, emax] 
        emin (array2d) in GeV (lower bound of integration =~ cut-off energy)
        theta (array2d) in rad
        out: array2d of shape (theta.shape[0], emin.shape[1]) 
        """
        outFlux = np.zeros(shape=(theta.shape[0], emin.shape[0])) #emin.shape[1]
        for i, t in enumerate(theta):
            for j, e in enumerate(emin):
                ##compute integral flux
                try : 
                    res, _ = self.get_transmitted_flux( emin=e, emax=emax, theta=t, model=model, sys_unc_df=sys_unc_df, **kwargs) 
                    outFlux[i,j]= res
                except : 
                    print("min, emax, theta = ", e, emax, t*180/np.pi)
        return outFlux
    
    
    def __call__(self, energy, theta, sys_unc:float=0.15, model:str="gaisser"):
        df = self.ComputeDiffFlux(energy, theta, model)
        unc_df = sys_unc*df
        return df, unc_df


    
    def number_expected_muons(self, dt:float, acceptance:np.ndarray, rho:float, thickness:np.ndarray, azimuth:np.ndarray, theta:np.ndarray, func_flux, *args, **kwargs ):
        """
        Inputs: 
        - dt (run duration [s])
        - acceptance  [cm^2.sr] : array(2Nxy-1, 2Nxy-1)
        - rho (density medium [g/cm^3])
        - thickness (apparent obstacle thickness [m]) : array(2Nxy-1, 2Nxy-1)
        - azimuth [rad] : array(2Nxy-1, 2Nxy-1)
        - theta [rad] : array(2Nxy-1, 2Nxy-1)
        - integrated flux func(theta, opacity) [1/(cm^2.sr.s)] computed from a given model (e.g Gaisser) : array(2Nxy-1, 2Nxy-1)
        Returns:
        - number of expected muons in each telescope pixel: array(2Nxy-1, 2Nxy-1)
        """
        expected_opacity = (thickness*1e2 * rho) *1e-2 #g/cm2 -> hg/cm2=mwe
        #integrand   = lambda theta, phi, opacity : func_flux(theta,opacity)[0] * np.sin(theta)
        #azimuth, theta = azimuth, theta*np.pi/180
        #phi_bounds  = np.array([ np.concatenate(([[phi[0]-dphi[0]/2, phi[0]+dphi[0]/2]], [[phi[i]-dphi[i-1]/2, phi[i]+dphi[i]/2] if i!=len(phi)-1 else [phi[i]-dphi[i-1]/2, phi[i]+dphi[i-1]/2]  for i in range(1,len(phi))]) ) for phi, dphi in zip( azimuth,np.diff(azimuth) )   ]) #
        #theta_bounds = np.array([ np.concatenate(([[theta[0]-dtheta[0]/2, theta[0]+dtheta[0]/2]], [[theta[i]-dtheta[i-1]/2, theta[i]+dtheta[i]/2] if i!=len(theta)-1 else [theta[i]-dtheta[i-1]/2, theta[i]+dtheta[i-1]/2]  for i in range(1,len(theta))]) ) for theta, dtheta in zip( theta.T,np.diff(theta.T) )   ]) #*np.pi/180
        #func_IntFlux_pixel = lambda i,j :  dblquad( integrand, theta_bounds[i,j,0]*np.pi/180, theta_bounds[i,j,1]*np.pi/180, phi_bounds[i,j,0]*np.pi/180, phi_bounds[i,j,1]*np.pi/180, args=(expected_opacity[i,j],) )[0]
        #func_Nexp = lambda i,j : dt * acceptance[i,j] * func_flux(theta[i,j],expected_opacity[i,j]) #func_IntFlux_pixel(i,j)
        
        Nexp,flux_exp = np.zeros_like(acceptance), np.zeros_like(acceptance)
        for i in range(acceptance.shape[0]): 
            for j in range(acceptance.shape[1]): 
                if np.isnan(thickness[i,j]): 
                    flux_exp[i,j], Nexp[i,j] = np.nan, np.nan
                    continue 
                flux_exp[i,j] = func_flux(theta[i,j],expected_opacity[i,j])
                Nexp[i,j] = dt * acceptance[i,j] * flux_exp[i,j]
        return flux_exp, Nexp



if __name__=="__main__":
    pass
   