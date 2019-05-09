# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

tur_eps = 1000. # turbulent flow epsilon
kB = 1.381e-23
flu_T = 130 + 273 # K
flu_niu = 2.5e-5 # fluid niu
flu_rho = 0.88
flu_P = 101325.
flu_miu = flu_niu * flu_rho
flu_moleculediam = 364. * 1e-12 # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
flu_lamda = kB * flu_T / (np.sqrt(2.) * np.pi * flu_moleculediam**2 * flu_P) # https://en.wikipedia.org/wiki/Mean_free_path
par_rho = 2250. # kg/m3

def setMatPlotLib():
    global fontsize_multiple, mpl_tick_size, figsize_dflt
    #matplotlib.rcParams['font.family'] = 'SimHei'
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    matplotlib.rcParams['mathtext.default'] = 'regular'
    #matplotlib.rcParams['font.family'] = 'SimHei' #'Times New Roman' 
    fontsize_multiple = 1.0
    matplotlib.rcParams['font.size'] = int(38 * fontsize_multiple)
    mpl_tick_size = int(34 * fontsize_multiple)
    matplotlib.rcParams['axes.labelsize'] = mpl_tick_size + 3
    matplotlib.rcParams['axes.titlesize'] = mpl_tick_size + 0
    matplotlib.rcParams['xtick.labelsize'] = mpl_tick_size
    matplotlib.rcParams['ytick.labelsize'] = mpl_tick_size
    matplotlib.rcParams['xtick.major.pad'] = 12
    matplotlib.rcParams['xtick.major.size'] = 12
    matplotlib.rcParams['xtick.minor.size'] = 6
    matplotlib.rcParams['ytick.major.pad'] = 12
    matplotlib.rcParams['ytick.major.size'] = matplotlib.rcParams['xtick.major.size'] 
    matplotlib.rcParams['ytick.minor.size'] = matplotlib.rcParams['ytick.minor.size']
    matplotlib.rcParams['legend.fontsize'] = mpl_tick_size-4 #'large' #'medium' #mpl_tick_size #'large'
    matplotlib.rcParams['legend.markerscale'] = 1.0 #mpl_tick_size #'large'
    matplotlib.rcParams['lines.markersize'] = 8
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    #matplotlib.rcParams['text.usetex'] = True
    figsize_dflt = (12,9) #(12, 9)

def drawBeta(beta_B, beta_S, beta_I, dp, ipar):
    plt.figure(figsize=figsize_dflt)
    plt.plot(dp*1e6, beta_B, 'r-', markersize = 1)
#    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
#    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('particle diamter, micron')
    plt.ylabel('beta (collision rate, m3/s)')
    plt.plot(dp*1e6, beta_S, 'b-', markersize = 1)
    plt.plot(dp*1e6, beta_I, 'g-', markersize = 1)
    plt.plot(dp*1e6, beta_B + beta_S + beta_I, 'm--', markersize = 1)
    plt.legend(['Brownian', 'Turbulent shearing', 'Turbulent inertia', 'Total']) #, fontsize=tick_size)
    foo = 1e9
    plt.title('Collision rate for ' + str(int(dp[ipar]*foo)/(foo/1e6)) + ' micron flyash particle at epislon = ' \
          + str(tur_eps) + ', 403K') 
    return


def calcBeta(tur_eps, dp, dp_tar, ldraw=True): # tur_eps:# turbulent flow epsilon, dp: particle diameter in m
    if dp is None:
        dp_logRng = np.linspace(-2,2,num=int(4/0.02))
        dp = np.power(np.ones_like(dp_logRng)*10, dp_logRng) * 1e-6 # in m
    else:
        dp = np.asarray(dp)        
    dp_N = dp.shape[0]
    assert dp_N >=2, 'dp_N must be >=2'
    assert dp_tar <= np.max(dp) and dp_tar >= np.min(dp), 'dp_tar is OUT of range!'
    #par_lamda = kB * flu_T / (np.sqrt(2.) * dp**2 * flu_P)
    par_Kn = flu_lamda / (dp/2)
    #par_Cs = np.ones_like(dp) # set to 1 for now, later will use Kn correction
    par_Cs = 1 + par_Kn * (1.27 + 0.4 * np.exp(-1.1 / par_Kn))
    par_m = np.pi/6 * (dp**3) * par_rho
    par_cbar = np.sqrt(8*kB*flu_T/(np.pi*par_m))
    par_D = par_Cs * kB * flu_T /(3 * np.pi * flu_miu * dp)
    par_l = 8 * par_D / (np.pi * par_cbar) # mean free path of particles
    par_g = ((dp + par_l)**3 - (dp**2 + par_l**2)**1.5) / (3 * dp *par_l) - dp
    #
    idx_big = np.argmax(dp>=dp_tar)
    if idx_big >= 1:
        ipar = idx_big if dp[idx_big]-dp_tar <= dp_tar - dp[idx_big-1] else idx_big - 1
    else:
        ipar = 0  
#    print 'ipar: ', ipar
    # --- Brownian diffusion ---
    beta_B = 2 * np.pi * (par_D[ipar] + par_D) * (dp[ipar] + dp) / \
     ((dp[ipar] + dp)/(dp[ipar] + dp + 2*np.sqrt(par_g[ipar]**2 + par_g**2)) \
      + 8*(par_D[ipar]+par_D)/np.sqrt(par_cbar[ipar]**2+par_cbar**2)/(dp[ipar]+dp))
    #--- Turbulent shearing ---
    beta_S = 1.3 * np.sqrt(tur_eps / flu_niu) * (dp[ipar]/2 + dp/2)**3
    beta_I = 5.7 * (dp[ipar]/2 + dp/2)**2 * np.abs((par_D[ipar]*par_m[ipar] - par_D*par_m)/kB/flu_T) \
    *(tur_eps**3/flu_niu)**0.25
    # --- Inception and impaction due induced by 'macroscopic' drift of bigger particles
    # considering particle ipar has a 'macroscopic' drift velocity relative to local fluid flow
    vel_drift = 2.0 # m/s
    deltau= vel_drift
    St_ipar = par_rho * deltau * dp**2. /(18. * flu_miu * dp[ipar])
    yita_p = (St_ipar/(St_ipar+0.65))**3.7
    Yc = np.sqrt(yita_p) * dp[ipar] / 2.
    beta_drift = 3.142 * Yc**2 * vel_drift
    beta_drift[dp > dp[ipar]*0.9] = 0. # only consider those particles that are smaller than dp[ipar]
    #    
#    beta_tot = beta_B + beta_S + beta_I + beta_drift
    if ldraw: drawBeta(beta_drift, beta_S, beta_I, dp, ipar)
    return beta_B, beta_S, beta_I, beta_drift


#%%
if __name__ == '__main__':
    setMatPlotLib()
    mode = 0
    if mode == 0:
        dp = np.asarray([2.5e-6, 20e-6], dtype=np.float)
        tur_eps_all = 10**(np.linspace(np.log10(10), np.log10(1000), 20))
        beta_B = np.asarray([], dtype=np.float)
        beta_S = np.asarray([], dtype=np.float)
        beta_I = np.asarray([], dtype=np.float)
        beta_drift = np.asarray([], dtype=np.float)
        for tur_eps in tur_eps_all:
            beta_B_i, beta_S_i, beta_I_i, beta_drift_i = calcBeta(tur_eps=tur_eps, dp=dp, dp_tar=dp[1], ldraw=False)
            beta_B = np.append(beta_B, beta_B_i[0])
            beta_S = np.append(beta_S, beta_S_i[0])
            beta_I = np.append(beta_I, beta_I_i[0])
            beta_drift = np.append(beta_drift, beta_drift_i[0])
        plt.figure(figsize=figsize_dflt)
        plt.plot(tur_eps_all, beta_B, 'r-', marker='s')
        plt.xticks() #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
        plt.yticks()
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\it\epsilon$ / $m^2$/$s^3$')#,fontsize=label_size)
        plt.ylabel(r'$\it\beta$ / $m^3$/s')#,fontsize=label_size)
        plt.plot(tur_eps_all, beta_S, 'b-', marker='o')
        plt.plot(tur_eps_all, beta_I, 'g-', marker='v')
        plt.plot(tur_eps_all, beta_drift, 'm-',marker='^')#, markersize = marker_size)
        plt.plot(tur_eps_all, beta_B + beta_S + beta_I + beta_drift, 'm--')#, markersize = marker_size)
        plt.legend(['Brownian', 'Turbulent shearing', 'Turbulent inertia', 'Drift', 'Total'], framealpha=0.5,
               loc=4) # 3/4 is left/right bottom
        foo = 1e9
#        title('Collision rate at variant turbulent energy dissipation rate', fontsize=label_size)
        plt.tight_layout()
    elif mode == 1:
#        kB = 1.381e-23
#        flu_T = 403 # K
#        flu_niu = 2.5e-5 # fluid niu
#        flu_rho = 0.88
#        flu_P = 101325.
#        flu_miu = flu_niu * flu_rho
#        flu_moleculediam = 364. * 1e-12 # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
#        flu_lamda = kB * flu_T / (np.sqrt(2.) * np.pi * flu_moleculediam**2 * flu_P) # https://en.wikipedia.org/wiki/Mean_free_path
        # tur_eps = 100. # turbulent flow epsilon
#        par_rho = 2200. # kg/m3
        dp_logRng = np.linspace(-2,2,num=int(4/0.02))
        dp = np.power(np.ones_like(dp_logRng)*10, dp_logRng) * 1e-6 # in m
        dp_N = dp.shape[0]
        #par_lamda = kB * flu_T / (np.sqrt(2.) * dp**2 * flu_P)
        par_Kn = flu_lamda / (dp/2)
        #par_Cs = np.ones_like(dp) # set to 1 for now, later will use Kn correction
        par_Cs = 1 + par_Kn * (1.27 + 0.4 * np.exp(-1.1 / par_Kn))
        par_m = np.pi/6 * (dp**3) * par_rho
        par_cbar = np.sqrt(8*kB*flu_T/(np.pi*par_m))
        par_D = par_Cs * kB * flu_T /(3 * np.pi * flu_miu * dp)
        par_l = 8 * par_D / (np.pi * par_cbar) # mean free path of particles
        par_g = ((dp + par_l)**3 - (dp**2 + par_l**2)**1.5) / (3 * dp *par_l) - dp
        ipar = 149 #149 for 9.88miu, 76 for 0.337miu, 111 for 1.70miu, 119 for 2.47 miu, 115 for 2miu, 50 for 0.1miu, 100 for 1.0 miu
        # --- Brownian diffusion ---
        beta_B = 2 * np.pi * (par_D[ipar] + par_D) * (dp[ipar] + dp) / \
         ((dp[ipar] + dp)/(dp[ipar] + dp + 2*np.sqrt(par_g[ipar]**2 + par_g**2)) \
          + 8*(par_D[ipar]+par_D)/np.sqrt(par_cbar[ipar]**2+par_cbar**2)/(dp[ipar]+dp))
        #%% --- Turbulent shearing ---
        beta_S = 1.3 * np.sqrt(tur_eps / flu_niu) * (dp[ipar]/2 + dp/2)**3
        beta_I = 5.7 * (dp[ipar]/2 + dp/2)**2 * np.abs((par_D[ipar]*par_m[ipar] - par_D*par_m)/kB/flu_T) \
        *(tur_eps**3/flu_niu)**0.25
        beta_tot = beta_B + beta_S + beta_I
        # considering particle ipar has a 'macroscopic' drift velocity relative to local fluid flow
        vel_drift = 2.0 # m/s
        deltau= vel_drift
        St_ipar = par_rho * deltau * dp**2. /(18. * flu_miu * dp[ipar])
        yita_p = (St_ipar/(St_ipar+0.65))**3.7
        Yc = np.sqrt(yita_p) * dp[ipar] / 2.
        beta_drift = 3.142 * Yc**2 * vel_drift
        beta_drift[dp > dp[ipar]*0.9] = 0.
        #%%
        drawBeta(beta_drift, beta_S, beta_I, dp, ipar)

