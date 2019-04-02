# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
#matplotlib.rcParams['font.family'] = 'SimHei'
matplotlib.rcParams['font.family'] = 'Times New Roman'

figsize_dflt = (8, 6)
label_size = 18
tick_size = 16
tur_eps = 10. # turbulent flow epsilon

def drawBeta(beta_B, beta_S, beta_I, dp, ipar):
    figure(figsize=figsize_dflt)
    plot(dp*1e6, beta_B, 'r-', markersize = 1)
    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('particle diamter, micron',fontsize=label_size)
    plt.ylabel('beta (collision rate, m3/s)',fontsize=label_size)
    plot(dp*1e6, beta_S, 'b-', markersize = 1)
    plot(dp*1e6, beta_I, 'g-', markersize = 1)
    plot(dp*1e6, beta_B + beta_S + beta_I, 'm--', markersize = 1)
    legend(['Brownian', 'Turbulent shearing', 'Turbulent inertia', 'Total'], fontsize=tick_size)
    foo = 1e9
    title('Collision rate for ' + str(int(dp[ipar]*foo)/(foo/1e6)) + ' micron flyash particle at epislon = ' \
          + str(tur_eps) + ', 403K',fontsize=label_size)
    return
#%%
if __name__ == '__main__':
    kB = 1.381e-23
    flu_T = 403 # K
    flu_niu = 2.5e-5 # fluid niu
    flu_rho = 0.88
    flu_P = 101325.
    flu_miu = flu_niu * flu_rho
    flu_moleculediam = 364. * 1e-12 # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
    flu_lamda = kB * flu_T / (np.sqrt(2.) * np.pi * flu_moleculediam**2 * flu_P) # https://en.wikipedia.org/wiki/Mean_free_path
#    tur_eps = 100. # turbulent flow epsilon
    par_rho = 2200. # kg/m3
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
    ipar = 76 #119 for 2.47 miu, 115 for 2miu, 50 for 0.1miu, 100 for 1.0 miu
    # --- Brownian diffusion ---
    beta_B = 2 * np.pi * (par_D[ipar] + par_D) * (dp[ipar] + dp) / \
     ((dp[ipar] + dp)/(dp[ipar] + dp + 2*np.sqrt(par_g[ipar]**2 + par_g**2)) \
      + 8*(par_D[ipar]+par_D)/np.sqrt(par_cbar[ipar]**2+par_cbar**2)/(dp[ipar]+dp))
    #%% --- Turbulent shearing ---
    beta_S = 1.3 * np.sqrt(tur_eps / flu_niu) * (dp[ipar]/2 + dp/2)**3
    beta_I = 5.7 * (dp[ipar]/2 + dp/2)**2 * np.abs((par_D[ipar]*par_m[ipar] - par_D*par_m)/kB/flu_T) \
    *(tur_eps**3/flu_niu)**0.25
    beta_tot = beta_B + beta_S + beta_I
    #%%
    drawBeta(beta_B, beta_S, beta_I, dp, ipar)

