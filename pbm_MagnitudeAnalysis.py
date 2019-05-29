# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#tur_eps = 1000. # turbulent flow epsilon
#kB = 1.381e-23
#flu_T = 130 + 273 # K
#flu_niu = 2.5e-5 # fluid niu
#flu_rho = 0.88
#flu_P = 101325.
#flu_miu = flu_niu * flu_rho
#flu_moleculediam = 364. * 1e-12 # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
#flu_lamda = kB * flu_T / (np.sqrt(2.) * np.pi * flu_moleculediam**2 * flu_P) # https://en.wikipedia.org/wiki/Mean_free_path
#par_rho = 2250. # kg/m3

def fluegas_rho(t=100.): # t is in oC
    fg_rho100C = 0.95; # 100℃下标准烟气的密度，《传热学》（第四版）560页
    rho = (100 + 273.)/(t + 273.) * fg_rho100C;
    return rho

def fluegas_niu(t=100.): # t is in ℃
    niu = 0.1194324e-4 + 0.892572e-7 * t + 0.830565e-10 * t**2 - 0.8619116e-14 * t**3 + 0.132968e-17 * t**4 # 标准烟气的运动粘度，红宝书，niu随温度的变化其实不能忽略
    return niu

def setParams(unittestcode=-1):
    params = {}
    if unittestcode == -1: # NOT a unittest if == -1
        params = dict(
            tur_eps = 745.7, # turbulent flow epsilon
            kB = 1.381e-23,
            flu_T = 130. + 273., # K
            flu_niu = 2.5e-5, # fluid niu
            flu_rho = 0.88,
            flu_P = 101325.,
            flu_moleculediam = 364. * 1e-12, # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
            par_rho = 2250., # kg/m3
            rel_vel = 0.,
        )
        params['flu_miu'] = params['flu_niu'] * params['flu_rho']
        params['flu_lamda'] = params['kB'] * params['flu_T'] / (np.sqrt(2.) * np.pi * params['flu_moleculediam']**2 * params['flu_P']) # https://en.wikipedia.org/wiki/Mean_free_path        
#        
    elif unittestcode == 0: # 'standard' test, do not modify
        flu_t = 55. # oC
        params = dict(
            tur_eps = 0.59, # turbulent flow epsilon
            kB = 1.381e-23,
            flu_T = flu_t + 273., # K
            flu_niu = fluegas_niu(t=flu_t), #2.5e-5, # fluid niu
            flu_rho = fluegas_rho(t=flu_t), #0.88,
            flu_P = 101325.,
            flu_moleculediam = 364. * 1e-12, # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
            par_rho = 2250., # kg/m3
            rel_vel = 0.,
            )
        params['flu_miu'] = params['flu_niu'] * params['flu_rho']
        params['flu_lamda'] = params['kB'] * params['flu_T'] / (np.sqrt(2.) * np.pi * params['flu_moleculediam']**2 * params['flu_P']) # https://en.wikipedia.org/wiki/Mean_free_path
#        
    elif unittestcode == 1: # not prepared yet
        pass
    return params
        
def setMatPlotLib():
    global fontsize_multiple, mpl_tick_size, figsize_dflt, total_markers
    matplotlib.rcParams['font.family'] = 'Times New Roman' # 'SimHei'
    matplotlib.rcParams['mathtext.default'] = 'regular'
    fontsize_multiple = 1.0
    matplotlib.rcParams['font.size'] = int(38 * fontsize_multiple)
    mpl_tick_size = int(34 * fontsize_multiple)
    matplotlib.rcParams['axes.labelsize'] = mpl_tick_size + 3
    matplotlib.rcParams['axes.titlesize'] = mpl_tick_size + 0
    matplotlib.rcParams['axes.titlepad'] = 20
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
    matplotlib.rcParams['lines.linewidth'] = 2
    matplotlib.rcParams['lines.markersize'] = 8
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
#    matplotlib.rcParams['text.usetex'] = True
    figsize_dflt = (12,9) #(12, 9)
    total_markers = 15
    return figsize_dflt, total_markers

def drawBeta(beta_B, beta_TS, beta_TI, beta_IC, beta_IM, beta_Total, dp, dp_tar):
    setMatPlotLib()
    plt.figure(figsize=figsize_dflt)
#    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
#    plt.yticks(fontsize=tick_size)
    markevery = int(len(beta_Total) / total_markers)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\it{d}_p$ ($\mu$m)')
    plt.ylabel(r'$\it{\beta}$ ($m^3$/s)')
    plt.plot(dp*1e6, beta_B, 'k-', marker='x', markevery=markevery)
    plt.plot(dp*1e6, beta_TS, 'b-', marker='v', markevery=markevery)
    plt.plot(dp*1e6, beta_TI, 'g-', marker='^', markevery=markevery)
    plt.plot(dp*1e6, beta_IC, 'c-', marker='o', markevery=markevery)
    plt.plot(dp*1e6, beta_IM, 'm-', marker='s', markevery=markevery)
    plt.plot(dp*1e6, beta_Total, 'r--', linewidth=1.5*matplotlib.rcParams['lines.linewidth'])
    plt.legend(['Brownian', 'Turbulent shearing', 'Turbulent inertia', 'Inertial interception', 'Inertial impaction' ,'Total'], framealpha=0.4) #, fontsize=tick_size)
    foo = 1e10
    plt.title('Collision rate for ' + str(int(dp_tar * foo)/(foo/1e6)) + ' micron flyash particle' ) #\
#          ' at epislon = ' + str(tur_eps) + ', 403K') 
    return


#def calcBeta(tur_eps, dp=[], dp_tars=[], rel_vel=[], dp_icAndim_ignore={}, ldraw=False): # tur_eps:# turbulent flow epsilon, dp: particle diameter in m
def calcBeta(params, dp=[], dp_tars=[], rel_vel=[], dp_icAndim_ignore={}, ldraw=False): # tur_eps:# turbulent flow epsilon, dp: particle diameter in m
    # this function calculate beta of kernel function values of "source" particles (dp) on "target particles" (dp_tar)
    # return: 1D (one-dimensional) numpy array, if dp_tar is a single value, array size is equal to dp
    #         2D numpy array (sort of a matrix), if dp_tar is a 1D numpy array or a 1D list, array size is len(dp_tar) rows X len(dp) columns
#    print "params['flu_T']: ", params['flu_T']
    if dp is None: # a simple test/case
        dp_logRng = np.linspace(-2,2,num=int(4/0.02))
        dp = np.power(np.ones_like(dp_logRng)*10, dp_logRng) * 1e-6 # in m
    else:
        dp = np.asarray(dp)        
    dp_N = dp.shape[0]
    assert dp_N >=2, 'dp_N must be >=2'
    dp_tars = np.asarray(dp_tars) 
    #par_lamda = kB * flu_T / (np.sqrt(2.) * dp**2 * flu_P)
    par_Kn = params['flu_lamda'] / (dp/2)
    #par_Cs = np.ones_like(dp) # set to 1 for now, later will use Kn correction
    par_Cs = 1 + par_Kn * (1.27 + 0.4 * np.exp(-1.1 / par_Kn))
    par_m = np.pi/6 * (dp**3) * params['par_rho']
    par_cbar = np.sqrt(8 * params['kB'] * params['flu_T'] / (np.pi * par_m))
    par_D = par_Cs * params['kB'] * params['flu_T'] /(3 * np.pi * params['flu_miu'] * dp)
    par_l = 8 * par_D / (np.pi * par_cbar) # mean free path of particles
    par_g = ((dp + par_l)**3 - (dp**2 + par_l**2)**1.5) / (3 * dp *par_l) - dp
    #    
    beta_B = np.zeros((len(dp_tars), len(dp)))
    beta_TS = np.zeros((len(dp_tars), len(dp)))
    beta_TI = np.zeros((len(dp_tars), len(dp)))
    beta_IC = np.zeros((len(dp_tars), len(dp)))
    beta_IM = np.zeros((len(dp_tars), len(dp)))
    beta_Total = np.zeros((len(dp_tars), len(dp)))
    # 
    if isinstance(params['rel_vel'], float) or isinstance(params['rel_vel'], int) :
        rel_vel = float(params['rel_vel']) * np.ones((len(dp_tars), len(dp))) # temporarily #############################
    #
    for i_tar, dp_tar in enumerate(dp_tars):
        assert dp_tar <= np.max(dp) and dp_tar >= np.min(dp), 'dp_tar is OUT of range!'
        idx_big = np.argmax(dp>=dp_tar)
        if idx_big >= 1:
            ipar = idx_big if dp[idx_big]-dp_tar <= dp_tar - dp[idx_big-1] else idx_big - 1
        else:
            ipar = 0  
#        print 'ipar: {}, dp(ipar):{}'.format(ipar, dp[ipar])
# --- Brownian diffusion ---
        beta_B[i_tar, :] = 2 * np.pi * (par_D[ipar] + par_D) * (dp[ipar] + dp) / \
                            ((dp[ipar] + dp)/(dp[ipar] + dp + 2*np.sqrt(par_g[ipar]**2 + par_g**2)) + \
                             8*(par_D[ipar]+par_D)/np.sqrt(par_cbar[ipar]**2+par_cbar**2)/(dp[ipar]+dp))
        # Brownian should be: 20190528 'text snapshot'. For future visual inspection and verification
                                    # --- Brownian diffusion ---
#        beta_B[i_tar, :] = 2 * np.pi * (par_D[ipar] + par_D) * (dp[ipar] + dp) / \
#                            ((dp[ipar] + dp)/(dp[ipar] + dp + 2*np.sqrt(par_g[ipar]**2 + par_g**2)) + \
#                             8*(par_D[ipar]+par_D)/np.sqrt(par_cbar[ipar]**2+par_cbar**2)/(dp[ipar]+dp))
# --- Turbulent shearing and intertial ---
        beta_TS[i_tar, :] = 1.3 * np.sqrt(params['tur_eps'] / params['flu_niu']) * (dp[ipar]/2. + dp/2.)**3
        beta_TI[i_tar, :] = 5.7 * (dp[ipar]/2. + dp/2.)**2 * np.abs((par_D[ipar]*par_m[ipar] - par_D*par_m)/params['kB']/params['flu_T']) * \
                            (params['tur_eps']**3/params['flu_niu'])**0.25
#                            
#       Turbulent shearing and intertial should be: 20190528 'text snapshot'. For future visual inspection and verification
        # --- Turbulent shearing and intertial ---
#        beta_TS[i_tar, :] = 1.3 * np.sqrt(params['tur_eps'] / params['flu_niu']) * (dp[ipar]/2. + dp/2.)**3
#        beta_TI[i_tar, :] = 5.7 * (dp[ipar]/2. + dp/2.)**2 * np.abs((par_D[ipar]*par_m[ipar] - par_D*par_m)/params['kB']/params['flu_T']) * \
#                            (params['tur_eps']**3/params['flu_niu'])**0.25
#
# --- InCeption and IMpaction capture induced by 'macroscopic' drift of bigger particles
        if dp[ipar] <= dp_icAndim_ignore['dp_tar_upbound']: # only consider dp[ipar] <= e.g. 1.67um(pm2.5)
            # IC refers Weber and Paddock, 1983, Interceptional and gravitational collision efficiencies for single collectors at intermediate Reynolds numbers.
            Reb = rel_vel[i_tar,:] * dp / params['flu_niu']
            yita_ic = 3./2. * np.power(dp[ipar]/dp, 2) * (1. + (3./16. * Reb) / (1 + 0.249 * np.power(Reb,0.56)))        
            beta_IC[i_tar, :] = yita_ic * np.pi/4 * np.power(dp, 2) * rel_vel[i_tar,:]
            beta_IC[i_tar, dp <= dp_icAndim_ignore['dp_src_lowbound']] = 0.  # only consider 'drift' caputure with particles bigger than a specified size
            # IM refers Ho and Sommerfeld 2002, modelling of micro-particle agglomeration in turbulent flows.
            St_ipar = params['par_rho'] * rel_vel[i_tar,:] * dp[ipar]**2. /(18. * params['flu_miu'] * dp)
            yita_p = (St_ipar/(St_ipar+0.65))**3.7
            Yc = np.sqrt(yita_p) * dp / 2.
            beta_IM[i_tar, :] = np.pi * Yc**2 * rel_vel[i_tar,:]
            beta_IM[i_tar, dp <= dp_icAndim_ignore['dp_src_lowbound']] = 0. # only consider 'drift' caputure with particles bigger than a specified size
        else:
            beta_IC[i_tar, :] = 0.
            beta_IM[i_tar, :] = 0.
        #    
    beta_Total = beta_B + beta_TS + beta_TI + beta_IC + beta_IM # beta_IM (if being a square matrix) is Non-Symmetric
    if ldraw and len(dp_tars) == 1:
        drawBeta(beta_B[0], beta_TS[0], beta_TI[0], beta_IC[0], beta_IM[0], beta_Total[0], dp, dp_tars[0])
    return beta_B, beta_TS, beta_TI, beta_IC, beta_IM, beta_Total


#%%
if __name__ == '__main__':
    setMatPlotLib()
    params = setParams()
    mode = 2
    if mode == 0: # in mode 0, we investigate variation of beta with >>epsilon<<, len(dp_tars) must be ONE
        dp = np.asarray([2.5e-6, 20e-6], dtype=np.float)
        tur_eps_all = 10**(np.linspace(np.log10(10), np.log10(1000), 20))
        beta_B = np.asarray([], dtype=np.float)
        beta_TS = np.asarray([], dtype=np.float)
        beta_TI = np.asarray([], dtype=np.float)
        beta_IC = np.asarray([], dtype=np.float)
        beta_IM = np.asarray([], dtype=np.float)
        i_tardp = 0 # on which we are investigating variation of beta with epsilon
        i_srcdp = 1 # i_tardp and i_srcdp can be the same. they are refering to 'dp'
        for tur_eps in tur_eps_all:
            params['tur_eps'] = tur_eps
            beta_B_i, beta_TS_i, beta_TI_i, beta_IC_i, beta_IM_i, beta_Total_ = calcBeta(params, dp=dp, dp_tars=[dp[i_tardp]], \
                                                         dp_icAndim_ignore={'dp_tar_upbound':2.5e-6,'dp_src_lowbound':10e-6},ldraw=False)
            beta_B = np.append(beta_B, beta_B_i[0][i_srcdp])
            beta_TS = np.append(beta_TS, beta_TS_i[0][i_srcdp])
            beta_TI = np.append(beta_TI, beta_TI_i[0][i_srcdp])
            beta_IC = np.append(beta_IC, beta_IC_i[0][i_srcdp])
            beta_IM = np.append(beta_IM, beta_IM_i[0][i_srcdp])
        plt.figure(figsize=figsize_dflt)
        plt.xticks() #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
        plt.yticks()
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$\it\epsilon$ / $m^2$/$s^3$')#,fontsize=label_size)
        plt.ylabel(r'$\it\beta$ / $m^3$/s')#,fontsize=label_size)
        plt.plot(tur_eps_all, beta_B, 'k-', marker='s')
        plt.plot(tur_eps_all, beta_TS, 'b-', marker='o')
        plt.plot(tur_eps_all, beta_TI, 'g-', marker='v')
        plt.plot(tur_eps_all, beta_IC, 'c-',marker='x')#, markersize = marker_size)
        plt.plot(tur_eps_all, beta_IM, 'm-',marker='^')#, markersize = marker_size)
        plt.plot(tur_eps_all, beta_B + beta_TS + beta_TI + beta_IM, 'm--')#, markersize = marker_size)
        plt.legend(['Brownian', 'Turbulent shearing', 'Turbulent inertia','Drift interception', 'Drift impaction', 'Total'], framealpha=0.5, loc=4) # 3/4 is left/right bottom
#        foo = 1e9
#        title('Collision rate at variant turbulent energy dissipation rate', fontsize=label_size)
        plt.tight_layout()
    elif mode == 2: # in mode 2, we investigate variation of beta with >>src particle sizes<< at fixed epsilon and len(dp_tars) must be ONE
#        dp = np.asarray([2.5e-6, 10.0e-6, 20e-6, 35e-6, 50e-6], dtype=np.float)
        tar_dp_ = 1.67e-6
#        dp_logRng = np.linspace(np.log10(tar_dp_ * 1.0e6), 2, num=int(4/0.02))
        dp_logRng = np.linspace(np.log10(0.01), 2, num=int(600))
        dp = np.power(np.ones_like(dp_logRng)*10, dp_logRng) * 1e-6 # in m
        dp_tar = dp[np.argmax(dp>=tar_dp_)] # dp[0]
#        dp = dp[dp>=tar_dp_]        
        beta_B, beta_TS, beta_TI, beta_IC, beta_IM, beta_Total = calcBeta(params, dp=dp, dp_tars=[dp_tar], \
                                          dp_icAndim_ignore={'dp_tar_upbound':2.5e-6,'dp_src_lowbound':2.0e-6}, ldraw=False)
        drawBeta(beta_B[0], beta_TS[0], beta_TI[0], beta_IC[0], beta_IM[0], beta_Total[0], dp, dp_tar)
#        
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
        beta_TS = 1.3 * np.sqrt(tur_eps / flu_niu) * (dp[ipar]/2 + dp/2)**3
        beta_TI = 5.7 * (dp[ipar]/2 + dp/2)**2 * np.abs((par_D[ipar]*par_m[ipar] - par_D*par_m)/kB/flu_T) \
        *(tur_eps**3/flu_niu)**0.25
        beta_tot = beta_B + beta_TS + beta_TI
        # considering particle ipar has a 'macroscopic' drift velocity relative to local fluid flow
        vel_drift = 2.0 # m/s
        deltau= vel_drift
        St_ipar = par_rho * deltau * dp**2. /(18. * flu_miu * dp[ipar])
        yita_p = (St_ipar/(St_ipar+0.65))**3.7
        Yc = np.sqrt(yita_p) * dp[ipar] / 2.
        beta_IM = 3.142 * Yc**2 * vel_drift
        beta_IM[dp > dp[ipar]*0.9] = 0.
        #%%
        drawBeta(beta_IM, beta_TS, beta_TI, dp, ipar)

