# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 10:56:07 2018

@author: gliu

reference: Direct simulation Monte Carlo method for particle coagulation and aggregation.Kruis 2000.pdf
"""
import numpy as npy
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from pbm_MagnitudeAnalysis import drawBeta
import copy
from scipy import interpolate
matplotlib.rcParams['font.family'] = 'Times New Roman'


figsize_dflt = (8, 6)
label_size = 18
tick_size = 16

kB = 1.381e-23
flu_T = 403 # K
flu_niu = 2.5e-5 # fluid niu
flu_rho = 0.88
flu_P = 101325.
flu_miu = flu_niu * flu_rho
flu_moleculediam = 364. * 1e-12 # using kinetic diameter, https://en.wikipedia.org/wiki/Kinetic_diameter
flu_lamda = kB * flu_T / (npy.sqrt(2.) * npy.pi * flu_moleculediam**2 * flu_P) # https://en.wikipedia.org/wiki/Mean_free_path
tur_eps = 10. # turbulent flow epsilon
par_rho = 2200. # kg/m3

def drawPSD(dps, np, figttl=''):
    figure(figsize=figsize_dflt)
    plot(dps, np, 'b-', markersize = 5)
    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('particle diamter, micron',fontsize=label_size)
#    plt.ylabel('particle number density(?)',fontsize=label_size)
#    legend(['Brownian', 'Turbulent shearing', 'Turbulent inertia', 'Total'], fontsize=tick_size)
#    foo = 1e7
    title(figttl, fontsize=label_size)

def calcBeta(dpin=None, linmiu=False): # linmiu=True means input dp is in miu m
    if dpin is None:
        dp_logRng = npy.linspace(-2,2,num=int(4/0.02))
        dp = npy.power(npy.ones_like(dp_logRng)*10, dp_logRng) * 1e-6 # in m
    else:
        dp = copy.deepcopy(dpin)        
    if linmiu:
        dp = dp * 1e-6 # in m
    dp_N = dp.shape[0]
    # following dp should be in unit of 'm' not 'miu m'
    #par_lamda = kB * flu_T / (npy.sqrt(2.) * dp**2 * flu_P)
    par_Kn = flu_lamda / (dp/2)
    #par_Cs = npy.ones_like(dp) # set to 1 for now, later will use Kn correction
    par_Cs = 1 + par_Kn * (1.27 + 0.4 * npy.exp(-1.1 / par_Kn))
    par_m = npy.pi/6 * (dp**3) * par_rho
    par_cbar = npy.sqrt(8*kB*flu_T/(npy.pi*par_m))
    par_D = par_Cs * kB * flu_T /(3 * npy.pi * flu_miu * dp)
    par_l = 8 * par_D / (npy.pi * par_cbar) # mean free path of particles
    par_g = ((dp + par_l)**3 - (dp**2 + par_l**2)**1.5) / (3 * dp *par_l) - dp
#    ipar = 100 #115 for 2miu, 50 for 0.1miu, 100 for 1.0 miu
    betaAll = npy.zeros((dp_N, dp_N))
    for ipar in range(dp_N):
        # --- Brownian diffusion ---
        beta_B = 2 * npy.pi * (par_D[ipar] + par_D) * (dp[ipar] + dp) / \
         ((dp[ipar] + dp)/(dp[ipar] + dp + 2*npy.sqrt(par_g[ipar]**2 + par_g**2)) \
          + 8*(par_D[ipar]+par_D)/npy.sqrt(par_cbar[ipar]**2+par_cbar**2)/(dp[ipar]+dp))
        #% --- Turbulent shearing ---
        beta_S = 1.3 * npy.sqrt(tur_eps / flu_niu) * (dp[ipar]/2 + dp/2)**3
        beta_I = 5.7 * (dp[ipar]/2 + dp/2)**2 * npy.abs((par_D[ipar]*par_m[ipar] - par_D*par_m)/kB/flu_T) \
        *(tur_eps**3/flu_niu)**0.25
        betaAll[ipar,:] = beta_B + beta_S + beta_I
#    drawBeta(beta_B, beta_S, beta_I, dp, ipar)
    return betaAll
    
#%%
def getnp(dps=[1.0], linmiu=True):
    if not linmiu:
        dps = dps * 1e6 # convert to in miu m
    Np = npy.array([5e14, 1e11, 1e9])
    dpg = npy.array([0.08, 2.0, 10.0]) # miu m
    sgmpg = npy.array([1.5, 2.0, 1.5])
    ln_sgmpg = npy.log(sgmpg)
    np = [npy.sum(Np / npy.sqrt(2*npy.pi)/ln_sgmpg * npy.exp(-(npy.log(dp/dpg))**2/2/ln_sgmpg**2) \
                  / dp) for dp in dps]
    delta_dp = dps[1:] - dps[:-1]
    np_bin = delta_dp * np[:-1]
    return np, np_bin

def getInitialPSD(dp=None):
    if dp is None:
        #dp=10**(npy.arange(-1, 2, 0.05))
        dp = 10**npy.linspace(npy.log10(0.08), npy.log10(20), 200) # miu m
    np, np_bin = getnp(dps=dp)
    #drawPSD(dp, np, 'np ~ dp')
    log_dp = npy.log10(dp*1e-6)
    par_m = npy.pi / 6 * (dp*1e-6)**3 * par_rho * 1e6 # mg
#    delta_dp=dp[1:]-dp[:-1]
#    np_bin = delta_dp*np[:-1]
    Np0 = npy.sum(np_bin)
    par_m_bin = par_m[:-1] * np_bin
    print 'Np0: ', Np0
    print 'total mass: ', npy.sum(par_m_bin)
    #
    #drawPSD(dp[:-1], np_bin/Np0, 'np_bin/Np0 ~ np')
    #delta_m = par_m[1:]- par_m[:-1]
    delta_log_dp = log_dp[1:]- log_dp[:-1]
    drawPSD(dp[:-1]*1e-6, par_m_bin/delta_log_dp, 'dM/dlog(dp) ~ dp')

#%%
laddOneBig = True # add One group of big particles etc 50um
mass_OB = 5e3 # mg
dpOB = 50e-6 # m
dp_norm0 = 1e-6 # in m
v_norm0 = npy.pi/6 * dp_norm0**3
#    
v_min = 0.08**3 #0.01 miu
v_max = 20**3 #100 miu
v_ratio = 1.08
Nv = int(npy.log10(v_max/v_min) / npy.log10(v_ratio))
ns = npy.arange(Nv+3)
Nv = len(ns) -1
vs = v_min * v_ratio ** ns
vs_tmp = (vs[:-1] + vs[1:])/2
vs_h = (vs[:-1] + vs[1:])/2
vs_l = vs_h[:-1].copy()
vs_l = npy.concatenate(([vs[0]],vs_l))
vs = vs[:-1]
vs_bounds = npy.column_stack((vs, vs_l, vs_h))
dps = npy.power((vs * v_norm0 * 6 / npy.pi), 1./3.)
B = calcBeta(dpin = dps, linmiu=False)
if npy.max(npy.abs(B-B.T)) > 1e-23: # B must be symmetric (B.T must = B)
    print 'potential error!'
log_dp = npy.log10(dps)
np_, Np_bin = getnp(dps=dps, linmiu=False)
if laddOneBig: npy.append(dps, dpOB)
par_m = npy.pi / 6 * dps**3 * par_rho * 1e6 # mg
if laddOneBig:
    Np_bin = npy.append(Np_bin, mass_OB/par_m[-1])
    par_m_bin = par_m * Np_bin
else:
    par_m_bin = par_m[:-1] * Np_bin
delta_log_dp = log_dp[1:]- log_dp[:-1]
#drawPSD(dps[:-1], par_m_bin/delta_log_dp, 'dM/dlog(dp) ~ dp')
#%%
lEvolv = True
if lEvolv:
    if not laddOneBig: Np_bin = npy.append(Np_bin, 0)
    dt = 0.001
    max_iters = 1000 # 不能过长，使得计算终点时刻的v_max尺寸的颗粒团数目仍然非常少，否则tot_count的统计值会不准，因为超出了粒径范围
#    par_counts = npy.zeros_like(vs)
    par_counts = Np_bin
    tot_counts = npy.zeros((max_iters+1,1))
    tot_counts[0] = npy.sum(par_counts)
    parN_t = npy.zeros((max_iters+1, Nv+1)) # first column is time
    parN_t[0, :] = npy.insert(par_counts, 0, 0.)
    for i_t in range(max_iters):
        t = (i_t+1) * dt
        dndt = npy.zeros_like(vs)
        for i in range(Nv):        
            ni = par_counts[i]
            vi = vs[i]
            for j in range(i, Nv):
                nj = par_counts[j]
                vj = vs[j]
                beta_ij = B[i,j] #A * (vi/v_min + vj/v_min) #((i+1) + (j+1))
                vij = vi + vj
                k = npy.argmax(vs_h > vij)
#                if i == Nv -1 and j == Nv -1 :
#                    print 'k:',k
#                    print vs_h
#                    print 'vij:',vij
                if k == 0:
#                    print i_t, i,j
                    k = Nv-1
                if i != j:
                    dndt[k] = dndt[k] + vij / vs[k] * ni * nj * beta_ij
                    dndt[i] = dndt[i] - ni * nj * beta_ij
                    dndt[j] = dndt[j] - ni * nj * beta_ij
                else:
                    dndt[k] = dndt[k] + vij / vs[k] * (0.5 * ni * nj * beta_ij)
                    dndt[i] = dndt[i] - (0.5 * ni * nj * beta_ij) * 2
        d_counts = dndt * dt
        par_counts = par_counts + d_counts
        parN_t[i_t+1, :] = npy.insert(par_counts, 0, t) 
        tot_counts[i_t+1] = npy.sum(par_counts)            
    ts = npy.arange(max_iters + 1) * dt
    primary_par_counts = vs / v_min
    count_percent = par_counts / tot_counts[-1]
    parM_t = copy.deepcopy(parN_t)
    parM_t[:,1:] = parN_t[:,1:] * par_m
    tot_mass = npy.sum(parM_t[:,1:], 1)
    yita = 100*(1 - parM_t[-1,1:-2]/parM_t[0,1:-2]) # percentage of removal
    rmvl = npy.column_stack((dps[:-2], yita))
    print 'N/N0 at end of simulation given by PBM: ', tot_counts[-1]/tot_counts[0]     
    print 'Ratio of second biggest particle counts between end and start of simulation: ', \
    parN_t[-1,-2]/parN_t[0,-2]
    #%%
    figure(figsize=figsize_dflt)
    plot(ts, tot_counts, 'b-', markersize = 5)
    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
    plt.yticks(fontsize=tick_size)
#    plt.xscale('')
#    plt.yscale('log')
    plt.xlabel('time (s)',fontsize=label_size)
    plt.ylabel('particle total number', fontsize=label_size)
    title('Evolution of particle total counts with time', fontsize=label_size)
    #%%
    figure(figsize=figsize_dflt)
    plot(dps[0:-1]*1e6, parM_t[0,1:-1], 'b-', markersize = 5)
    t1_ = 0.5
    it1_ = npy.argmax(ts>=t1_)
    plot(dps[0:-1]*1e6, parM_t[it1_, 1:-1], 'm-', markersize = 5)
    plot(dps[0:-1]*1e6, parM_t[-1, 1:-1], 'r-', markersize = 5)
    legend(['t=0 s','t=' + str(ts[it1_]) + 's', 't=' + str(ts[-1]) + 's'], fontsize=tick_size)
    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('dp (micron)',fontsize=label_size)
    plt.ylabel('dM/dlog(dp), mg/m3', fontsize=label_size)
    title('dM/dlog(dp) ~ dp', fontsize=label_size)
    #%%
    figure(figsize=figsize_dflt)
    plot(dps[:-2]*1e6, yita, 'xr-', markersize = 5)
    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
#    plt.yscale('log')
    plt.xlabel('dp (micron)',fontsize=label_size)
    plt.ylabel("Percentage of 'removal', %", fontsize=label_size)
    title("Percentage of 'removal' ~ dp", fontsize=label_size)
    #%% interpolate to certain dp points
    dp_tar = 10**(npy.arange(-1, 1, 0.01) - 6)
    f_tar = interpolate.interp1d(dps[:-2], yita, kind='slinear')
    yita_tar = f_tar(dp_tar)
    rmvl_tar = npy.column_stack((dp_tar, yita_tar))
    t2_ = 0.5
    it2_ = npy.argmax(ts>=t2_)
    yita2 = 100*(1 - parM_t[it2_,1:-2]/parM_t[0,1:-2]) 
    f_tar2 = interpolate.interp1d(dps[:-2], yita2, kind='slinear')
    yita_tar2 = f_tar2(dp_tar)
    rmvl_tar2 = npy.column_stack((dp_tar, yita_tar2))
    figure(figsize=figsize_dflt)
    plot(dp_tar*1e6, yita_tar2, 'm-', markersize = 5)
    plot(dp_tar*1e6, yita_tar, 'r-', markersize = 5)
    legend(['t=' + str(ts[it2_]) + 's', 't=' + str(ts[-1]) + 's'], fontsize=tick_size)
    plt.xticks(fontsize=tick_size) #, rotation=90) # plot.tick_params(axis='both', which='major', labelsize=10)
    plt.yticks(fontsize=tick_size)
    plt.xscale('log')
#    plt.yscale('log')
    plt.xlabel('dp (micron)',fontsize=label_size)
    plt.ylabel("Percentage of 'removal', %", fontsize=label_size)
    title("Percentage of 'removal' ~ dp", fontsize=label_size)