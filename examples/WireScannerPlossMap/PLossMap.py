"""
This example shows rescale the power loss from simulation for 
300-nm NEG coated absorber in Wire Scanner by
fitting impedance mode by mode and compute the power loss.

It plots the power loss map using package pyvista.

@date: 09/05/2025
@author: J. Li, E. de la Fuente
"""

import sys
sys.path.append(r'../../') # if bihc is not pip installed
import bihc

import iddefix
import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv


'''Impedance comparison between origin and fitting'''
'''Fitting resonator parameters are calculated by IDDEFIX'''
Z_origin = bihc.Impedance()
Z_origin.getImpedanceFromCST(r'.\Zz.txt')
fmax = np.max(Z_origin.f)

fresamp_n = 10000
frequencies = np.linspace(0, fmax, fresamp_n, endpoint=True)
pars = {
    1: [1.07e+01, 12.24, 1.497e+09],
    2: [3.06e+01, 6.67, 1.920e+09],
    3: [1.10e+01, 6.59, 2.683e+09],
    4: [4.13e+00, 8.44, 3.297e+09],
    5: [3.99e+00, 66.55, 3.797e+09],
    6: [4.20e+00, 10.03, 4.059e+09],
    }
Z_fitting = iddefix.Impedances.n_Resonator_longitudinal_imp(frequencies, pars)

'''plot setting'''
fontsize = 30
fig_size = (10,5.5)

'''plot Zr'''
plt.figure(figsize=fig_size)
plt.plot(Z_origin.f/1e9, np.real(Z_origin.Zr), lw=5, c='black', label='Origin', alpha=0.9)
plt.plot(frequencies/1e9, Z_fitting.real, lw=3, c='red', label='Fitting', alpha=0.9)

ax = plt.gca()
ax_color = 'black'
ax.spines['bottom'].set_linewidth(2)
ax.spines['top'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
plt.tick_params(width=2)
plt.legend(loc='best',fontsize=fontsize-5)
plt.xlabel('frequency [GHz]', fontsize=fontsize)
plt.ylabel('$\\Re(Z_{\\rm L})$ $\\rm [\\Omega]$', fontsize=fontsize, color=ax_color)
plt.xticks(fontsize=fontsize);plt.yticks(fontsize=fontsize, color=ax_color)
ax.yaxis.get_offset_text().set(size=fontsize-2)
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.tight_layout()
plt.grid(linestyle = '--', linewidth=2)
plt.show()


'''Create beam object'''
bunch_length = 1e-9
N_density = 2.3e11
q_value = 0.6
shift = 20E6

from bihc.fillingschemes import fillingSchemeLHC_standard
fillingScheme = fillingSchemeLHC_standard(ninj=10)

'''Computing power of every mode'''
def calModePower(frequency, par, index:int):
    
    print(f'Mode {index} power is calculating...')
    
    Resonator_par = pars[index]
    Z = bihc.Impedance(f=Z_origin.f) # specify the frequency when create the object
    Z.getResonatorImpedance(Rs=Resonator_par[0], Qr=Resonator_par[1], fr=Resonator_par[2]) # input the resonator parameters from IDDEFIX
    fmax = np.max(Z.f)
    
    beam_qgauss = bihc.Beam(bunchLength=bunch_length, Np=N_density, 
                            bunchShape='q-GAUSSIAN', qvalue=q_value,
                            fillingScheme=fillingScheme, machine='LHC', spectrum='numeric', fmax=fmax, 
                            verbose=False)

    # get power
    Z_p_scan = beam_qgauss.getShiftedPloss(Z, shift=shift)[1]
    # Z_p_mini = np.min(Z_p_scan)
    Z_p_maxi = np.max(Z_p_scan)
    # Z_p_aver = np.average(Z_p_scan)
    
    # return Z_p_mini, Z_p_maxi, Z_p_aver
    return Z_p_maxi


mode_number = 6
ModePLoss_Max = np.zeros(mode_number)
for i in np.arange(mode_number):
    mode_index = i+1
    ModePLoss_Max[i] = calModePower(frequency=Z_origin.f, par=pars, index=mode_index)




'''load CST result of loss'''
data = np.loadtxt('./Loss_in_Metals.txt', skiprows=2)
MetalLoss = data[:,1]
DiLoss = np.loadtxt('./Loss_in_Dielectrics.txt', skiprows=2)[:,1]
PLossTotal = MetalLoss+DiLoss


'''Rescal the power loss map exported from CST'''
# load CST result of ploss map
path = './CST_PLossMap/'
Pmonitors = sorted(glob.glob(path+'mode_*'))

P = None
for i in range(mode_number):

    print(f'Monitor {i+1}')

    # x [mm] y [mm] z [mm]  P [W/m^3]
    d = pd.read_csv(Pmonitors[i], skiprows=[0,1], sep=r'\s+', header=None, low_memory=False)
    d[d=='-nan(ind)'] = 0
    x, y, z, PowerMap = np.array(d[0]), np.array(d[1]), np.array(d[2]), np.array(d[3],dtype='float')

    if P is None:
        P = np.zeros_like(PowerMap)

    P += PowerMap*ModePLoss_Max[i]/PLossTotal[i] #[W/m3]


# Using pyvista polydata point cloud
surf = pv.read('./All.stl')
surf_ab = pv.read('./absorber.stl')

# using voxelize to cteate a voxel model with more mesh cells
# It is very important for those model with flat surfaces
surf_ab = pv.voxelize(surf_ab, density=surf_ab.length / 200, check_surface=False)


pdata = pv.PolyData(np.vstack((x,y,z)).transpose())
pdata['Ploss Density Map [W/m^3]'] = P
Radius = 8
fieldonsurf = surf_ab.interpolate(pdata, radius=Radius, sharpness=2,)


pl = pv.Plotter()
pl.add_mesh(fieldonsurf, cmap='jet', opacity=1.0,
    # show_edges=True,
    interpolate_before_map=True,
    )
pl.add_mesh(surf, opacity=0.1, color='white', show_edges=False)
pl.add_axes()
pl.enable_3_lights()
pl.camera_position = 'zy'
pl.camera.azimuth += 60
pl.camera.elevation += 10
pl.show(full_screen=False)


