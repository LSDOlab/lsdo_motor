import numpy as np
import csdl
from python_csdl_backend import Simulator

from lsdo_motor.motor_submodels.TC1_magnet_mec_model import MagnetMECModel
from lsdo_motor.permeability.mu_fitting import permeability_fitting

# SETUP PERMEABILITY FITTING
file_name = 'Magnetic_alloy_silicon_core_iron_C.tab'
order=10

mu_fitting = permeability_fitting(
    file_name=file_name,
    test=False,
    order=order
)

fit_coeff_dep_H = mu_fitting[0]
fit_coeff_dep_B = mu_fitting[1]

# CSDL MODEL RUN
m = MagnetMECModel(
    fit_coeff_dep_H=fit_coeff_dep_H,
    fit_coeff_dep_B=fit_coeff_dep_B,
)
rep = csdl.GraphRepresentation(m)

sim = Simulator(rep)
# MANUALLY UPDATE VARIABLES BASED ON ZEYU'S MATLAB CODE
sim['tooth_pitch'] = 0.0325
sim['tooth_width'] = 0.0171
sim['slot_height'] = 0.0239
sim['alpha_i'] = 0.8104
sim['pole_pitch'] = 0.0975
sim['l_ef'] = 0.2755
sim['height_yoke_stator'] = 0.0226
sim['L_j1'] = 0.058
sim['air_gap_depth'] = 0.0013
sim['K_theta'] = 1.0835
sim['A_f2'] = 0.0011
sim['bm'] = 0.0751
sim['phi_r'] = 0.0232
sim['lambda_m'] = 3.4210e-06

print(sim['B_delta'])
#sim.visualize_implementation()
# ONES NEEDED IN FUTURE COMPUTATIONS
print('B_delta: ', sim['B_delta'])
print('phi_air: ', sim['phi_air'])
print('F_total: ', sim['F_total'])
print('F_delta: ', sim['F_delta'])

# ONES NEEDED JUST FOR CHECKING
print('H_y: ', sim['H_y'])
print('phi_f: ', sim['phi_f'])
print('phi_s: ', sim['phi_s'])
print('phi_mag: ', sim['phi_mag'])

