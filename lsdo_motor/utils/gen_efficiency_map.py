import numpy as np 
import csdl
from python_csdl_backend import Simulator
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

from lsdo_motor.utils.efficiency_map_models import ExtractUpperLimitTorqueModel, EfficiencyMapAnalysisModel

def gen_efficiency_map(L, D, V_lim=400, I_rated=115, plot=True, flux_weakening=True,
                       p=6, m=3, Z=36, RPM_limit=15000, num_RPM_step=30, num_torque_step=80, ax_plot=None):
    '''
    This function generates the data and a figure for the efficiency map of a PMSM.
    Inputs:
        - L: motor length (m)
        - D: motor diameter (of the inner stator) (m)
        - V_lim: Voltage limit (V); default is 400
        - I_rated: rated current of the machine (A); default is 123
        - plot: boolean flag for if the user wants to directly generate a plot
        - flux_weakening: boolean flag to (de)activate flux weakening

        - p: pole pairs (default is 6)
        - m: number of phases (default is 3)
        - Z: number of slots (default is 36)
        - RPM_limit: upper bound on MOTOR RPM of efficiency map (default is 15000)
        - RPM_step: step size in discretization of RPM (default is 30)
        - torque_step: step size in discretization of torque (default is 50)
    '''
    fit_coeff_dep_H = np.array([1.12832651, 1., 1., 1.])
    fit_coeff_dep_B = np.array([2.29955130e-03, 1.03633806e+01, -3.16820279e+00])

    RPM_range = np.linspace(500,RPM_limit, num_RPM_step) # NOTE: MAY NEED TO CHANGE LOWER BOUND
    omega_range = RPM_range * 2.*np.pi/60. * p
    num_active_nodes = len(omega_range)
    # region Evaluation of upper limit torque curve
    torque_lim_model = ExtractUpperLimitTorqueModel(
        component_name='torque_lim',
        pole_pairs=p,
        phases=m,
        num_slots=Z,
        V_lim=V_lim,
        rated_current=I_rated,
        fit_coeff_dep_H=fit_coeff_dep_H,
        fit_coeff_dep_B=fit_coeff_dep_B,
        num_active_nodes=num_active_nodes,
        flux_weakening=flux_weakening
    )
    sim = Simulator(torque_lim_model)
    sim['motor_diameter'] = D
    sim['motor_length'] = L
    if flux_weakening:
        sim['omega'] = omega_range
    sim.run()

    if flux_weakening:
        torque_voltage_limit = sim['T_lim']
    else:
        torque_voltage_limit = sim['T_upper_lim_curve']
    upper_torque_curve = sim['T_upper_lim_curve']
    R = sim['Rdc']
    L_d = sim['L_d']
    L_q = sim['L_q']
    PsiF = sim['PsiF']
    B_delta = sim['B_delta']
    motor_variables = sim['motor_parameters']
    # endregion
    
    # if debug:
    #     fig = plt.figure(1)
    #     plt.plot(omega_range*60/p/2/np.pi, upper_torque_curve, 'k', linewidth=5)
    #     plt.xlabel('RPM')
    #     plt.ylabel('EM Torque')
    #     plt.show()
    #     exit()

    torque_grid_step = num_torque_step # 80 TO ADJUST NUMBER OF STEPS BETWEEN LOWER (NONZERO) AND UPPER LIMIT TORQUE
    omega_grid = np.resize(omega_range, (torque_grid_step, len(omega_range)))
    torque_grid = np.linspace(10, upper_torque_curve, torque_grid_step)
    torque_voltage_limit_grid = np.resize(torque_voltage_limit, (torque_grid_step, len(omega_range)))

    omega_grid_vector = np.reshape(omega_grid, (torque_grid_step * len(omega_range), ))
    torque_grid_vector = np.reshape(torque_grid, (torque_grid_step * len(omega_range), ))
    torque_voltage_limit_grid_vector = np.reshape(torque_voltage_limit_grid, (torque_grid_step * len(omega_range), ))

    # region generating data for efficiency map
    efficiency_map_analysis_model = EfficiencyMapAnalysisModel(
        pole_pairs=p,
        phases=m,
        num_slots=Z,
        V_lim=V_lim,
        rated_current=I_rated,
        num_active_nodes=len(omega_grid_vector),
        flux_weakening=flux_weakening
    )
    sim_eff_map = Simulator(efficiency_map_analysis_model)

    sim_eff_map['omega'] = omega_grid_vector
    sim_eff_map['T_em'] = torque_grid_vector
    sim_eff_map['Rdc'] = R
    sim_eff_map['L_d'] = L_d
    sim_eff_map['L_q'] = L_q 
    sim_eff_map['PsiF'] = PsiF 
    sim_eff_map['T_lim'] = torque_voltage_limit_grid_vector
    sim_eff_map['motor_diameter'] = D
    sim_eff_map['B_delta'] = B_delta
    sim_eff_map['motor_parameters'] = motor_variables
    sim_eff_map.run()

    efficiency_active = sim_eff_map['efficiency_active']
    # endregion

    # region Plotting
    efficiency_grid = np.zeros((torque_grid_step + 1,len(omega_range) + 1))
    efficiency_grid[1:, 1:] = np.reshape(efficiency_active, (torque_grid_step, len(omega_range)))
    omega_grid_plot = np.zeros_like(efficiency_grid)
    omega_grid_plot[1:, 1:] = np.reshape(omega_grid, (torque_grid_step, len(omega_range)))
    omega_grid_plot[0, 1:] = omega_range
    torque_grid_plot = np.zeros_like(efficiency_grid)
    torque_grid_plot[1:, 1:] = np.reshape(torque_grid, (torque_grid_step, len(omega_range)))
    torque_grid_plot[1:, 0] = torque_grid_plot[1:,1]

    levels_f = np.linspace(0,1,31)
    levels = np.array([30, 80, 88, 90, 92, 93, 94, 95, 96, 97, 98])*0.01
    # levels = np.array([0.8, 0.9, 0.93, 0.96, 0.97, 0.975, 0.977, 0.9784])
    
    if ax_plot:
        ax_plot.contourf(omega_grid_plot*60/2/p/np.pi, torque_grid_plot, efficiency_grid, cmap='RdGy_r', levels=levels_f)
        # plt.contourf(omega_grid_plot*60/2/np.pi, torque_grid_plot, efficiency_grid, cmap='jet', levels=levels_f)
        # plt.co
        # lorbar()
        # ax_plot.clim(0,1)
        contours = ax_plot.contour(omega_grid_plot*60/2/p/np.pi, torque_grid_plot, efficiency_grid, colors='black', levels=levels)
        # contours = plt.contour(omega_grid_plot*60/2/np.pi, torque_grid_plot, efficiency_grid, colors='black', levels=levels)
        ax_plot.clabel(contours, inline=True, fontsize=12)
        ax_plot.plot(omega_range*60/2/p/np.pi, upper_torque_curve, 'k', linewidth=3)
        ax_plot.set_xlabel('Speed [RPM]', fontsize=12)
        ax_plot.set_ylabel('Torque [Nm]', fontsize=12)
        # ax_plot.set_xticks(ticks=np.arange(0,15001,3000), fontsize=12)
        # ax_plot.set_yticks(fontsize=12)
        ax_plot.grid(False)
        

    else:
        sns.set_style("ticks")
        plt.figure(2)
        plt.contourf(omega_grid_plot*60/2/p/np.pi, torque_grid_plot, efficiency_grid, cmap='RdGy_r', levels=levels_f)
        # plt.contourf(omega_grid_plot*60/2/np.pi, torque_grid_plot, efficiency_grid, cmap='jet', levels=levels_f)
        # plt.colorbar()
        plt.clim(0,1)
        contours = plt.contour(omega_grid_plot*60/2/p/np.pi, torque_grid_plot, efficiency_grid, colors='black', levels=levels)
        # contours = plt.contour(omega_grid_plot*60/2/np.pi, torque_grid_plot, efficiency_grid, colors='black', levels=levels)
        plt.clabel(contours, inline=True, fontsize=12)
        plt.plot(omega_range*60/2/p/np.pi, upper_torque_curve, 'k', linewidth=3)
        plt.xlabel('Speed [RPM]', fontsize=12)
        plt.ylabel('Torque [Nm]', fontsize=12)
        plt.xticks(ticks=np.arange(0,15001,3000), fontsize=12)
        plt.yticks(fontsize=12)
        plt.grid(False)
        plt.show()
    # endregion

    return sim_eff_map, ax_plot

if __name__ == '__main__':
    sim = gen_efficiency_map(L=0.086, D=0.182, V_lim=400, num_RPM_step=20, num_torque_step=20, flux_weakening=True)
