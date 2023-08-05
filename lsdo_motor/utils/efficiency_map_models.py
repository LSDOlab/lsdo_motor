import numpy as np 
import csdl
from python_csdl_backend import Simulator

from lsdo_motor.core.TC1_motor_sizing_model import TC1MotorSizingModel
from lsdo_motor.core.motor_submodels.TC1_magnet_mec_model import MagnetMECModel
from lsdo_motor.core.motor_submodels.TC1_inductance_mec_model import InductanceModel
from lsdo_motor.core.motor_submodels.TC1_torque_limit_model import TorqueLimitModel
from lsdo_motor.core.motor_submodels.TC1_flux_weakening_model import FluxWeakeningBracketModel
from lsdo_motor.core.motor_submodels.TC1_implicit_em_torque_model import EMTorqueModel

class ExtractUpperLimitTorqueModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('component_name')
        self.parameters.declare('pole_pairs') # 6
        self.parameters.declare('phases') # 3
        self.parameters.declare('num_slots') # 36
        self.parameters.declare('V_lim')
        self.parameters.declare('rated_current')
        self.parameters.declare('fit_coeff_dep_H') # FITTING COEFFICIENTS (X = H, B = f(H))
        self.parameters.declare('fit_coeff_dep_B') # FITTING COEFFICIENTS (X = B, H = g(B))
        self.parameters.declare('num_active_nodes')
        self.parameters.declare('flux_weakening')

        self.motor_variable_names = [
            'outer_stator_radius', 'pole_pitch', 'tooth_pitch', 'air_gap_depth', 'l_ef',
            'rotor_radius', 'turns_per_phase', 'Acu',  'tooth_width', 'height_yoke_stator',
            'slot_bottom_width', 'slot_height', 'slot_width_inner', 'Tau_y', 'L_j1', 'Kdp1',
            'bm', 'Am_r', 'phi_r', 'lambda_m', 'alpha_i', 'Kf', 'K_phi', 'K_theta', 'A_f2',
        ]

    def define(self):
        m = self.parameters['phases']
        p = self.parameters['pole_pairs']
        Z = self.parameters['num_slots']
        V_lim = self.parameters['V_lim']
        rated_current = self.parameters['rated_current']
        fit_coeff_dep_H = self.parameters['fit_coeff_dep_H']
        fit_coeff_dep_B = self.parameters['fit_coeff_dep_B']
        num_active_nodes = self.parameters['num_active_nodes']

        component_name=self.parameters['component_name']

        flux_weakening=self.parameters['flux_weakening']

        D_i = self.declare_variable('motor_diameter') # inner radius of stator
        L = self.declare_variable('motor_length') # effective length of motor

        self.add(
            TC1MotorSizingModel(
                pole_pairs=p,
                phases=m,
                num_slots=Z,
                rated_current=rated_current,
                component_name=component_name
            ),
            'TC1_motor_sizing_model',
        )
        # outputs: resistance, mass, motor_variables
        motor_variables = self.declare_variable('motor_parameters', shape=(27,)) # array of motor sizing outputs
        for i in range(motor_variables.shape[0]-2):
            self.register_output(self.motor_variable_names[i], motor_variables[i])
        self.declare_variable('T_em_max')
        Rdc = self.declare_variable('Rdc') # DC resistance

        # gear_ratio = 4.
        # omega_rotor = self.declare_variable('rpm', shape=num_active_nodes)
        # omega = self.register_output('omega', omega_rotor * gear_ratio * 2*np.pi/60)'
        omega = self.declare_variable('omega', shape=(num_active_nodes))

        self.add(
            MagnetMECModel(
                fit_coeff_dep_H=fit_coeff_dep_H,
                fit_coeff_dep_B=fit_coeff_dep_B,
            ),
            'magnet_MEC_model',
        )

        self.add(
            InductanceModel(
                pole_pairs=p,
                phases=m,
                num_slots=Z,
                rated_current=rated_current,
                fit_coeff_dep_H=fit_coeff_dep_H,
                fit_coeff_dep_B=fit_coeff_dep_B,
            ),
            'inductance_MEC_model',
        )
        L_d = self.declare_variable('L_d')
        L_q = self.declare_variable('L_q')
        phi_air = self.declare_variable('phi_air')
        W_1 = self.declare_variable('turns_per_phase')
        PsiF = self.register_output('PsiF', W_1 * phi_air)

        R_expanded = self.register_output('R_expanded', csdl.expand(Rdc, (num_active_nodes,)))
        L_d_expanded = self.register_output('L_d_expanded', csdl.expand(L_d, (num_active_nodes,)))
        L_q_expanded = self.register_output('L_q_expanded', csdl.expand(L_q, (num_active_nodes,)))
        PsiF_expanded = self.register_output('PsiF_expanded', csdl.expand(PsiF, (num_active_nodes,)))

        T_em_max = self.declare_variable('T_em_max')

        if flux_weakening:
            self.add(
                    TorqueLimitModel(
                        pole_pairs=p,
                        V_lim=V_lim,
                        num_nodes=num_active_nodes
                    ),
                    'torque_limit_model'
            )
            # OUTPUT OF torque_limit_model GIVES THE UPPER LIMIT CURVE BASED ON FLUX WEAKENING
            T_lim = self.declare_variable('T_lim', shape=(num_active_nodes,))

            T_upper_lim_curve = self.register_output(
                    'T_upper_lim_curve',
                    # -csdl.log(csdl.exp(-T_lim) + csdl.exp(-csdl.expand(T_em_max, (num_active_nodes,))))
                    csdl.min(csdl.reshape(T_lim, new_shape=(num_active_nodes, )),
                                csdl.expand(T_em_max, (num_active_nodes, )))
                )
        else:
            T_upper_lim_curve = self.register_output(
                    'T_upper_lim_curve',
                    csdl.expand(T_em_max, (num_active_nodes, ))
                )

class EfficiencyMapAnalysisModel(csdl.Model):
    '''
    INPUTS TO THIS MODEL:
        - tau_max (maximum torque)
        - rotor diameter
        - base speed of motor
    '''

    def initialize(self):
        self.parameters.declare('pole_pairs') # 6
        self.parameters.declare('phases') # 3
        self.parameters.declare('num_slots') # 36
        self.parameters.declare('V_lim')
        self.parameters.declare('rated_current')
        self.parameters.declare('num_active_nodes')
        self.parameters.declare('flux_weakening')

        self.motor_variable_names = [
            'outer_stator_radius', 'pole_pitch', 'tooth_pitch', 'air_gap_depth', 'l_ef',
            'rotor_radius', 'turns_per_phase', 'Acu',  'tooth_width', 'height_yoke_stator',
            'slot_bottom_width', 'slot_height', 'slot_width_inner', 'Tau_y', 'L_j1', 'Kdp1',
            'bm', 'Am_r', 'phi_r', 'lambda_m', 'alpha_i', 'Kf', 'K_phi', 'K_theta', 'A_f2',
            'Rdc', 'T_em_max' # LAST TWO ARE SPECIAL VARIABLES
        ]

    def define(self):

        m = self.parameters['phases']
        p = self.parameters['pole_pairs']
        Z = self.parameters['num_slots']
        V_lim = self.parameters['V_lim']
        rated_current = self.parameters['rated_current']
        num_active_nodes = self.parameters['num_active_nodes']
        flux_weakening = self.parameters['flux_weakening']

        D_i = self.declare_variable('D_i') # inner radius of stator

        Rdc = self.declare_variable('Rdc')
        L_d = self.declare_variable('L_d')
        L_q = self.declare_variable('L_q')
        PsiF = self.declare_variable('PsiF')

        omega = self.declare_variable('omega', shape=(num_active_nodes,))
        T_em = self.declare_variable('T_em', shape=(num_active_nodes,))

        T_lim = self.declare_variable('T_lim', shape=(num_active_nodes,))
        T_lower_lim = self.declare_variable('T_lower_lim', val=0., shape=(num_active_nodes,))

        R_expanded = self.register_output('R_expanded', csdl.expand(Rdc, (num_active_nodes,)))
        L_d_expanded = self.register_output('L_d_expanded', csdl.expand(L_d, (num_active_nodes,)))
        L_q_expanded = self.register_output('L_q_expanded', csdl.expand(L_q, (num_active_nodes,)))
        PsiF_expanded = self.register_output('PsiF_expanded', csdl.expand(PsiF, (num_active_nodes,)))

        D = (3*p*(L_d_expanded-L_q_expanded))
        if flux_weakening:
            # FLUX WEAKENING BRACKETING IMPLICIT MODEL (TO FIND BRACKET LIMIT)
            a_bracket = self.register_output(
                'a_bracket',
                D**2 * ((omega*L_q_expanded)**2 + R_expanded**2)
            )
            c_bracket = self.register_output(
                'c_bracket',
                (3*p*PsiF_expanded)**2 * (R_expanded**2 + (omega*L_q_expanded)**2) + \
                12*p*omega*R_expanded*T_lim*(L_d_expanded-L_q_expanded)**2 - (V_lim*D)**2
            )
            d_bracket = self.register_output(
                'd_bracket',
                -12*p*PsiF_expanded*T_lim*(R_expanded**2 + omega**2*L_d_expanded*L_q_expanded)
            )
            e_bracket = self.register_output(
                'e_bracket',
                4*T_lim**2*(R_expanded**2 + (omega*L_d_expanded)**2)
            )

            self.add(
                FluxWeakeningBracketModel(
                    pole_pairs=p,
                    num_nodes=num_active_nodes
                ),
                'flux_weakening_bracket_method'
            )

            Iq_fw_bracket = self.declare_variable('Iq_fw_bracket', shape=(num_active_nodes, ))
            Id_fw_bracket = self.declare_variable('Id_fw_bracket', shape=(num_active_nodes, ))

        self.add(
            EMTorqueModel(
                pole_pairs=p,
                V_lim=V_lim,
                num_nodes=num_active_nodes,
                rated_current=rated_current,
                phases=m,
                motor_variable_names=self.motor_variable_names,
                mode='efficiency_map',
                flux_weakening=flux_weakening
            ),
            'efficiency_map_model'
        )