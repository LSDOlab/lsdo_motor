import numpy as np
from csdl import Model, ScipyKrylov, NewtonSolver
from python_csdl_backend import Simulator
import csdl

from TC1_motor_model.motor_submodels.TC1_magnet_mec_model import MagnetMECModel
from TC1_motor_model.motor_submodels.TC1_inductance_mec_model import InductanceModel
from TC1_motor_model.motor_submodels.TC1_torque_limit_model import TorqueLimitModel
from TC1_motor_model.motor_submodels.TC1_implicit_em_torque_model import EMTorqueModel
from TC1_motor_model.motor_submodels.TC1_flux_weakening_model import FluxWeakeningModel, FluxWeakeningBracketModel
from TC1_motor_model.motor_submodels.TC1_mtpa_model import MTPAModel
from TC1_motor_model.motor_submodels.TC1_post_processing_model import PostProcessingModel


'''
THIS MODEL CONTAINS THE ENTIRE MOTOR ANALYSIS MODEL.
THE GOAL OF THIS MODEL IS TO IMPLICITLY SOLVE FOR POWER AND EFFICIENCY
GIVEN A RESIDUAL CONSIDERING TORQUE AND EFFICIENCY.

FEOC = FOR EACH OPERATING CONDITION
INPUTS:
    - motor geometry (from sizing model)
    - RPM and load torque (FEOC)
    - voltage limit
    - voltage (FEOC); fed back in from the battery model

OUTPUTS:
    - power and efficiency (FEOC); fed into other models
    - EM torque (FEOC)
'''

class ParseActiveOperatingConditions(csdl.CustomExplicitOperation):
    def initialize(self):
        self.parameters.declare('num_nodes')
        self.parameters.declare('num_active_nodes') # COMING FROM CADDEE

    def define(self):

        self.num_nodes = self.parameters['num_nodes']
        self.num_active_nodes = self.parameters['num_active_nodes']

        self.add_input('omega_rotor', shape=(self.num_nodes,))
        self.add_input('load_torque_rotor', shape=(self.num_nodes,))

        self.add_output('omega_rotor_active', shape=(self.num_active_nodes,))
        self.add_output('load_torque_rotor_active', shape=(self.num_active_nodes,))
        self.add_output('selection_indices', shape=(self.num_nodes,self.num_active_nodes))

        self.declare_derivatives('omega_rotor_active', 'omega_rotor')
        self.declare_derivatives('load_torque_rotor_active', 'load_torque_rotor')

    def compute(self, inputs, outputs):

        omega_rotor = inputs['omega_rotor']
        load_torque_rotor = inputs['load_torque_rotor']

        omega_rotor_active = []
        load_torque_rotor_active = []
        self.selection_indices = np.zeros((self.num_nodes, self.num_active_nodes))

        active_counter = 0
        for i in range(self.num_nodes):
            if omega_rotor[i] != 0 and load_torque_rotor[i] != 0:
                self.selection_indices[i, active_counter] = 1
                omega_rotor_active.append(omega_rotor[i])
                load_torque_rotor_active.append(load_torque_rotor[i])
                active_counter += 1
        
        outputs['omega_rotor_active'] = np.array(omega_rotor_active)
        outputs['load_torque_rotor_active'] = np.array(load_torque_rotor_active)
        outputs['selection_indices'] = self.selection_indices
        
    def compute_derivatives(self, inputs, derivatives):
        derivatives['omega_rotor_active', 'omega_rotor'] = np.transpose(self.selection_indices)
        derivatives['load_torque_rotor_active', 'load_torque_rotor'] = np.transpose(self.selection_indices)


class TC1MotorAnalysisModel(Model):
    def initialize(self, model_test=False):
        self.parameters.declare('pole_pairs') # 6
        self.parameters.declare('phases') # 3
        self.parameters.declare('num_slots') # 36
        self.parameters.declare('V_lim')
        self.parameters.declare('rated_current')
        self.parameters.declare('fit_coeff_dep_H') # FITTING COEFFICIENTS (X = H, B = f(H))
        self.parameters.declare('fit_coeff_dep_B') # FITTING COEFFICIENTS (X = B, H = g(B))
        self.parameters.declare('num_nodes')
        self.parameters.declare('num_active_nodes')
        self.parameters.declare('model_test', default=False)

        self.motor_variable_names = [
            'outer_stator_radius', 'pole_pitch', 'tooth_pitch', 'air_gap_depth', 'l_ef',
            'rotor_radius', 'turns_per_phase', 'Acu',  'tooth_width', 'height_yoke_stator',
            'slot_bottom_width', 'slot_height', 'slot_width_inner', 'Tau_y', 'L_j1', 'Kdp1',
            'bm', 'Am_r', 'phi_r', 'lambda_m', 'alpha_i', 'Kf', 'K_phi', 'K_theta', 'A_f2',
        ]

    def define(self):
        # INPUT PARAMETERS
        m = self.parameters['phases']
        p = self.parameters['pole_pairs']
        Z = self.parameters['num_slots']
        V_lim = self.parameters['V_lim']
        rated_current = self.parameters['rated_current']
        fit_coeff_dep_H = self.parameters['fit_coeff_dep_H']
        fit_coeff_dep_B = self.parameters['fit_coeff_dep_B']
        num_nodes = self.parameters['num_nodes']
        num_active_nodes = self.parameters['num_active_nodes']
        model_test=self.parameters['model_test']

        # DECLARE VARIABLES FROM SIZING & UPSTREAM MODELS
        D_i = self.declare_variable('D_i') # Diameter (DV or input)
        T_em_max = self.declare_variable('T_em_max') # Max torque (structural)
        Rdc = self.declare_variable('Rdc') # Resistance
        motor_variables = self.declare_variable('motor_variables', shape=(25,)) # array of motor sizing outputs
        for i in range(motor_variables.shape[0]):
            self.register_output(self.motor_variable_names[i], motor_variables[i])

        # ========================= MAGNET MEC =========================
        self.add(
            MagnetMECModel(
                fit_coeff_dep_H=fit_coeff_dep_H,
                fit_coeff_dep_B=fit_coeff_dep_B,
            ),
            'magnet_MEC_model',
        )

        # ========================= INDUCTANCE MEC =========================
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

        # ========================= TORQUE AND RPM INPUT ADJUSTMENT =========================
        omega_rotor = self.declare_variable('omega_rotor', shape=(num_nodes,))
        load_torque_rotor = self.declare_variable('load_torque_rotor', shape=(num_nodes,))

        omega_rotor_active, load_torque_rotor_active, selection_indices = csdl.custom(
            omega_rotor, load_torque_rotor,
            op=ParseActiveOperatingConditions(num_nodes=num_nodes, num_active_nodes=num_active_nodes)
        )

        omega_rotor_active = self.register_output('omega_rotor_active', omega_rotor_active)
        load_torque_rotor_active = self.register_output('load_torque_rotor_active', load_torque_rotor_active)

        R_expanded = self.register_output('R_expanded', csdl.expand(Rdc, (num_active_nodes,)))
        L_d_expanded = self.register_output('L_d_expanded', csdl.expand(L_d, (num_active_nodes,)))
        L_q_expanded = self.register_output('L_q_expanded', csdl.expand(L_q, (num_active_nodes,)))
        PsiF_expanded = self.register_output('PsiF_expanded', csdl.expand(PsiF, (num_active_nodes,)))

        # ========================= GEARBOX =========================
        gear_ratio = 4.
        omega = self.register_output('omega', omega_rotor_active * gear_ratio * 2*np.pi/60)
        load_torque = self.register_output('load_torque', load_torque_rotor_active/gear_ratio)

        # ========================= FINDING UPPER TORQUE LIMIT (DISCRIMINANT = 0 CASE) =========================
        self.add(
                TorqueLimitModel(
                    pole_pairs=p,
                    V_lim=V_lim,
                    num_nodes=num_active_nodes
                ),
                'torque_limit_model'
        )
        T_lim = self.declare_variable('T_lim', shape=(num_active_nodes,))
        T_lower_lim = self.declare_variable('T_lower_lim', val=0., shape=(num_active_nodes,))
        
        D = (3*p*(L_d_expanded-L_q_expanded))

        # ========================= FLUX WEAKENING BRACKETING IMPLICIT MODEL (TO FIND BRACKET LIMIT) =========================
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

        if model_test == False:
            self.add(
                EMTorqueModel(
                    pole_pairs=p,
                    V_lim=V_lim,
                    num_nodes=num_active_nodes,
                    rated_current=rated_current,
                    phases=m,
                    motor_variable_names=self.motor_variable_names
                ),
                'implicit_em_torque_model'
            )

            input_power_active = self.declare_variable('input_power_active', shape=(num_active_nodes,))
            input_power = csdl.matvec(selection_indices, input_power_active)
            self.register_output('input_power', input_power)

            '''
            NOTE TO SELF:
            Make this active calculation available for all relevant outputs:
                - efficiency, input power, output power, load torque, EM torque (potentially others)
            '''
            
            # CALCULATING UPPER LIMIT TORQUE CURVE  
            T_upper_lim_curve = self.register_output(
                'T_upper_lim_curve',
                # -csdl.log(csdl.exp(-T_lim) + csdl.exp(-csdl.expand(T_em_max, (num_active_nodes,))))
                csdl.min(csdl.reshape(T_lim, new_shape=(num_active_nodes, )),
                            csdl.expand(T_em_max, (num_active_nodes, )))
            )

            max_torque_constraint = self.register_output(name='max_torque_constraint',
                                                            var=T_upper_lim_curve-load_torque)


        else:
    # ========================= IMPLICIT T_EM MODEL =========================
            T_em = self.declare_variable('T_em', shape=(num_nodes,)) # STATE of implicit model

            # FLUX WEAKENING MODEL
            self.add(
                FluxWeakeningModel(
                    pole_pairs=p,
                    V_lim=V_lim,
                    num_nodes=num_nodes
                ),
                'flux_weakening_model',
            )
            
            # MTPA MODEL
            self.add(
                MTPAModel(
                    pole_pairs=p,
                    num_nodes=num_nodes
                ),
                'mtpa_model',
            )
            
            # SMOOTHING

            # psi_d = L_d_expanded*I_d*np.sqrt(2) + PsiF_expanded # NOT USED
            # psi_q = L_q_expanded*I_d*np.sqrt(2) # NOT USED

            I_q_rated = self.declare_variable('I_q_temp')
            I_q_rated_expanded = csdl.expand(I_q_rated, shape=(num_nodes,))

            f_i = 5000*p/60 # rated omega from sizing model = 3000
            U_d = -R_expanded*rated_current*np.sin(0.6283) - 2*np.pi*f_i*L_q_expanded*I_q_rated_expanded
            U_q = R_expanded*rated_current*np.sin(0.6283) + 2*np.pi*f_i*(PsiF_expanded - L_d_expanded*I_q_rated_expanded)

            U_rated = self.register_output(
                'voltage_amplitude',
                (U_d**2 + U_q**2)**(1/2)
            )

            Iq_fw = self.declare_variable('Iq_fw', shape=(num_nodes,))
            Iq_MTPA = self.declare_variable('Iq_MTPA', shape=(num_nodes,)) # CHECK NAMING SCHEME FOR VARIABLE
            k = 1 # ASK SHUOFENG WHAT THIS IS

            I_q = (csdl.exp(k*(U_rated - V_lim))*Iq_fw + Iq_MTPA) / (csdl.exp(k*(U_rated - V_lim)) + 1.0)
            I_d = (T_em / (1.5*p*I_q) - PsiF_expanded) / (L_d_expanded-L_q_expanded) # CHECK SIZE OF COMPUTATIONS HERE
            
            current_amplitude = self.register_output(
                'current_amplitude',
                (I_q**2 + I_d**2)**0.5
            )

            
            ''' POWER LOSS CALCULATIONS '''
            # load power
            # eq of the form P0 = speed * torque
            P0 = load_torque * omega
            self.register_output('output_power', P0)
            frequency = omega*p/60

            # copper loss
            R_expanded = csdl.expand(Rdc, (num_nodes,))
            P_copper = m*R_expanded*current_amplitude**2 # NEED TO CHECK IF THIS IS ELEMENT-WISE SQUARING
            
            # eddy_loss
            a = 0.00055 # lamination thickness in m
            sigma_c = 2e6 # bulk conductivity (2000000 S/m)
            l_ef = self.declare_variable('l_ef')
            D1 = self.declare_variable('outer_stator_radius')
            D_i = self.declare_variable('D_i')
            Acu = self.declare_variable('Acu')
            B_delta = self.declare_variable('B_delta')
            B_delta_expanded = csdl.expand(B_delta, (num_nodes,))
            
            K_e = (a*np.pi)**2 * sigma_c/6
            V_s = csdl.expand((np.pi*l_ef*(D1-D_i)**2)-36*l_ef*Acu, (num_nodes,)); # volume of stator
            K_c = 0.822;
            P_eddy = K_e*V_s*(B_delta_expanded*frequency)**2; # eddy loss

            # hysteresis loss
            K_h = 100
            n_h = 2
            P_h = K_h*V_s*frequency*B_delta_expanded**n_h

            # stress loss
            P_stress = 0.01*P0

            # windage & friction loss
            k_r = 4 # roughness coefficient
            fr = 3e-3 # friction coefficient
            rho_air = 1.225 # air density kg/m^3
            D2 = self.declare_variable('rotor_radius')
            l_ef_expanded = csdl.expand(l_ef, (num_nodes,))
            D2_expanded = csdl.expand(D2, (num_nodes,))
            P_wo = k_r*np.pi*fr*rho_air*(2*np.pi*frequency)**2*l_ef_expanded*D2_expanded**4

            # total losses
            P_loss = P_copper + P_eddy + P_h + P_stress + P_wo
            input_power = self.register_output('input_power', P0 + P_loss)
            efficiency = self.register_output('efficiency', P0/input_power)
                
                # residual = model.register_output(
                #     'residual',
                #     load_torque - efficiency*T_em
                # )
                # bracket between 0 and smoothing of torque 

                # FROM HERE, THE OUTPUT POWER AND EFFICIENCY ARE PROPOGATED TO THE
                # BATTERY ANALYSIS MODELS

        