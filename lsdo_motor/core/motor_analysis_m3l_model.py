import m3l
import numpy as np
from lsdo_motor.core.TC1_motor_analysis_model import TC1MotorAnalysisModel

def_coeff_H = np.array([1.92052530e-03,  1.03633835e+01, -2.98809161e+00])
def_coeff_B = np.array([1.12832651, 1., 1., 1.])

class MotorAnalysis(m3l.ExplicitOperation):
    def initialize(self, kwargs):
        self.parameters.declare('component_name')
        self.parameters.declare('pole_pairs', default=6)
        self.parameters.declare('phases', default=3)
        self.parameters.declare('num_slots', default=36)
        self.parameters.declare('rated_current', default=120.)
        self.parameters.declare('num_nodes', default=1)
        self.parameters.declare('V_lim')
        self.parameters.declare('fit_coeff_dep_H', default=def_coeff_H) # FITTING COEFFICIENTS (X = H, B = f(H))
        self.parameters.declare('fit_coeff_dep_B', default=def_coeff_B) # FITTING COEFFICIENTS (X = B, H = g(B))

    def compute(self):
        model = TC1MotorAnalysisModel(
            pole_pairs=self.pole_pairs,
            phases=self.phases,
            num_slots=self.num_slots,
            num_nodes=self.num_nodes,
            rated_current=self.I_rated,
            V_lim=self.V_lim,
            fit_coeff_dep_H=self.fit_coeff_dep_H,
            fit_coeff_dep_B=self.fit_coeff_dep_B,
        )
        return model
    
    def evaluate(self, 
                 torque: m3l.Variable=None,
                 motor_parameters: m3l.Variable=None,
                 design_condition=None) -> m3l.Variable:
        self.pole_pairs = self.parameters['pole_pairs']
        self.phases = self.parameters['phases']
        self.num_slots = self.parameters['num_slots']
        self.I_rated = self.parameters['rated_current']
        self.num_nodes = self.parameters['num_nodes']
        self.V_lim = self.parameters['V_lim']
        self.fit_coeff_dep_H = self.parameters['fit_coeff_dep_H']
        self.fit_coeff_dep_B = self.parameters['fit_coeff_dep_B']

        component_name = self.parameters['component_name']

        if design_condition:
            dc_name = design_condition.parameters['name']
            self.name = f'{dc_name}_{component_name}_analysis_model'
        else:
            self.name = f'{component_name}_analysis_model'

        self.arguments={}
        self.arguments['load_torque_rotor'] = torque
        self.arguments['motor_parameters'] = motor_parameters

        output_power = m3l.Variable(
            name='output_power',
            shape=(self.num_nodes,),
            operation=self
        )
        # efficiency = m3l.Variable(
        #     name='efficiency',
        #     shape=(self.num_nodes,),
        #     operation=self
        # )

        return output_power
    # efficiency