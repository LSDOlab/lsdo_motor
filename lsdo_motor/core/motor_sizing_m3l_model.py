import m3l
import numpy as np
from lsdo_motor.core.TC1_motor_sizing_model import TC1MotorSizingModel

class MotorSizing(m3l.ExplicitOperation):
    def initialize(self):
        self.parameters.declare('component_name')
        self.parameters.declare('pole_pairs', default=6)
        self.parameters.declare('phases', default=3)
        self.parameters.declare('num_slots', default=36)
        self.parameters.declare('rated_current', default=120.)
        self.parameters.declare('num_nodes', default=1)

    def compute(self):
        model = TC1MotorSizingModel(
            pole_pairs=self.pole_pairs,
            phases=self.phases,
            num_slots=self.num_slots,
            rated_current=self.I_rated,
        )
        return model
    
    def evaluate(self, 
                 motor_diameter: m3l.Variable=None,
                 motor_length: m3l.Variable=None,
                 design_condition=None) -> m3l.Variable:
        self.pole_pairs = self.parameters['pole_pairs']
        self.phases = self.parameters['phases']
        self.num_slots = self.parameters['num_slots']
        self.I_rated = self.parameters['rated_current']
        self.num_nodes = self.parameters['num_nodes']
        component_name = self.parameters['component_name']

        if design_condition:
            dc_name = design_condition.parameters['name']
            self.name = f'{dc_name}_{component_name}_motor_sizing_model'
        else:
            self.name = f'{component_name}_motor_sizing_model'

        self.arguments = {}
        self.arguments['D_i'] = motor_diameter
        self.arguments['L'] = motor_length

        motor_mass = m3l.Variable(
            name='motor_mass',
            shape=(self.num_nodes,),
            operation=self
        )
        motor_parameters = m3l.Variable(
            name='motor_parameters',
            shape=(self.num_nodes,),
            operation=self
        ) # motor_parameters includes resistance and max torque

        return motor_mass, motor_parameters
    
    # NOTE: # ADD PREFIX TO NAMES OF INPUTS AND OUTPUTS