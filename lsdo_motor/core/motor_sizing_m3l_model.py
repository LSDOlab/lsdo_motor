import m3l
import numpy as np
from lsdo_motor.core.TC1_motor_sizing_model import TC1MotorSizingModel

class MotorSizing(m3l.ExplicitOperation):
    def initialize(self, kwargs):
        self.parameters.declare('rotor_component', default='lsdolab')
        self.parameters.declare('pole_pairs', default=6)
        self.parameters.declare('phases', default=3)
        self.parameters.declare('num_slots', default=36)
        self.parameters.declare('rated_current', default=120.)

    def compute(self):
        model = TC1MotorSizingModel(
            module=self,
            pole_pairs=self.pole_pairs,
            phases=self.phases,
            num_slots=self.num_slots,
            rated_current=self.I_rated,
            component_name=self.component_name,
        )
        return model
    
    def evaluate(self) -> m3l.Variable:
        self.pole_pairs = self.parameters['pole_pairs']
        self.phases = self.parameters['phases']
        self.num_slots = self.parameters['num_slots']
        self.I_rated = self.parameters['rated_current']
        component = self.parameters['rotor_component']

        if isinstance(component, str):
            self.component_name = component
        else:
            self.component_name = component.parameters['name']

        self.name = f'{self.component_name}_motor_sizing_model'
        self.arguments = {}

        motor_mass = m3l.Variable(
            name='mass',
            shape=(1,),
            operation=self
        )

        motor_cg = m3l.Variable(
            name='cg_vector',
            shape=(3,),
            operation=self
        )

        motor_inertia = m3l.Variable(
            name='inertia_tensor',
            shape=(3, 3),
            operation=self
        )

        motor_parameters = m3l.Variable(
            name='motor_parameters',
            shape=(27,),
            operation=self
        ) # motor_parameters includes resistance and max torque

        return motor_mass, motor_cg, motor_inertia, motor_parameters
    
    # NOTE: # ADD PREFIX TO NAMES OF INPUTS AND OUTPUTS