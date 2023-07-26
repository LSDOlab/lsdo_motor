import numpy as np
import csdl
from python_csdl_backend import Simulator

from lsdo_motor import TC1MotorAnalysisModel, TC1MotorSizingModel

def_coeff_H = np.array([1.92052530e-03,  1.03633835e+01, -2.98809161e+00])
def_coeff_B = np.array([1.12832651, 1., 1., 1.])

class FullMotorModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('num_nodes', default=1)
        self.parameters.declare('num_active_nodes', default=None)
        self.parameters.declare('pole_pairs', default=6)
        self.parameters.declare('phases', default=3)
        self.parameters.declare('num_slots', default=36)
        self.parameters.declare('V_lim', default=400)
        self.parameters.declare('rated_current', default=120.)
        self.parameters.declare('fit_coeff_dep_H', default=def_coeff_H) # FITTING COEFFICIENTS (X = H, B = f(H))
        self.parameters.declare('fit_coeff_dep_B', default=def_coeff_B) # FITTING COEFFICIENTS (X = B, H = g(B))

    def define(self):
        num_nodes = self.parameters['num_nodes']
        num_active_nodes = self.parameters['num_active_nodes']
        pole_pairs = self.parameters['pole_pairs']
        phases = self.parameters['phases']
        num_slots = self.parameters['num_slots']
        V_lim = self.parameters['V_lim']
        I_rated = self.parameters['rated_current']
        fit_coeff_dep_H = self.parameters['fit_coeff_dep_H']
        fit_coeff_dep_B = self.parameters['fit_coeff_dep_B']

        self.add(
            TC1MotorSizingModel(
                pole_pairs=pole_pairs,
                phases=phases,
                num_slots=num_slots,
                rated_current=I_rated
            ),
            'sizing_model'
        )

        self.add(
            TC1MotorAnalysisModel(
                pole_pairs=pole_pairs,
                phases=phases,
                num_slots=num_slots,
                V_lim=V_lim,
                rated_current=I_rated,
                fit_coeff_dep_B=fit_coeff_dep_B,
                fit_coeff_dep_H=fit_coeff_dep_H,
                num_nodes=num_nodes,
                num_active_nodes=num_active_nodes,
                use_caddee=False
            ),
            'analysis_model'
        )


D_i = 0.182
L = 0.086

load_torque_rotor = np.array([172.317]) # CL, CR, 0, H, H
omega_rotor = np.array([3375.10]) # CL, CR, 0, H, H

num_nodes=len(load_torque_rotor)

m = FullMotorModel(
    num_nodes=num_nodes,
    V_lim=800
)

sim = Simulator(m)
sim['D_i'] = D_i
sim['L'] = L

sim['omega_rotor'] = omega_rotor
sim['load_torque_rotor'] = load_torque_rotor

sim.run()