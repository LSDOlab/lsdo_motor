__version__ = '0.0.1'

from lsdo_motor.core.motor_analysis_m3l_model import MotorAnalysis, evaluate_multiple_motor_analysis_models
from lsdo_motor.core.motor_sizing_m3l_model import MotorSizing, evaluate_multiple_motor_sizing_models

from lsdo_motor.core.TC1_motor_analysis_model import TC1MotorAnalysisModel
from lsdo_motor.core.TC1_motor_sizing_model import TC1MotorSizingModel

from pathlib import Path
ROOT = Path(__file__).parents[0]
GEOMETRY_PATH = ROOT / 'core' / 'sample_geometry' / 'geometry'
IMPORTS_PATH = ROOT / 'core' / 'sample_geometry' / 'imports'