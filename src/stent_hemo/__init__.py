"""Hemodynamics of a Stenosed Artery and the Effect of a Stent.

Pulsatile, shear-thinning blood flow in a stenosed artery (50 % and 70 %
area stenosis), a parametric lattice stent, and the change in wall shear
stress, pressure drop and recirculation after stenting. Analytic core
(Poiseuille, Womersley, Young-Tsai, GCI) always runs; OpenFOAM, SU2 and
CalculiX are auto-detected wrappers marked TBD when absent.
"""

from .core import BloodProperties, FlowWaveform, StentParams, VesselGeometry
from .gci import grid_convergence_index
from .hand_calc import flow_hand_calc
from .womersley import womersley_profile, womersley_wss
from .young_tsai import stenosis_pressure_drop

__version__ = "1.0.0"

__all__ = [
    "BloodProperties",
    "FlowWaveform",
    "StentParams",
    "VesselGeometry",
    "flow_hand_calc",
    "womersley_profile",
    "womersley_wss",
    "stenosis_pressure_drop",
    "grid_convergence_index",
    "__version__",
]
