"""Property-based tests using Hypothesis.

Tests physics invariants over random valid inputs for the stenosis hemodynamics model:
- Young-Tsai pressure drop increases monotonically with stenosis severity
- Hand calculations produce positive physical quantities
- Vessel radius stays positive for all valid severities
"""

from __future__ import annotations

from hypothesis import given, settings
from hypothesis import strategies as st

from stent_hemo.core import BloodProperties, FlowWaveform, VesselGeometry
from stent_hemo.hand_calc import flow_hand_calc
from stent_hemo.young_tsai import stenosis_pressure_drop

# Strategies for valid physical inputs
severity_strategy = st.floats(min_value=0.0, max_value=0.9, allow_nan=False)
length_strategy = st.floats(min_value=10e-3, max_value=200e-3, allow_nan=False)
diameter_strategy = st.floats(min_value=1e-3, max_value=10e-3, allow_nan=False)


@settings(max_examples=50, deadline=None)
@given(severity=severity_strategy, length=length_strategy, d0=diameter_strategy)
def test_young_tsai_monotonic_in_severity(severity, length, d0):
    """Young-Tsai pressure drop increases monotonically with stenosis severity."""
    blood = BloodProperties()
    wf = FlowWaveform()
    q = wf.mean_flow()

    v_low = VesselGeometry(severity=max(0.0, severity - 0.1), length=length, d0=d0)
    v_high = VesselGeometry(severity=severity, length=length, d0=d0)

    drop_low = stenosis_pressure_drop(v_low, blood, q).dp_total
    drop_high = stenosis_pressure_drop(v_high, blood, q).dp_total

    assert drop_high >= drop_low * 0.99  # monotonic within numerical tolerance


@settings(max_examples=50, deadline=None)
@given(severity=severity_strategy, length=length_strategy, d0=diameter_strategy)
def test_hand_calc_positive(severity, length, d0):
    """Hand calculations produce positive physical quantities."""
    blood = BloodProperties()
    wf = FlowWaveform()
    v = VesselGeometry(severity=severity, length=length, d0=d0)
    hc = flow_hand_calc(v, blood, wf)

    assert hc.u_mean >= 0
    assert hc.reynolds > 0
    assert hc.womersley > 0
    assert hc.tau_poiseuille >= 0
    assert hc.dp_vessel >= 0


@settings(max_examples=50, deadline=None)
@given(severity=severity_strategy, length=length_strategy, d0=diameter_strategy)
def test_radius_stays_positive(severity, length, d0):
    """Vessel radius stays positive for all valid severities."""
    import numpy as np

    v = VesselGeometry(severity=severity, length=length, d0=d0)
    z = np.linspace(0.0, v.length, 50)
    r = v.radius_at(z)
    assert np.all(r > 0)
    assert np.all(r <= v.r0 + 1e-12)
