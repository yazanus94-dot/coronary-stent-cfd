"""Shared fixtures and configuration for the Coronary-Stent-CFD test suite."""

from __future__ import annotations

import numpy as np
import pytest

from stent_hemo.core import BloodProperties, FlowWaveform, VesselGeometry
from stent_hemo.struct_grid import StructGrid


@pytest.fixture
def blood() -> BloodProperties:
    """Default blood properties."""
    return BloodProperties()


@pytest.fixture
def wf() -> FlowWaveform:
    """Default coronary flow waveform."""
    return FlowWaveform()


@pytest.fixture
def vessel_straight() -> VesselGeometry:
    """Unstenosed straight vessel, 30 mm length."""
    return VesselGeometry(severity=0.0, length=30e-3)


@pytest.fixture
def vessel_stenosis_50() -> VesselGeometry:
    """50% area stenosis vessel, 90 mm length."""
    return VesselGeometry(severity=0.5, length=90e-3)


@pytest.fixture
def vessel_stenosis_70() -> VesselGeometry:
    """70% area stenosis vessel, 90 mm length."""
    return VesselGeometry(severity=0.7, length=90e-3)


@pytest.fixture
def grid_coarse(vessel_straight: VesselGeometry) -> StructGrid:
    """Coarse structured grid (20x8) for fast tests."""
    return StructGrid(vessel_straight, 20, 8)


@pytest.fixture
def grid_medium(vessel_straight: VesselGeometry) -> StructGrid:
    """Medium structured grid (40x14)."""
    return StructGrid(vessel_straight, 40, 14)


@pytest.fixture
def grid_stenosis(vessel_stenosis_50: VesselGeometry) -> StructGrid:
    """Structured grid for 50% stenosis (80x24)."""
    return StructGrid(vessel_stenosis_50, 80, 24)


# Tolerance constants used across test modules
POISEUILLE_DP_TOL = 0.15  # 15% relative tolerance on straight-tube dp
REYNOLDS_LAMINAR_MAX = 500  # coronary regime upper bound
GCI_ACCEPT_PCT = 5.0  # max acceptable GCI on fine grid
