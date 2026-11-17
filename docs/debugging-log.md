# Debugging Log

Real issues encountered while building the in-house axisymmetric FV solver,
with root causes and regression guards. Each entry lists the symptom, the
root cause as actually diagnosed, the fix, and the guard that now prevents
silent recurrence.

## BUG-1: Patankar relaxation formula error (root cause of wrong converged state)

**Symptom:** Steady solves converged (tiny mass residual) but the straight-tube
profile was not Poiseuille and the centerline velocity was systematically off.

**Root cause:** The under-relaxation source term was implemented as
`b += (1 - relax)/relax * diag * u_guess * relax`. The trailing `* relax`
cancels part of the relaxation denominator and shifts the fixed point of
every momentum solve, so the converged field satisfied a *different* equation
than the PDE.

**Fix:** `b += (1 - relax)/relax * diag * u_guess` in both `_u_momentum`
and `_v_momentum`.

**Regression guard:** `test_steady_poiseuille_small_grid` and the
three-level Poiseuille grid study in `validation.py`: after the fix the
dp error drops to 1.4 % on the fine grid with a clean first-order trend.

## BUG-2: Double-assembled matrix coefficients

**Symptom:** Interior diffusion/convection coefficients came out twice as
large as hand-calculated values; solutions were over-damped.

**Root cause:** The matrix-assembly helper was visited once for each ordered
pair of neighbouring cells, but it assembled *both* sides of the pair each
time. SciPy COO duplicate entries are summed, so every interior coefficient
was doubled.

**Fix:** The helper now adds only the calling cell's side per visit
(diagonal + cross-entry); the partner cell contributes its own side when the
loop reaches it.

**Regression guard:** Unit tests on a straight tube reproduce the exact
Poiseuille dp against the closed-form value.

## BUG-3: Inlet/wall Dirichlet double counting

**Symptom:** Mass flux at the inlet did not equal the prescribed waveform
flow; the first axial cells showed spurious negative velocities.

**Root cause:** For a Dirichlet boundary face, the convective contribution
was added to both the matrix diagonal and the RHS. Only the RHS belongs
there (the boundary value is known).

**Fix:** Boundary convection goes only into `b`; boundary diffusion uses the
half-cell distance (`df = 2*mu*A/dz` at the inlet, `(mu+mu)*A/dru` at the
wall).

**Regression guard:** `test_transient_inlet_harmonic_periodicity` (exact
inlet profile reproduction) plus the mass-residual monitor in every solve.

## BUG-4: v-control-volume geometry used integer ring centers

**Symptom:** The radial-velocity system was inconsistent with the mass
operator; the checkerboard persisted at coarse grids.

**Root cause:** The radial (v) control volumes were constructed spanning
integer ring centers `f_{k-1}` to `f_{k+1}` instead of the half-integer
ring centers `(f_{k-1}+f_k)/2` to `(f_k+f_{k+1})/2` where the staggered v
values live.

**Fix:** All v-CV quantities (`vv`, `arv_in/out`, `azv_lo/hi`) rebuilt from
the half-integer mid-fractions.

**Regression guard:** Transient straight-tube run shows smooth (non-checkerboard)
radial velocity; `test_transient_no_divergence_short`.

## BUG-5: Corrections taken from total pressure instead of the correction field

**Symptom:** Steady solutions depended on the initial pressure guess; the
SIMPLE iteration did not converge to a unique field.

**Root cause:** Velocity corrections used the total pressure `p` instead of
the Poisson solution `pp`. After `p += relax * pp`, correcting with `p`
double-applies previous corrections.

**Fix:** `u[j] -= d_u * rp * (pp[j] - pp[j-1]) / dz` (and the v analogue).

**Regression guard:** Repeated runs from different initial guesses converge
to the same dp within 1e-6 (validated in the GCI study).

## BUG-6: Missing density in the transient term

**Symptom:** The first transient attempt diverged to NaN within a few steps.

**Root cause:** The transient (pseudo-time) source term `V/dt * u_old` was
missing the `rho` factor, making the transient matrix inconsistent with the
convection/diffusion scaling.

**Fix:** `coeff = rho * V / dt` in both momentum systems.

**Regression guard:** `test_transient_stays_finite` (integration) runs
transient cycles on stenosis geometry and asserts finite history.

## BUG-7: Transient Courant instability on fine grids (Carreau feedback loop)

**Symptom:** Transient Womersley runs on 120x24 diverged at dt = 2 ms while
the same dt was stable on 60x12 and on straight-tube geometry.

**Root cause:** Courant number scales with dr, and at the throat the
shear-thinning feedback accelerates divergence: overshoot -> higher shear ->
viscosity drops to mu_inf -> less damping -> larger overshoot.

**Fix:** dt <= 1 ms on the production grid (verified stable, ~7x slower);
post-processing defaults to 60x12/dt=2 ms for the WSS metrics and raises a
loud `RuntimeError` if NaN appears.

**Regression guard:** `test_transient_no_divergence_short` plus the NaN
guard in `postprocess.run_postprocess`.

## BUG-8: Richardson-extrapolation sign error in the GCI

**Symptom:** The GCI module returned `f_exact` farther from the analytic
value than the fine-grid value itself.

**Root cause:** `f_exact = f_fine + (f_medium - f_fine)/(r^p - 1)` used the
wrong sign; the correct Richardson form is
`f_exact = f_fine - (f_medium - f_fine)/(r^p - 1)`.

**Fix:** Sign corrected in `gci.py`.

**Regression guard:** `test_gci_on_exact_quadratic` recovers p = 2 and the
exact solution from an analytically constructed family; the validation GCI
now reproduces the analytic Poiseuille dp to 0.06 %.

## Documented blocker: Gmsh 3-D tetrahedral meshing

**Symptom:** OpenCASCADE STEP solids fail HXT 3-D meshing (reproducible
across gmsh versions tried).

**Resolution:** Not fixed in this project; the structured (z, r) wedge grid
is the primary mesh path and is what the in-house solver, OpenFOAM and SU2
cases all consume. `meshing.py` is kept as a documented fallback that raises
a clear, actionable error instead of failing silently.
