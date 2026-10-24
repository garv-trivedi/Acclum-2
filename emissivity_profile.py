"""Region-wise spectral surface-flux prescriptions for a truncated SS disk.

Equations are transcribed from the user's
'Emissivity_functions_of_regions_of_truncated_accretion_disk.pdf'.
The cgs expressions are converted to SI at the public function boundary.

Region (a): Wien-limit Planck surface flux (normalization matched to Eq. 9).
Region (b): homogeneous or exponentially stratified prescription.
Region (c): Planck surface flux.

These are emergent spectral surface-flux functions F_nu, not volume emissivity j_nu.
"""

import numpy as np

from Temp_Profile import (
    temp_region_a,
    temp_region_b,
    density_region_b,
    temp_region_c,
)

# CGS constants required by the supplied emission equations
H_CGS = 6.62607015e-27             # erg s
K_B_CGS = 1.380649e-16             # erg K^-1
C_CGS = 2.99792458e10              # cm s^-1
SIGMA_T_CGS = 6.6524587321e-25     # cm^2
M_P_CGS = 1.67262192369e-24        # g

# 1 erg s^-1 cm^-2 Hz^-1 = 1e-3 W m^-2 Hz^-1
CGS_FLUX_TO_SI = 1.0e-3


def _dimensionless_energy(nu_hz, temperature_k):
    """Return x = h nu / (k_B T); nu is Hz and T is K."""
    nu = float(nu_hz)
    T = np.asarray(temperature_k, dtype=float)

    if not np.isfinite(nu) or nu <= 0.0:
        raise ValueError("Frequency nu_hz must be finite and > 0 Hz.")
    if np.any(~np.isfinite(T)) or np.any(T <= 0.0):
        raise ValueError("Temperature must be finite and > 0 K.")

    return H_CGS * nu / (K_B_CGS * T)


def planck_surface_flux_nu_cgs(nu_hz, temperature_k):
    r"""Eq. (9): B(x) = 2*pi*h*nu^3/c^2 /(exp(x)-1), in cgs F_nu units.

    This is the Planck spectral surface flux used as B(x) in the supplied PDF,
    i.e. pi times Planck specific intensity.
    """
    nu = float(nu_hz)
    x = _dimensionless_energy(nu, temperature_k)
    prefactor = 2.0 * np.pi * H_CGS * nu**3 / C_CGS**2

    x_safe = np.minimum(x, 700.0)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        result = prefactor / np.expm1(x_safe)
    result = np.where(x >= 700.0, 0.0, result)
    return result


def free_free_chi_nu(nu_hz, temperature_k):
    r"""Eq. (8): chi_nu = 4.1e-23 (1-exp(-x))/(T^(7/2) x^3).

    The numeric coefficient and constants follow the supplied PDF's cgs convention.
    """
    T = np.asarray(temperature_k, dtype=float)
    x = _dimensionless_energy(nu_hz, T)
    numerator = -np.expm1(-x)  # stable evaluation of 1 - exp(-x)
    return 4.1e-23 * numerator / (T**(7.0 / 2.0) * x**3)


def scale_height_region_b_cm(r_rs, alpha, m, mdot):
    r"""Region-(b) vertical scale height H0 = z0,b, in cm.

    Uses the region-(b) structure relation in the supplied corrected SS
    region summary: z0,b = 1.2e4 alpha^(-1/10) mdot^(1/5) m^(9/10)
                    r^(21/20) (1-r^(-1/2))^(1/5) cm,
    with paper radius r = R/(3 R_S) = r_rs/3.

    If your current Temp_Profile.py uses a different normalization for z0,b,
    keep that convention consistent here too.
    """
    r_paper = np.asarray(r_rs, dtype=float) / 3.0
    if np.any(r_paper <= 1.0):
        raise ValueError("Region-(b) scale height requires R/R_S > 3.")
    boundary = 1.0 - r_paper**(-0.5)
    return (
        1.2e4
        * alpha**(-1.0 / 10.0)
        * mdot**(1.0 / 5.0)
        * m**(9.0 / 10.0)
        * r_paper**(21.0 / 20.0)
        * boundary**(1.0 / 5.0)
    )


def emission_flux_nu_region_a(nu_hz, r_rs, alpha, m):
    r"""Region (a), Eq. (2): F_a(x) proportional to x^3 exp(-x).

    The PDF gives only proportionality. To allow an absolute SED to be summed,
    this implementation uses the Wien limit of the PDF's Planck surface-flux
    Eq. (9): F_nu,a = (2*pi*h*nu^3/c^2) exp(-x).

    Returns W m^-2 Hz^-1. This normalization choice should be documented if
    comparing absolute region-a fluxes with another prescription.
    """
    T_a = temp_region_a(r_rs, alpha, m)
    x = _dimensionless_energy(nu_hz, T_a)
    prefactor = 2.0 * np.pi * H_CGS * float(nu_hz)**3 / C_CGS**2
    F_cgs = prefactor * np.exp(-np.minimum(x, 745.0))
    return F_cgs * CGS_FLUX_TO_SI


def emission_flux_nu_region_b_homogeneous(
    nu_hz, r_rs, alpha, m, mdot, f1=1.0, eta_E=0.06
):
    r"""Region (b), Eq. (3): homogeneous-atmosphere spectral surface flux.

    F_b,hom = [3 chi_nu n / (sigma_T m_p)]^(1/2) B(x).
    Density is expected in cm^-3; constants and chi follow cgs convention.
    Returns W m^-2 Hz^-1.
    """
    T_b = temp_region_b(r_rs, alpha, m, mdot, f1, eta_E)
    n_b = density_region_b(r_rs, alpha, m, mdot, f1, eta_E)
    chi = free_free_chi_nu(nu_hz, T_b)
    B_cgs = planck_surface_flux_nu_cgs(nu_hz, T_b)

    factor_squared = 3.0 * chi * n_b / (SIGMA_T_CGS * M_P_CGS)
    if np.any(~np.isfinite(factor_squared)) or np.any(factor_squared < 0.0):
        raise ValueError("Invalid homogeneous Region-(b) emission prefactor.")

    return np.sqrt(factor_squared) * B_cgs * CGS_FLUX_TO_SI


def emission_flux_nu_region_b_exponential(
    nu_hz, r_rs, alpha, m, mdot, f1=1.0, eta_E=0.06, H0_cm=None
):
    r"""Region (b), Eq. (6): exponentially stratified-atmosphere flux.

    F_b,exp = [3 chi_nu / (sigma_T^2 m_p^2 H0)]^(1/3) B(x).

    If H0_cm is omitted, use the Region-(b) SS vertical scale-height relation
    implemented in scale_height_region_b_cm(). A supplied H0_cm must be in cm.
    Returns W m^-2 Hz^-1.
    """
    T_b = temp_region_b(r_rs, alpha, m, mdot, f1, eta_E)
    chi = free_free_chi_nu(nu_hz, T_b)
    B_cgs = planck_surface_flux_nu_cgs(nu_hz, T_b)

    if H0_cm is None:
        H0 = scale_height_region_b_cm(r_rs, alpha, m, mdot)
    else:
        H0 = np.asarray(H0_cm, dtype=float)
        if np.any(~np.isfinite(H0)) or np.any(H0 <= 0.0):
            raise ValueError("H0_cm must be finite and > 0 cm.")

    factor_cubed = (
        3.0 * chi /
        (SIGMA_T_CGS**2 * M_P_CGS**2 * H0)
    )
    if np.any(~np.isfinite(factor_cubed)) or np.any(factor_cubed < 0.0):
        raise ValueError("Invalid exponential Region-(b) emission prefactor.")

    return np.cbrt(factor_cubed) * B_cgs * CGS_FLUX_TO_SI


def emission_flux_nu_region_b(
    nu_hz, r_rs, alpha, m, mdot,
    f1=1.0, eta_E=0.06, model="homogeneous", H0_cm=None
):
    """Select one of the two Region-(b) emission prescriptions."""
    model = model.strip().lower()
    if model in {"homogeneous", "homogeneous medium", "hom"}:
        return emission_flux_nu_region_b_homogeneous(
            nu_hz, r_rs, alpha, m, mdot, f1, eta_E
        )
    if model in {"exponential", "exponentially stratified", "stratified"}:
        return emission_flux_nu_region_b_exponential(
            nu_hz, r_rs, alpha, m, mdot, f1, eta_E, H0_cm
        )
    raise ValueError("Region-b model must be 'homogeneous' or 'exponential'.")


def emission_flux_nu_region_c(nu_hz, r_rs, alpha, m, mdot, f1=1.0, eta_E=0.06):
    r"""Region (c), Eq. (9): Planck blackbody spectral surface flux.

    Returns W m^-2 Hz^-1.
    """
    T_c = temp_region_c(r_rs, alpha, m, mdot, f1, eta_E)
    return planck_surface_flux_nu_cgs(nu_hz, T_c) * CGS_FLUX_TO_SI
