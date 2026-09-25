import numpy as np


# ============================================================
# REGION (a)
# ============================================================

def temp_region_a(r_rs, alpha, m):
    """
    Temperature in Region (a).

    r_rs = R / R_S
    m    = M_BH / M_sun

    Returns K.
    """

    r_rs = np.asarray(r_rs, dtype=float)

    return (
        2.3e7
        * 3.0 ** (3.0 / 4.0)
        * (alpha * m) ** (-1.0 / 4.0)
        * r_rs ** (-3.0 / 4.0)
    )


def density_region_a(
    r_rs,
    alpha,
    m,
    mdot,
    f1=1.0,
    eta_E=0.06
):
    """
    Number density in Region (a).

    r_rs = R / R_S
    mdot = Shakura-Sunyaev dimensionless mdot

    Returns cm^-3.
    """

    r_rs = np.asarray(r_rs, dtype=float)

    return (
        4.3e17
        * 3.0 ** (-3.0 / 2.0)
        * (
            f1 * (0.06 / eta_E)
        ) ** (-2.0)
        * alpha ** (-1.0)
        * m ** (-1.0)
        * mdot ** (-2.0)
        * r_rs ** (3.0 / 2.0)
        * (
            1.0 - np.sqrt(3.0 / r_rs)
        ) ** (-2.0)
    )


# ============================================================
# REGION (b)
# ============================================================

def temp_region_b(
    r_rs,
    alpha,
    m,
    mdot,
    f1=1.0,
    eta_E=0.06
):
    """
    Temperature in Region (b).

    Returns K.
    """

    r_rs = np.asarray(r_rs, dtype=float)

    return (
        3.1e8
        * 3.0 ** (9.0 / 10.0)
        * (
            f1 * (0.06 / eta_E)
        ) ** (2.0 / 5.0)
        * alpha ** (-1.0 / 5.0)
        * m ** (-1.0 / 5.0)
        * mdot ** (2.0 / 5.0)
        * r_rs ** (-9.0 / 10.0)
        * (
            1.0 - np.sqrt(3.0 / r_rs)
        ) ** (2.0 / 5.0)
    )


def density_region_b(
    r_rs,
    alpha,
    m,
    mdot,
    f1=1.0,
    eta_E=0.06
):
    """
    Number density in Region (b).

    Returns cm^-3.
    """

    r_rs = np.asarray(r_rs, dtype=float)

    return (
        4.2e24
        * 3.0 ** (33.0 / 20.0)
        * (
            f1 * (0.06 / eta_E)
        ) ** (2.0 / 5.0)
        * alpha ** (-7.0 / 10.0)
        * m ** (-7.0 / 10.0)
        * mdot ** (2.0 / 5.0)
        * r_rs ** (-33.0 / 20.0)
        * (
            1.0 - np.sqrt(3.0 / r_rs)
        ) ** (2.0 / 5.0)
    )


# ============================================================
# REGION (c)
# ============================================================

def temp_region_c(
    r_rs,
    alpha,
    m,
    mdot,
    f1=1.0,
    eta_E=0.06
):
    """
    Temperature in Region (c).

    Returns K.
    """

    r_rs = np.asarray(r_rs, dtype=float)

    return (
        8.6e7
        * 3.0 ** (3.0 / 4.0)
        * (
            f1 * (0.06 / eta_E)
        ) ** (3.0 / 10.0)
        * alpha ** (-1.0 / 5.0)
        * m ** (-1.0 / 5.0)
        * mdot ** (3.0 / 10.0)
        * r_rs ** (-3.0 / 4.0)
        * (
            1.0 - np.sqrt(3.0 / r_rs)
        ) ** (3.0 / 10.0)
    )


def density_region_c(
    r_rs,
    alpha,
    m,
    mdot,
    f1=1.0,
    eta_E=0.06
):
    """
    Number density in Region (c).

    Returns cm^-3.
    """

    r_rs = np.asarray(r_rs, dtype=float)

    return (
        3.0e25
        * 3.0 ** (15.0 / 8.0)
        * (
            f1 * (0.06 / eta_E)
        ) ** (11.0 / 12.0)
        * alpha ** (-7.0 / 10.0)
        * m ** (-7.0 / 10.0)
        * mdot ** (11.0 / 12.0)
        * r_rs ** (-15.0 / 8.0)
        * (
            1.0 - np.sqrt(3.0 / r_rs)
        ) ** (11.0 / 20.0))
