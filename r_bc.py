from scipy.optimize import brentq


def r_bc_equation(r, mdot):
    """
    Shakura-Sunyaev equation for the b/c boundary.

    r = R / (3 R_S)
    """

    RHS = (
        6.3e3
        * mdot ** (2.0 / 3.0)
        * (1.0 - r ** (-0.5)) ** (2.0 / 3.0)
    )

    return r - RHS


def calculate_r_bc(mdot):
    """
    Return both mathematical roots, converted to R/R_S.
    """

    if mdot <= 0:
        raise ValueError("mdot must be greater than zero.")

    A = 6.3e3 * mdot ** (2.0 / 3.0)

    root1 = None
    root2 = None

    try:
        root1 = brentq(
            lambda r: r_bc_equation(r, mdot),
            1.000001,
            2.0
        )
    except ValueError:
        pass

    try:
        root2 = brentq(
            lambda r: r_bc_equation(r, mdot),
            2.0,
            A
        )
    except ValueError:
        pass

    if root1 is None and root2 is None:
        raise ValueError(
            "No positive roots were found for the given mdot."
        )

    root1_Rs = None if root1 is None else 3.0 * root1
    root2_Rs = None if root2 is None else 3.0 * root2

    return root1_Rs, root2_Rs

##############################################################################################################################################################################################################################################################################
