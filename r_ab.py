from scipy.optimize import brentq


def r_ab_equation(r, alpha, m, mdot):
    """
    Shakura-Sunyaev equation for the a/b boundary.

    r is the Shakura-Sunyaev dimensionless radius:
        r = R / (3 R_S)
    """

    RHS = (
        150.0
        * (alpha * m) ** (2.0 / 21.0)
        * mdot ** (16.0 / 21.0)
    )

    LHS = (
        r / (1.0 - r ** (-0.5)) ** (16.0 / 21.0)
    )

    return LHS - RHS


def calculate_r_ab(alpha, m, mdot):
    """
    Return both mathematical roots, converted to R/R_S.

    Returns
    -------
    root1_Rs, root2_Rs
        root1 = inner mathematical root
        root2 = outer/physical a-b boundary
    """

    if alpha <= 0:
        raise ValueError("alpha must be greater than zero.")

    if m <= 0:
        raise ValueError("m must be greater than zero.")

    if mdot <= 0:
        raise ValueError("mdot must be greater than zero.")

    r_min = (29.0 / 21.0) ** 2

    root1 = brentq(
        lambda r: r_ab_equation(r, alpha, m, mdot),
        1.000001,
        r_min
    )

    root2 = brentq(
        lambda r: r_ab_equation(r, alpha, m, mdot),
        r_min,
        1.0e8
    )

    # SS coordinate r = R/(3 Rs) → app coordinate R/Rs
    root1_Rs = 3.0 * root1
    root2_Rs = 3.0 * root2

    return root1_Rs, root2_Rs    
##############################################################################################################################################################################################################################################################################


