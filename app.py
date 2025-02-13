import math

def forced_convection_flat_plate(rho, u, L, mu, cp, k, heat_flux=False):
    """
    Calculates the heat transfer coefficient for forced convection over a flat plate.
    Determines laminar or turbulent flow based on the Reynolds number. For this geometry,
    the critical Reynolds number is around 500,000. In the case of mixed boundary layer conditions,
    that is, flow is in the transition zone close to the critical Reynolds number. This function will
    assume that transition zone to be within 1% of the critical Reynolds number.

    Assumptions:
        1. Isothermal conditions. If you want to instead consider a flat plate with uniform heat flux, set the heat_flux to True when calling the function.
        2. dp/dx = 0 
        3. Constant fluid properties
        4. Viscous dissipation is negligible

    Args:
        rho: Fluid density (kg/m^3)
        u: Fluid velocity (m/s)
        L: Characteristic length (m)
        mu: Dynamic viscosity (Pa.s or kg/(m.s))
        cp: Specific heat at constant pressure (J/(kg.K))
        k: Thermal conductivity (W/(m.K))

    Returns:
        h: Heat transfer coefficient (W/(m^2.K)) or None if the Reynolds number is in a transitional range.
    """

    Re = (rho * u * L) / mu
    Pr = (cp * mu) / k
    Rec = 500_000
    if heat_flux==False:
        if Re <= Rec*0.99 and Pr > 0.6:  # Laminar
            Nu = 0.664 * Re**0.5 * Pr**(1/3) # Local Nusselt number. The average Nusselt for these conditions is twice this amount
        elif Re >= Rec*1.0 and 0.6 < Pr < 60:
            Nu = 0.023 * Re**0.8 * Pr**0.4
        elif Rec*0.99 < Re < Rec*1.01:
            Nu = Pr**(1/3)*(0.664*Re**(1/2) + 0.037*(Re**(4/5)))
        else:
            print("Encountered an error.")
            return None
    
    if heat_flux==True:
        if Re < Rec and Pr > 0.6:
            Nu = 0.453*Re**(1/2)*Pr*(1/3)
        elif Re >= Rec and 0.6 < Pr < 60:
            Nu = 0.0308*Re**(4/5)*Pr**(1/3)


    h = (Nu * k) / L
    return h


# Example usage:
rho = 1.2  # Air density at room temperature (kg/m^3)
u = 5.0  # Velocity (m/s)
L = 1.0  # Length (m)
mu = 1.8e-5  # Dynamic viscosity of air at room temperature (Pa.s)
cp = 1005.0  # Specific heat of air at constant pressure (J/(kg.K))
k = 0.026  # Thermal conductivity of air at room temperature (W/(m.K))

h = forced_convection_flat_plate(rho, u, L, mu, cp, k)
print(f"Heat transfer coefficient: {round(h,4)} W/(m^2.K)")
