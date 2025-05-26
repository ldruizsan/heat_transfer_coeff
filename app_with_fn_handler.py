import math

# --- Fluid Property Data (Example - could be from a library or more complex structure) ---
# For simplicity, we'll pass them directly as in your original example.

# --- Constants ---
# General constants
PI = math.pi

# Geometry-specific critical Reynolds numbers (examples)
REC_CRIT_FLAT_PLATE = 500_000
REC_CRIT_CYLINDER_LAMINAR_MAX = 200_000 # Example, varies with source

# --- Nusselt Number Correlations for Specific Geometries/Conditions ---

def _nusselt_flat_plate_forced_detailed(Re, Pr, L_char, thermal_condition="isothermal", rec_crit=REC_CRIT_FLAT_PLATE, **kwargs):
    """
    Calculates average Nusselt number for forced convection over a flat plate
    and returns detailed information about the correlation used.
    """
    result = {
        "Nu_value": None,
        "correlation_name": "N/A",
        "correlation_formula_str": "N/A",
        "conditions_applied": [
            f"Geometry: Flat Plate, Length = {L_char}m",
            f"Reynolds Number (Re): {Re:.2e}",
            f"Prandtl Number (Pr): {Pr:.2f}"
        ],
        "explanation_notes": [],
        "warnings": []
    }

    re_laminar_max = rec_crit * 0.99
    re_turbulent_min = rec_crit * 1.01
    result["conditions_applied"].append(f"Critical Re for transition assumed around: {rec_crit:.2e} (Transition zone: {re_laminar_max:.2e} - {re_turbulent_min:.2e})")

    if thermal_condition == "isothermal":
        result["conditions_applied"].append("Thermal Condition: Isothermal Plate")
        if Re < re_laminar_max:  # Laminar
            result["conditions_applied"].append(f"Flow Regime: Laminar (Re < {re_laminar_max:.2e})")
            if Pr >= 0.6:
                result["Nu_value"] = 0.664 * Re**0.5 * Pr**(1/3)
                result["correlation_name"] = "Average Nusselt Number for Laminar Flow over Isothermal Flat Plate"
                result["correlation_formula_str"] = "Nu_L = 0.664 * Re_L^(1/2) * Pr^(1/3)"
                result["conditions_applied"].append("Pr >= 0.6")
                result["explanation_notes"].append("This is a standard correlation for average heat transfer over an isothermal flat plate in laminar flow.")
                result["explanation_notes"].append("Assumes constant fluid properties and negligible viscous dissipation.")
            else:
                result["warnings"].append(f"Prandtl number {Pr:.2f} is out of range (Pr >= 0.6) for this laminar isothermal correlation.")
        elif re_laminar_max <= Re <= re_turbulent_min:  # Transition (user's specific)
            result["conditions_applied"].append(f"Flow Regime: Transitional ({re_laminar_max:.2e} <= Re <= {re_turbulent_min:.2e})")
            if Pr >= 0.6:
                result["Nu_value"] = Pr**(1/3) * (0.664 * Re**0.5 + 0.037 * Re**(4/5)) # User's original formula
                result["correlation_name"] = "User-defined Transition Flow Correlation for Isothermal Flat Plate"
                result["correlation_formula_str"] = "Nu_L = Pr^(1/3) * (0.664 * Re_L^(1/2) + 0.037 * Re_L^(4/5))"
                result["conditions_applied"].append("Pr >= 0.6 (assumed)")
                result["explanation_notes"].append("This correlation is based on the user's provided formula for the defined narrow transition zone.")
                result["explanation_notes"].append("Note: Standard mixed-flow correlations often cover a broader range after the critical Reynolds number.")
            else:
                result["warnings"].append(f"Prandtl number {Pr:.2f} is out of range (Pr >= 0.6 assumed) for this transitional isothermal correlation.")
        elif Re > re_turbulent_min:  # Turbulent
            result["conditions_applied"].append(f"Flow Regime: Turbulent (Re > {re_turbulent_min:.2e})")
            if 0.6 <= Pr <= 60:
                result["Nu_value"] = 0.023 * Re**0.8 * Pr**0.4 # User's correlation
                result["correlation_name"] = "User-defined Turbulent Flow Correlation for Isothermal Flat Plate (Dittus-Boelter like)"
                result["correlation_formula_str"] = "Nu_L = 0.023 * Re_L^(0.8) * Pr^(0.4)"
                result["conditions_applied"].append("0.6 <= Pr <= 60")
                result["explanation_notes"].append("This correlation is based on the user's provided formula, similar in form to Dittus-Boelter (often for internal flow).")
                result["explanation_notes"].append("For external turbulent flow over a flat plate, Nu_L = 0.0296 * Re_L^0.8 * Pr^(1/3) is also common.")
            else:
                result["warnings"].append(f"Prandtl number {Pr:.2f} is out of range (0.6 <= Pr <= 60) for this turbulent isothermal correlation.")
    elif thermal_condition == "constant_heat_flux":
        result["conditions_applied"].append("Thermal Condition: Uniform Heat Flux Plate")
        if Re < re_laminar_max:  # Laminar
            result["conditions_applied"].append(f"Flow Regime: Laminar (Re < {re_laminar_max:.2e})")
            if Pr >= 0.6:
                result["Nu_value"] = 0.680 * Re**0.5 * Pr**(1/3)
                result["correlation_name"] = "Average Nusselt Number for Laminar Flow over Flat Plate with Uniform Heat Flux"
                result["correlation_formula_str"] = "Nu_L = 0.680 * Re_L^(1/2) * Pr^(1/3)"
                result["conditions_applied"].append("Pr >= 0.6")
                result["explanation_notes"].append("This is a standard correlation for average heat transfer over a flat plate with uniform heat flux in laminar flow.")
            else:
                result["warnings"].append(f"Prandtl number {Pr:.2f} is out of range (Pr >= 0.6) for this laminar constant heat flux correlation.")
        elif re_laminar_max <= Re <= re_turbulent_min: # Transition
            result["conditions_applied"].append(f"Flow Regime: Transitional ({re_laminar_max:.2e} <= Re <= {re_turbulent_min:.2e})")
            result["warnings"].append("No specific correlation is implemented for transitional flow with uniform heat flux in this narrow zone.")
            result["explanation_notes"].append("A mixed boundary layer correlation would typically be used if the flow starts laminar and transitions to turbulent over the plate length.")
        elif Re > re_turbulent_min:  # Turbulent
            result["conditions_applied"].append(f"Flow Regime: Turbulent (Re > {re_turbulent_min:.2e})")
            if 0.6 <= Pr <= 60:
                result["Nu_value"] = 0.0308 * Re**(4/5) * Pr**(1/3)
                result["correlation_name"] = "Average Nusselt Number for Turbulent Flow over Flat Plate with Uniform Heat Flux"
                result["correlation_formula_str"] = "Nu_L = 0.0308 * Re_L^(4/5) * Pr^(1/3)"
                result["conditions_applied"].append("0.6 <= Pr <= 60")
                result["explanation_notes"].append("This correlation is used for turbulent flow over a flat plate with uniform heat flux.")
            else:
                result["warnings"].append(f"Prandtl number {Pr:.2f} is out of range (0.6 <= Pr <= 60) for this turbulent constant heat flux correlation.")

    if result["Nu_value"] is None and not result["warnings"]: # If no Nu and no specific warning yet
        result["warnings"].append(f"Could not determine Nusselt number. Conditions might be outside implemented ranges or a specific correlation might be missing for Re={Re:.2e}, Pr={Pr:.2f}, thermal_condition='{thermal_condition}'.")

    # Clean up conditions_applied if Nu is None to avoid confusion
    if result["Nu_value"] is None:
        result["correlation_name"] = "N/A"
        result["correlation_formula_str"] = "N/A"
        # Keep basic conditions like Re, Pr, Geometry, Thermal Condition
        base_conditions = [c for c in result["conditions_applied"] if any(keyword in c for keyword in ["Geometry:", "Reynolds Number", "Prandtl Number", "Thermal Condition:", "Critical Re"])]
        # Add flow regime if determined
        flow_regime_condition = next((c for c in result["conditions_applied"] if "Flow Regime:" in c), None)
        if flow_regime_condition:
            base_conditions.append(flow_regime_condition)
        result["conditions_applied"] = base_conditions

    return result


def _nusselt_cylinder_crossflow_forced_detailed(Re, Pr, L_char, thermal_condition="isothermal", **kwargs):
    """
    Calculates average Nusselt number for forced convection over a cylinder in crossflow.
    L_char here is the cylinder diameter.
    Returns detailed information.
    """
    result = {
        "Nu_value": None, "correlation_name": "N/A", "correlation_formula_str": "N/A",
        "conditions_applied": [
            f"Geometry: Cylinder in Crossflow, Diameter = {L_char}m",
            f"Reynolds Number (Re): {Re:.2e}",
            f"Prandtl Number (Pr): {Pr:.2f}",
            f"Thermal Condition: {thermal_condition}" # Assuming isothermal for this example
        ],
        "explanation_notes": [], "warnings": []
    }

    # Example: Churchill and Bernstein correlation (valid for wide Re, Pr range)
    # This correlation is generally applicable for Pe = Re*Pr > 0.2
    if Re * Pr > 0.2:
        term1_numerator = 0.62 * (Re**0.5) * (Pr**(1/3))
        term1_denominator = (1 + (0.4 / Pr)**(2/3))**0.25
        term1 = term1_numerator / term1_denominator
        term2_factor = (1 + (Re / 282000)**(5/8))**(4/5)
        
        result["Nu_value"] = 0.3 + term1 * term2_factor
        result["correlation_name"] = "Churchill and Bernstein Correlation"
        result["correlation_formula_str"] = "Nu_D = 0.3 + [0.62*Re_D^(1/2)*Pr^(1/3) / (1+(0.4/Pr)^(2/3))^(1/4)] * [1+(Re_D/282000)^(5/8)]^(4/5)"
        result["conditions_applied"].append("Pe = Re*Pr > 0.2")
        result["explanation_notes"].append("The Churchill and Bernstein correlation is a comprehensive equation valid over a wide range of Re and Pr numbers for flow over a cylinder.")
        result["explanation_notes"].append("It's particularly useful as it covers many flow regimes smoothly.")
    else:
        result["warnings"].append(f"Pe = Re*Pr ({Re*Pr:.2e}) is not > 0.2, Churchill and Bernstein correlation may not be accurate. Consider other low Pe correlations if available.")
        # Placeholder for Hilpert or other low Re correlations if you want to add them as fallbacks
        # For now, we just warn.

    if result["Nu_value"] is None and not result["warnings"]:
         result["warnings"].append(f"Could not determine Nusselt number for Cylinder, Re={Re:.2e}, Pr={Pr:.2f}")

    return result

# --- Main Heat Transfer Coefficient Calculator ---

# Using string constants or Enums for types is good practice
GEOMETRY_FLAT_PLATE = "flat_plate"
GEOMETRY_CYLINDER_CROSSFLOW = "cylinder_crossflow"
# ... other geometries

THERMAL_ISOTHERMAL = "isothermal"
THERMAL_CONSTANT_HEAT_FLUX = "constant_heat_flux"

NUSSELT_HANDLERS = {
    (GEOMETRY_FLAT_PLATE, "forced"): _nusselt_flat_plate_forced_detailed,
    (GEOMETRY_CYLINDER_CROSSFLOW, "forced"): _nusselt_cylinder_crossflow_forced_detailed,
    # Add more handlers here:
    # (GEOMETRY_SPHERE, "forced"): _nusselt_sphere_forced,
    # (GEOMETRY_FLAT_PLATE, "natural"): _nusselt_flat_plate_natural, # Natural convection would use Grashof number
}

def calculate_heat_transfer_coefficient(
    rho, u, L_char, mu, cp, k,
    geometry_type, # String constant e.g. GEOMETRY_FLAT_PLATE
    flow_type="forced", # e.g., "forced", "natural"
    thermal_condition="isothermal", # e.g., "isothermal", "constant_heat_flux"
    **kwargs # For additional parameters like surface temperature if needed by some correlations
):
    """
    Calculates the average heat transfer coefficient for various geometries and flow conditions.
    """
    # Input validation for basic properties
    if not all(isinstance(val, (int, float)) and val > 0 for val in [rho, L_char, mu, cp, k]):
        if u <=0 and flow_type == "forced": # u can be 0 for natural convection, but not for forced
             pass # allow u=0 if not forced, Re will be 0
        else:
            return {
                "h_value": None, "Nu_details": None,
                "error": "Invalid input: rho, u (for forced), L_char, mu, cp, k must be positive numbers."
            }

    Re = (rho * u * L_char) / mu if u > 0 and mu > 0 else 0 # Reynolds number
    Pr = (cp * mu) / k  # Prandtl number

    handler_key = (geometry_type, flow_type)
    nusselt_calculator = NUSSELT_HANDLERS.get(handler_key)

    nu_details = None
    h_value = None
    error_message = None

    if nusselt_calculator:
        if flow_type == "forced":
            nu_details = nusselt_calculator(Re, Pr, L_char, thermal_condition=thermal_condition, **kwargs)
        # elif flow_type == "natural":
        #     g = 9.81 # gravitational acceleration
        #     # beta = 1/T_film # Volumetric thermal expansion coefficient for ideal gas
        #     # delta_T = abs(kwargs.get('temp_surface', T_bulk) - T_bulk) # Temperature difference
        #     # Gr = (g * beta * delta_T * L_char**3) / (mu/rho)**2 # Grashof number
        #     # nu_details = nusselt_calculator(Gr, Pr, L_char, thermal_condition=thermal_condition, **kwargs)
        else:
            error_message = f"Flow type '{flow_type}' not currently supported."
    else:
        error_message = f"No handler found for geometry '{geometry_type}' and flow type '{flow_type}'."

    if nu_details and nu_details.get("Nu_value") is not None:
        h_value = (nu_details["Nu_value"] * k) / L_char

    return {
        "h_value": h_value,
        "Nu_details": nu_details, # This will contain all the explanatory text
        "error": error_message
    }

# --- Example Usage of the Expanded Structure ---
if __name__ == "__main__":
    # Fluid properties (Air at approx 27C/300K)
    rho_air = 1.177  # kg/m^3
    mu_air = 1.846e-5 # Pa.s
    cp_air = 1007    # J/(kg.K)
    k_air = 0.0262   # W/(m.K)

    # Scenario 1: Flat Plate (re-using your example values)
    print("\n--- Flat Plate Example ---")
    L_plate = 1.0  # m
    u_plate_flow = 5.0  # m/s
    result_plate_iso = calculate_heat_transfer_coefficient(
        rho_air, u_plate, L_plate, mu_air, cp_air, k_air,
        geometry_type=GEOMETRY_FLAT_PLATE,
        thermal_condition=THERMAL_ISOTHERMAL
    )
    if result_plate_iso["error"]:
        print(f"Error: {result_plate_iso['error']}")
    elif result_plate_iso["h_value"] is not None:
        print(f"Flat Plate (Isothermal): h = {result_plate_iso['h_value']:.4f} W/(m^2.K)")
        print(f"  Correlation: {result_plate_iso['Nu_details']['correlation_name']}")
        print(f"  Formula: {result_plate_iso['Nu_details']['correlation_formula_str']}")
        print(f"  Explanation: {'; '.join(result_plate_iso['Nu_details']['explanation_notes'])}")
        if result_plate_iso['Nu_details']['warnings']:
            print(f"  Warnings: {'; '.join(result_plate_iso['Nu_details']['warnings'])}")
    else:
        print(f"Could not calculate h for Isothermal Flat Plate.")
        if result_plate_iso['Nu_details'] and result_plate_iso['Nu_details']['warnings']:
            print(f"  Reason/Warnings: {'; '.join(result_plate_iso['Nu_details']['warnings'])}")

    result_plate_q = calculate_heat_transfer_coefficient(
        rho_air, u_plate, L_plate, mu_air, cp_air, k_air,
        geometry_type=GEOMETRY_FLAT_PLATE,
        thermal_condition=THERMAL_CONSTANT_HEAT_FLUX
    )
    # ... (similar printing logic for result_plate_q)

    # Scenario 2: Cylinder in Crossflow (Example)
    print("\n--- Cylinder in Crossflow Example ---")
    D_cylinder = 0.05  # m (characteristic length is diameter)
    u_cylinder = 10.0  # m/s
    
    result_cylinder = calculate_heat_transfer_coefficient(
        rho_air, u_cylinder, D_cylinder, mu_air, cp_air, k_air,
        geometry_type=GEOMETRY_CYLINDER_CROSSFLOW,
        thermal_condition=THERMAL_ISOTHERMAL
        # kwargs like temp_surface might be needed for some correlations
    )
    # ... (similar printing logic for result_cylinder)
