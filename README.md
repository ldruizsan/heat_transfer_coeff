# Heat Transfer Coefficient Calculator

This Python script provides a flexible framework for calculating the average heat transfer coefficient (`h`) for various geometries and flow conditions, primarily focusing on convection. It uses a modular design where specific Nusselt number correlations are handled by dedicated functions, making it extensible for new scenarios.

## Core Functionality

The primary function for users is `calculate_heat_transfer_coefficient`. This function acts as a dispatcher, selecting the appropriate Nusselt number correlation based on the specified geometry, flow type, and thermal boundary condition.

## Main Function: `calculate_heat_transfer_coefficient`

This function calculates the average heat transfer coefficient.
```python
def calculate_heat_transfer_coefficient(
    rho, u, L_char, mu, cp, k,
    geometry_type,
    flow_type="forced",
    thermal_condition="isothermal",
    **kwargs
    ):
# ... implementation ...
```

### Parameters:

*   `rho` (float): Fluid density in kg/m³.
*   `u` (float): Fluid velocity in m/s. For natural convection, this might be set to 0 or not directly used by the Nusselt correlation, which would rely on Grashof/Rayleigh numbers.
*   `L_char` (float): Characteristic length in meters. The definition of `L_char` depends on the geometry:
    *   For a flat plate: Length of the plate in the direction of flow.
    *   For a cylinder in crossflow: Diameter of the cylinder.
    *   For internal pipe flow: Hydraulic diameter of the pipe.
*   `mu` (float): Dynamic viscosity of the fluid in Pa.s (or kg/(m.s)).
*   `cp` (float): Specific heat capacity of the fluid at constant pressure in J/(kg.K).
*   `k` (float): Thermal conductivity of the fluid in W/(m.K).
*   `geometry_type` (str): A string specifying the geometry. Supported values are defined as constants (see Constants section below).
    *   Example: `GEOMETRY_FLAT_PLATE`, `GEOMETRY_CYLINDER_CROSSFLOW`.
*   `flow_type` (str, optional): Specifies the type of convection. Defaults to `"forced"`.
    *   `"forced"`: Forced convection (driven by external means like a fan or pump).
    *   `"natural"`: Natural convection (driven by buoyancy forces due to temperature differences). (Note: Full implementation for natural convection requires Grashof number calculation and relevant correlations).
*   `thermal_condition` (str, optional): Specifies the thermal boundary condition at the surface. Defaults to `"isothermal"`. Supported values are defined as constants (see Constants section below).
    *   `"isothermal"`: The surface is maintained at a constant temperature.
    *   `"constant_heat_flux"`: A uniform heat flux is applied to the surface.
*   `**kwargs` (dict, optional): Additional keyword arguments that might be required by specific Nusselt number correlations. Examples:
    *   `rec_crit` (float): Critical Reynolds number for transition on a flat plate. Defaults to `500_000` if not provided to `_nusselt_flat_plate_forced` via this mechanism.
    *   For natural convection (if implemented): `temp_surface` (surface temperature), `temp_bulk` (bulk fluid temperature), `beta` (volumetric thermal expansion coefficient).

### Returns:

*   `float`: The calculated average heat transfer coefficient (`h`) in W/(m².K).
*   `None`: If the Nusselt number cannot be determined (e.g., unsupported conditions, parameters out of correlation range, or missing handler).

## Internal Nusselt Handler Functions

The script uses internal "handler" functions to calculate the Nusselt number (`Nu`) for specific scenarios. These are typically not called directly by the end-user but are invoked by `calculate_heat_transfer_coefficient`.

Examples of such functions include:

*   `_nusselt_flat_plate_forced(Re, Pr, L, thermal_condition, rec_crit)`: Calculates `Nu` for forced convection over a flat plate.
*   `_nusselt_cylinder_crossflow_forced(Re, Pr, L_char, thermal_condition)`: Calculates `Nu` for forced convection over a cylinder in crossflow.

These functions take Reynolds number (`Re`), Prandtl number (`Pr`), characteristic length (`L` or `L_char`), and the `thermal_condition` as primary arguments, along with any other specific parameters (like `rec_crit`). They return the calculated average Nusselt number or `None`.

## Constants

The script defines string constants for specifying geometry, flow type, and thermal conditions to improve code readability and prevent errors from typos.

### Geometry Types:
*   `GEOMETRY_FLAT_PLATE = "flat_plate"`
*   `GEOMETRY_CYLINDER_CROSSFLOW = "cylinder_crossflow"`
    *(Add more as they are implemented)*

### Thermal Conditions:
*   `THERMAL_ISOTHERMAL = "isothermal"`
*   `THERMAL_CONSTANT_HEAT_FLUX = "constant_heat_flux"`

## Sample Scenarios & Usage

Below are examples demonstrating how to use the `calculate_heat_transfer_coefficient` function.

```python
import math # Only if used directly in user script, not needed for app_with_fn_handler.py usage

# Assuming app_with_fn_handler.py is in the same directory or accessible in PYTHONPATH
from app_with_fn_handler import (
    calculate_heat_transfer_coefficient,
    GEOMETRY_FLAT_PLATE,
    GEOMETRY_CYLINDER_CROSSFLOW,
    THERMAL_ISOTHERMAL,
    THERMAL_CONSTANT_HEAT_FLUX
    )

# Fluid properties for Air at approximately 27°C (300K)
rho_air = 1.177    # Density (kg/m^3)
mu_air = 1.846e-5  # Dynamic viscosity (Pa.s)
cp_air = 1007      # Specific heat capacity (J/(kg.K))
k_air = 0.0262     # Thermal conductivity (W/(m.K))

print("--- Heat Transfer Coefficient Calculator Examples ---")

# --- Scenario 1: Forced Convection over an Isothermal Flat Plate ---
print("\n--- Scenario 1: Flat Plate (Isothermal) ---")
L_plate = 1.0  # Characteristic length: length of the plate (m)
u_plate_flow = 5.0  # Fluid velocity (m/s)

h_plate_isothermal = calculate_heat_transfer_coefficient(
    rho=rho_air,
    u=u_plate_flow,
    L_char=L_plate,
    mu=mu_air,
    cp=cp_air,
    k=k_air,
    geometry_type=GEOMETRY_FLAT_PLATE,
    flow_type="forced",  # Default, can be omitted
    thermal_condition=THERMAL_ISOTHERMAL # Default, can be omitted
    # For flat plate, can also pass rec_crit if a different value is desired:
    # rec_crit=600000
    )

if h_plate_isothermal is not None:
    print(f"Flat Plate (Isothermal): Average h = {h_plate_isothermal:.4f} W/(m^2.K)")
else:
    print("Flat Plate (Isothermal): Could not calculate h.")

# --- Scenario 2: Forced Convection over a Flat Plate with Uniform Heat Flux ---
print("\n--- Scenario 2: Flat Plate (Uniform Heat Flux) ---")
# Using same fluid properties, plate length, and velocity as Scenario 1

h_plate_heatflux = calculate_heat_transfer_coefficient(
    rho=rho_air,
    u=u_plate_flow,
    L_char=L_plate,
    mu=mu_air,
    cp=cp_air,
    k=k_air,
    geometry_type=GEOMETRY_FLAT_PLATE,
    thermal_condition=THERMAL_CONSTANT_HEAT_FLUX
    )

if h_plate_heatflux is not None:
    print(f"Flat Plate (Uniform Heat Flux): Average h = {h_plate_heatflux:.4f} W/(m^2.K)")
else:
    print("Flat Plate (Uniform Heat Flux): Could not calculate h.")

# --- Scenario 3: Forced Convection over an Isothermal Cylinder in Crossflow ---
print("\n--- Scenario 3: Cylinder in Crossflow (Isothermal) ---")
D_cylinder = 0.05  # Characteristic length: diameter of the cylinder (m)
u_cylinder_flow = 10.0  # Fluid velocity (m/s)

h_cylinder_isothermal = calculate_heat_transfer_coefficient(
    rho=rho_air,
    u=u_cylinder_flow,
    L_char=D_cylinder,
    mu=mu_air,
    cp=cp_air,
    k=k_air,
    geometry_type=GEOMETRY_CYLINDER_CROSSFLOW,
    thermal_condition=THERMAL_ISOTHERMAL
    )

if h_cylinder_isothermal is not None:
    print(f"Cylinder (Isothermal, Crossflow): Average h = {h_cylinder_isothermal:.4f} W/(m^2.K)")
else:
    # The _nusselt_cylinder_crossflow_forced in the example code is a placeholder
    # and might return None if specific Re/Pr ranges are not covered.
    print("Cylinder (Isothermal, Crossflow): h could not be calculated (correlation might be missing for these specific Re/Pr conditions).")

```

## How to Expand

The calculator is designed to be extensible:

1.  **Define a New Nusselt Handler Function:**
    Create a new Python function (e.g., `_nusselt_sphere_forced(Re, Pr, L_char, thermal_condition, **kwargs)`) that implements the Nusselt number correlations for the new geometry/condition. This function should calculate and return the average Nusselt number.

2.  **Add to `NUSSELT_HANDLERS` Dictionary:**
    In `app_with_fn_handler.py`, add an entry to the `NUSSELT_HANDLERS` dictionary. The key should be a tuple `(GEOMETRY_CONSTANT, "flow_type_string")`, and the value should be the reference to your new handler function.

    ```python
    # Example:
    GEOMETRY_SPHERE = "sphere" # Define new constant

    NUSSELT_HANDLERS = {
        (GEOMETRY_FLAT_PLATE, "forced"): _nusselt_flat_plate_forced,
        (GEOMETRY_CYLINDER_CROSSFLOW, "forced"): _nusselt_cylinder_crossflow_forced,
        (GEOMETRY_SPHERE, "forced"): _nusselt_sphere_forced, # Add new entry
        # ... other handlers
    }
    ```

3.  **Implement Correlations:**
    Ensure your new handler function correctly implements well-established Nusselt number correlations for various flow regimes (laminar, transitional, turbulent) and Prandtl number ranges applicable to the new scenario. Cite sources for correlations if possible within comments.

## Prerequisites

*   Python 3.x
*   The standard Python `math` module is used internally. No external libraries are required for the core functionality as presented.

## Notes on Correlations

*   The accuracy of the calculated heat transfer coefficient heavily depends on the chosen Nusselt number correlations and their range of applicability (Reynolds number, Prandtl number, Grashof number for natural convection, etc.).
*   Always refer to heat transfer textbooks (e.g., Incropera & DeWitt, Cengel & Ghajar, Holman) for appropriate correlations and their limitations.
*   Fluid properties (`rho`, `mu`, `cp`, `k`, `Pr`) are often evaluated at the film temperature (`T_film = (T_surface + T_bulk) / 2`). For high accuracy, especially when properties vary significantly with temperature, this might require an iterative solution if the surface temperature itself depends on `h`. The current script assumes properties are provided at a relevant mean temperature.



This README should provide a good starting point for users to understand and use your `app_with_fn_handler.py` script. Remember to replace placeholder comments or extend sections like "Supported Geometries and Conditions" as you further develop the script.
