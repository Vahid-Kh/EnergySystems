YIELD_STRENGTHS = {
    "SS304": 215,  # MPa
    "SS316": 205,
    "Aluminum 6061": 276,
    "Copper": 70,
    "Polyethylene": 20,
    "Polypropylene": 25,
    "Nylon 6": 45,
    "PVC": 52,
}

def fetch_yield_strength(material):
    """Fetch the yield strength for the given material."""
    return YIELD_STRENGTHS.get(material)

def calculate_required_thickness(material, pressure, phase, outer_diameter_mm):
    """
    Calculate the required wall thickness based on material, pressure, phase, and outer diameter.

    Parameters:
    - material: str, the material of the pipe
    - pressure: float, the internal pressure in MPa
    - phase: str, either 'single' or 'two'
    - outer_diameter_mm: float, the outer diameter in millimeters
    Returns:
    - required_thickness_mm: float, the required wall thickness in millimeters
    """
    allowable_stress = fetch_yield_strength(material)
    if allowable_stress is None:
        raise ValueError("Material not found or yield strength not available.")

    safety_factor = 3 if phase == 'single' else 5

    # Barlow's formula: t = (P * D) / (2 * S / SF)
    required_thickness_mm = (pressure * outer_diameter_mm) / (2 * (allowable_stress / safety_factor))

    return required_thickness_mm

# Example usage:
try:
    material = "SS304"
    pressure = 50.0  # MPa
    phase = "single"
    outer_diameter_mm = 1000  # mm

    required_thickness = calculate_required_thickness(material, pressure, phase, outer_diameter_mm)
    print(f"Required wall thickness for {material} at {pressure} MPa: {required_thickness:.3f} mm")
except ValueError as e:
    print(e)
