"""bASIC

Boundary Conditions

FWS/CW Supply temperature: 15–25 [°C]
FWS/CW Return temperature: 25–40 [°C]
TCS/CW Supply, cold plates: 30–45 [°C]
TCS/CW Return, cold plates: 40–60 [°C]

FWS Design Pressure: PN 6, 10 or 16 [Bar] – 10 Bar being common.
TCS Design Pressure: PN 2–3 [Bar] range is common.

FWS fluid media, type: 25–30% Polypropylene- or ethylene glycol
TCS fluid media, type: Deionized water + additive or single-phase coolants like: Shell XS5, Lubrizol IM20, etc.

"""


from tabulate import tabulate
import math

# Try to import CoolProp
try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    PropsSI = None

# Yield strengths for materials
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

Piping_Costing = {
    "material": 62.5,
    "installation": 17.5,
    "flush&balance":17,
    "other":3
    }
PipeLenghtPerMW = ((15*2+16+5)*2)/2 # Pipe length forward and return per MW, in drawing we have 2.2 MW roughly

def fetch_yield_strength(material):
    return YIELD_STRENGTHS.get(material)

def calculate_required_thickness(material, pressure_bar, phase, outer_diameter_mm):
    allowable_stress = fetch_yield_strength(material)
    if allowable_stress is None:
        return "Material not found"
    safety_factor = 3 if phase == '1Ø' else 5
    pressure_mpa = pressure_bar / 10  # Convert bar to MPa
    required_thickness_mm = (pressure_mpa * outer_diameter_mm) / (2 * (allowable_stress / safety_factor))
    return round(required_thickness_mm, 3)

class Refrigerant:
    def __init__(self, name, gwp, phase, pfas, Vlimit, material="SS304"):
        self.name = name
        self.gwp = gwp
        self.phase = phase
        self.pfas = pfas
        self.Vlimit = Vlimit
        self.material = material
        self.psat_25C = self.calculate_psat_25C()
        self.pressure_used_bar = self.determine_pressure_used()
        self.kW_per_m2 = self.calc_kW_per_m2()
        self.area_for_1MW_m2 = self.calculate_area_for_1MW()
        self.pipe_diameter_mm = self.calculate_pipe_diameter()
        self.pipe_thickness_mm = self.calculate_pipe_thickness()

    def calculate_psat_25C(self):
        if PropsSI is None:
            return "CoolProp not available"
        try:
            pressure_pa = PropsSI("P", "T", 25 + 273.15, "Q", 0, self.name)
            return round(pressure_pa / 1e5, 2)
        except:
            return "N/A"

    def determine_pressure_used(self):
        if self.phase == "1Ø":
            return 3.0  # bar
        return self.psat_25C if isinstance(self.psat_25C, (int, float)) else 1.0

    def calc_kW_per_m2(self):
        if PropsSI is None:
            return "CoolProp not available"
        try:
            T = 25 + 273.15
            dT = 10
            Q = 0.75
            P_amb = 1e5
            area = 1
            if self.phase == "1Ø" and self.name in ["INCOMP::MEG-25%", "Water"]:
                cp = PropsSI("C", "T", T + dT / 2, "P", P_amb, self.name)
                rho = PropsSI("D", "T", T, "P", P_amb, self.name)
                heat_flow = area * self.Vlimit * rho * cp * dT
                return round(heat_flow / 1000, 3)

            if self.phase == "2Ø":
                rho = PropsSI("D", "T", T, "Q", Q, self.name)
                volumetric_flow = self.Vlimit * area
                mass_flow = rho * volumetric_flow
                delta_h = PropsSI("H", "T", T, "Q", Q, self.name) - PropsSI("H", "T", T, "Q", 0.01, self.name)
                cooling_power = mass_flow * delta_h
                return round(cooling_power / 1000, 3)

        except Exception as e:
            return f"Error: {e}"

    def calculate_area_for_1MW(self):
        try:
            if isinstance(self.kW_per_m2, (int, float)) and self.kW_per_m2 > 0:
                return round(1000 / self.kW_per_m2, 3)
            else:
                return "N/A"
        except:
            return "N/A"

    def calculate_pipe_diameter(self):
        try:
            if isinstance(self.area_for_1MW_m2, (int, float)):
                diameter_m = math.sqrt(4 * self.area_for_1MW_m2 / math.pi)
                return round(diameter_m * 1000, 2)  # Convert to mm
            else:
                return "N/A"
        except:
            return "N/A"

    def calculate_pipe_thickness(self):
        try:
            if isinstance(self.pipe_diameter_mm, (int, float)):
                return calculate_required_thickness(self.material, self.pressure_used_bar, self.phase, self.pipe_diameter_mm)
            else:
                return "N/A"
        except Exception as e:
            return f"Error: {e}"

# Create refrigerant objects
refrigerants = [
    Refrigerant("CO2", 1, "2Ø", "No", 20.0),
    Refrigerant("Ammonia", 1.37, "2Ø", "No", 25.0),
    Refrigerant("R290", 0.02, "2Ø", "No", 18.0),
    Refrigerant("R152a", 124, "2Ø", "Yes", 12.0),
    Refrigerant("R1234ze(E)", 1.37, "2Ø", "Yes", 10.0),
    Refrigerant("Water", 0, "2Ø", "No", 30.0),
    Refrigerant('INCOMP::MEG-25%', 0, "1Ø", "No", 2.0),
    Refrigerant("Water", 0, "1Ø", "No", 2.0)
]

# Prepare data for tabulation
table_data = [[
    r.name, r.gwp, r.phase, r.pfas, r.psat_25C, r.Vlimit, r.kW_per_m2,
    r.pressure_used_bar, r.area_for_1MW_m2, r.pipe_diameter_mm, r.pipe_thickness_mm
] for r in refrigerants]

headers = [
    "Media", "t CO2 eq.", "2Ø / 1Ø", "PFAS", "Psat @25°C (bar)",
    "Vlimit", "kW/m²", "Pressure Used (bar)", "Area for 1MW (m²)",
    "Pipe Diameter (mm)", "Pipe Thickness (mm)"
]

# Print the table
print(tabulate(table_data, headers=headers, tablefmt="grid", colalign=("center",)*len(headers)))
