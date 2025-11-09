from tabulate import tabulate

# Try to import CoolProp
try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    PropsSI = None

# Define a Refrigerant class
class Refrigerant:
    def __init__(self, name, gwp, phase, pfas, Vlimit):
        self.name = name
        self.gwp = gwp
        self.phase = phase
        self.pfas = pfas
        self.psat_25C = self.calculate_psat_25C()
        self.Vlimit = Vlimit
        self.kW_LPM = self.calc_kW_LPM()

    def calculate_psat_25C(self):
        P_amb = 1e5  # Ambient pressure in Pa
        if self.name == "1Ø":
            return round(P_amb / 1e5, 2)  # Return ambient pressure in bar

        if PropsSI is None:
            return "CoolProp not available"
        try:
            pressure_pa = PropsSI("P", "T", 25 + 273.15, "Q", 0, self.name)
            return round(pressure_pa / 1e5, 2)  # Convert to bar
        except:
            return "N/A"

    def calc_kW_LPM(self):
        import math
        if PropsSI is None:
            return "CoolProp not available"
        try:
            # Get density at 25°C and quality = 0.75
            T = 25 + 273.15  # Temperature in Kelvin
            dT = 10  #Kelvin
            Q = 0.75  # Quality (75% vapor)
            P_amb = 1e5  # Ambient pressure in Pa
            area = 1

            if self.phase == "1Ø" and self.name in ["INCOMP::MEG-25%", "Water"]:
                cp_water = PropsSI("C", "T", T + dT / 2, "P", P_amb, self.name)
                rho_water = PropsSI("D", "T", T, "P", P_amb, self.name)
                heat_flow_per_area = area * self.Vlimit * rho_water * cp_water * dT
                return round(heat_flow_per_area/1000 , 3)

            if self.phase == "2Ø":
                diameter_m = math.sqrt(4 * area / math.pi)
                rho = PropsSI("D", "T", T, "Q", Q, self.name)  # Density in kg/m³

                # Convert Vlimit (m/s) and density to volumetric flow rate in L/min
                # Assume 1 cm² = 1e-4 m² cross-sectional area
                volumetric_flow_m3_s = self.Vlimit * area  # m³/s

                # Mass flow rate = density * volumetric flow
                mass_flow_kg_s = rho * volumetric_flow_m3_s

                # Cooling capacity = mass flow * enthalpy difference (assume Δh = 200 kJ/kg as a placeholder)
                delta_h = PropsSI("H", "T", T, "Q", Q, self.name)-PropsSI("H", "T", T, "Q", 0.01, self.name)  # J/kg (placeholder, ideally should be calculated)
                cooling_power_W = mass_flow_kg_s * delta_h
                cooling_power_kW = cooling_power_W / 1000

                # Return cooling power per LPM
                return round(cooling_power_kW , 3)

        except Exception as e:
            return f"Error: {e}"




# Create refrigerant objects with PFAS classification

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
table_data = [[r.name, r.gwp, r.phase, r.pfas, r.psat_25C, r.Vlimit, r.kW_LPM] for r in refrigerants]
headers = ["Media", "t CO2 eq.", "2Ø / 1Ø", "PFAS", "Psat @25°C (bar)", "Vlimit", "kW/1m^2", ]

# Print the table
print(tabulate(table_data, headers=headers, tablefmt="grid", colalign=("center", "center", "center", "center", "center", "center", "center")))
