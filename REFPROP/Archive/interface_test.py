from refprop_interface import RefpropFluid

fluid = "R1234zee"
rp = RefpropFluid(fluid, rpprefix=r"C:\Program Files (x86)\REFPROP")
rp.TQ(50, 0.0)

print("rho =", rp.rho, "kg/m^3")
print("a   =", rp.a, "m/s")
print("mu  =", rp.mu, "Pa·s")
print(rp.state())   # now includes Z (float in single-phase, None in two-phase)
