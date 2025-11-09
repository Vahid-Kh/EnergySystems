import tkinter as tk
from tkinter import ttk, messagebox

# SI material allowable stresses in MPa
MATERIALS = {
    'SS 304 (Welded Pipe)': 110,      # MPa
    'SS 304 (Plate/Seamless)': 129,   # MPa
    'SS 316 (Welded Pipe)': 110,      # MPa
    'SS 316 (Plate/Seamless)': 129,   # MPa
    # Add more as needed
}

def solve_barlows(known, vals, sf):
    P, S, T, D = [vals.get(x) for x in ['P', 'S', 'T', 'D']]
    if sf is None or sf == 0: sf = 1
    try:
        # All values in MPa, mm
        if P is None:
            return (2 * S * T) / (D * sf)
        if S is None:
            return (P * D * sf) / (2 * T)
        if T is None:
            return (P * D * sf) / (2 * S)
        if D is None:
            return (2 * S * T) / (P * sf)
    except (ZeroDivisionError, TypeError):
        return None

def calculate():
    missing = [v.get() == '' for v in [P_var, S_var, T_var, D_var]]
    if sum(missing) != 1:
        messagebox.showerror("Input Error", "Leave exactly one variable blank.")
        return
    vals = {
        'P': None if P_var.get() == '' else float(P_var.get()),
        'S': None if S_var.get() == '' else float(S_var.get()),
        'T': None if T_var.get() == '' else float(T_var.get()),
        'D': None if D_var.get() == '' else float(D_var.get()),
    }
    sf = None if SF_var.get() == '' else float(SF_var.get())

    # Fill in the allowable stress from the material, unless S is being solved for
    if vals['S'] is None and material_var.get() in MATERIALS:
        vals['S'] = MATERIALS[material_var.get()]
    elif vals['S'] is not None and material_var.get() in MATERIALS:
        vals['S'] = MATERIALS[material_var.get()]
    result = solve_barlows(missing.index(True), vals, sf)
    if result is not None:
        # Decide format of result
        idx = missing.index(True)
        units = ["MPa", "MPa", "mm", "mm"]
        result_var.set(f"{result:.3f} {units[idx]}")
    else:
        result_var.set("Error in calculation.")

root = tk.Tk()
root.title("Barlow’s Formula Calculator (SI Units)")

mainframe = ttk.Frame(root, padding="10")
mainframe.grid(row=0, column=0, sticky='NSEW')

labels = [
    "P (pressure, MPa)",
    "S (allowable stress, MPa)",
    "t (wall thickness, mm)",
    "D (outside diameter, mm)",
    "Safety Factor (optional)"
]
vars = [tk.StringVar() for _ in range(5)]
P_var, S_var, T_var, D_var, SF_var = vars

for i, (label, var) in enumerate(zip(labels, vars)):
    ttk.Label(mainframe, text=label).grid(row=i, column=0, sticky='W', pady=2)
    ttk.Entry(mainframe, textvariable=var).grid(row=i, column=1, pady=2)

material_var = tk.StringVar(value="SS 304 (Welded Pipe)")
ttk.Label(mainframe, text="Material:").grid(row=5, column=0, sticky='W', pady=2)
ttk.Combobox(mainframe, textvariable=material_var,
             values=list(MATERIALS.keys()),
             state='readonly').grid(row=5, column=1, pady=2)

ttk.Label(mainframe, text="Result:").grid(row=6, column=0, sticky='W', pady=5)
result_var = tk.StringVar()
ttk.Label(mainframe, textvariable=result_var, foreground='blue').grid(row=6, column=1, pady=5)

ttk.Button(mainframe, text="Calculate", command=calculate).grid(row=7, column=0, columnspan=2, pady=10)

root.mainloop()
