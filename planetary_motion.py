import plm_3_1
import numpy as np

def curve2order(x, y):
    A = np.column_stack([x**2, x*y, y**2, x, y, np.ones_like(x)])
    _, _, Vt = np.linalg.svd(A)
    coeffs = Vt[-1, :]
    
    A, B, C, D, E, F = coeffs
    x0 = (C*D - B*E) / (B**2 - A*C)
    y0 = (A*E - B*D) / (B**2 - A*C)
    
    beta = 0.5 * np.arctan2(B, A - C)
    cos_beta = np.cos(beta)
    
    num = 2 * (A*E**2 + C*D**2 - B*D*E + (B**2 - A*C) * F)
    den1 = (B**2 - A*C) * ((A + C) + np.sqrt((A - C)**2 + B**2))
    den2 = (B**2 - A*C) * ((A + C) - np.sqrt((A - C)**2 + B**2))
    
    a = np.sqrt(abs(num / den1)) if den1 != 0 else np.nan
    b = np.sqrt(abs(num / den2)) if den2 != 0 else np.nan
    
    curve_type = 1 if A*C > 0 else -1
    
    return a, b, x0, y0, cos_beta, curve_type

x0, y0 = 1.0e15, 0
v0x, v0y = 0, 2e7
Mc, m = 1e40, 2e39
dt = 60*60*24
N = 3000 # make sure N is enough to close your ellipse, otherwise the calculations may be wrong
x, y, t = plm_3_1.planet2D(x0, y0, v0x, v0y, Mc, m, dt, N)

a, b, x0, y0, cos_beta, curve_type = curve2order(x, y)
print(f"a = {a:.3e}, b = {b:.3e}, x0 = {x0:.3e}, y0 = {y0:.3e}, cos(beta) = {cos_beta:.3e}, type = {curve_type}")