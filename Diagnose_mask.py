#!/usr/bin/env python3
"""Compare the solid mask baked into two _final_velocity.dat files at one x column."""
import numpy as np

NX, NY = 100, 20
X_SLICE = 90

def load_mask_and_vx(case_name):
    data = np.loadtxt(case_name + "_final_velocity.dat", comments="#")
    vx = np.zeros((NY, NX))
    solid = np.zeros((NY, NX), dtype=bool)
    for row in data:
        x, y = int(row[0]), int(row[1])
        vx[y, x] = row[2]
        solid[y, x] = row[4] > 0.5
    return vx[:, X_SLICE], solid[:, X_SLICE]

vx_n, solid_n = load_mask_and_vx("finite_vol_4_syringe_SIMPLE")
vx_pl, solid_pl = load_mask_and_vx("finite_vol_6_syringe_SIMPLE_shearthinning")

print(f"{'y':>3}  {'solid(newt)':>12}  {'solid(pl)':>10}  {'agree':>6}  {'vx(newt)':>10}  {'vx(pl)':>10}")
for y in range(NY):
    agree = "yes" if solid_n[y] == solid_pl[y] else "*** NO ***"
    print(f"{y:>3}  {str(solid_n[y]):>12}  {str(solid_pl[y]):>10}  {agree:>6}  {vx_n[y]:>10.4f}  {vx_pl[y]:>10.4f}")