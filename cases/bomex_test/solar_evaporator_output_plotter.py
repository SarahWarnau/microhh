import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc

import matplotlib.pyplot as plt

def open_default(case, time="0000000"):
    fname = f"{case}.default.{time}.nc"
    groups = [''] + list(nc.Dataset(fname).groups.keys())
    prof_dict = {}
    for group in groups:
        prof_dict[group] = xr.open_dataset(fname, group=group, decode_times=False)
        prof = xr.merge(prof_dict.values())
        prof['time'] = np.round(prof.time.values)
    return prof

ds = open_default("bomex")

# print(ds.variables)
dz = ds.isel(z=0).z.values
# dz = ds.isel(zh=1).zh.values
LE = ds.isel(zh=0).qt_flux*ds.isel(zh=0).rhoh*2.5e6/dz
H = ds.isel(zh=0).thl_flux*ds.isel(zh=0).rhoh*1005/dz
T_atm = ds.isel(z=0)["T"]

thl_atm0 = ds.isel(z=0)["thl"]
thl_atm1 = ds.isel(z=1)["thl"]
thl_atm2 = ds.isel(z=2)["thl"]

thl_tech = ds["thl_bot"]

qt_atm0 = ds.isel(z=0)["qt"]
qt_atm1 = ds.isel(z=1)["qt"]
qt_atm2 = ds.isel(z=2)["qt"]

qt_tech = ds["qt_bot"]

RH0 = ds.isel(z=0)["rh"]
RH1 = ds.isel(z=1)["rh"]
RH2 = ds.isel(z=2)["rh"]
time = ds.time

z0 = np.round(ds.isel(z=0).z.values)
z1 = np.round(ds.isel(z=1).z.values)
z2 = np.round(ds.isel(z=2).z.values)

fig, ax = plt.subplots(2,2, figsize=(12,8))

ax[0,0].set_title("Surface fluxes")
ax[0,0].plot(time, LE, label="LE")
ax[0,0].plot(time, H, label="H")
ax[0,0].plot(time, LE+H, label="LE+H")
ax[0,0].axhline(1000, label="Rnet", c='tab:green', ls=':')
ax[0,0].legend(ncol=4)
ax[0,0].set_ylabel("[W m-2]")

ax[0,1].set_title("Liquid water pot. temp.")
ax[0,1].plot(time, thl_atm0, label=f"thl at z={z0} m", c='k', alpha=0.5)
ax[0,1].plot(time, thl_atm1, label=f"thl at z={z1} m", c='k', alpha=0.3)
ax[0,1].plot(time, thl_atm2, label=f"thl at z={z2} m", c='k', alpha=0.1)
ax[0,1].plot(time, thl_tech, label="thl_tech", c='k')
ax[0,1].legend(ncol=2)
ax[0,1].set_ylabel("[K]")

ax[1,0].set_title("Relative humidity")
ax[1,0].plot(time, RH0*100, label=f'RH at z={z0} m', c='k', alpha=0.5)
ax[1,0].plot(time, RH1*100, label=f'RH at z={z1} m', c='k', alpha=0.3)
ax[1,0].plot(time, RH2*100, label=f'RH at z={z2} m', c='k', alpha=0.1)
ax[1,0].legend()
ax[1,0].set_ylabel("[%]")

ax[1,1].set_title("Total specific humidity")
ax[1,1].plot(time, qt_atm0*1000, label=f"qt at z={z0} m", c='k', alpha=0.5)
ax[1,1].plot(time, qt_atm1*1000, label=f"qt at z={z1} m", c='k', alpha=0.3)
ax[1,1].plot(time, qt_atm2*1000, label=f"qt at z={z2} m", c='k', alpha=0.1)
ax[1,1].plot(time, qt_tech*1000, label="qt_tech", c='k')
ax[1,1].legend(ncol=2)
ax[1,1].set_ylabel(["g kg-1"])

# ds.rh.plot()
for axs in ax.flatten():
    axs.set_xlabel('time [s]')
    
fig.tight_layout()
fig.savefig("solar_evap_test_plot.png")
plt.show()