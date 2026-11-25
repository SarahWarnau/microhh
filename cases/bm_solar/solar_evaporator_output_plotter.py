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

datasets = [open_default("bm_solar", time="ref"), open_default("bm_solar")]
linestyles = [":", "-"]
names = ["ref", "solar"]

fig, axs = plt.subplots(3,3, figsize=(16,12))
ax = axs.flatten()
for i, ds in enumerate(datasets):
    name = names[i]
    ls = linestyles[i]
    
    # print(ds.variables)
    dz = ds.isel(z=0).z.values
    # dz = ds.isel(zh=1).zh.values
    LE = ds.isel(zh=0).qt_flux*ds.isel(zh=0).rhoh*2.5e6/dz
    H = ds.isel(zh=0).thl_flux*ds.isel(zh=0).rhoh*1005/dz
    
    EF = LE/(LE+H)
    
    thl_atm0 = ds.sel(z=0, method="nearest")["T"]
    thl_atm1 = ds.sel(z=10, method="nearest")["T"]
    thl_atm2 = ds.sel(z=50, method="nearest")["T"]
    
    thl_tech = ds["thl_bot"]
    
    qt_atm0 = ds.sel(z=0, method="nearest")["qt"]
    qt_atm1 = ds.sel(z=10, method="nearest")["qt"]
    qt_atm2 = ds.sel(z=50, method="nearest")["qt"]
    
    qt_tech = ds["qt_bot"]
    
    RH0 = ds.sel(z=0, method="nearest")["rh"]
    RH1 = ds.sel(z=10, method="nearest")["rh"]
    RH2 = ds.sel(z=50, method="nearest")["rh"]
    
    thl_prof = ds.isel(time=-1)["thl"]
    qt_prof = ds.isel(time=-1)["qt"]*1000
    RH_prof = ds.isel(time=-1)["rh"]*100
    blh = ds.isel(time=-1)["zi"]
    z = ds.z.values
    
    time = ds.time
    
    z0 = np.round(ds.sel(z=0, method="nearest").z.values)
    z1 = np.round(ds.sel(z=10, method="nearest").z.values)
    z2 = np.round(ds.sel(z=50, method="nearest").z.values)
    
    j=0
    ax[j].set_title('Total water')
    ax[j].plot(time, ds['qt_path'], ls=ls, c='k', label=name)
    ax[j].legend()
    ax[j].set_ylabel("[kg m-2]")
    
    j+=1
    ax[j].set_title("Surface fluxes")
    ax[j].plot(time, LE, label=f"LE {name}", c='tab:blue', ls=ls)
    ax[j].plot(time, H, label=f"H {name}", c='tab:red', ls=ls)
    ax[j].plot(time, LE+H, label=f"LE+H {name}", c='tab:green', ls=ls)
    if i == 1:
        ax[j].axhline(1000, label=f"Rnet {name}", c='k', ls='-')
    ax[j].legend(ncol=4)
    ax[j].set_ylabel("[W m-2]")
    
    j+=1
    ax[j].set_title("Temperature")
    ax[j].plot(time, thl_atm0, label=f"={z0} m {name}", c='k', alpha=0.5, ls=ls)
    # ax[0,1].plot(time, thl_atm1, label=f"thl at z={z1} m", c='k', alpha=0.3)
    ax[j].plot(time, thl_atm2, label=f"z={z2} m {name}", c='k', alpha=0.1, ls=ls)
    ax[j].plot(time, thl_tech, label=f"surf. {name}", c='k', ls=ls)
    ax[j].legend(ncol=2)
    ax[j].set_ylabel("[K]")
    
    j+=1
    ax[j].set_title("Relative humidity")
    ax[j].plot(time, RH0*100, label=f'z={z0} m {name}', c='k', alpha=0.5, ls=ls)
    # ax[1,0].plot(time, RH1*100, label=f'RH at z={z1} m', c='k', alpha=0.3)
    ax[j].plot(time, RH2*100, label=f'z={z2} m {name}', c='k', alpha=0.1, ls=ls)
    ax[j].legend()
    ax[j].set_ylabel("[%]")
    
    j+=1
    ax[j].set_title("Total specific humidity")
    ax[j].plot(time, qt_atm0*1000, label=f"z={z0} m {name}", c='k', alpha=0.5, ls=ls)
    # ax[1,1].plot(time, qt_atm1*1000, label=f"qt at z={z1} m", c='k', alpha=0.3)
    ax[j].plot(time, qt_atm2*1000, label=f"z={z2} m {name}", c='k', alpha=0.1, ls=ls)
    ax[j].plot(time, qt_tech*1000, label=f"surf. {name}", c='k', ls=ls)
    ax[j].legend(ncol=2)
    ax[j].set_ylabel("[g kg-1]")
    
    j+=1
    ax[j].set_title("Evaporative fraction")
    ax[j].plot(time, EF, label=name, c='k', alpha=0.5, ls=ls)
    if i == 1:
        ax[j].plot(time, LE/1000, c='r', ls=ls, label="Tech. efficiency = LE/Rnet")
    ax[j].legend()
    ax[j].set_ylabel("[-]")
    ax[j].set_ylim(0,1)
    
    
    for k, prof in enumerate([thl_prof, qt_prof, RH_prof]):
        profname = ["thl", "qt", "RH"][k]
        xlabel = ["[K]", "[g/kg]", "[%]"][k]
        j+=1
        ax[j].set_title(f"{profname} profile")
        ax[j].axhline(blh, ls=ls, label=f"ABL height {name}")
        ax[j].plot(prof, z, label=f"{profname} {name}", c='k', ls=ls)
        ax[j].legend()
        ax[j].set_xlabel(xlabel)
        ax[j].set_ylabel("[m]")
    

    
# ds.rh.plot()
for axss in ax[:-3]:
    axss.set_xlabel('time [s]')
    
fig.tight_layout()
fig.savefig("solar_evap_test_plot.png")
plt.show()