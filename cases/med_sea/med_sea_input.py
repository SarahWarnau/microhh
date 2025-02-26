import numpy as np
import netCDF4 as nc
import pandas as pd
import os
    
float_type = 'f8'

case = 'med_sea'
profiles = 'med_sea_JJA'

# Get number of vertical levels and size from .ini file
with open(f'{case}.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)

df_ini = pd.read_csv(f"{profiles}_ini.csv").iloc[1:]

thl   = np.interp(z, df_ini.z, df_ini.theta)
qt    = np.interp(z, df_ini.z, df_ini.q)
u     = np.interp(z, df_ini.z, df_ini.ws)
ugeo  = u
v     = np.zeros(z.size)
vgeo  = np.zeros(z.size)


# write the data to a file
nc_file = nc.Dataset(f"{case}_input.nc", mode="w", datamodel="NETCDF4", clobber=True)
nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
nc_ugeo  = nc_group_init.createVariable("u_geo" , float_type, ("z"))
nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
nc_vgeo  = nc_group_init.createVariable("v_geo" , float_type, ("z"))


nc_z    [:] = z    [:]
nc_thl  [:] = thl  [:]
nc_qt   [:] = qt   [:]
nc_u    [:] = u    [:]
nc_ugeo [:] = ugeo [:]
nc_v    [:] = v    [:]
nc_vgeo [:] = vgeo [:]

nc_file.close()

# Surface settings
df_surf = pd.read_csv(f"{profiles}_ini.csv").iloc[0]
print("pbot =", df_surf['p'].item())
print("sbot[thl] =", df_surf['theta'].item())
print("sbot[qt]", df_surf['q'].item())