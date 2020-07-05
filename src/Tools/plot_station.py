import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from down_data import sta_list

mpl.rcParams['font.sans-serif'] = ['Helvetical']
mpl.rcParams['axes.unicode_minus'] = False
mpl.rc('xtick', labelsize=9)
mpl.rc('ytick', labelsize=9)
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
plt.rcParams['savefig.dpi'] = 300

# Load the coordinate of IGS Core & MGEX sites, The CSV files are
# exported from: http://www.igs.org/network
igs = np.recfromcsv('igs.csv', names=True, encoding='utf-8')

fig = plt.figure(figsize=(4.8, 2.97))
# Set projection
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
# Add ocean and land
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE, linewidth=0.1)

# sites = ['ABMF', 'DJIG', 'FAA1', 'HARB',
#          'JFNG', 'KIRU', 'MAYG', 'MGUE',
#          'NNOR', 'NYA2', 'PADO', 'SGOC', 'SGPO', 'TOW2', 'URUM']
sites = sta_list

lists = []
for i in igs:
    for j in sites:
        if i['site'] == j:
            lists.append(i)

plot_sites = pd.array(lists)

# Add MGEX & IGS core sites
ax.plot(plot_sites['long'], plot_sites['lat'], 'o', color='tomato', transform=ccrs.Geodetic(), ms=1.0)

for j in plot_sites:
    plt.text(j['long'], j['lat'], j['site'], transform=ccrs.Geodetic(), FontSize=10)

# Plot gridlines
# ax.gridlines(linestyle='--', LineWidth=0.05)
ax.set_global()
ax.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

fig.savefig('sites.tiff', bbox_inches='tight', dpi=300)
plt.show()
