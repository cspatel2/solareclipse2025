# %%
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_body,SkyCoord, BaseGeodeticRepresentation
from astropy.time import Time
import astropy.units as u
from matplotlib.patches import Rectangle, Circle
from datetime import datetime, timedelta
from pytz import timezone, UTC
import sunpy
import sunpy.coordinates
import matplotlib as mpl
usetex = False

if not usetex:
    # computer modern math text
    mpl.rcParams.update({'mathtext.fontset': 'cm'})
mpl.rc('font', **{'family': 'serif',
       'serif': ['Times' if usetex else 'Times New Roman']})
# for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
mpl.rc('text', usetex=usetex)
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=10) 

#%%
def solar_eclipse(dt:datetime, location:EarthLocation):
    """calculates solar eclipse at the observer location.

    Args:
        dt (datetime): time in UTC. 
        location (EarthLocation): Observer location in geodetic coordinates (longitude, latitude, height above a reference ellipsoid).Longitudes are measured increasing to the east, so west longitudes are negative.

    Returns:
        sol (astropy.coordinates.sky_coordinate.SkyCoord): Postion of sun in Altitude-Azimuth coordinate system at observer location.\n
        moon (astropy.coordinates.sky_coordinate.SkyCoord): Postion of moon in Altitude-Azimuth coordinate system at observer location.\n
        ecp (float): fraction of eclipse with values [0(No eclipse), 1(Totality)]. 

    """    
    time = Time(dt,scale = 'utc')
    ALTAZ = AltAz(obstime=time, location=location)
    sol = get_sun(time).transform_to(ALTAZ)
    moon = get_body('moon', time).transform_to(ALTAZ)
    ecp = sunpy.coordinates.sun.eclipse_amount(location.get_itrs(time),moon_radius='minimum') #faction of eclipse
    return sol, moon, float(ecp)

#%%
# Location: North Hero, Vermont
location = EarthLocation(lat=44.83*u.deg, lon=-73.27*u.deg, height=43*u.m)

#date/time: both edt and utc centered around eclipse peak at location.
eastern = timezone('US/Eastern')
peak_dt = datetime(2024, 4, 8, 15,27,50)
offset_minutes = [21, 48, 82]
peak_dt = eastern.localize(peak_dt, is_dst=True) #day time savings
res = [peak_dt]
for m in offset_minutes:
    delta = timedelta(minutes = m)
    res.append(peak_dt + delta)
    res.append(peak_dt - delta)
res.sort()

times = dict(edt = res, utc = [t.astimezone(UTC) for t in res])


#%%
#Set up plot
TIME = 'utc'

fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=300, tight_layout=True, squeeze=True)
ax.set_xlabel("Azimuth (°)")
ax.set_ylabel("Elevation (°)")
ax.grid(True, linestyle='--', alpha=0.6, lw = .1)
ax.set_aspect('equal')
# ax.scatter(solaz,solalt)
fontsize = 8
ax.set_title(f'Eclipse Geometry \n North Hero, Vermont ({location.lon.degree:0.1f} W, {location.lat.degree:0.1f} N)')

#plot slit
#center of slit (x,y) = slit_az, slit_alt
slit_alt = 60 #alt, deg
slit_az = 230 #az, deg
slit_width, slit_height = 0.5, 16
bottom = slit_alt - slit_height/2
left = slit_az - slit_width/2
slit = Rectangle((left,bottom), slit_width, slit_height, facecolor = 'grey', edgecolor = 'black', lw = .3)
ax.add_patch(slit)
ax.text(slit_az + 1.5, slit_alt, 'HiT&MIS FOV', rotation = 270, 
                         horizontalalignment='center', verticalalignment='center', fontsize=fontsize)


#initalize 
solalt,solaz = [], []
moonalt,moonaz = [], []
eclipse = []
solr, moonr = .5/2, .5/2 #angluar radius

#find solar and lunar postions and eclipse %
for i,tutc in enumerate(times[TIME.lower()]):
    sol,moon,ecp = solar_eclipse(tutc,location)
    solalt.append(sol.alt.degree)
    solaz.append(sol.az.degree)
    moonalt.append(moon.alt.degree)
    moonaz.append(moon.az.degree)
    eclipse.append(ecp*100)


    #plotsun and moon
    s = Circle((solaz[i],solalt[i]), solr, color = 'red')
    m = Circle((moonaz[i],moonalt[i]), moonr, color = 'black')
    ax.add_patch(s)
    ax.add_patch(m)
    hhmm = tutc.strftime(f"%H:%M {TIME.upper()}")
    print(hhmm)
    ax.text(solaz[i]-.25, solalt[i] - 3.75, f'{hhmm} \n{eclipse[i]:.0f} %', ha='center', fontsize=fontsize)



ax.set_xlim(np.min(solaz+moonaz)-5,np.max(solaz+moonaz)+5 )
ax.set_ylim(np.min(solalt+moonalt)-5,np.max(solalt+moonalt+[slit_alt+slit_height/2])+2 )

plt.savefig(f'Eclipse_geometry_{TIME.upper()}.png') 
# %%
