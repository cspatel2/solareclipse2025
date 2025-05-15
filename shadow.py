# %%
from typing import Tuple
from datetime import datetime,timedelta
from pytz import timezone, UTC
from astropy.coordinates import get_body, get_sun, ITRS, EarthLocation, CartesianRepresentation, AltAz, earth,  GCRS
from scipy.optimize import root_scalar
import astropy.units as u
from astropy.time import Time
import numpy as np
import sys
from geopy.distance import geodesic
from pyproj import Geod
import matplotlib.pyplot as plt
import matplotlib as mpl
from eclipse_geometry import solar_eclipse
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
# %%


def get_umbra(dt: datetime, location: EarthLocation=None, printbool:bool=False):# -> Tuple[earth.GeodeticLocation, float, float, float]:
    
    """ Calulates the location, size, and inclination of the umbra during an eclipse given an observer location and datetime.

    Args:
        dt (datetime): datetime object with tz = pytz.UTC.
        location (EarthLocation): observer location in geodetic coordinates (longitude, latitude, height above a reference ellipsoid). When using the geodetic forms, Longitudes are measured increasing to the east, so west longitudes are negative. Defaults to None.
        printbool (bool, optional): if True, all return values will print. Defaults to False.

    Returns:
        umbra_loc (earth.GeodeticLocation): lat, lon of the umbra center in deg.\n
        Ru (float): Radius of the umbra on the ground in km.\n
        phi (float): Angle the umbra makes with the zenith (of observer location) in deg. Nan if observer location is None.\n\n

        None if there is no Eclipse.
    """
    Re = 6367.607 * u.km  # average radius of earth (this is from a latitude of the north hero,VT)
    Rm = 1737.4 * u.km  # average radius of moon
    Rs = 696340 * u.km  # average radius of sun

    time = Time(dt, scale='utc')  # astropy Time object
    moon = get_body('moon', time).transform_to(ITRS(obstime=time))
    sol = get_sun(time).transform_to(ITRS(obstime=time))

    # moon vector
    mvec = moon.cartesian.xyz.to(u.km).value  # geocentric moon vector
    M = np.sqrt(np.dot(mvec, mvec))  # magnitude of geocentric moon vector
    mhat = mvec/M  # unit vector in the firextion of geocentric moon vector

    # sun nector
    svec = sol.cartesian.xyz.to(u.km).value  # geocentric solar vector

    # shadow (umbra) vector
    gvec = svec-mvec  # shadow vector from moon -> sun direct
    # magnitude of shadow vector from moon -> sun direction (distance between sun and moon)
    G = np.sqrt(np.dot(gvec, gvec))
    ghat = gvec / G  # unit vector of shadow from moon -> sun direction

    # the center of the umbra (shadow) vector is ~u = U(- g^)
    # calc the magnitude (U) of umbra vector (~u)
    a = 1
    b = -2*M * np.dot(mhat, ghat)
    c = np.dot(mvec, mvec) - Re.value**2
    sub = b**2 - 4*a*c
    if sub < 0:
        print('No eclipse at this time.')
        return None
    U = (-b - np.sqrt(sub)) / (2*a)
    uvec = mvec - U*ghat  # umbra vector from moon -> earth direction where it intersects earth
   

    # umbra location in lat lon
    # umbra center in geocentric coordinate system
    umbra = ITRS(CartesianRepresentation(*(uvec), unit=u.km), obstime=time)
    # umbra in geodetic coordinate system (lat,lon) in degrees.
    umbra_loc = umbra.earth_location.to_geodetic()
    ALTAZ = AltAz(location=location, obstime=time)
    umbra_azel = umbra.transform_to(ALTAZ)
    
    # for the radius of umbra (Theta)
    # distance between umbra tip to moon
    L = Rm.value * G / (Rs.value - Rm.value)
    lvec = L*u.km*ghat  # vector from umbra tip -> moon center
    # theta = np.arcsin(Rm.value/L)  # radians
    theta = np.arcsin((Rm / L).value)
    Ru = (L-(M-Re.value))*np.tan(theta)

    if location == None:
        phi = np.nan
    else:
        # zenith vector at observer location
        obs = location.get_itrs(obstime=time)
        # geocentric observer zenith vector
        ovec = obs.cartesian.xyz.to(u.km).value
        # magnitude of geocentric oberver zenith vector
        O = np.sqrt(np.dot(ovec, ovec))
        ohat = ovec/O  # unit vector of geocentric observer zenith vector

        # angle between umbra vector and zenith vector (phi)
        phi = np.degrees(np.arccos(np.dot(ghat, ohat)))
        # check step: phi should be equal to the solar zenith angle at location
        ALTAZ = AltAz(obstime=time, location=location)
        sol = get_sun(time).transform_to(ALTAZ)
        za = 90 - sol.alt.degree
        # print(f'phi = {phi}, za = {za}, diff = {np.abs(phi-za)}')

    if printbool:
        print(
            f'Umbra location: LAT {umbra_loc.lat.value:0.2f}, LON {umbra_loc.lon.value:0.2f}')
        print(f'radius of umbra: {Ru:0.2f}')
        print(f'width of path of totality: {2*Ru:0.2f}')

        print(
            f'Angle of umbra from zenith(observer) = {phi:0.2f} ± {np.abs(phi-za):0.2f} deg')
        print(f'Solar Location: ALT {sol.alt.value:0.1f}, AZ {sol.az.value:0.1f}')

    return umbra_loc, Ru, phi


def get_azimuth_totality_path(start_dt:datetime,end_dt:datetime,obs_loc:EarthLocation)->float:
    """calculates the direction of path of totality. 0 deg -> North.

    Args:
        start_dt (datetime):  start time. start_time < end_time
        end_dt (datetime): end time.
        obs_loc (EarthLocation): observer location in geodetic coordinates (longitude, latitude, height above a reference ellipsoid). When using the geodetic forms, Longitudes are measured increasing to the east, so west longitudes are negative.

    Returns:
        float: average direction for path of totality. angle in  degrees.
    """     
    loc1,_,_ = get_umbra(start_dt, obs_loc)
    loc2,_,_ = get_umbra(end_dt, obs_loc)
    geod = Geod(ellps="WGS84")
    az,_,_ = geod.inv(loc1.lon.value,loc1.lat.value,loc2.lon.value,loc2.lat.value) 
    return az

def get_signed_dist_to_umbra(dt:datetime, az_totality_path:float, obs_loc:EarthLocation)->float:
    """returns signed distance of umbra center  relative to a observer location on earth.

    Args:
        dt (datetime): current time.
        az_totality_path (float): average direction of path of totality. angle in degrees with 0 deg -> North.
        obs_loc (EarthLocation): observer location.

    Returns:
        float: distance in km. 
    """    
    umbra_loc,_,_ = get_umbra(dt,obs_loc)

    geod = Geod(ellps="WGS84")
    #distance from observer location to umbra center
    az, _,dist = geod.inv(obs_loc.lon.value,obs_loc.lat.value, umbra_loc.lon.value, umbra_loc.lat.value) 

    #angle diff between current az and path of totality az
    azdif = (az-az_totality_path+180)%360 - 180 #between [-180, 180]
    
    sign = np.sign(np.cos(np.radians(azdif)))

    return sign*dist/1000 #m -> km

def horizontal(altitude: float, za: float, hoffset: float = 0) -> float:
    """Calculates the horizontal distance in km for a given zenith angle and altitude.

    Args:
        altitude (float): altitude in km.
        za (float):  zenith angle in degrees.
        hoffset (float, optional): horizontal offset in km.

    Returns:
        float: horizontal distance in km. 
    """
    return altitude/np.cos(np.deg2rad(za)) + hoffset

def get_umbra_speed(location:EarthLocation, start_dt:datetime, end_dt:datetime, t_interval_s:int):
    """Prints Speed of the umbra at a given time.

    Args:
        location (EarthLocation): observer location in geodetic coordinates (longitude, latitude, height above a reference ellipsoid). When using the geodetic forms, Longitudes are measured increasing to the east, so west longitudes are negative.
        start_dt (datetime): start time in UTC.
        end_dt (datetime): end time in UTC.
        t_interval_s (int): time interval in seconds.
    """    
    interval = timedelta(seconds=t_interval_s)
    geod = Geod(ellps="WGS84")
    
    umbra_positions = []
    timestamps = []

    dt = start_dt
    while dt <= end_dt:
        result = get_umbra(dt, location)
        if result:
            loc, _, _ = result
            umbra_positions.append((loc.lat.value, loc.lon.value))
        else:
            umbra_positions.append(None)
        timestamps.append(dt)
        dt += interval

    # Calculate speeds between each time step
    print("\nTime      | Speed (km/min) | Distance (km)")
    for i in range(1, len(umbra_positions)):
        prev = umbra_positions[i - 1]
        curr = umbra_positions[i]
        if prev and curr:
            az, _, dist = geod.inv(prev[1], prev[0], curr[1], curr[0])  # lon, lat order
            speed = dist/interval.seconds *60 / 1000  # in km/min
            print(f"{timestamps[i].strftime('%H:%M:%S')} | {speed:12.2f} | {dist/1000:10.2f}")
        else:
            print(f"{timestamps[i].strftime('%H:%M:%S')} |       --       |     --")
#%%

#%%
if __name__ == '__main__':
    location = EarthLocation(lat=44.77*u.deg, lon=-73.29*u.deg) #North Hero, VT
    peak_dt = datetime(2024,4,8,19,27,46, tzinfo = UTC)
    tdel = 10
    start_dt =  peak_dt - timedelta(minutes = tdel)
    end_dt = peak_dt + timedelta(minutes = tdel)

    #umbra location, radius, and zenith angle of umbra
    uloc,ru,phi = get_umbra(peak_dt,location,True)
    # direction of eclipse path described by azimuth
    eclipse_path = get_azimuth_totality_path(start_dt,end_dt,location)
    # distace along the path of totality wrt the observer location
    dist = get_signed_dist_to_umbra(peak_dt,eclipse_path,location)
    print(f'distance of umbra from location: {dist:0.2f} km')
    #umbra speed over the time interval at dt = 120s
    get_umbra_speed(location,start_dt,end_dt,120)

    plot = True
    if plot:
        #time stamps and corresponding transparency for plotting
        offset_minutes = [5,10,20]
        max_alpha = 1
        tau = 5  # control how fast it fades
        alphas = [max_alpha] # transparency
        obstime_utc = [peak_dt] 

        for minutes in offset_minutes: #creates tstamps at each ofset_min on eitherside of the center dt.
            alpha = max_alpha * np.exp(-minutes / tau)
            alphas.append(alpha)
            alphas.append(alpha)
            delta = timedelta(minutes=minutes)
            obstime_utc.append(peak_dt - delta)  # time before
            obstime_utc.append(peak_dt + delta)  # time after
        sorted_pairs = sorted(zip(obstime_utc, alphas), key=lambda x: x[0])
        obstime_utc, alphas = zip(*sorted_pairs)

        #direction of path of totality
        eclipse_path = get_azimuth_totality_path(start_dt,end_dt,location)

        #PLOTTING
        fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8), dpi=300,tight_layout=True, squeeze=True)
        altitude = np.linspace(0, 500, 10)  # km

        #plot slit
        slit_el, slit_az = 60, 230  # deg
        fov = [horizontal(a,90-slit_el) for a in -altitude]
        ax.plot(fov, altitude, color='k', lw=1)

        #plot umbra at time stamps chosen above
        for i, t in enumerate(obstime_utc):
            uloc,ru,phi = get_umbra(t, location,False)
            d = get_signed_dist_to_umbra(t,eclipse_path,location)
            line1 = [horizontal(a,phi,d-ru) for a in altitude]
            line2 = [horizontal(a,phi,d+ru) for a in altitude]
            ax.plot(line1, altitude, color='orange', lw=0.5)
            ax.plot(line2, altitude, color='orange', lw=0.5)
            _,_,ecp = solar_eclipse(t,location)
            ax.fill_betweenx(altitude,line1,line2,alpha=.5, color='orange')

            
            hhmm = t.strftime(f"%H:%M {UTC}")
            x = (line2[-1] + line2[0])/2
            y = 300
            ax.text(x+40, y, f'{hhmm} - {ecp*100:.0f} %',rotation=90-phi-10,ha='center', va='bottom',rotation_mode='default', transform_rotates_text=True, fontsize = 8)

        ax.grid(color = 'grey', ls = '--', lw = 0.5, alpha = 0.3)
        # ax.set_xlim(-200,200)
        ax.set_ylim(0,400)
        ax.set_xlabel('Horizontal Distance along Path of Totality (km)')
        ax.set_ylabel('Altitude (km)')






