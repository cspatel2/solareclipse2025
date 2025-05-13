#%%
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj, Transformer

# Constants
EARTH_RADIUS_KM = 6378.137  # Earth's equatorial radius in kilometers

# North Hero, VT coordinates
lat_center = 44.8062  # degrees
lon_center = -73.2732  # degrees

# Besselian elements at a specific time (example values)
# These should be obtained from NASA's data for the exact time of interest
x = 0.1234  # in Earth radii
y = -0.5678  # in Earth radii
l2 = 0.0056  # umbral radius in Earth radii

# Convert x, y, l2 to kilometers
x_km = x * EARTH_RADIUS_KM
y_km = y * EARTH_RADIUS_KM
l2_km = l2 * EARTH_RADIUS_KM

# Define a local projection centered on North Hero, VT
proj = Proj(proj='aeqd', lat_0=lat_center, lon_0=lon_center, datum='WGS84')
transformer = Transformer.from_proj(proj, proj, always_xy=True)

# Compute the position of the umbra center relative to North Hero
umbra_x, umbra_y = transformer.transform(lon_center + x_km / (EARTH_RADIUS_KM * np.pi / 180),
                                         lat_center + y_km / (EARTH_RADIUS_KM * np.pi / 180))

# Plotting
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_title('Umbra Projection on April 8, 2024')
ax.set_xlabel('East-West Distance (km)')
ax.set_ylabel('North-South Distance (km)')

# Plot the umbra
umbra_circle = plt.Circle((umbra_x, umbra_y), l2_km, color='black', alpha=0.5)
ax.add_artist(umbra_circle)

# Plot North Hero, VT
ax.plot(0, 0, 'ro', label='North Hero, VT')

# Set limits
ax.set_xlim(-200, 200)
ax.set_ylim(-200, 200)
ax.legend()
ax.grid(True)
plt.axis('equal')
plt.show()

# %%
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_body
from astropy.coordinates import solar_system_ephemeris
import astropy.units as u

# Step 1: Define the observer location (North Hero, VT)
north_hero = EarthLocation(lat=44.8062*u.deg, lon=-73.2732*u.deg, height=30*u.m)

# Step 2: Set the time of maximum eclipse in UTC
eclipse_time = Time("2024-04-08 18:17:00")  # Approximate time of max eclipse at North Hero

# Step 3: Compute Sun and Moon positions
with solar_system_ephemeris.set('de432s'):
    sun_coord = get_sun(eclipse_time).transform_to(AltAz(obstime=eclipse_time, location=north_hero))
    moon_coord = get_body('moon',eclipse_time).transform_to(AltAz(obstime=eclipse_time, location=north_hero))

# Step 4: Compute the shadow vector in local frame (Alt-Az)
shadow_alt = moon_coord.alt - sun_coord.alt
shadow_az = moon_coord.az - sun_coord.az

# Step 5: Estimate umbral cone parameters
# Sun and Moon angular sizes at eclipse time
sun_radius_deg = (0.53 / 2) * u.deg  # approximate
moon_radius_deg = (0.52 / 2) * u.deg  # approximate

# Shadow cone angle (approximate small angle assumption)
cone_angle = sun_radius_deg - moon_radius_deg  # negative means umbral shadow converges

# Shadow length L from Moon to tip of umbra
moon_distance_km = 384400  # average Earth-Moon distance
sun_distance_km = 149600000  # average Earth-Sun distance

# Shadow radius at Earth's surface (simplified)
umbra_radius_km = moon_distance_km * (sun_radius_deg - moon_radius_deg).to(u.rad).value

# For plotting: X-axis spans ~100 km around North Hero, Y is altitude up to 100 km
x = np.linspace(-50, 50, 500)
ground = np.zeros_like(x)

# Shadow cone projected: we'll assume it's a triangle for the umbra
umbra_top = 100  # km altitude where cone starts
umbra_bottom = 0  # intersects the Earth
umbra_width_km = 2 * abs(umbra_radius_km)  # full width at base

# Simple triangle cone with base centered on 0
umbra_left = -umbra_width_km / 2
umbra_right = umbra_width_km / 2

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot topography (flat for now)
ax.plot(x, ground, color='green', label='Ellipsoid')

# Draw the umbra as a shaded triangle
umbra_x = [umbra_left, 0, umbra_right]
umbra_y = [0, umbra_top, 0]
ax.fill(umbra_x, umbra_y, color='gray', alpha=0.5, label='Umbra')

# Annotate the point of interest
ax.plot(0, 0, 'go', label='North Hero, VT (Ellipsoid)')
ax.plot(0, 0.01, 'ro', label='North Hero, VT (Topography)')  # Slight offset for clarity

ax.set_xlabel('Position relative to North Hero, VT (km)')
ax.set_ylabel('Altitude (km)')
ax.set_title('Eclipse Umbra Projection at North Hero, VT\n2024-04-08 18:17 UTC')
ax.set_ylim(0, 110)
ax.set_xlim(-50, 50)
ax.legend()
ax.grid(True)

plt.show()

# %%
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

# Observer location and eclipse time
north_hero = EarthLocation(lat=44.8062*u.deg, lon=-73.2732*u.deg, height=30*u.m)
eclipse_time = Time("2024-04-08 18:17:00")

# Compute the Sun's altitude at the given time and location
altaz = AltAz(obstime=eclipse_time, location=north_hero)
sun_altitude = get_sun(eclipse_time).transform_to(altaz).alt
angle_of_incidence = 90 * u.deg - sun_altitude  # angle from horizontal

# Define plot extent (in km)
x_range = 100  # km to each side of center
ground = np.array([-x_range, x_range])

# Length of umbra shadow (approximate based on average umbra length)
umbra_length = 100 * u.km  # you can refine this with Besselian elements later

# Compute umbra tip height using the incidence angle
height_km = np.tan(angle_of_incidence.to(u.rad)) * x_range

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

# Ground line
ax.plot(ground, [0, 0], 'k', lw=2)

# Umbra cone from left and right edges to tip
umbra_left = [-x_range, 0]
umbra_right = [x_range, 0]
umbra_tip = [0, height_km.value]

ax.plot([umbra_left[0], umbra_tip[0]], [umbra_left[1], umbra_tip[1]], 'gray', lw=2)
ax.plot([umbra_right[0], umbra_tip[0]], [umbra_right[1], umbra_tip[1]], 'gray', lw=2)

# Fill umbra region
ax.fill([umbra_left[0], umbra_tip[0], umbra_right[0]], 
        [umbra_left[1], umbra_tip[1], umbra_right[1]], 
        'gray', alpha=0.3, label="Umbra")

# Labels and annotations
ax.set_xlabel("Horizontal Distance from North Hero, VT (km)")
ax.set_ylabel("Altitude (km)")
ax.set_title("2024 Eclipse Umbra Projection at North Hero, VT")
ax.grid(True)
ax.axvline(0, color='blue', linestyle='--', label='North Hero Center')
ax.legend()

plt.tight_layout()
plt.show()

# %%
from astropy.coordinates import EarthLocation
import astropy.units as u
import numpy as np

def earth_radius_astropy(lat_deg, lon_deg):
    # Create EarthLocation at zero height (sea level)
    loc = EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=0*u.m)

    # Convert to Cartesian (ITRS = Earth-centered frame)
    x = loc.x.to(u.km).value
    y = loc.y.to(u.km).value
    z = loc.z.to(u.km).value

    # Radius from Earth's center
    r = np.sqrt(x**2 + y**2 + z**2)
    return r

# Example: North Hero, VT (~44.81° N, -73.27° E)
radius_km = earth_radius_astropy(44.81, -73.27)
print(f"Earth radius at North Hero, VT is approximately {radius_km:.3f} km")
# %%
from astropy.coordinates import get_sun, get_body, ITRS, CartesianRepresentation
from astropy.time import Time
from astropy import units as u
import numpy as np

def umbra_center_position(time_utc):
    time = Time(time_utc)

    # Get positions in ITRS (Earth-centered frame)
    moon = get_body('moon',time).transform_to(ITRS(obstime=time))
    sun = get_sun(time).transform_to(ITRS(obstime=time))

    # Extract Cartesian vectors in km
    r_moon = moon.cartesian.xyz.to(u.km).value
    r_sun = sun.cartesian.xyz.to(u.km).value

    # Shadow direction (from Moon toward Sun)
    d = r_sun - r_moon
    d_unit = d / np.linalg.norm(d)

    # Radius of Earth (mean)
    R_earth = 6371.0  # km

    # Ray from Moon to Earth: r = r_moon + t*(-d_unit)
    # Solve for intersection with sphere: ||r|| = R_earth
    # => ||r_moon - t*d_unit||^2 = R_earth^2

    a = 1
    b = -2 * np.dot(r_moon, d_unit)
    c = np.dot(r_moon, r_moon) - R_earth**2

    # Solve quadratic for t
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return None  # no intersection

    t = (-b - np.sqrt(discriminant)) / (2*a)

    # Umbra center position
    r_umbra = r_moon - t * d_unit

    # Convert to geodetic coordinates (lat, lon, height)
    loc = ITRS(CartesianRepresentation(r_umbra * u.km), obstime=time)
    geodetic = loc.earth_location.to_geodetic()

    return {
        "latitude": geodetic.lat.deg,
        "longitude": geodetic.lon.deg,
        "height_km": geodetic.height.to(u.km).value,
        "x_km": r_umbra[0],
        "y_km": r_umbra[1],
        "z_km": r_umbra[2],
    }

# Example usage
result = umbra_center_position("2024-04-08 18:17:00")
print(result)

# %%
from datetime import datetime
from pytz import UTC
from astropy.coordinates import get_body, get_sun, ITRS
from astropy.time import Time

dt = datetime(2024,4,8,19,25, tzinfo = UTC) #input object, datetime

time = Time(dt,scale = 'utc') #astropy time object
moon = get_body('moon',time).transform_to(ITRS(obstime=time))
sun = get_sun(time).transform_to(ITRS(obstime=time))
# %%
