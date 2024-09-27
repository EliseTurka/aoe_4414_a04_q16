# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  Converts ECEF vector components to SEZ
#  See "Fundamentals of Astrodynamics and Applications, Fourth Edition" by
#  David A. Vallado, pages 172-173
# Parameters:
#  o_x_km: ECEF origin of SEZ frame
#  o_y_km: ECEF origin of SEZ frame
#  o_z_km: ECEF origin of SEZ frame
#  x_km: ECEF x-position
#  y_km: ECEF y-position
#  x_km: ECEF x-position
#
# Output:
#  Prints the s, e, and z (SEZ) components of the ground station
#
# Written by Elise Turka
# Other contributors: None
#
# This work is licensed under CC BY-SA 4.0

# import Python modules
import math # math module
import sys  # argv

# "constants"
R_E_KM = 6378.137
E_E    = 0.081819221456

# helper functions


# initialize script arguments
o_x_km = float('nan') # ECEF x in SEZ frame
o_y_km = float('nan') # ECEF y in SEZ frame
o_z_km = float('nan') # ECEF z in SEZ frame
x_km = float('nan') # ECEF x-position
y_km = float('nan') # ECEF y-position
z_km = float('nan') # ECEF z-position

# parse script arguments
if len(sys.argv)==7:
  o_x_km = float(sys.argv[1])
  o_y_km = float(sys.argv[2])
  o_z_km = float(sys.argv[3])
  x_km = float(sys.argv[4])
  y_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
   'Usage: '\
   'python3 ecef_to_llh.py o_x_km o_y_km o_z_km x_km y_km z_km'\
  )
  exit()

# setting up vectors

# calculating ECEF FROM STATION TO OBJECT
ecef_x = o_x_km - x_km
ecef_y = o_y_km - y_km
ecef_z = o_z_km - z_km

# calculating LAT, LON, HAE OF GROUND STATION
lon_rad = math.atan2(y_km,x_km)
lon_deg = lon_rad*180.0/math.pi

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(z_km/math.sqrt(x_km**2+y_km**2+z_km**2))
r_lon_km = math.sqrt(x_km**2+y_km**2)
prev_lat_rad = float('nan')

def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0-math.pow(ecc,2.0) * math.pow(math.sin(lat_rad),2.0))

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1
  
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

# calculating INVERSE ROTATIONS
s_km = ecef_x*math.sin(lat_rad)*math.cos(lon_rad)+ecef_y*math.sin(lat_rad)*math.sin(lon_rad)-ecef_z*math.cos(lat_rad)
e_km = ecef_y*math.cos(lon_rad)-ecef_x*math.sin(lon_rad)
z_km = math.cos(lon_rad)*math.cos(lat_rad)*ecef_x+math.cos(lat_rad)*math.sin(lon_rad)*ecef_y+math.sin(lat_rad)*ecef_z

# print statements
print(s_km)
print(e_km)
print(z_km)