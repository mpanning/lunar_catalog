"""
Calculate geodetic distance between moonquake epicenters and Apollo landing
elements
"""
from geopy import distance, units
import pandas as pd

DMQlocsfile = 'Nakamura2005Locs.csv'
Apollolocsfile = 'LanderLocs.csv'

# Lunar ellipsoid defined as tuple of equatorial radius, polar radius
# and flattening
LunarEllipsoid = (1738.1, 1736.0, 0.0012)
LunarMeanRadius = 1737.4

DMQnumber = 33
lander = 'A14'
element = 'ALSEP'

locs_df = pd.read_csv(DMQlocsfile, index_col=0, encoding='utf-8')
lander_df = pd.read_csv(Apollolocsfile, index_col=0, encoding='utf-8')

# print(locs_df)
# print(lander_df)

# print(locs_df.dtypes)
# locs_df = locs_df.astype({'Latitude': 'float64'})
# print(locs_df.dtypes)


slat = locs_df.loc[locs_df['Number'] == DMQnumber]['Latitude'].iloc[0]
slon = locs_df.loc[locs_df['Number'] == DMQnumber]['Longitude'].iloc[0]

# print('slat', slat)
# print('slon', slon)

rlat = lander_df.loc[(lander_df['Lander'] == lander)
                     & (lander_df['Element'] == element)]['Latitude'].iloc[0]
rlon = lander_df.loc[(lander_df['Lander'] == lander)
                     & (lander_df['Element'] == element)]['Longitude'].iloc[0]

# print('rlat', rlat)
# print('rlon', rlon)

geodesicdist = distance.distance((slat, slon), (rlat, rlon),
                                 ellipsoid=LunarEllipsoid)
gcdist = distance.great_circle((slat, slon), (rlat, rlon),
                               radius=LunarMeanRadius)
gcdegrees = units.degrees(gcdist.km/LunarMeanRadius)

print("The distance from A{} to {} {} is:".format(DMQnumber, lander, element))
print("{} km (geodesic distance)".format(geodesicdist.km))
print("{} km (great circle distance)".format(gcdist.km))
print("{} degrees (great circle distance)".format(gcdegrees))



