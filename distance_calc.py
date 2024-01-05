"""
Calculate geodetic distance between moonquake epicenters and Apollo landing
elements
"""
from geopy import distance
import pandas as pd

DMQlocsfile = 'Nakamura2005Locs.csv'
Apollolocsfile = 'LanderLocs.csv'

# Lunar ellipsoid defined as tuple of equatorial radius, polar radius
# and flattening
LunarEllipsoid = (1738.1, 1736.0, 0.0012)

DMQnumber = 33
Apollo = 'A14'
element = 'ALSEP'

locs_df = pd.read_csv(DMQlocsfile, index_col=0, encoding='utf-8')
lander_df = pd.read_csv(Apollolocsfile, index_col=0, encoding='utf-8')

print(locs_df)
print(lander_df)

print(locs_df.dtypes)
# locs_df = locs_df.astype({'Latitude': 'float64'})
# print(locs_df.dtypes)


lat1 = locs_df.loc[locs_df['Number'] == DMQnumber]['Latitude'].iloc[0]
lon1 = locs_df.loc[locs_df['Number'] == DMQnumber]['Longitude'].iloc[0]

print('lat1', lat1)
print('lon1', lon1)
