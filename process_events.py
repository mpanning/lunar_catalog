"""
Script to loop through events in catalog and extract only those with
identified nests and amplitude estimates at more than one station and
then calculate distances
"""
import pandas as pd
import numpy as np
from geopy.distance import distance
import matplotlib.pyplot as plt
from statistics import linear_regression
import math

catcsvfile = "levent.csv"

df = pd.read_csv(catcsvfile, index_col=0, encoding='utf-8')

# First down select to only those with identified updated nest
nest_df = df.loc[~df['2004Class'].isna()]

# Pick the ones that also have more than one amplitude pick
select_df = nest_df.loc[((nest_df[['A1112Amp', 'A14Amp', 'A15Amp',
                                   'A16Amp']] > 0).sum(axis=1)) > 1]
# print(select_df)

# Reduce to only needed columns and add in columns for distances
reduce_df = select_df[['A1112Amp', 'A14Amp', 'A15Amp', 'A16Amp', '2004Class']]
reduce_df['2004Class'] = reduce_df['2004Class'].str.slice(1, 4)
reduce_df.insert(1, 'A1112Dist', np.nan)
reduce_df.insert(3, 'A14Dist', np.nan)
reduce_df.insert(5, 'A15Dist', np.nan)
reduce_df.insert(7, 'A16Dist', np.nan)
# print(reduce_df)

# Loop over these events and calculate distances
Apollolocsfile = 'LanderLocs.csv'
DMQlocsfile = 'Nakamura2005Locs.csv'

lander_df = pd.read_csv(Apollolocsfile, index_col=0, encoding='utf-8')
A12lat = lander_df.loc[(lander_df['Lander'] == 'A12') &
                       (lander_df['Element'] == 'ALSEP')]['Latitude'].iloc[0]
A12lon = lander_df.loc[(lander_df['Lander'] == 'A12') &
                       (lander_df['Element'] == 'ALSEP')]['Longitude'].iloc[0]
A14lat = lander_df.loc[(lander_df['Lander'] == 'A14') &
                       (lander_df['Element'] == 'ALSEP')]['Latitude'].iloc[0]
A14lon = lander_df.loc[(lander_df['Lander'] == 'A14') &
                       (lander_df['Element'] == 'ALSEP')]['Longitude'].iloc[0]
A15lat = lander_df.loc[(lander_df['Lander'] == 'A15') &
                       (lander_df['Element'] == 'ALSEP')]['Latitude'].iloc[0]
A15lon = lander_df.loc[(lander_df['Lander'] == 'A15') &
                       (lander_df['Element'] == 'ALSEP')]['Longitude'].iloc[0]
A16lat = lander_df.loc[(lander_df['Lander'] == 'A16') &
                       (lander_df['Element'] == 'ALSEP')]['Latitude'].iloc[0]
A16lon = lander_df.loc[(lander_df['Lander'] == 'A16') &
                       (lander_df['Element'] == 'ALSEP')]['Longitude'].iloc[0]

# Lunar ellipsoid defined as tuple of equatorial radius, polar radius
# and flattening
LunarEllipsoid = (1738.1, 1736.0, 0.0012)
LunarMeanRadius = 1737.4
locs_df = pd.read_csv(DMQlocsfile, index_col=0, encoding='utf-8')

count = 0
drops = []
for i in range(len(reduce_df)):
    DMQind = int(reduce_df.iloc[i, reduce_df.columns.get_loc('2004Class')])
    # Get source latitude and longitude if it's a DMQ in Nakamura (2005)
    try:
        slat = locs_df.loc[locs_df['Number'] == DMQind]['Latitude'].iloc[0]
    except:
        drops.append(i)
        continue
    slon = locs_df.loc[locs_df['Number'] == DMQind]['Longitude'].iloc[0]
    # print(ind, DMQind, slat, slon)

    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A1112Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A1112Dist')] = distance((slat, slon), (A12lat, A12lon), ellipsoid=LunarEllipsoid).km
    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A14Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A14Dist')] = distance((slat, slon), (A14lat, A14lon), ellipsoid=LunarEllipsoid).km
    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A15Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A15Dist')] = distance((slat, slon), (A15lat, A15lon), ellipsoid=LunarEllipsoid).km
    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A16Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A16Dist')] = distance((slat, slon), (A16lat, A16lon), ellipsoid=LunarEllipsoid).km
    count = count + 1
    
    
print("Found {} events with more than one amplitude and a DMQ location".format(count))
# Remove the rows without Nakamura (2005) locations
reduce_df.drop(reduce_df.index[drops], inplace=True)

csvoutfile = 'AmpDistance.csv'
reduce_df.to_csv(csvoutfile, encoding='utf-8')

# Loop again and calculate amplitude ratios within events and append to arrays of distance
# difference and amplitude ratio (assumes simple exponential decay by epicentral distance)
ddiff = []
ampratio = []
for i in range(len(reduce_df)):
    row = reduce_df.iloc[i].dropna() # strips out the missing stations for each event
    namps = int((len(row)-1)/2)
    if namps == 2:
        ddiff.append(row.iloc[3]-row.iloc[1])
        ampratio.append(row.iloc[2]/row.iloc[0])
    elif namps == 3:
        ddiff.append(row.iloc[3]-row.iloc[1])
        ampratio.append(row.iloc[2]/row.iloc[0])
        ddiff.append(row.iloc[5]-row.iloc[1])
        ampratio.append(row.iloc[4]/row.iloc[0])
        ddiff.append(row.iloc[5]-row.iloc[3])
        ampratio.append(row.iloc[4]/row.iloc[2])
    elif namps == 4:
        ddiff.append(row.iloc[3]-row.iloc[1])
        ampratio.append(row.iloc[2]/row.iloc[0])
        ddiff.append(row.iloc[5]-row.iloc[1])
        ampratio.append(row.iloc[4]/row.iloc[0])
        ddiff.append(row.iloc[7]-row.iloc[1])
        ampratio.append(row.iloc[6]/row.iloc[0])
        ddiff.append(row.iloc[5]-row.iloc[3])
        ampratio.append(row.iloc[4]/row.iloc[2])
        ddiff.append(row.iloc[7]-row.iloc[3])
        ampratio.append(row.iloc[6]/row.iloc[2])
        ddiff.append(row.iloc[7]-row.iloc[5])
        ampratio.append(row.iloc[6]/row.iloc[4])
    else:
        raise RuntimeError("Unexpected number of amplitudes")

# Flip the order of any negative ddiff values
for i in range(len(ddiff)):
    if ddiff[i] < 0:
        ddiff[i] = -1*ddiff[i]
        ampratio[i] = 1./ampratio[i]
        
plt.scatter(ddiff, ampratio)
plt.yscale("log")

# Calculate best fit line between ddiff and log(ampratio) constrained
# to pass through origin

x = np.array(ddiff)
x = x[:, np.newaxis]
y = np.log10(np.array(ampratio))
a, _, _, _ = np.linalg.lstsq(x, y)

# xmin = np.amin(x)
xmax = np.amax(x)
ypred = [1.0, math.pow(10, a[0]*xmax)]

# print(a, xmin, xmax)
plt.plot([0, xmax], ypred, color='red')
plt.savefig("Ampratio.png")
        
    
