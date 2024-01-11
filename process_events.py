"""
Script to loop through events in catalog and extract only those with
identified nests and amplitude estimates at more than one station and
then calculate distances
"""
import pandas as pd
import numpy as np
from geopy.distance import distance
import matplotlib.pyplot as plt
import math

def spheuclid(lat1, lon1, r1, lat2, lon2, r2):
    theta1 = math.radians(90.0-lat1)
    phi1 = math.radians(lon1)
    theta2 = math.radians(90.0-lat2)
    phi2 = math.radians(lon2)
    distance = math.sqrt(r1**2 + r2**2
                         - 2.0*r1*r2*(math.sin(theta1)*math.sin(theta2)
                                      *math.cos(phi1-phi2)
                                      + math.cos(theta1)*math.cos(theta2)))
    return distance
    
catcsvfile = "levent.csv"
ifNormalize = True
pltxlims = [0, 1250]
pltylims = [0.02, 30]
euclidxlims = [0, 650]
euclidylims = [0.03, 30]

df = pd.read_csv(catcsvfile, index_col=0, encoding='utf-8')

# First down select to only those with identified updated nest
nest_df = df.loc[~df['2004Class'].isna()]

# Pick the ones that also have more than one amplitude pick
select_df = nest_df.loc[((nest_df[['A1112Amp', 'A14Amp', 'A15Amp',
                                   'A16Amp']] > 0).sum(axis=1)) > 1]
# print(select_df)
StationMeans = select_df.loc[:,['A1112Amp', 'A14Amp', 'A15Amp', 'A16Amp']].mean()
print("Station mean amplitude values")
print(StationMeans)

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

if ifNormalize:
    reduce_df["A1112Amp"] = reduce_df["A1112Amp"].div(StationMeans.loc["A1112Amp"])
    reduce_df["A14Amp"] = reduce_df["A14Amp"].div(StationMeans.loc["A14Amp"])
    reduce_df["A15Amp"] = reduce_df["A15Amp"].div(StationMeans.loc["A15Amp"])
    reduce_df["A16Amp"] = reduce_df["A16Amp"].div(StationMeans.loc["A16Amp"])
    csvoutfile = 'AmpDistanceNorm.csv'
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
print("Proportionality factor: {}".format(a[0]))

# xmin = np.amin(x)
xmax = np.amax(x)
ypred = [1.0, math.pow(10, a[0]*xmax)]
plt.plot([0, xmax], ypred, color='red')
plt.xlim(pltxlims[0], pltxlims[1])
plt.ylim(pltylims[0], pltylims[1])

# print(a, xmin, xmax)
if ifNormalize:
    pngfile = 'AmpratioNorm.png'
else:
    pngfile = 'Ampratio.png'
plt.savefig(pngfile)
        
# Try again with euclidean distance for those with depths    
count = 0
drops = []
for i in range(len(reduce_df)):
    DMQind = int(reduce_df.iloc[i, reduce_df.columns.get_loc('2004Class')])
    # Get source latitude and longitude if it's a DMQ in Nakamura (2005)
    try:
        depth = locs_df.loc[locs_df['Number'] == DMQind]['Depth'].iloc[0]
    except:
        drops.append(i)
        continue
    if math.isnan(depth):
        drops.append(i)
        continue
    slat = locs_df.loc[locs_df['Number'] == DMQind]['Latitude'].iloc[0]
    slon = locs_df.loc[locs_df['Number'] == DMQind]['Longitude'].iloc[0]
    srad = LunarMeanRadius - depth

    
    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A1112Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A1112Dist')] = spheuclid(slat, slon, srad, A12lat, A12lon, LunarMeanRadius)
    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A14Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A14Dist')] = spheuclid(slat, slon, srad, A14lat, A14lon, LunarMeanRadius)
    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A15Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A15Dist')] = spheuclid(slat, slon, srad, A15lat, A15lon, LunarMeanRadius)
    if (reduce_df.iloc[i, reduce_df.columns.get_loc('A16Amp')] > 0):
        reduce_df.iloc[i, reduce_df.columns.get_loc('A16Dist')] = spheuclid(slat, slon, srad, A16lat, A16lon, LunarMeanRadius)
    count = count + 1

print("Found {} events with more than one amplitude and a DMQ depth".format(count))
# Remove the rows without Nakamura (2005) depths
reduce_df.drop(reduce_df.index[drops], inplace=True)

csvoutfile = 'AmpDistanceEuclid.csv'
reduce_df.to_csv(csvoutfile, encoding='utf-8')

if ifNormalize:
    reduce_df["A1112Amp"] = reduce_df["A1112Amp"].div(StationMeans.loc["A1112Amp"])
    reduce_df["A14Amp"] = reduce_df["A14Amp"].div(StationMeans.loc["A14Amp"])
    reduce_df["A15Amp"] = reduce_df["A15Amp"].div(StationMeans.loc["A15Amp"])
    reduce_df["A16Amp"] = reduce_df["A16Amp"].div(StationMeans.loc["A16Amp"])
    csvoutfile = 'AmpDistanceNormEuclid.csv'
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
        print(row)
        raise RuntimeError("Unexpected number of amplitudes")

# Flip the order of any negative ddiff values
for i in range(len(ddiff)):
    if ddiff[i] < 0:
        ddiff[i] = -1*ddiff[i]
        ampratio[i] = 1./ampratio[i]
        
plt.clf()
plt.scatter(ddiff, ampratio, label="Nakamura amplitude ratios")
plt.yscale("log")

# Calculate best fit line between ddiff and log(ampratio) constrained
# to pass through origin

x = np.array(ddiff)
x = x[:, np.newaxis]
y = np.log10(np.array(ampratio))
a, _, _, _ = np.linalg.lstsq(x, y)
print("Proportionality factor: {}".format(a[0]))

# xmin = np.amin(x)
xmax = np.amax(x)
ypred = [1.0, math.pow(10, a[0]*xmax)]
plt.plot([0, xmax], ypred, color='red', label="Best fit exponential function")
plt.legend()
plt.xlabel("Difference in source to station distance (km)")
plt.ylabel("Amplitude ratio between stations")
plt.xlim(euclidxlims[0], euclidxlims[1])
plt.ylim(euclidylims[0], euclidylims[1])

# print(a, xmin, xmax)
if ifNormalize:
    pngfile = 'AmpratioNormEuclid.png'
else:
    pngfile = 'AmpratioEuclid.png'
plt.savefig(pngfile, dpi=300)
