"""
Code to calculate the expected reduction in event detections at the poles
relative to the detections within the Apollo footprint
"""

import pandas as pd
import numpy as np
import math

# function to convert lat lon on unit sphere to x,y,z
def latlon2xyz(lat, lon):
    theta = math.radians(90.0-lat)
    phi = math.radians(lon)
    x = math.sin(theta)*math.cos(phi)
    y = math.sin(theta)*math.sin(phi)
    z = math.cos(theta)
    return (x, y, z)

# function to get Euclidean distance between two points in spherical coords
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


Apollolocsfile = 'LanderLocs.csv'
DMQlocsfile = 'Nakamura2005Locs.csv'
LunarMeanRadius = 1737.4

# Values for detection reduction calculation
C_est = -0.0005398 # Slope from log-linear fit
bvalue = 1.0


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

# Find the center of the triangle defined by A12, A15, and A16
A12xyz = np.array(latlon2xyz(A12lat, A12lon))
A15xyz = np.array(latlon2xyz(A15lat, A15lon))
A16xyz = np.array(latlon2xyz(A16lat, A16lon))

centerxyz = np.cross((A12xyz-A15xyz), (A16xyz-A15xyz))
centerlat = math.degrees(math.asin(centerxyz[2]/math.sqrt(centerxyz[0]**2 +
                                                          centerxyz[1]**2 +
                                                          centerxyz[2]**2)))
centerlon = math.degrees(math.asin(centerxyz[1]/math.sqrt(centerxyz[0]**2 +
                                                          centerxyz[1]**2)))
# print(centerlat, centerlon)

# Find the mean depth of the events
locs_df = pd.read_csv(DMQlocsfile, index_col=0, encoding='utf-8')
meanDepth = locs_df.loc[:,['Depth']].mean().iloc[0]

# print(meanDepth)

# Find average distance to Apollo stations from center of network at DMQ depth
A12tocenter = spheuclid(A12lat, A12lon, LunarMeanRadius, centerlat, centerlon,
                        LunarMeanRadius-meanDepth)
A14tocenter = spheuclid(A14lat, A14lon, LunarMeanRadius, centerlat, centerlon,
                        LunarMeanRadius-meanDepth)
A15tocenter = spheuclid(A15lat, A15lon, LunarMeanRadius, centerlat, centerlon,
                        LunarMeanRadius-meanDepth)
A16tocenter = spheuclid(A16lat, A16lon, LunarMeanRadius, centerlat, centerlon,
                        LunarMeanRadius-meanDepth)
NetDistance = 0.25*(A12tocenter + A14tocenter + A15tocenter + A16tocenter)
print("Average distance to center of network: {}".format(NetDistance))

# Find average distance to Apollo stations from north pole at DMQ depth
A12tonorthpole = spheuclid(A12lat, A12lon, LunarMeanRadius, 90., 0.,
                           LunarMeanRadius-meanDepth)
A14tonorthpole = spheuclid(A14lat, A14lon, LunarMeanRadius, 90., 0.,
                           LunarMeanRadius-meanDepth)
A15tonorthpole = spheuclid(A15lat, A15lon, LunarMeanRadius, 90., 0.,
                           LunarMeanRadius-meanDepth)
A16tonorthpole = spheuclid(A16lat, A16lon, LunarMeanRadius, 90., 0.,
                           LunarMeanRadius-meanDepth)
NpoleDistance = 0.25*(A12tonorthpole + A14tonorthpole + A15tonorthpole
                      + A16tonorthpole)
print("Average distance to north pole: {}".format(NpoleDistance))

# Find average distance to Apollo stations from south pole at DMQ depth
A12tosouthpole = spheuclid(A12lat, A12lon, LunarMeanRadius, -90., 0.,
                           LunarMeanRadius-meanDepth)
A14tosouthpole = spheuclid(A14lat, A14lon, LunarMeanRadius, -90., 0.,
                           LunarMeanRadius-meanDepth)
A15tosouthpole = spheuclid(A15lat, A15lon, LunarMeanRadius, -90., 0.,
                           LunarMeanRadius-meanDepth)
A16tosouthpole = spheuclid(A16lat, A16lon, LunarMeanRadius, -90., 0.,
                           LunarMeanRadius-meanDepth)
SpoleDistance = 0.25*(A12tosouthpole + A14tosouthpole + A15tosouthpole
                      + A16tosouthpole)
print("Average distance to south pole: {}".format(SpoleDistance))

R_npole = math.exp(C_est*(NpoleDistance - NetDistance))
R_spole = math.exp(C_est*(SpoleDistance - NetDistance))

print("Relative amplitude reduction: {} (north) {} (south)".format(R_npole,
                                                                   R_spole))

Rpole = 0.5*(R_npole + R_spole)
DetectionFactor = 1./math.pow(Rpole,-2.*bvalue/3.)

print("Reduction factor for pole detections: {}".format(DetectionFactor))
