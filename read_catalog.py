"""
Read Yosio's catalog into a pandas dataframe and output it into a csv file for
easier future processing
"""
import pandas as pd

catfile = "UTIGTR_0018/levent.1008.dat"
catcsvfile = "levent.csv"
colspecs = [(2,4), (5,8), (9,13), (14,18), (19,23), (23,27), (27,31), (31,35),
            (36,40), (41,45), (46,76), (76,77), (77,80), (81, 85)]
names = ['Year', 'DOY', 'StartTime', 'StopTime', 'A1112Amp', 'A14Amp',
         'A15Amp', 'A16Amp', 'Availability', 'Quality', 'Comments',
         'EventType', 'DeepClass', '2004Class']

events_df = pd.read_fwf(catfile, colspecs=colspecs, header=None, names=names)

print(events_df)

events_df.to_csv(catcsvfile, encoding='utf-8')

# Select only rows corresponding to events in known nests
# dfA = df.loc[df['EventType'] == 'A']

# print(dfA.to_string())

# Test read csv file

# df2 = pd.read_csv(csvfile, index_col=0, encoding='utf-8')

# print(df2)

# Read copy and paste version of Nakamura (2005) locations tables

locscsvfile = 'Nakamura2005Locs_cp.csv'
locs_df = pd.read_csv(locscsvfile, index_col=False, encoding='utf-8')

# print(locs_df)

# Move uncertainty values to own columns
locs_df[['Latitude', 'LatSigma']] = locs_df['Latitude'].str.split(' ± ',
                                                                  expand=True)
locs_df[['Longitude', 'LonSigma']] = locs_df['Longitude'].str.split(' ± ',
                                                                  expand=True)
locs_df[['Depth', 'DepSigma']] = locs_df['Depth'].str.split(' ± ',
                                                                  expand=True)
# Reorder to put uncertainties after values
locs_df = locs_df[['Number', 'Latitude', 'LatSigma', 'Longitude', 'LonSigma',
                   'Depth', 'DepSigma']]
# Fix negative sign formatting
locs_df['Latitude'] = locs_df['Latitude'].str.replace('−','-')
locs_df['Longitude'] = locs_df['Longitude'].str.replace('−','-')
print(locs_df)

# Output to csv file
locscsvout = 'Nakamura2005Locs.csv'
locs_df.to_csv(locscsvout, encoding='utf-8')
