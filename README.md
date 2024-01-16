# lunar_catalog
Codes for using the Nakamura lunar event catalog

The catalog itself is in the subdirectory `UTIGTR_0018`. This is simply the unzipped version of the zipfile of the long period event catalog available as a [technical report from the University of Texas Institute of Geophysics](https://repositories.lib.utexas.edu/items/3c162b9a-ce04-4900-8f6e-63147eb7d4db). Details are available through that link and in the file `UTIGTR_0018/UTIGTR_0018.pdf`.

Full citation for the catalog is
> Nakamura, Y., Latham, G.V., Dorman, H.J., and Harris, J.E. "Passive Seismic Experiment, Long Period Event Catalog, Final Version (1969 Day 202 - 1977 Day 273, ALSEP Stations 11, 12, 13, 14, 15, and 16)." University of Texas Institute for Geophysics Technical Report No. 18 (originally published 19 June 1981, revised 2 October 2008), 2p plus appendices.

The repository also contains deep moonquake nest locations directly copied and pasted from [Nakamura (2005, Tables 2 and 3)](https://doi.org/10.1029/2004JE002332) in the csv file `Nakamura2005Locs_cp.csv`. The full citation for those locations are
> Nakamura, Y. (2005), Farside deep moonquakes and deep interior of the Moon, J. Geophys. Res., 110, E01001, doi:10.1029/2004JE002332.

Finally, the repository also includes coordinates of the Apollo lander elements in `LanderLocs.csv`. These are sourced from the [NASA Space Science Data Coordinated Archive (NSSDCA)](https://nssdc.gsfc.nasa.gov/planetary/lunar/lunar_sites.html). Currently, the codes only used the ALSEP locations from Apollos 12, 14, 15, and 16.

The code `read_catalog.py` reads the original fixed width format catalog file (`UTIGTR_0018/levent.1008.dat`) and outputs it into a csv file (`levent.csv`) for easier use in other python processing. It also reads the copy and pasted location file, `Nakamura2005Locs_cp.csv` and strips out some formatting to put it in an easier format as `Nakamura2005Locs.csv`.

The code `distance_calc.py` is just a simple code to calculate geodetic and great-circle distances between a deep moonquake (DMQ) nest and a specific Apollo ALSEP station and output it to the screen. The output is not used in other codes at this time. 
