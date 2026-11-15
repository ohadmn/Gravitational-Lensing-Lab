# Data

Put the photometry files here, one per event, named as in the table at the
top of `main.m` (e.g. `SMC001.dat`, `BLG002.dat`).

The expected format is the OGLE Early Warning System `phot.dat` file:
whitespace-separated columns, of which the first three are used

| column | content                         |
|--------|---------------------------------|
| 1      | time (HJD, days)                |
| 2      | I-band magnitude                |
| 3      | magnitude error                 |

The events analysed in the lab are from the OGLE-III EWS 2008 season
(<https://ogle.astrouw.edu.pl/ogle3/ews/ews.html>). The data files are not
redistributed in this repository.
