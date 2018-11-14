#!/usr/bin/env python
# ------------------------------------------------------------------------------
# Python 3 example on how to query by position the ESO SSA service
# and retrieve the spectra science files with SNR higher than a given threshold.
#
# Name: eso_ssa.py
# Version: 2018-04-09
#
# In case of questions, please send an email to: usd-help@eso.org
# with the following subject: ASG programmatic access (eso_ssa.py)
#
# Useful to read: http://www.astropy.org/astropy-tutorials/Coordinates.html
# ------------------------------------------------------------------------------

import sys

import pyvo as vo
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

# --------------------------------------------------
# Define the end point and the SSA service to query:
# --------------------------------------------------

ssap_endpoint = "http://archive.eso.org/ssap"

ssap_service = vo.dal.SSAService(ssap_endpoint)

# A typical SSA invocation: ssap_endpoint + '?' + 'REQUEST=queryData&POS=197.44888,-23.38333&SIZE=0.5'
# where
#   197.44888 is the right ascension J2000 (decimal degrees),
#   -23.38333 is the declination J2000 (decimal degrees),
#   0.5       is the diameter of the cone to search in (decimal degrees).
# Within astroquery this is achieved by the following lines:
#   mytarget = SkyCoord(197.44888, -23.38333, unit='deg')
#   mysize   = Quantity(0.5, unit="deg")
#   ssap_resultset = ssap_service.search(pos=mytarget, diameter=mysize)
# Really the coordinates are ICRS, but for all practical effects when searching, there is no difference with J2000.


# --------------------------------------------------
# Prepare search parameters
# searching by cone around NGC 4993,
# with a diameter or 0.5 deg (radius=0.25 deg).
# --------------------------------------------------
target = "NGC 4993"
diameter = 0.5

print()
print("Looking for spectra around target %s in a cone of diameter %f deg."   %(target, diameter))
print("Querying the ESO SSAP service at %s" %(ssap_endpoint))

# --------------------------------------------------
# The actual position of the selected target
# is queried by the from_name() function,
# which queries the CDS SESAME service
# (http://cdsweb.u-strasbg.fr/cgi-bin/Sesame).
# --------------------------------------------------

print("The provided target is being resolved by SESAME...")
pos = SkyCoord.from_name(target)
size = Quantity(diameter, unit="deg")
print("SESAME coordinates for %s: %s" % (target, pos.to_string()))

# see: http://docs.astropy.org/en/stable/coordinates/skycoord.html
# In case you know better the coordinates, then use:
# my_icrs_pos = SkyCoord(197.44888, -23.38333, unit='deg')
#
# Or in case you know the galactic coordinates instead:
# my_gal_pos=SkyCoord(308.37745107, 39.29423547, frame='galactic', unit='deg')
#    in which case you will have to use:
# ssap_service.search(pos=my_gal_pos.fk5, diameter=size)
#    given that my_fk5_pos = my_gal_pos.fk5

# --------------------------------------------------
# Query in that defined cone (pos, size):
# --------------------------------------------------
print("Performing a Simple Spectral Access query...")
ssap_resultset = ssap_service.search(pos=pos.fk5, diameter=size)

# NOTE: The ESO coordinate system is: FK5. You would not be off by more than 20mas by querying by pos==pos.icrs (or simply pos=pos) instead.

# --------------------------------------------------
# define the output fields you are interested in;
# uppercase fields are the one defined by the SSAP standard as valid input fields
# --------------------------------------------------
fields = ["COLLECTION", "TARGETNAME", "s_ra", "s_dec", "APERTURE",
          "em_min", "em_max", "SPECRP", "SNR", "t_min", "t_max",
          "CREATORDID", "access_url"]

# --------------------------------------------------
# Print the blank-separated list of fields
# (one line for each of the spectra)
# with the following formatting rules:
# - Do not show the CREATORDID in the stdout:
# - In Python 3, to pretty-print a 'bytes' value, this must be decoded
# - Wavelengths are expressed in meters, for display they are converted to nanometers
# Also, count how many spectra have a SNR > min_SNR
# --------------------------------------------------
min_SNR = 70
count_high_SNR_files=0
separator=' '
for row in ssap_resultset:
   if row["SNR"] > min_SNR:
       count_high_SNR_files += 1
   for field in fields:
       if field == "CREATORDID":
          continue
       value = row[field]
       if isinstance(value, bytes):
          print(value.decode().rjust(16), end=separator)
       elif isinstance(value, str):
          print(value.rjust(16), end=separator)
       else:
          if (field == "em_min" or field == "em_max"):
             value *= 1E9
          print('%16.10f' % (value), end=separator)
   print()

print()

# --------------------------------------------------
# Download those spectra that have SNR > min_SNR
# The name of the file on disk will be the file name
# defined by the creator of such file (field: CREATORDID).
# --------------------------------------------------
import urllib
prompt = "Of the above spectra, download the (%d) ones with SNR > %d, [y|n]:" % (count_high_SNR_files, min_SNR)
shall_I = input(prompt) # Python 3

if shall_I != "y":
    print("Stopping here, without downloading any file")
    sys.exit(0)

print("Downloading files with SNR > %d (if any)" % (min_SNR))

id={}
for row in ssap_resultset:
   if row["SNR"] > min_SNR:
      dp_id = row["dp_id"].decode()
      origfile = row["CREATORDID"].decode()[23:]
      id[origfile] = dp_id
      # The downloaded file is saved with the name provided by the creator of the file: origfile.
      # Though, care should be taken, because reduced products
      # generated by external users might have colliding CREATORDID!
      # This "demo" script does not take into consideration this risk.
      print("Fetching file with SNR=%f: %s.fits renamed to %s" %(row["SNR"], dp_id, origfile))
      urllib.request.urlretrieve(row["access_url"].decode(), row["CREATORDID"].decode()[23:])

print("End of execution")
