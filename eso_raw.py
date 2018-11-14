#!/usr/bin/env python

import sys
import math
from datetime import date

from pyvo.dal import tap
import urllib

from astropy.coordinates import SkyCoord
from astropy.units import Quantity

# schema: http://archive.eso.org/tap_obs/tables
ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
tapobs = tap.TAPService(ESO_TAP_OBS)
top = "TOP %d" % 1000
exp = 300


def query(date):
    query = """SELECT %s object, ra, dec, tpl_start, prog_id, filter_path, dp_id, access_url, dp_type, exp_start, exposure
    from dbo.raw
    where dp_type='DARK' and exp_start >= '%sT00:00:00Z' and exp_start <= '%sT23:59:59Z' and exposure > %f """ % (top, date, date, exp)

    res = tapobs.search(query=query)
    return res


def download(rows):
    for row in rows:
        if row["access_url"].decode() != "proprietary":
            fn = row["dp_id"].decode() + '.fits.Z'
            print("Download %s file..." % fn)
            urllib.request.urlretrieve(row["access_url"].decode(), fn)


def from_date_to_date(start_date, end_date):
    from datetime import timedelta

    def daterange(start_date, end_date):
        for n in range(int((end_date - start_date).days)):
            yield start_date + timedelta(n)

    for single_date in daterange(start_date, end_date):
        #print (single_date.strftime("%Y-%m-%d"))
        rows = query(single_date)
        download(rows)


start_date = date(2013, 1, 2)
end_date = date(2015, 6, 2)
from_date_to_date(start_date, end_date)
