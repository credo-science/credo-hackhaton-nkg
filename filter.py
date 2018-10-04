import matplotlib
import numpy as np

# Set up matplotlib
import matplotlib.pyplot as plt
#matplotlib inline

from astropy.io import fits


def get_fits_data(image_file):
    hdu_list = fits.open(image_file)
    image_data = hdu_list[0].data
    hdu_list.close()
    return image_data


def filter(out, in1, in2):

    data1 = get_fits_data(in1)
    data2 = get_fits_data(in2)

    mean1 = np.mean(data1)
    mean2 = np.mean(data2)

    median1 = np.median(data1)
    median2 = np.median(data2)

    min1 = np.min(data1)
    min2 = np.min(data2)

    p1_25 = np.percentile(data1, 0.1)
    p1_75 = np.percentile(data1, 96.5)

    p2_25 = np.percentile(data2, 0.1)
    p2_75 = np.percentile(data2, 96.5)

    max1 = np.max(data1)
    max2 = np.max(data2)

    rs1 = mean1 - min1
    rs2 = mean2 - min2

    for y in range(0, data1.shape[0]):
        for x in range(0, data1.shape[1]):
            v1 = data1[y][x]
            #v1p = data1[y][x+1]
            #v1d = float(v1) - float(v1p)

            v2 = data2[y][x]
            #v2p = data2[y][x+1]
            #v2d = float(v2) - float(v2p)

            d1 = float(v1) - p1_25
            d2 = float(v2) - p2_25

            dn1 = p1_75 - p1_25
            dn2 = p2_75 - p2_25
            rs = float(dn1) / float(dn2)

            norm = (v2 - p2_25) * rs

            preview = norm + p1_25

            data2[y][x] = max(0.0, float(v1) - float(norm))
            #if v1d > 30:
            #    data2[y][x] = v1p
            #else:
            #    data2[y][x] = data1[y][x]
        print('Done line: %d' % y)

    #data2[0][0] = max1

    hdu = fits.open(in1)
    hdu[0].data = data2
    #hdu.header['telescop'] = 'CREDO'
    hdu.writeto(out, overwrite=True)


image_file2 = '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/SC79282.fits'
image_file = '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/SC79283.fits'
