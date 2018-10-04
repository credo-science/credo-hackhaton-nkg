import os
import pathlib
from math import degrees

import ccdproc
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from astropy import units as u
from astropy.convolution import Gaussian2DKernel

from astropy.nddata import CCDData

from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import ImageNormalize, SqrtStretch
from photutils import detect_sources, detect_threshold, source_properties, deblend_sources, EllipticalAperture
from scipy.misc import imrotate

image_file2 = '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/SC79282.fits'
image_file = '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/SC79283.fits'


def save_as_png(img, thr1, thr2, png):
    d = thr2 - thr1
    if d == 0:
        factor = 256 / thr1
    else:
        factor = 256 / d

    i8 = img - thr1
    i8 = i8 * factor
    i8[i8 < 0] = 0
    i8[i8 >= 256*128] = 0
    i8[i8 >= 255] = 255

    img = Image.fromarray(i8.astype(np.uint8))
    img.save(png)


def filter(out_filter, out_marked, out_dir, in1, in2, desc):

    image = CCDData.read(in1, unit="adu")
    #save_as_png(image.data, 800, 1200, '/tmp/a.png')

    dark = CCDData.read(in2, unit="adu")

    dark_sub = ccdproc.subtract_dark(image, dark, dark_exposure=240*u.second, data_exposure=240*u.second)

    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    threshold = detect_threshold(dark_sub, snr=2.)
    segm = detect_sources(dark_sub, threshold, npixels=5, filter_kernel=kernel)

    hdu = fits.open(in1)
    hdu[0].data = dark_sub
    #hdu.header['telescop'] = 'CREDO'
    hdu.writeto(out_filter, overwrite=True)
    hdu.close()

    hdu = fits.open(in1)
    hdu[0].data = segm
    #hdu.header['telescop'] = 'CREDO'
    hdu.writeto(out_marked, overwrite=True)
    hdu.close()

    thr1 = np.percentile(dark_sub, 25)
    #thr2 = np.max(dark_sub)
    thr2 = np.percentile(dark_sub, 99.999)

    save_as_png(dark_sub, thr1, thr2, out_filter + ".png")

    npixels = 5
    segm_deblend = deblend_sources(dark_sub, segm, npixels=npixels,
                    filter_kernel=kernel, nlevels=32,
                    contrast=0.001)

    cat = source_properties(segm, segm_deblend)
    r = 3.  # approximate isophotal extent

    dots_count = 0
    comet_count = 0
    worm_count = 0
    group_count = 0

    dots_area = 0
    comet_area = 0
    worm_area = 0
    group_area = 0


    for obj in cat:
        position = (obj.xcentroid.value, obj.ycentroid.value)
        x, y = position
        cropx, cropy = 60, 60
        startx = int(x - (cropx // 2))
        starty = int(y - (cropy // 2))
        #crop = dark_sub[10:60,10:60]
        crop = dark_sub[max(0, starty):min(starty + cropy, dark_sub.data.shape[1]), max(startx, 0):min(startx + cropx, dark_sub.data.shape[0])]
        area = obj.area.value

        fn = '%s/%04d-%04dx%04d' % (out_dir, int(area), int(position[0]), int(position[1]))
        fnc = '%s/class/%04d-%04dx%04d' % (out_dir, int(area), int(position[0]), int(position[1]))

        hdu = fits.open(in1)
        hdu[0].data = crop
        hdu.writeto(fn + '.fits', overwrite=True)
        hdu.close()

        clsasse = 'dot'
        group = False
        if obj.ellipticity.value >= 0.6:
            clsasse = 'comet'
            comet_area += obj.area.value
            comet_count += 1
        elif obj.ellipticity.value >= 0.2:
            clsasse = 'worm'
            worm_area += obj.area.value
            worm_count += 1
        else:
            dots_area += obj.area.value
            dots_count += 1

        for i in cat:
            xx = i.xcentroid.value
            yy = i.ycentroid.value
            if pow(xx - x, 2) + pow(yy - y, 2) < pow(30, 2) and i.id != obj.id:
                group = True
                group_area += obj.area.value
                group_count += 1
                break

        suffix = "-%s-%s-%.3f-%.3f-%.3f.png" % (clsasse, 'g' if group else 'n', obj.ellipticity.value, obj.elongation.value, obj.eccentricity.value)

        save_as_png(crop, thr1, thr2, fn + suffix)

        #cl = segm.data[int(obj.ymin.value):int(obj.ymax.value), int(obj.xmin.value):int(obj.xmax.value)]
        save_as_png(obj.data_cutout, np.min(obj.data_cutout), np.max(obj.data_cutout), fnc + suffix)
        #rotated = imrotate(obj.data_cutout, degrees(obj.orientation.value))
        #save_as_png(rotated, np.min(rotated), np.max(rotated), fnc + "-rotated.png")

        #if obj.area.value >= 39:
        #    print('big')

        #a = obj.semimajor_axis_sigma.value * r
        #b = obj.semiminor_axis_sigma.value * r
        #theta = obj.orientation.value
        #print("x: %d, y: %d" % (int(position[0]), int(position[1])))
        #apertures.append(EllipticalAperture(position, a, b, theta=theta))
    count = dots_count + comet_count + worm_count
    area = dots_area + comet_area + worm_area
    print('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (desc, dots_count, dots_area, comet_count, comet_area, worm_count, worm_area, count, area, group_count, group_area))


files = [
['SC79282.fits',
'SC79283.fits',
'SC79284.fits',
'SC79285.fits',
'SC79286.fits',
'SC79287.fits',
'SC79288.fits',
'SC79289.fits',
'SC79288.fits'],

# ['SC79297.fits',
# 'SC79298.fits',
# 'SC79299.fits',
# 'SC79300.fits',
# 'SC79301.fits',
# 'SC79302.fits',
# 'SC79303.fits',
# 'SC79304.fits',
# 'SC79305.fits',
# 'SC79306.fits',
# 'SC79305.fits'],

# ['SC79401.fits',
# 'SC79402.fits',
# 'SC79403.fits',
# 'SC79404.fits',
# 'SC79405.fits',
# 'SC79406.fits',
# 'SC79405.fits'],

['SC79423.fits',
'SC79424.fits',
'SC79425.fits',
'SC79426.fits',
'SC79427.fits',
'SC79428.fits',
'SC79429.fits',
'SC79430.fits',
'SC79431.fits',
'SC79432.fits',
'SC79433.fits',
'SC79432.fits'],
#
# ['SC79494.fits',
# 'SC79495.fits',
# 'SC79496.fits',
# 'SC79497.fits',
# 'SC79498.fits',
# 'SC79499.fits',
# 'SC79498.fits'],

['SC79583.fits',
'SC79584.fits',
'SC79585.fits',
'SC79586.fits',
'SC79587.fits',
'SC79588.fits',
'SC79589.fits',
'SC79590.fits',
'SC79591.fits',
'SC79590.fits']
]

#image_file2 = '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/SC79282.fits'
#image_file = '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/SC79283.fits'
image_dir = '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/'

print('#file\tdots count\tdots area\tcomet count\tcomet area\tworms count\tworms area\tsum count\tsum area\tgroup count\tgroup area')

pathlib.Path('/tmp/credo/detections/').mkdir(parents=True, exist_ok=True)
for g in files:
    for i in range(0, len(g) - 1):
        fn = g[i]
        dir = '/tmp/credo/detections/' + fn + "/"
        pathlib.Path(dir + "class/").mkdir(parents=True, exist_ok=True)
        filter('/tmp/credo/cleared-' + fn, '/tmp/credo/marked-' + fn, dir, image_dir + fn, image_dir + g[i+1], fn)

#filter('/tmp/credo/cleared.fits', '/tmp/credo/marked.fits', '/tmp/credo/detections', image_file, image_file2)
