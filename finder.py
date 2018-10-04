import matplotlib
import numpy as np
from scipy import ndimage

# Set up matplotlib
import matplotlib.pyplot as plt
#matplotlib inline

from astropy.io import fits


def find_outlier_pixels(data,tolerance=3,worry_about_edges=True):
    #This function finds the hot or dead pixels in a 2D dataset.
    #tolerance is the number of standard deviations used to cutoff the hot pixels
    #If you want to ignore the edges and greatly speed up the code, then set
    #worry_about_edges to False.
    #
    #The function returns a list of hot pixels and also an image with with hot pixels removed

    from scipy.ndimage import median_filter
    #blurred = median_filter(Z, size=2)
    #difference = data - blurred
    threshold = 10*np.std(data)

    #find the hot pixels, but ignore the edges
    #hot_pixels = np.nonzero((np.abs(difference[1:-1,1:-1])>threshold) )
    #hot_pixels = np.array(hot_pixels) + 1 #because we ignored the first row and first column

    fixed_image = np.copy(data) #This is the image with the hot pixels removed
    #for y,x in zip(hot_pixels[0],hot_pixels[1]):
    #    fixed_image[y,x]=blurred[y,x]

    if worry_about_edges == True:
        height,width = np.shape(data)

        ###Now get the pixels on the edges (but not the corners)###

        #left and right sides
        for index in range(1,height-1):
            #left side:
            med  = np.median(data[index-1:index+2,0:2])
            diff = np.abs(data[index,0] - med)
            if diff>threshold:
                #hot_pixels = np.hstack(( hot_pixels, [[index],[0]]  ))
                fixed_image[index,0] = med

            #right side:
            med  = np.median(data[index-1:index+2,-2:])
            diff = np.abs(data[index,-1] - med)
            if diff>threshold:
                #hot_pixels = np.hstack(( hot_pixels, [[index],[width-1]]  ))
                fixed_image[index,-1] = med

        #Then the top and bottom
        for index in range(1,width-1):
            #bottom:
            med  = np.median(data[0:2,index-1:index+2])
            diff = np.abs(data[0,index] - med)
            if diff>threshold:
                #hot_pixels = np.hstack(( hot_pixels, [[0],[index]]  ))
                fixed_image[0,index] = med

            #top:
            med  = np.median(data[-2:,index-1:index+2])
            diff = np.abs(data[-1,index] - med)
            if diff>threshold:
                #hot_pixels = np.hstack(( hot_pixels, [[height-1],[index]]  ))
                fixed_image[-1,index] = med

        ###Then the corners###

        #bottom left
        med  = np.median(data[0:2,0:2])
        diff = np.abs(data[0,0] - med)
        if diff>threshold:
            #hot_pixels = np.hstack(( hot_pixels, [[0],[0]]  ))
            fixed_image[0,0] = med

        #bottom right
        med  = np.median(data[0:2,-2:])
        diff = np.abs(data[0,-1] - med)
        if diff>threshold:
            #hot_pixels = np.hstack(( hot_pixels, [[0],[width-1]]  ))
            fixed_image[0,-1] = med

        #top left
        med  = np.median(data[-2:,0:2])
        diff = np.abs(data[-1,0] - med)
        if diff>threshold:
            #hot_pixels = np.hstack(( hot_pixels, [[height-1],[0]]  ))
            fixed_image[-1,0] = med

        #top right
        med  = np.median(data[-2:,-2:])
        diff = np.abs(data[-1,-1] - med)
        if diff>threshold:
            #hot_pixels = np.hstack(( hot_pixels, [[height-1],[width-1]]  ))
            fixed_image[-1,-1] = med

    return fixed_image


def finder(out_dir, image):
    hdu_list = fits.open(image)
    image_data = hdu_list[0].data
    #hdu_list.close()

    image_data = find_outlier_pixels(image_data) #ndimage.gaussian_filter(image_data, 2)

    sx = ndimage.sobel(image_data, axis=0, mode='constant')
    sy = ndimage.sobel(image_data, axis=1, mode='constant')
    sob = np.hypot(sx, sy)

    hdu_list[0].data = image_data
    hdu_list.writeto(out_dir, overwrite=True)

finder('/tmp/credo/sobel.fits', '/home/nkg/no_backup/projects/pk/credo/hackhaton/Dark frames/SC79283.fits')
