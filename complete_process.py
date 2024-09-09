'''
PLEASE MAKE SURE THAT YOU HAVE COPIED THE SCRIPT TO YOUR GROUP'S DIRECTORY
AND ARE WORKING ON THE COPIED VERSION

This is the incomplete version of the calibration script that you will be
using to process all of your observatory data. You will need to finish
this script and submit it to the TA's as a group. Individually, you are
required to write an ipython notebook that calls functions from this script
to reduce your week 1 observations as described in Lab 1.1. While developing
the functions in this script you may want to run them on real data to see if
the output products make sense.

N.B. We have retained the implementation that allows for your darks, flats,
and science images to have different exposure times assuming you do the
correct scaling. While we will aim to only input observations of the same
exposure time, this will allow your code to be generalizable, both in the
case of traditional CCD astronomy (where darks can easily be scaled) as well
as in ht event that we need to improvise.
'''
# Here are the libraries we need. Note that wherever you see np, that
# stands for Numpy. 'import' loads an external library.

import numpy as np
from astropy.io import fits
import astropy.stats as stats


def average_bias(biasfiles):
    """ Produce an average bias image from a list of individual bias exposures.

    Args:
        biasfiles ([str,...]): A list of string specifying the path to the fits
            files with the bias frames.

    Returns:
        (np.array): The average bias frame. Units of counts.
    """

    # Open each file using fitsio. Fill this in!
    biasdata= [fits.getdata(filepath) for filepath in biasfiles]
    biascube= np.asarray(biasdata) # Use numpy to turn the biasdata list into an array

    averagebias = np.mean(biascube, axis=0) # Use numpy to average over the first axis. What kind
    # of average should we use? [WE CHOSE MEAN!]

    return averagebias


def average_dark(darkfiles,averagebias):
    """ Produce an average dark image from a list of individual dark exposures.

    Args:
        darkfiles ([str,...]): A list of string specifying the path to the fits
            files with the dark frames.
        averagebias (np.array): The average bias array returned from
            average_bias.

    Returns:
        (np.array): The average dark frame. Units of counts per second.
    """

    # Read in the darks (very similar to the bias version)
    darkdata = [fits.getdata(filepath) for filepath in darkfiles]

    # Make a list of exposure times (hint: you'll need to use the header)
    darkexpo = [fits.getheader(i)['EXPTIME'] for i in darkfiles]

    # Turn both lists into numpy arrays
    darkcube = np.asarray(darkdata)
    darkexpocube = np.asarray(darkexpo)

    darklist=[]
    # The zip command here loops over two lists at the same time. Iterating
    # over dark cube should return individual darks and iterating over
    # darkexpocube should return times.
    for (image,time) in zip(darkcube,darkexpocube):

        # Correct the dark for the bias. How do you do this?
        cleandark = image - averagebias
        # Normalize the dark for the exposure time. How do you do this?
        normdark = cleandark / time

        # We'll append this normalized dark to our list.
        darklist.append(normdark)

    # As before, we're going to want to turn our list into an array
    cleandarkcube = np.asarray(darklist)

    # How do we want to average these darks?
    averagedark = np.mean(cleandarkcube, axis=0)
  
    return averagedark


def average_flat(flatfiles,averagebias,averagedark):
    """ Produce an average flat image from a list of individual dark exposures.

    Args:
        flatfiles ([str,...]): A list of string specifying the path to the fits
            files with the flat frames.
        averagebias (np.array): The average bias array returned from
            average_bias.
        averagedark (np.array): The average dark array returned from
            average_dark.

    Returns:
        (np.array): The average flat frame. No units.
"""

    # As before, read in the flats from the fits
    flatdata=[fits.getdata(i) for i in flatfiles]

    # We'll also need the exposure time for each flat.
    flatexpo=[fits.getheader(i)['EXPTIME'] for i in flatfiles]

    # Make them into numpy arrays.
    flatcube=np.asarray(flatdata)
    flatexpocube=np.asarray(flatexpo)

    # Now we'll iterate through the flats like before and make our list
    # of cleaned and normalized flats.
    cleanflatlist=[]

    for (image,time) in zip(flatcube,flatexpocube):
      # First we need to correct the flat for the bias and dark. How?
      cleanflat = image - averagedark*time - averagebias
      # Now we need to normalize the flat. What normalization is this?
      normflat = cleanflat / stats.sigma_clipped_stats(cleanflat,sigma=3)[0] #TODO: meep why is sigma clip weird?
      #normflat= cleanflat / np.mean(cleanflat, axis=0)

      # Append this to our list
      cleanflatlist.append(normflat)

    # Turn our list into a cube
    cleancube = np.asarray(cleanflatlist)
    # What kind of average should we take? (This is different than Lab0! You
    # may need to explore more robust algorithms for cleaning/clipping data)
    averageflat = np.mean(cleancube, axis=0)
    #averageflat = stats.sigma_clipped_stats(cleancube, sigma=3)[0] #TODO: this somehow won't work...
    print(averageflat)

    return np.asarray(averageflat)


def science_exposure(rawscidata,rawsciheader,averagebias,averagedark,
    averageflat):
    """ Produce an science frame by combining all the calibration steps.

    Args:
        rawscidata (np.array): The raw science data
        rawsciheader (FITS header; Dictionary): A fits header file for the
            raw science image.
        averagebias (np.array): The average bias array returned from
            average_bias.
        averagedark (np.array): The normalized, average dark array returned
            from average_dark.
        averageflat (np.array): The normalized, average flat array returned
            from average_flat.

    Returns:
        (np.array): The science image.
    """

    # Grab the exposure time from the header
    expotime = rawsciheader['EXPTIME']

    # Fill in the formula from image calibration from the parts you have
    scienceimage = (rawscidata - (averagebias + averagedark*expotime)) / averageflat
    return scienceimage

# Copied from work in notebook lab 1.3!

def create_coadd_image(sciencefilelist, ref_image, averagebias, averagedark,
    averageflat):
    """ Produce a single coadded image from a list of exposures.

    Args:
        sciencefilelist ([str,...]): A list of strings containing the
            filenames for each of the exposure to coadd.
        ref_image (np.array): An image to use as the reference image when
            aligning the coadds.
        averagebias (np.array): The average bias array returned from
            average_bias.
        averagedark (np.array): The normalized, average dark array returned
            from average_dark.
        averageflat (np.array): The normalized, average flat array returned
            from average_flat.

    Returns:
        (np.array): The coadded image.

    Notes:
        Not all the images you want to coadd will always have the same exposure
        time. What is the best way to deal with this? Also, choose your
        reference image wisely. The first or last image in a series of exposures
        is probably not a good choice.
    """

    # Get x and y dimensions of the reference image
    xsz, ysz = ref_image.shape
    
    # Number of images to coadd
    calct = len(sciencefilelist)
    
    # Create an empty stack of images
    imstack = np.zeros((xsz, ysz, calct))
    
    # exptime variable
    tot_exptime = 0

    for i, filename in enumerate(sciencefilelist):
        # Read image
        raw_hdu = fits.open(filename)
        raw_data =raw_hdu[0].data
        sci_header = raw_hdu[0].header
        final_image = science_exposure(raw_data,sci_header,averagebias,averagedark,
    averageflat)
        
        # Reproject image to match the reference frame
        im_match, footprint = reproject_interp(final_image, ref_image.header)
        
        # Get the exposure time
        exptime = sci_header['EXPTIME']
        
        # Update the total exposure time
        tot_exptime += exptime
        
        # Put it in the stack, assuming they are aligned
        imstack[:, :, i] = im_match

        # Print progress 
        if i % 5 == 0:
            print(f"Processed {i} out of {calct}.")

    # Take median along the third axis to obtain the final coadded image
    coadd = np.median(imstack, axis=2)
    
    # Update header information # Q: Do we need to do this??
   # headref = ref_image.header
   # headref['EXPTIME'] = tot_exptime
   # headref['EXPOSURE'] = tot_exptime

    return coadd