#!/usr/bin/env python
# -*- coding: utf-8 -*-
# wersja: 09 VII 2022

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import argcomplete

parser = ArgumentParser(description='Convert FITS to png', epilog='Tab completion is supported', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", "--fname",  type=str,                            help="Low brightness limit",               required=True)
parser.add_argument("-c", "--center", type=float,                 nargs=2, help="Cutout center [RA [deg] Dec [deg]]", required=True)
parser.add_argument("--vmin",         type=float, default=-0.001,          help="Low brightness limit")
parser.add_argument("--vmax",         type=float, default=0.01,            help="High brightness limit")
parser.add_argument("-r", "--ra",     type=float, default=0.5,             help="RA  cutout size [deg]. 0 here and in --dec to keep original size")
parser.add_argument("-d", "--dec",    type=float, default=0.5,             help="Dec cutout size [deg]. 0 here and in --ra  to keep original size")
argcomplete.autocomplete(parser)
args = parser.parse_args()


from astropy.io import fits
import numpy as np
from imageio import imwrite
hdul = fits.open(args.fname)
image_data = hdul[0].data
image_data = np.squeeze(image_data)

fNameSize = len(args.fname)
fNamePart = args.fname[:fNameSize - 5]
pngFname = fNamePart + "_cutout.png"

# bez kadrowania (ra, dec = 0)
if (args.ra==0) and (args.dec==0):
    image_data[image_data > args.vmax] = args.vmax
    image_data[image_data < args.vmin] = args.vmin
    # Scale image_data to range [0, 1] 
    image_data = (image_data - args.vmin)/(args.vmax - args.vmin)
    # Convert to 16-bit integer  
    image_data = ((2**16-1)*image_data).astype(np.uint16)
    # Invert y axis
    image_data = image_data[::-1, :]
    imwrite(pngFname, image_data, compression=0)

# z kadrowaniem
else:
    from astropy.nddata import Cutout2D
    from astropy import wcs, coordinates, units as u
    import sys
    w = wcs.WCS(hdul[0].header)
    centerOfCutout = coordinates.SkyCoord(args.center[0]*u.deg, args.center[1]*u.deg, frame='fk5')
    size=[args.ra, args.dec]*u.deg
    try:
        cutout = Cutout2D(image_data, centerOfCutout, size, wcs=w.celestial)
    except ValueError as e:
        print(str(e))
        sys.exit('Zadany wycinek prawdopodobnie wychodzi poza obszar mapy.')

    fitsFname = fNamePart + "_cutout.fits"
    hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
    hdu.writeto(fitsFname, overwrite=True)

    image_data=cutout.data
    image_data[image_data > args.vmax] = args.vmax
    image_data[image_data < args.vmin] = args.vmin
    image_data = (image_data - args.vmin)/(args.vmax - args.vmin) 
    image_data = ((2**16-1)*image_data).astype(np.uint16)
    image_data = image_data[::-1, :]
    imwrite(pngFname, image_data, compression=0)


try:
    import matplotlib.pyplot as plt
    plt.imshow(image_data)
    pltFname = fNamePart + "_cutout_matplot.png"
    plt.imsave(pltFname, image_data)
    plt.show()
