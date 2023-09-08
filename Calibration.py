'''
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SALT RSS Calibration

This program performs the flux calibration and heliocentric velocity correction
for a FITS file reduced by the SALT RSS pipeline.

The target FITS file is typically a 2D flux spectrum of a galaxy, but may but may be other astronomical objects.

*NOTE*: This requires a pre-existing sensitivity or sens FITS file.

The SALT RSS Data Reduction procedure is described in:
http://mips.as.arizona.edu/~khainline/salt_redux.html

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
'''

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Import Python Libraries
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

import os  # For bash commands
from pyraf import iraf  # For IRAF commands in Python
import astropy.io.fits as fits  # For FITS file handling

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Load IRAF Libraries
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

iraf.noao()
iraf.noao.twod()
iraf.noao.twod.longslit()
iraf.imutil()

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Flux Calibration
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''


def flux(file, err_file, basename):
    '''
    Applies the flux calibration for the astronomical object's FITS file using a pre-existing sensitivity or sens file.
    
    :param file [String]: Filename (without extension) of the FITS file.
    :param err_file [String]: Filename (without extension) of the error (1 sigma) FITS file.
    :param basename [String]: Basename of the FITS file. 
    
    :return: fc_file [String]: Filename (without extension) of the *flux corrected* FITS file.
    :return fc_err_file [String]: Filename (without extension) of the *flux corrected* error (1 sigma) FITS file.
    '''

    print('\n------------------------\nFlux Calibration\n-------------------------\n')

    # Load sensitivity FITS file.
    sens = raw_input('Please enter the complete filepath (with extension) of the sensitivity file:  ')

    # Set names of *flux corrected* FITS files.
    fc_file = 'fc_{}.fits'.format(basename)
    fc_err_file = 'fc_err_{}.fits'.format(basename)

    # Apply FLUX calibration
    iraf.noao.twod.longslit.calibrate(file, output=fc_file, sensitivity=sens, fnu='No')
    iraf.noao.twod.longslit.calibrate(err_file, output=fc_err_file, sensitivity=sens, fnu='No')

    # Extract header information
    hdu = fits.open(file)
    hdr = hdu[0].header
    RA = hdr['RA']
    DEC = hdr['DEC']
    UT = hdr['UTC-OBS']
    DATE = hdr['DATE-OBS']
    EPOCH = hdr['EPOCH']
    OBS = hdr['OBSERVAT']

    # Update header information of newly created files
    iraf.imutil.hedit(images=fc_file, fields='RA', value=RA)
    iraf.imutil.hedit(images=fc_file, fields='DEC', value=DEC)
    iraf.imutil.hedit(images=fc_file, fields='UT', value=UT)
    iraf.imutil.hedit(images=fc_file, fields='DATE-OBS', value=DATE)
    iraf.imutil.hedit(images=fc_file, fields='EPOCH', value=EPOCH)
    iraf.imutil.hedit(images=fc_file, fields='OBSERVAT', value=OBS)

    iraf.imutil.hedit(images=fc_err_file, fields='RA', value=RA)
    iraf.imutil.hedit(images=fc_err_file, fields='DEC', value=DEC)
    iraf.imutil.hedit(images=fc_err_file, fields='UT', value=UT)
    iraf.imutil.hedit(images=fc_err_file, fields='DATE-OBS', value=DATE)
    iraf.imutil.hedit(images=fc_err_file, fields='EPOCH', value=EPOCH)
    iraf.imutil.hedit(images=fc_err_file, fields='OBSERVAT', value=OBS)

    return fc_file, fc_err_file


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Heliocentric Velocity Correction
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''


def velcor(fc_file, fc_err_file, basename):
    '''
    
    Generates the bash code that can be used to apply the Heliocentric velocity correction to the *flux calibrated* FITS file.

    *Note*: The code is generated since the velocity correction command needs manual input through the command terminal.
    This just speeds up the process by generating some code you can initially copy and paste.
    
    :param: fc_file [String]: Filename (without extension) of the *flux corrected* FITS file.
    :param fc_err_file [String]: Filename (without extension) of the *flux corrected* error (1 sigma) FITS file.
    :param basename [String]: Basename of the FITS file.

    :return: None
    '''

    # Set names of *heliocentric velocity corrected* FITS files.
    hc_file = 'hc_{}.fits'.format(basename)
    hc_err_file = 'hc_err_{}.fits'.format(basename)
    
    print(
        '\n---------------------------------------\nHeliocentric Velocity Correction\n---------------------------------------\n')

    print(
        '\n\n>> Please complete the Heliocentric Velocity Correction & aperture attraction manually through PyRAF in bash:\n\n')

    # Print out the code to be used in PyRAF through bash.
    print('rvcorrect images={}'.format(fc_file))
    print('dopcor {} {} isvel+'.format(fc_file, hc_file))
    print('rvcorrect images={}'.format(fc_err_file))
    print('dopcor {} {} isvel+'.format(fc_err_file, hc_err_file))

    print(
        '\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\nData Reduction Completed for {}!\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n'.format(
            basename))

    print(
        '\n\n>> Remember to apply Galactic Extinction Correction!:\n\n')

    return None


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Apply Corrections
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''


def correct(basename):
    '''

    Automatically applies the flux calibration and generates the code to be used in PyRAF (through bash) to apply the heliocentric velocity correction.

    :param basename [String]: Basename of the FITS file.

    :return: None
    '''

    # Set file names.
    file = '{}.fits'.format(basename)
    err_file = 'err_{}.fits'.format(basename)

    # Apply flux calibration.
    fc_file, fc_err_file = flux(file, err_file, basename)
    # Generate code for heliocentric velocity correction.
    velcor(fc_file, fc_err_file, basename)

    return None


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Run Program
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

if __name__ == '__main__':

    # Set file path of data.
    path = raw_input('Please enter the path of the working directory:  ')
    # Change to path.
    os.chdir(path)
    
    # Specify which file you're working with. This program assumes the main file and 1 sigma error file have the same basename.
    # I.e:
    # basename.fits and err_basename.fits
    # where basename is the name of the astronomical object.
    basename = raw_input('Please enter filename (without extension) of the galaxy\'s FITS file:  ')
    
    # Run the main program.
    correct(basename)
