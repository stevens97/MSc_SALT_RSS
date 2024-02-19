rajin_code
==========
This piece of code is for the spectroscopic reduction of long-slit data from SALT (Southern African Large Telescope).
It is by no means the only way to reduce the data. What it provides though is convenience. Instead of running a few 
IRAF task by hand and on single or list of images, this piece of code, does the reduction at one go.

====================================================================================================================

The Requirements:
This code has been tested to run on IRAF 2.14, Pyraf 1.10 and Python 2.6
There is a know issue with Python 2.7 that IRAF pop-up windows fails to load due to some issue with Tkinter. In some 
cases though it could run.
It would be good to load the following packages: noao, onedspec, twodspec, longslit, apextract, stsdas, imred, ccdred
in your IRAF default login.cl file found in your IRAf home folder.
Also there are some files (ARC line-list (has extension '.txt'), LA-Cosmic file (lacos_spec.cl)) that you need to paste in your IRAF home folder for the code
to work

====================================================================================================================

Input:
How to input data so that the code run smoothly?
The code can only calibrate 1 ARC spectra when running.
Therefore it would be good idea to select all images (science, flats, etc) to be calibrated by that ARC and paste 
all of them (along with the ARC) in a folder.
The code NEEDS to have an ARC file inside the folder where it is running to run. In the eventuality that you want to 
use an already calibrated ARC file for your data, paste both the ARC fits file along with the iraf database folder
containing the caliubrated information. 
Also put the code inside the directory where the images are, to run it.
The code has the ability to create a master flat and use it to to flat-fielding of the data but for safety, if there 
are flats, the code will create 2 sets of data: 1 where flat-fielding is done and the other where no flat fielding is
applied

====================================================================================================================

Output:
The output of the code is a set of data which has been wavelength calibrated, background subtracted, cosmic ray removal
ccdgap interpolation and  tilt/aperture extraction
As the code output new fits file during the reduction, some of them will be moved automatically to the folder history
The fits file will be renamed to Object names, arc but be aware that to keep track of the data, the last 3 digits at 
the end of the name, are the frame numbers recorded from salt
Eg
If the file of the object NGC1010 has a nice name like mbxgpP201305120035.fits at the beginning, after the renaming 
process, the name will be NGC1010_035.fits
There on prefixes will be added with respect to the task that has been applied on the data. Below is a table for the 
prefix used and there meaning

prefix    -    task

'c'      -       ccdgap filling with an interpolation value,

'la'     -       cosmic ray removal through LA cosmic,

'fl'     -       Flat fielding applied,

't'      -       transformed frame -> wavelength calibration applied,

'b'      -       background subtraction applied,

'til'    -       tilt correction applied through Apall function, 

'err'    -       error file,

'terr'   -       transformed error file,

'dterr'  -       divided error file (tranformed error file is divided by corresponding background subtracted science file),

'tildterr'  -    tilt corrected error file (the frame has gone through the 3 steps above),

'ARC'     -      the Arc file,

'FLAT'    -      Flat field images,

'com_flat'  -    combined image of the flat fields,

'ilcom_flat'  -  illumination corrected combined flat field,

'master_flat' -  The master flat image which is a normalised image of the previous step,

====================================================================================================================
contact
In case you use the code and encounter any error, please contact me at: rajin250@yahoo.com. Please do copy paste the 
error message inside the email
