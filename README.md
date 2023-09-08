# SALT RSS Data Reduction Code

Collection of Python scripts to calibrate and process data from the Southern African Large Telescope Robert Stobie Spectrograph (SALT RSS).

<img src="https://mcdonaldobservatory.org/sites/default/files/images/news/gallery/salt.startrails.jpg" width="50%" height="50%">
Photo credit: https://mcdonaldobservatory.org/sites/default/files/images/news/gallery/salt.startrails.jpg

SALT RSS Documentation:
========================================

The SALT RSS data reduction procedure is well documented.

See the data reduction procedure here:
http://mips.as.arizona.edu/~khainline/salt_redux.html

Also see the data reduction FAQ:
https://astronomers.salt.ac.za/data/data-reduction-faq/

Python Scripts Included:
========================================

- Covolution.py: To convolve 2 images from the SALT RSS with each other.
- Calibration.py: To calibrate data from the SALT RSS.
- Sensitivity.py: For the generation of SALT RSS sensitivity files.
- Extract_SNR.py: To extract spectra up to a specific signal-to-noise ratio.

Prerequisites for using this program:
========================================

This program requires a programming environment with Python 3.6 installed.
With the following external libraries:
- PyRAF, 
- Astropy
- pathlib
- see: https://faculty1.coloradocollege.edu/~sburns/courses/18-19/pc362/Anaconda_IRAF_install.html for installation instructions.

