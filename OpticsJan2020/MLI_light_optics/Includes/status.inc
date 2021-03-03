c status and error flags
c ieig46 is a flag set by eig4 and eig6, and reset by
c eig4, eig6, fit, and opt.
c A value of 0 means everything is ok, and a value of 1 means
c that eigenvalues in eig4 or eig6 were found to be off the unit
c circle.
      common /status/ imbad
