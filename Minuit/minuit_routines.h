//-*-mode:c; mode:font-lock;-*-

#ifndef MYMINUIT_H
#define MYMINUIT_H

extern "C" {
  //! Initialize Minuit with I/O unit numbers for in, out, save
  void mninit_(const int*, const int*, const int*);
  //! Define a parameter, assigning values and bounds
  void mnparm_(int *  num, const char * chnam, double * stval,
                double * step,  double * bnd1 ,
                double * bnd2, int * ierror, int stringlen);
  //! Prototype of the function to be minimized.
  typedef void (mfcn)(int * npar, double * grad, double * fval,
                      double * xval, int * iflag, void * futil);
  //! Execute a Minuit command specified as a character string
  void mncomd_(mfcn * fcn, const char * chstr, int * ierr, void * futil,
               int stringlen);
  //! Execute a Minuit command
  void mnexcm_(mfcn * fcn, char * chcom, double * arglis, int * narg,
               int * ierflg, void * futil, int strln);
  //! Set I/O unit numbers
  void mintio_(const int * iread, const int * iwrite, const int * isave);
  //! Get current value of a parameter
  void mnpout_(int * num, char * chnam, double * val, double * error,
               double * bnd1, double * bnd2, int * ivarbl, int strln);
  //! Get current status of minimization
  void mnstat_(double * fmin, double * fedm, double * errdef, int * npari,
               int * nparx, int * istat);
  //! Specify a title for a problem
  void mnseti_(char * ctitle, int strln);
  //! Define a parameter, assigning values and bounds from variables
  void mnpars_(char * chstr, int * icondn, int strln);
  //! Get current value of covariance matrix
  void mnemat_(double * emat, int * ndim);
  //! Access current parameter errors
  void mnerrs_(int * num, double * eplus, double * eminus, double * eparab,
               double * globcc);
  //! Find a function contour with the MNContour method
  void mncont_(mfcn * fcn, int * num1, int * num2, int * npt, double * xpt,
               double * ypt, int * nfound, void * futil);
  //! Utility function used by Minuit: interactive or batch mode
  long intrac_(double *);
}

#endif
