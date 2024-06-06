#ifndef IFFT_H_INCLUDED
#define IFFT_H_INCLUDED

struct IDFT {
  double *idft;
  unsigned long size;
};

struct IDFT calculate_ifft(double *xkReal, double *xkImag, unsigned long *M);

#endif
