#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED

struct DFT {
  double *magnitude;
  double *phase;
  unsigned long size;
};

struct DFT calculate_fft(double *xn, unsigned long *M);

#endif
