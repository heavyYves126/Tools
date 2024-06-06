#ifndef STFT_H_INCLUDED
#define STFT_H_INCLUDED

struct STFT {
  unsigned long k;
  unsigned long m;
  double **dftMagnitude;
  double **dftPhase;
};

struct STFT calculate_stft(
  double *sn,
  unsigned long *N,
  unsigned long *L,
  unsigned long *overlap,
  unsigned short type
);

#endif
