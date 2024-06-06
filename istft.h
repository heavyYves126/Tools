#ifndef ISTFT_H_INCLUDED
#define ISTFT_H_INCLUDED

struct ISTFT {
  double *istft;
  unsigned long size;
};

struct ISTFT calculate_istft(
	struct STFT *stft,
	unsigned long *N,
	unsigned long *overlap,
	unsigned short type
);

#endif
