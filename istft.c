#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stft.h"
#include "istft.h"
#include "ifft.h"

#define PI 3.14159265358979323846

double ** create_matrix(
  unsigned long *rows,
  unsigned long *cols
);
void calculate_ifft_per_window(
  struct STFT *stft,
  double **istftW
);
void get_current_window(
  double *currentWindow,
  struct STFT *stft,
  unsigned long *colum,
  unsigned short type
);
void add_rectangular_windows(
  double *xnp,
  double **istftW,
  unsigned long *k,
  unsigned long *m
);
void add_hamming_windows(
  double *xnp,
  double **istftW,
  unsigned long *k,
  unsigned long *m
);

/**
 * Calculate the Inverse Short Time Fourier Transform of a given STFT.
 *
 * @param *stft      The STFT.
 * @param *L       The size of the signal
 * @param *overlap Number of samples to overlap between windows.
 * @param  type    Type of window to be used. 0 -> Hamming; 1 -> Rectangular
 * @return         The Short Time Fourier Transform
 */
struct ISTFT calculate_istft(
  struct STFT *stft,
  unsigned long *N,
  unsigned long *overlap,
  unsigned short type
)
{
  struct ISTFT istft;
  double **istftW;
  double *xnp;
  unsigned long i;

  xnp = calloc(*N, sizeof(double));

  if (NULL == xnp) {
    printf("Unable to allocate a variable for ISTFT\n");
    exit(0);
  }

  istftW = create_matrix(&stft->k, &stft->m);
  calculate_ifft_per_window(stft, istftW);

  switch (type) {
    case 0:
      add_hamming_windows(xnp, istftW, &stft->k, &stft->m);
      break;
    case 1:
      add_rectangular_windows(xnp, istftW, &stft->k, &stft->m);
      break;
  }

  for (i = 0; i < stft->k; i++){
    free(istftW[i]);
  }

  free(istftW);

  istft.istft = xnp;
  istft.size = *N;

  return istft;
}

/**
 * Add the rectangular windowing. Windows need 0% overlap
 * @param *sn      The output signal.
 * @param *stft    The STFT struct.
 * @param *k       Window length.
 * @param *m       Total windows.
 */
void add_rectangular_windows(
    double *xnp,
    double **istftW,
    unsigned long *k,
    unsigned long *m
)
{
  unsigned long i = 0, j = 0, n = 0;

  for (i = 0; i < *m; i++) {
    for (j = 0; j < *k; j++) {
      xnp[n] = istftW[j][i];
      n++;
    }
  }
}

/**
 * Add the Hamming windowing. Windows need 50% overlap
 * @param *sn      The output signal.
 * @param *stft    The STFT struct.
 * @param *k       Window length.
 * @param *m       Total windows.
 */
void add_hamming_windows(
    double *xnp,
    double **istftW,
    unsigned long *k,
    unsigned long *m
)
{
  unsigned long i = 0, j = 0, n = 0;
  unsigned short hop = 0;

  hop = (unsigned short) (*k / 2);

  for (j = 0; j < hop; j++) {
    xnp[n] = istftW[j][0];
    n++;
  }

  for (i = 0; i < *m - 1 ; i++) {
    for (j = hop ; j < *k; j++) {
      xnp[n] = istftW[j][i] + istftW[j - hop][i + 1];
      n++;
    }
  }

  for (j = hop; j < *k; j++) {
    xnp[n] = istftW[j][*m - 1];
    n++;
  }
}

/**
 * Calculate the Inverse Fast Fourier Transform of each window.
 *
 * @param *stft   The STFT struct.
 * @param *istftW A matrix to stored the IFFT.
 */
void calculate_ifft_per_window(struct STFT *stft, double **istftW)
{
  struct IDFT xifft;
  double *currentMagnitudeWindow, *currentPhaseWindow;;
  unsigned long i, j;

  currentMagnitudeWindow = calloc(stft->k, sizeof(double));
  currentPhaseWindow = calloc(stft->k, sizeof(double));

  if (NULL == currentMagnitudeWindow || NULL == currentPhaseWindow) {
    printf("Unable to allocate the memory ISTFT\n");
    exit(0);
  }

  for (i = 0; i < stft->m; i++) {
    get_current_window(currentMagnitudeWindow, stft, &i, 0);
    get_current_window(currentPhaseWindow, stft, &i, 1);

    xifft = calculate_ifft(currentMagnitudeWindow, currentPhaseWindow, &stft->k);
    for (j = 0; j < stft->k; j++) {
      istftW[j][i] = xifft.idft[j];
    }
    free(xifft.idft);
  }

  free(currentMagnitudeWindow);
  free(currentPhaseWindow);
}
