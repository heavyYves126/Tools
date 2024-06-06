#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stft.h"
#include "fft.h"

#define PI 3.14159265358979323846

unsigned long calculate_total_windows(
  unsigned long *L,
  unsigned long *overlap,
  unsigned long *N
);
double ** create_matrix(
  unsigned long *rows,
  unsigned long *cols
);
void perform_rectangular_windowing(
  double *sn,
  struct STFT *stft,
  unsigned long *overlap
);
void perform_hamming_windowing(
  double *sn,
  struct STFT *stft,
  unsigned long *overlap
);
void calculate_fft_per_window(
  struct STFT *stft
);
void get_current_window(
  double *currentWindow,
  struct STFT *stft,
  unsigned long *colum,
  unsigned short type
);

/**
 * Calculate the Short Time Fourier Transform of a fiven signal.
 *
 * @param *sn      The signal.
 * @param *N       The size of the signal.
 * @param *L       The size of the signal
 * @param *overlap Number of samples to overlap between windows.
 * @param  type    Type of window to be used. 0 -> Hamming; 1 -> Rectangular
 * @return         The Short Time Fourier Transform
 */
struct STFT calculate_stft(
  double *sn,
  unsigned long *N,
  unsigned long *L,
  unsigned long *overlap,
  unsigned short type
)
{
  struct STFT stft;
  unsigned long nowind = 0;

  nowind = calculate_total_windows(L, overlap, N);
  printf("# Windows %lu\n", nowind);

  stft.k = *L;
  stft.m = nowind;
  stft.dftMagnitude = create_matrix(L, &nowind);
  stft.dftPhase = create_matrix(L, &nowind);

  switch (type) {
    case 0:
      perform_hamming_windowing(sn, &stft, overlap);
      break;
    case 1:
      perform_rectangular_windowing(sn, &stft, overlap);
      break;
    default :
      perform_hamming_windowing(sn, &stft, overlap);
  }
  calculate_fft_per_window(&stft);

  return stft;
}

/**
 * Perform the windowing using a rectangular window.
 * @param *sn      The signal.
 * @param *stft    The STFT struct.
 * @param *L       length of the window.
 * @param *overlap Overlap between windows.
 */
void perform_rectangular_windowing(
  double *sn,
  struct STFT *stft,
  unsigned long *overlap
)
{
  unsigned long i = 0, j = 0, start = 0;

  for (i = 0; i < stft->m; i++) {
    for (j = 0; j < stft->k; j++) {
      stft->dftMagnitude[j][i] = sn[j + start];
    }
    start += (stft->k - *overlap);
  }
}

/**
 * Perform the windowing using a Hamming window.
 * @param *sn      The signal.
 * @param *stft    The STFT struct.
 * @param *L       length of the window.
 * @param *overlap Overlap between windows.
 */
void perform_hamming_windowing(
  double *sn,
  struct STFT *stft,
  unsigned long *overlap
)
{
  unsigned long i = 0, j = 0, start = 0;
  double w = 0.0;

  for (i = 0; i < stft->m; i++) {
    for (j = 0; j < stft->k; j++) {
      w = 0.54 - 0.46 * cos(2 * PI * (j + 1) / (stft->k - 1));
      stft->dftMagnitude[j][i] = sn[j + start] * w;
    }
    start += (stft->k - *overlap);
  }
}

/**
 * Calculate the total of windows according to the length (l) and the ovelap.
 * @param *L       Window length.
 * @param *overlap Number of samples to overlap between windows.
 * @param *N       The size of the signal.
 * @return         Total of windows.
 */
unsigned long calculate_total_windows(
  unsigned long *L,
  unsigned long *overlap,
  unsigned long *N
)
{
  unsigned long start = 0, total = 1;

  if (*L == *overlap) {
    return total;
  }

  if (*overlap > *L) {
    printf("Error: The overlap is greater than the length of the window\n");
    exit(0);
  }

  do {
    start += (*L - *overlap);
    ++total;
  } while ((start + *L) <= *N);

  total -= 1;

  return total;
}

/**
 * Create a matrix of specified rows and columns.
 * @param *rows Matrix rows.
 * @param *cols Matrix columns.
 * @return      Matrix.
 */
double ** create_matrix(unsigned long *rows, unsigned long *cols) {
  double ** matrix;
  unsigned long i;
  matrix = (double **) calloc(*rows, sizeof(double));

  if (NULL == matrix) {
    printf("Unable to allocate the matrix STFT\n");
    exit(0);
  }

  for (i = 0; i < *rows; i++) {
    matrix[i] = (double*) calloc(*cols, sizeof(double));

    if (NULL == matrix[i]) {
      printf("Unable to allocate the column %lu STFT\n", i);
      exit(0);
    }
  }

  return matrix;
}

/**
 * Calculate the Fast Fourier Transform of each window.
 *
 * @param *stft The STFT struct.
 */
void calculate_fft_per_window(struct STFT *stft)
{
  struct DFT xwfft;
  double *currentWindow;
  unsigned long i, j;

  currentWindow = calloc(stft->k, sizeof(double));

  if (NULL == currentWindow) {
    printf("Unable to allocate the matrix STFT\n");
    exit(0);
  }

  for (i = 0; i < stft->m; i++) {
    get_current_window(currentWindow, stft, &i, 0);
    xwfft = calculate_fft(currentWindow, &stft->k);
    for (j = 0; j < stft->k; j++) {
      stft->dftMagnitude[j][i] = xwfft.magnitude[j];
      stft->dftPhase[j][i] = xwfft.phase[j];
    }

    free(xwfft.magnitude);
    free(xwfft.phase);
  }

  free(currentWindow);
}

/**
 * Get the current window (column) of the matrix.
 *
 * @param *currentWindow  Current windows array.
 * @param *stft           The STFT struct.
 * @param *colum          Number of the column
 * @param *type           Type of data to get.
 */
void get_current_window(
  double *currentWindow,
  struct STFT *stft,
  unsigned long *colum,
  unsigned short type
)
{
  unsigned long i;

  for (i = 0; i < stft->k; i++) {
    switch (type) {
      case 0:
        currentWindow[i] = stft->dftMagnitude[i][*colum];
        break;
      case 1:
        currentWindow[i] = stft->dftPhase[i][*colum];
        break;
    }
  }
}
