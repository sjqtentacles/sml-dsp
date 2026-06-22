(* fft.sig

   The discrete Fourier transform over `Complex.t` arrays, plus a couple of
   real-input conveniences. All functions are pure and allocate fresh output
   arrays; the input is never mutated.

   Conventions (unnormalized forward transform):

     fft  x : X[k] = sum_{j=0}^{n-1} x[j] * exp(-2*pi*i*j*k/n)
     ifft X : x[j] = (1/n) sum_{k=0}^{n-1} X[k] * exp(+2*pi*i*j*k/n)

   so `ifft (fft x) = x` up to floating-point rounding. Power-of-two lengths
   use an iterative radix-2 Cooley-Tukey transform; other lengths use
   Bluestein's chirp-z algorithm, so every length is handled in O(n log n).
   An empty array transforms to an empty array. *)

signature FFT =
sig
  (* Forward DFT. *)
  val fft : Complex.t array -> Complex.t array

  (* Inverse DFT (the 1/n-normalized conjugate transform). *)
  val ifft : Complex.t array -> Complex.t array

  (* Forward DFT of a purely real signal (each sample embedded as x + 0i). *)
  val rfft : real array -> Complex.t array

  (* Inverse real FFT: the real part of `ifft`, discarding the (rounding-sized)
     imaginary component. For a spectrum produced from a real signal this
     recovers that signal, i.e. `irfft (rfft x) = x` up to rounding. *)
  val irfft : Complex.t array -> real array

  (* Separable 2D forward / inverse DFT of a row-major matrix (an array of
     equal-length rows): transform every row, then every column. `ifft2` is
     1/(rows*cols)-normalized, so `ifft2 (fft2 m) = m` up to rounding.
     Rectangular matrices are supported; an empty matrix (no rows) maps to an
     empty matrix, and a matrix of empty rows maps to the same shape. *)
  val fft2  : Complex.t array array -> Complex.t array array
  val ifft2 : Complex.t array array -> Complex.t array array

  (* Discrete cosine transform, type II (unnormalized):

       dct  x : X[k] = sum_{j=0}^{n-1} x[j] * cos(pi*(2j+1)*k/(2n))

     and its exact inverse, the type-III transform scaled by 1/n:

       idct X : x[j] = (1/n) * (X[0] + 2 * sum_{k=1}^{n-1}
                                          X[k] * cos(pi*(2j+1)*k/(2n)))

     so `idct (dct x) = x` up to rounding. Both map an empty array to an empty
     array. *)
  val dct  : real array -> real array
  val idct : real array -> real array

  (* Linear convolution of two real sequences via the FFT. The result has
     length `length a + length b - 1` (empty if either input is empty), and
     equals the direct sum-of-products convolution up to rounding. *)
  val convolve : real array * real array -> real array
end
