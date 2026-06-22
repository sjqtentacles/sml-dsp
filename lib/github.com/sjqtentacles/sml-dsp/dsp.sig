(* dsp.sig

   A frequency-domain-aware digital-signal-processing toolkit in pure Standard
   ML, built on the vendored `Complex` and `Fft` libraries.

   Everything is pure and deterministic: no FFI, threads, clock or randomness.
   Reals are compared with an epsilon by callers; the only transcendentals are
   `Math.sin`/`Math.cos`/`Math.sqrt` (both compilers defer to the same `libm`),
   so identical inputs yield identical outputs under MLton and Poly/ML.

   Conventions: signals are `real array`; frequencies and sample rates are in
   Hz; window/FIR lengths are sample counts. Filter transfer functions follow
   the RBJ audio-EQ cookbook with the standard biquad form

       H(z) = (b0 + b1 z^-1 + b2 z^-2) / (a0 + a1 z^-1 + a2 z^-2).

   `convolve`/`applyFir` return the full linear convolution (length
   `len a + len b - 1`); `applyIir`/`applyBiquad` return a signal the same
   length as the input (zero initial state). *)

signature DSP =
sig
  (* ---- Window functions (real array of length n) ---- *)
  structure Window :
  sig
    datatype kind = Rectangular | Hann | Hamming | Blackman | Bartlett

    (* `make k n` builds the length-n symmetric window of kind k. n <= 1
       returns an all-ones window (the degenerate case). Standard (N-1)-
       denominator periodic-symmetric formulas. *)
    val make : kind -> int -> real array
  end

  (* ---- Biquad filter design (RBJ cookbook) + evaluation ---- *)
  structure Biquad :
  sig
    (* Unnormalized transfer-function coefficients (a0 is not forced to 1). *)
    type coeffs = { b0 : real, b1 : real, b2 : real,
                    a0 : real, a1 : real, a2 : real }

    (* Design parameters. `cutoff` is the corner/center frequency in Hz, `q`
       the quality factor (0.707 ~ Butterworth), `gainDb` the peak/shelf gain
       in decibels (used only by peakingEq/lowShelf/highShelf), `rate` the
       sample rate in Hz. *)
    type spec = { cutoff : real, q : real, gainDb : real, rate : int }

    val lowpass   : spec -> coeffs
    val highpass  : spec -> coeffs
    val bandpass  : spec -> coeffs   (* constant 0 dB peak gain *)
    val notch     : spec -> coeffs
    val peakingEq : spec -> coeffs
    val lowShelf  : spec -> coeffs
    val highShelf : spec -> coeffs

    (* Direct-Form-II-transposed application over a signal, zero initial
       state; output length = input length. *)
    val apply : coeffs -> real array -> real array

    (* Evaluate H(e^jw) at `freq` Hz: |H| and arg H (radians) where
       w = 2*pi*freq/rate. *)
    val freqResponse : coeffs -> { freq : real, rate : int }
                              -> { magnitude : real, phase : real }
  end

  (* ---- FIR design (windowed sinc) ---- *)
  structure Fir :
  sig
    (* Windowed-sinc designs returning `taps` coefficients. lowpass is
       normalized to unity DC gain. *)
    val lowpass  : { cutoff : real, taps : int, window : Window.kind, rate : int }
                   -> real array
    val highpass : { cutoff : real, taps : int, window : Window.kind, rate : int }
                   -> real array
    val bandpass : { lo : real, hi : real, taps : int,
                     window : Window.kind, rate : int } -> real array

    (* Apply FIR taps by direct convolution (full output). *)
    val apply : real array -> real array -> real array
  end

  (* ---- Short-time Fourier transform ---- *)
  structure Stft :
  sig
    type params = { frame : int, hop : int, window : Window.kind }

    (* Per-frame windowed FFT. Frames start at 0, hop, 2*hop, ... covering the
       whole signal (last frames zero-padded). Each result is a length-`frame`
       Complex spectrum. *)
    val analyze : params -> real array -> Complex.t array list

    (* Weighted overlap-add reconstruction (windowed synthesis normalized by
       the summed squared window). Output length =
       (#frames - 1) * hop + frame. Inverts `analyze` up to rounding on the
       interior. *)
    val synthesize : params -> Complex.t array list -> real array

    (* |X[k]| for one frame's spectrum. *)
    val magnitude : Complex.t array -> real array

    (* Magnitude spectrum of every frame (= map magnitude o analyze). *)
    val spectrogram : params -> real array -> real array list
  end

  (* ---- Top-level convenience entry points ---- *)

  (* Direct-Form-II-transposed biquad (alias of Biquad.apply). *)
  val applyBiquad  : Biquad.coeffs -> real array -> real array

  (* |H| and phase of a biquad at a frequency (alias of Biquad.freqResponse). *)
  val freqResponse : Biquad.coeffs -> { freq : real, rate : int }
                                   -> { magnitude : real, phase : real }

  (* Apply FIR taps by direct convolution (alias of Fir.apply). *)
  val applyFir : real array -> real array -> real array

  (* Generic direct-form IIR: y is the same length as x, zero initial state,
     coefficients given as lists b = [b0,b1,...], a = [a0,a1,...] (a0 <> 0). *)
  val applyIir : { b : real list, a : real list } -> real array -> real array

  (* Direct (sum-of-products) linear convolution, length len a + len b - 1. *)
  val convolve : real array * real array -> real array

  (* FFT-based linear convolution via the vendored `Fft.convolve`; agrees with
     `convolve` up to rounding. *)
  val fftConvolve : real array * real array -> real array
end
