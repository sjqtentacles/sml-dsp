# sml-dsp

[![CI](https://github.com/sjqtentacles/sml-dsp/actions/workflows/ci.yml/badge.svg)](https://github.com/sjqtentacles/sml-dsp/actions/workflows/ci.yml)

A frequency-domain-aware **digital-signal-processing** toolkit in pure Standard
ML: window functions, RBJ-cookbook biquad design and evaluation, windowed-sinc
FIR design, generic IIR, direct & FFT-based convolution, and a windowed
STFT/ISTFT pair — all built on the vendored
[`sml-complex`](https://github.com/sjqtentacles/sml-complex) and
[`sml-fft`](https://github.com/sjqtentacles/sml-fft) libraries.

No FFI, no threads, no clock, no randomness: the same inputs always produce the
same outputs under **MLton** and **Poly/ML**. Reals are compared with an
epsilon and printed through a forced-decimal `fmtReal`; the only transcendentals
are `Math.sin`/`cos`/`sqrt`, which both compilers defer to the same `libm`, so
results agree to the last ulp.

This goes well beyond the tiny oscillator/biquad toolkit bundled with
`sml-wav`: it is a genuine, spectrum-aware DSP library (FIR design, FFT
convolution, short-time Fourier analysis/synthesis).

- **`Dsp.Window`** — `Rectangular`, `Hann`, `Hamming`, `Blackman`, `Bartlett`.
- **`Dsp.Biquad`** — `lowpass`/`highpass`/`bandpass`/`notch`/`peakingEq`/
  `lowShelf`/`highShelf` design returning transfer-function coefficients, a
  Direct-Form-II-transposed `apply`, and `freqResponse` (|H| and phase).
- **`Dsp.Fir`** — windowed-sinc `lowpass`/`highpass`/`bandpass` + `apply`.
- **`Dsp.Stft`** — `analyze`/`synthesize` (weighted overlap-add),
  `magnitude`, `spectrogram`.
- Top-level `applyBiquad`, `freqResponse`, `applyFir`, `applyIir`, `convolve`,
  `fftConvolve`.

## API

```sml
structure Dsp : sig
  structure Window : sig
    datatype kind = Rectangular | Hann | Hamming | Blackman | Bartlett
    val make : kind -> int -> real array
  end
  structure Biquad : sig
    type coeffs = { b0:real, b1:real, b2:real, a0:real, a1:real, a2:real }
    type spec   = { cutoff:real, q:real, gainDb:real, rate:int }
    val lowpass   : spec -> coeffs
    val highpass  : spec -> coeffs
    val bandpass  : spec -> coeffs
    val notch     : spec -> coeffs
    val peakingEq : spec -> coeffs
    val lowShelf  : spec -> coeffs
    val highShelf : spec -> coeffs
    val apply        : coeffs -> real array -> real array
    val freqResponse : coeffs -> { freq:real, rate:int }
                              -> { magnitude:real, phase:real }
  end
  structure Fir : sig
    val lowpass  : { cutoff:real, taps:int, window:Window.kind, rate:int } -> real array
    val highpass : { cutoff:real, taps:int, window:Window.kind, rate:int } -> real array
    val bandpass : { lo:real, hi:real, taps:int, window:Window.kind, rate:int } -> real array
    val apply    : real array -> real array -> real array
  end
  structure Stft : sig
    type params = { frame:int, hop:int, window:Window.kind }
    val analyze     : params -> real array -> Complex.t array list
    val synthesize  : params -> Complex.t array list -> real array
    val magnitude   : Complex.t array -> real array
    val spectrogram : params -> real array -> real array list
  end
  val applyBiquad  : Biquad.coeffs -> real array -> real array
  val freqResponse : Biquad.coeffs -> { freq:real, rate:int }
                                   -> { magnitude:real, phase:real }
  val applyFir     : real array -> real array -> real array
  val applyIir     : { b:real list, a:real list } -> real array -> real array
  val convolve     : real array * real array -> real array
  val fftConvolve  : real array * real array -> real array
end
```

The transfer function is
`H(z) = (b0 + b1 z^-1 + b2 z^-2) / (a0 + a1 z^-1 + a2 z^-2)`. `convolve`,
`fftConvolve` and `applyFir` return the full linear convolution (length
`len a + len b - 1`); `applyIir`/`applyBiquad` return a signal the same length
as the input with zero initial state.

## Example

```sml
val rate = 48000
val lp = Dsp.Biquad.lowpass { cutoff = 1000.0, q = 0.707, gainDb = 0.0, rate = rate }
val { magnitude, ... } = Dsp.freqResponse lp { freq = 0.0, rate = rate }  (* ~ 1.0 *)
val taps = Dsp.Fir.lowpass { cutoff = 800.0, taps = 11, window = Dsp.Window.Hamming, rate = 8000 }
val y = Dsp.fftConvolve (Array.fromList [1.0, ~2.0, 3.0, 0.5],
                         Array.fromList [0.5, 1.0, ~0.5])
```

Running [`examples/demo.sml`](examples/demo.sml) with `make example` prints:

```
Biquad lowpass (cutoff 1000 Hz, Q 0.707, rate 48000):
  coeffs b = [0.004278, 0.008555, 0.004278]
  coeffs a = [1.092310, -1.982890, 0.907690]
  |H(f)| magnitude response:
    f =      0 Hz  |H| = 1.0000
    f =    500 Hz  |H| = 0.9702
    f =   1000 Hz  |H| = 0.7070
    f =   2000 Hz  |H| = 0.2406
    f =   8000 Hz  |H| = 0.0129
    f =  20000 Hz  |H| = 0.0003

Windowed-sinc FIR lowpass (cutoff 800 Hz, 11 taps, Hamming, rate 8000):
  taps = [
    h[0] = 0.000000
    h[1] = 0.009304
    h[2] = 0.047578
    h[3] = 0.122364
    h[4] = 0.202247
    h[5] = 0.237016
    h[6] = 0.202247
    h[7] = 0.122364
    h[8] = 0.047578
    h[9] = 0.009304
    h[10] = 0.000000
  ]
  sum of taps (DC gain) = 1.000000

Direct vs FFT convolution of [1,-2,3,0.5] * [0.5,1,-0.5]:
  direct      = [0.5000, 0.0000, -1.0000, 4.2500, -1.0000, -0.2500]
  fft-based   = [0.5000, 0.0000, -1.0000, 4.2500, -1.0000, -0.2500]
```

## Build & test

Requires [MLton](http://mlton.org/) and/or [Poly/ML](https://polyml.org/).

```sh
make test        # build + run the suite under MLton
make test-poly   # run the suite under Poly/ML
make all-tests   # both
make example     # build + run the demo
make clean
```

## Installing with smlpkg

```sh
smlpkg add github.com/sjqtentacles/sml-dsp
smlpkg sync
```

Reference `lib/github.com/sjqtentacles/sml-dsp/dsp.mlb` from your own `.mlb`
(MLton / MLKit), or feed `sources.mlb` to `tools/polybuild` (Poly/ML).

## Vendored dependencies

The libraries `sml-complex` and `sml-fft` are vendored under
`lib/github.com/sjqtentacles/` and listed (in dependency order) in
`sources.mlb` ahead of `dsp.sig`/`dsp.sml`. The copied `.sig`/`.sml` files are
byte-identical to their upstream canonical sources.

## Layout

```
sml.pkg                                       smlpkg manifest (requires complex, fft)
Makefile                                      MLton + Poly/ML targets
.github/workflows/ci.yml                      CI: MLton + Poly/ML 5.9.1
lib/github.com/sjqtentacles/
  sml-complex/  complex.sig complex.sml       (vendored)
  sml-fft/      fft.sig fft.sml               (vendored)
  sml-dsp/
    dsp.sig    DSP signature
    dsp.sml    Dsp implementation
    sources.mlb  vendored deps, then own sources
    dsp.mlb      public basis
examples/
  demo.sml       biquad + FIR + convolution walkthrough
test/
  harness.sml    shared assertion harness
  test.sml       canonical window/biquad/FIR/STFT vectors (42 checks)
  entry.sml / main.sml
tools/polybuild  Poly/ML build wrapper
```

## Tests

42 deterministic checks: hand-computed window vectors (Hann/Hamming/Blackman/
Bartlett at length 5, symmetry); biquad |H| at DC and Nyquist proven exactly
(lowpass 1/0, highpass 0/1, bandpass 0/0, notch 1/1, peakingEq unity-DC and
+6 dB ≈ 2× at center); DF2T `apply` cross-checked against the generic
direct-form IIR; hand-computed `convolve`/`applyIir` vectors; `fftConvolve`
agreeing with direct `convolve`; FIR symmetry and unity DC gain; a
DCT round-trip through the vendored `Fft`; and an STFT→ISTFT round-trip
reconstructing a sine on the interior. `make all-tests` verifies identical
output under both compilers.

## License

MIT. See [LICENSE](LICENSE).
