(* Tests for sml-dsp. Reference vectors are computed by hand (window formulas,
   biquad transfer-function values at DC/Nyquist) or cross-checked against
   closed forms and the vendored Fft. Reals are compared with an epsilon. *)

structure Tests =
struct
  open Harness
  structure W = Dsp.Window
  structure BQ = Dsp.Biquad

  val eps = 1e~9
  fun approx e (a, b) = Real.abs (a - b) <= e
  fun checkAppx name e (expected, actual) =
    check name (approx e (expected, actual))

  fun arrApprox e (xs, a) =
    length xs = Array.length a
    andalso List.all (fn i => approx e (List.nth (xs, i), Array.sub (a, i)))
                     (List.tabulate (length xs, fn i => i))
  fun checkArr name e (xs, a) = check name (arrApprox e (xs, a))

  fun arrEqArr e (a, b) =
    Array.length a = Array.length b
    andalso List.all (fn i => approx e (Array.sub (a, i), Array.sub (b, i)))
                     (List.tabulate (Array.length a, fn i => i))

  fun fmtReal n r =
    let val s = Real.fmt (StringCvt.FIX (SOME n)) r
    in if String.isPrefix "~" s then "-" ^ String.extract (s, 1, NONE) else s end

  fun runAll () =
    let
      (* ----------------------- windows ----------------------- *)
      val () = section "Window: canonical length-5 vectors"
      val () = checkArr "rectangular 4" eps ([1.0,1.0,1.0,1.0], W.make W.Rectangular 4)
      val () = checkArr "hann 5" eps ([0.0,0.5,1.0,0.5,0.0], W.make W.Hann 5)
      val () = checkArr "hamming 5" eps ([0.08,0.54,1.0,0.54,0.08], W.make W.Hamming 5)
      val () = checkArr "blackman 5" eps ([0.0,0.34,1.0,0.34,0.0], W.make W.Blackman 5)
      val () = checkArr "bartlett 5" eps ([0.0,0.5,1.0,0.5,0.0], W.make W.Bartlett 5)
      val () = checkInt "hann length 32" (32, Array.length (W.make W.Hann 32))
      val () = checkArr "degenerate n=1" eps ([1.0], W.make W.Hann 1)
      val () = check "hann symmetric (64)"
                 (let val w = W.make W.Hann 64
                  in List.all (fn i => approx eps (Array.sub (w, i), Array.sub (w, 63 - i)))
                              (List.tabulate (64, fn i => i)) end)

      (* ----------------------- biquad ------------------------ *)
      val rate = 48000
      val lp = BQ.lowpass   { cutoff = 1000.0, q = 0.707, gainDb = 0.0, rate = rate }
      val hp = BQ.highpass  { cutoff = 1000.0, q = 0.707, gainDb = 0.0, rate = rate }
      val bp = BQ.bandpass  { cutoff = 1000.0, q = 1.0,   gainDb = 0.0, rate = rate }
      val nt = BQ.notch     { cutoff = 1000.0, q = 1.0,   gainDb = 0.0, rate = rate }
      val pk = BQ.peakingEq { cutoff = 1000.0, q = 1.0,   gainDb = 6.0, rate = rate }
      fun mag c f = #magnitude (BQ.freqResponse c { freq = f, rate = rate })
      val nyq = real rate / 2.0

      val () = section "Biquad: lowpass frequency response"
      val () = checkAppx "lowpass |H| at DC ~ 1" 1e~9 (1.0, mag lp 0.0)
      val () = checkAppx "lowpass |H| at Nyquist ~ 0" 1e~9 (0.0, mag lp nyq)
      val () = check "lowpass attenuates 20 kHz" (mag lp 20000.0 < 0.05)
      val () = check "lowpass passes DC vs stopband" (mag lp 0.0 > mag lp 20000.0)

      val () = section "Biquad: highpass / bandpass / notch / peaking"
      val () = checkAppx "highpass |H| at DC ~ 0" 1e~9 (0.0, mag hp 0.0)
      val () = checkAppx "highpass |H| at Nyquist ~ 1" 1e~9 (1.0, mag hp nyq)
      val () = checkAppx "bandpass |H| at DC ~ 0" 1e~9 (0.0, mag bp 0.0)
      val () = checkAppx "bandpass |H| at Nyquist ~ 0" 1e~9 (0.0, mag bp nyq)
      val () = checkAppx "notch |H| at DC ~ 1" 1e~9 (1.0, mag nt 0.0)
      val () = checkAppx "notch |H| at Nyquist ~ 1" 1e~9 (1.0, mag nt nyq)
      val () = check "notch deep at center" (mag nt 1000.0 < 0.05)
      val () = checkAppx "peakingEq |H| at DC ~ 1" 1e~9 (1.0, mag pk 0.0)
      val () = check "peakingEq boosts center (+6 dB ~ 2x)"
                 (approx 0.02 (2.0, mag pk 1000.0))

      val () = section "Biquad: apply == generic IIR"
      val xsig = Array.fromList [1.0, 2.0, ~1.0, 0.5, 3.0, ~2.0, 0.0, 1.0, 0.25, ~0.75]
      val viaBiquad = Dsp.applyBiquad lp xsig
      val viaIir = Dsp.applyIir
                     { b = [#b0 lp, #b1 lp, #b2 lp], a = [#a0 lp, #a1 lp, #a2 lp] } xsig
      val () = check "DF2T biquad matches direct-form IIR" (arrEqArr 1e~9 (viaBiquad, viaIir))
      val () = check "impulse response is finite/defined"
                 (let val imp = Array.fromList [1.0,0.0,0.0,0.0,0.0]
                      val r = Dsp.applyBiquad lp imp
                  in approx 1e~9 (#b0 lp / #a0 lp, Array.sub (r, 0)) end)

      (* --------------------- convolution --------------------- *)
      val () = section "convolve (direct)"
      val () = checkArr "[1,2,3]*[1,1]" eps ([1.0,3.0,5.0,3.0],
                 Dsp.convolve (Array.fromList [1.0,2.0,3.0], Array.fromList [1.0,1.0]))
      val () = checkArr "[1,2]*[3,4]" eps ([3.0,10.0,8.0],
                 Dsp.convolve (Array.fromList [1.0,2.0], Array.fromList [3.0,4.0]))
      val () = check "empty convolve" (Array.length (Dsp.convolve (Array.fromList [], Array.fromList [1.0])) = 0)

      val () = section "fftConvolve agrees with direct convolve"
      val ca = Array.fromList [1.0, ~2.0, 3.0, 0.5, ~1.5, 2.0]
      val cb = Array.fromList [0.5, 1.0, ~0.5, 2.0]
      val () = check "fftConvolve == convolve (small)"
                 (arrEqArr 1e~9 (Dsp.convolve (ca, cb), Dsp.fftConvolve (ca, cb)))
      val () = checkArr "fftConvolve [1,2,3]*[1,1]" 1e~9 ([1.0,3.0,5.0,3.0],
                 Dsp.fftConvolve (Array.fromList [1.0,2.0,3.0], Array.fromList [1.0,1.0]))

      (* ------------------------ IIR -------------------------- *)
      val () = section "applyIir (generic direct form)"
      val () = checkArr "moving average b=[.5,.5]" eps ([0.5,1.5,2.5,3.5],
                 Dsp.applyIir { b = [0.5,0.5], a = [1.0] }
                              (Array.fromList [1.0,2.0,3.0,4.0]))
      val () = checkArr "one-pole y=x+0.5y[-1]" 1e~12 ([1.0,0.5,0.25,0.125],
                 Dsp.applyIir { b = [1.0], a = [1.0, ~0.5] }
                              (Array.fromList [1.0,0.0,0.0,0.0]))

      (* ------------------------ FIR -------------------------- *)
      val () = section "Fir: windowed-sinc lowpass"
      val flp = Dsp.Fir.lowpass { cutoff = 1000.0, taps = 21, window = W.Hamming, rate = 8000 }
      val () = checkInt "lowpass tap count" (21, Array.length flp)
      val () = check "lowpass symmetric"
                 (List.all (fn i => approx 1e~12 (Array.sub (flp, i), Array.sub (flp, 20 - i)))
                           (List.tabulate (21, fn i => i)))
      val () = check "lowpass unity DC gain"
                 (approx 1e~9 (1.0, Array.foldl op+ 0.0 flp))
      val fhp = Dsp.Fir.highpass { cutoff = 1000.0, taps = 21, window = W.Hamming, rate = 8000 }
      val () = checkInt "highpass tap count" (21, Array.length fhp)
      val () = check "highpass symmetric"
                 (List.all (fn i => approx 1e~12 (Array.sub (fhp, i), Array.sub (fhp, 20 - i)))
                           (List.tabulate (21, fn i => i)))
      val () = check "applyFir == convolve(taps, sig)"
                 (let val sg = Array.fromList [1.0,2.0,3.0,4.0,5.0]
                  in arrEqArr 1e~12 (Dsp.applyFir flp sg, Dsp.convolve (flp, sg)) end)

      (* ------------- DCT round-trip via vendored Fft --------- *)
      val () = section "Vendored Fft: DCT round-trip"
      val dx = Array.fromList [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]
      val () = check "idct (dct x) = x"
                 (arrEqArr 1e~9 (dx, Fft.idct (Fft.dct dx)))

      (* ------------------ STFT / ISTFT ----------------------- *)
      val () = section "Stft: analyze / synthesize round-trip"
      val sigLen = 64
      val xsine = Array.tabulate (sigLen, fn i =>
                    Math.sin (2.0 * Math.pi * 5.0 * real i / real sigLen))
      val params = { frame = 16, hop = 4, window = W.Hann }
      val spectra = Dsp.Stft.analyze params xsine
      val recon = Dsp.Stft.synthesize params spectra
      val () = checkInt "frame count" (16, List.length spectra)
      val () = check "each spectrum length = frame"
                 (List.all (fn s => Array.length s = 16) spectra)
      val () = check "interior reconstruction matches sine (1e~6)"
                 (List.all (fn i => approx 1e~6 (Array.sub (xsine, i), Array.sub (recon, i)))
                           (List.tabulate (40, fn i => i + 12)))
      val () = checkInt "spectrogram frame count" (16, List.length (Dsp.Stft.spectrogram params xsine))
      val () = check "magnitude nonnegative"
                 (List.all (fn s => List.all (fn i => Array.sub (s, i) >= 0.0)
                                             (List.tabulate (Array.length s, fn i => i)))
                           (Dsp.Stft.spectrogram params xsine))
    in
      Harness.run ()
    end

  val run = runAll
end
