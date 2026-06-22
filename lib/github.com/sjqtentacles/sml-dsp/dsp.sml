(* dsp.sml

   Implementation of DSP, sealed behind the signature. Built on the vendored
   `Complex` and `Fft` structures. Pure and deterministic. *)

structure Dsp :> DSP =
struct
  structure C = Complex

  val pi = Math.pi

  (* ------------------------------------------------------------------ *)
  (* small array helpers                                                *)
  (* ------------------------------------------------------------------ *)

  fun tab (n, f) = Array.tabulate (n, f)
  fun len a = Array.length a
  fun sub (a, i) = Array.sub (a, i)

  (* ------------------------------------------------------------------ *)
  (* windows                                                            *)
  (* ------------------------------------------------------------------ *)

  structure Window =
  struct
    datatype kind = Rectangular | Hann | Hamming | Blackman | Bartlett

    fun make k n =
      if n <= 0 then Array.fromList []
      else if n = 1 then tab (1, fn _ => 1.0)
      else
        let
          val m = real (n - 1)
          fun w i =
            let val x = real i
            in
              case k of
                Rectangular => 1.0
              | Hann => 0.5 - 0.5 * Math.cos (2.0 * pi * x / m)
              | Hamming => 0.54 - 0.46 * Math.cos (2.0 * pi * x / m)
              | Blackman =>
                  0.42 - 0.5 * Math.cos (2.0 * pi * x / m)
                       + 0.08 * Math.cos (4.0 * pi * x / m)
              | Bartlett =>
                  1.0 - Real.abs ((x - m / 2.0) / (m / 2.0))
            end
        in
          tab (n, w)
        end
  end

  (* ------------------------------------------------------------------ *)
  (* generic direct-form IIR (Direct Form I)                            *)
  (* ------------------------------------------------------------------ *)

  fun applyIir { b, a } x =
    let
      val bA = Array.fromList b
      val aA = Array.fromList a
      val nb = len bA
      val na = len aA
      val () = if na = 0 then raise Domain else ()
      val a0 = sub (aA, 0)
      val n = len x
      val y = Array.array (n, 0.0)
      fun loop i =
        if i >= n then ()
        else
          let
            (* feed-forward: sum b[k] x[i-k] *)
            fun ff (k, acc) =
              if k >= nb then acc
              else if i - k < 0 then ff (k + 1, acc)
              else ff (k + 1, acc + sub (bA, k) * sub (x, i - k))
            (* feed-back: sum a[k] y[i-k], k >= 1 *)
            fun fb (k, acc) =
              if k >= na then acc
              else if i - k < 0 then fb (k + 1, acc)
              else fb (k + 1, acc + sub (aA, k) * sub (y, i - k))
            val v = (ff (0, 0.0) - fb (1, 0.0)) / a0
          in
            Array.update (y, i, v);
            loop (i + 1)
          end
    in
      loop 0; y
    end

  (* ------------------------------------------------------------------ *)
  (* biquad design + evaluation                                         *)
  (* ------------------------------------------------------------------ *)

  structure Biquad =
  struct
    type coeffs = { b0 : real, b1 : real, b2 : real,
                    a0 : real, a1 : real, a2 : real }
    type spec = { cutoff : real, q : real, gainDb : real, rate : int }

    (* common intermediate quantities *)
    fun omega ({ cutoff, rate, ... } : spec) = 2.0 * pi * cutoff / real rate
    fun alphaOf (sp as { q, ... } : spec) = Math.sin (omega sp) / (2.0 * q)
    fun ampOf ({ gainDb, ... } : spec) = Math.pow (10.0, gainDb / 40.0)

    fun lowpass sp =
      let val w0 = omega sp val c = Math.cos w0 val al = alphaOf sp
      in { b0 = (1.0 - c) / 2.0, b1 = 1.0 - c, b2 = (1.0 - c) / 2.0,
           a0 = 1.0 + al, a1 = ~2.0 * c, a2 = 1.0 - al } end

    fun highpass sp =
      let val w0 = omega sp val c = Math.cos w0 val al = alphaOf sp
      in { b0 = (1.0 + c) / 2.0, b1 = ~(1.0 + c), b2 = (1.0 + c) / 2.0,
           a0 = 1.0 + al, a1 = ~2.0 * c, a2 = 1.0 - al } end

    fun bandpass sp =
      let val w0 = omega sp val c = Math.cos w0 val al = alphaOf sp
      in { b0 = al, b1 = 0.0, b2 = ~al,
           a0 = 1.0 + al, a1 = ~2.0 * c, a2 = 1.0 - al } end

    fun notch sp =
      let val w0 = omega sp val c = Math.cos w0 val al = alphaOf sp
      in { b0 = 1.0, b1 = ~2.0 * c, b2 = 1.0,
           a0 = 1.0 + al, a1 = ~2.0 * c, a2 = 1.0 - al } end

    fun peakingEq sp =
      let val w0 = omega sp val c = Math.cos w0 val al = alphaOf sp
          val amp = ampOf sp
      in { b0 = 1.0 + al * amp, b1 = ~2.0 * c, b2 = 1.0 - al * amp,
           a0 = 1.0 + al / amp, a1 = ~2.0 * c, a2 = 1.0 - al / amp } end

    fun lowShelf sp =
      let val w0 = omega sp val c = Math.cos w0 val al = alphaOf sp
          val amp = ampOf sp val sq = 2.0 * Math.sqrt amp * al
      in { b0 = amp * ((amp + 1.0) - (amp - 1.0) * c + sq),
           b1 = 2.0 * amp * ((amp - 1.0) - (amp + 1.0) * c),
           b2 = amp * ((amp + 1.0) - (amp - 1.0) * c - sq),
           a0 = (amp + 1.0) + (amp - 1.0) * c + sq,
           a1 = ~2.0 * ((amp - 1.0) + (amp + 1.0) * c),
           a2 = (amp + 1.0) + (amp - 1.0) * c - sq } end

    fun highShelf sp =
      let val w0 = omega sp val c = Math.cos w0 val al = alphaOf sp
          val amp = ampOf sp val sq = 2.0 * Math.sqrt amp * al
      in { b0 = amp * ((amp + 1.0) + (amp - 1.0) * c + sq),
           b1 = ~2.0 * amp * ((amp - 1.0) + (amp + 1.0) * c),
           b2 = amp * ((amp + 1.0) + (amp - 1.0) * c - sq),
           a0 = (amp + 1.0) - (amp - 1.0) * c + sq,
           a1 = 2.0 * ((amp - 1.0) - (amp + 1.0) * c),
           a2 = (amp + 1.0) - (amp - 1.0) * c - sq } end

    (* Direct-Form-II-transposed, zero initial state. *)
    fun apply ({ b0, b1, b2, a0, a1, a2 } : coeffs) x =
      let
        val nb0 = b0 / a0 val nb1 = b1 / a0 val nb2 = b2 / a0
        val na1 = a1 / a0 val na2 = a2 / a0
        val n = len x
        val y = Array.array (n, 0.0)
        fun loop (i, z1, z2) =
          if i >= n then ()
          else
            let
              val xn = sub (x, i)
              val yn = nb0 * xn + z1
              val z1' = nb1 * xn - na1 * yn + z2
              val z2' = nb2 * xn - na2 * yn
            in
              Array.update (y, i, yn);
              loop (i + 1, z1', z2')
            end
      in
        loop (0, 0.0, 0.0); y
      end

    fun freqResponse ({ b0, b1, b2, a0, a1, a2 } : coeffs) { freq, rate } =
      let
        val w = 2.0 * pi * freq / real rate
        val e1 = C.complex (Math.cos (~w), Math.sin (~w))
        val e2 = C.complex (Math.cos (~2.0 * w), Math.sin (~2.0 * w))
        val num = C.add (C.complex (b0, 0.0),
                         C.add (C.scale (b1, e1), C.scale (b2, e2)))
        val den = C.add (C.complex (a0, 0.0),
                         C.add (C.scale (a1, e1), C.scale (a2, e2)))
        val h = C.divide (num, den)
      in
        { magnitude = C.abs h, phase = C.arg h }
      end
  end

  fun applyBiquad c x = Biquad.apply c x
  fun freqResponse c r = Biquad.freqResponse c r

  (* ------------------------------------------------------------------ *)
  (* convolution                                                        *)
  (* ------------------------------------------------------------------ *)

  fun convolve (a, b) =
    let
      val la = len a val lb = len b
    in
      if la = 0 orelse lb = 0 then Array.fromList []
      else
        let
          val outLen = la + lb - 1
          fun out n =
            let
              (* sum_{k} a[k] * b[n-k], valid k in [max(0,n-lb+1) .. min(n,la-1)] *)
              val lo = Int.max (0, n - lb + 1)
              val hi = Int.min (n, la - 1)
              fun loop (k, acc) =
                if k > hi then acc
                else loop (k + 1, acc + sub (a, k) * sub (b, n - k))
            in
              loop (lo, 0.0)
            end
        in
          tab (outLen, out)
        end
    end

  fun fftConvolve (a, b) = Fft.convolve (a, b)

  (* ------------------------------------------------------------------ *)
  (* FIR design                                                         *)
  (* ------------------------------------------------------------------ *)

  structure Fir =
  struct
    (* normalized sinc: sin(pi x)/(pi x), sinc 0 = 1 *)
    fun sinc x =
      if Real.abs x < 1e~12 then 1.0
      else Math.sin (pi * x) / (pi * x)

    (* ideal windowed lowpass (not yet normalized), fcn = cutoff/rate *)
    fun rawLowpass (fcn, taps, wk) =
      let
        val w = Window.make wk taps
        val center = real (taps - 1) / 2.0
        fun h i =
          let val m = real i - center
          in 2.0 * fcn * sinc (2.0 * fcn * m) * sub (w, i) end
      in
        tab (taps, h)
      end

    fun sumA a =
      let
        fun loop (i, acc) = if i >= len a then acc else loop (i + 1, acc + sub (a, i))
      in loop (0, 0.0) end

    fun scaleA (k, a) = tab (len a, fn i => k * sub (a, i))

    fun lowpass { cutoff, taps, window, rate } =
      let
        val fcn = cutoff / real rate
        val h = rawLowpass (fcn, taps, window)
        val s = sumA h
      in
        if Real.abs s < 1e~12 then h else scaleA (1.0 / s, h)
      end

    (* highpass by spectral inversion of a unity-DC lowpass *)
    fun highpass { cutoff, taps, window, rate } =
      let
        val lp = lowpass { cutoff = cutoff, taps = taps, window = window, rate = rate }
        val center = (taps - 1) div 2
        fun h i = (if i = center then 1.0 else 0.0) - sub (lp, i)
      in
        tab (taps, h)
      end

    (* bandpass = lowpass(hi) - lowpass(lo) *)
    fun bandpass { lo, hi, taps, window, rate } =
      let
        val lpHi = lowpass { cutoff = hi, taps = taps, window = window, rate = rate }
        val lpLo = lowpass { cutoff = lo, taps = taps, window = window, rate = rate }
      in
        tab (taps, fn i => sub (lpHi, i) - sub (lpLo, i))
      end

    fun apply taps x = convolve (taps, x)
  end

  fun applyFir taps x = Fir.apply taps x

  (* ------------------------------------------------------------------ *)
  (* STFT / ISTFT                                                       *)
  (* ------------------------------------------------------------------ *)

  structure Stft =
  struct
    type params = { frame : int, hop : int, window : Window.kind }

    fun numFrames (sigLen, frame, hop) =
      if sigLen <= 0 orelse frame <= 0 orelse hop <= 0 then 0
      else (sigLen + hop - 1) div hop

    fun analyze ({ frame, hop, window } : params) x =
      let
        val n = len x
        val w = Window.make window frame
        val nf = numFrames (n, frame, hop)
        fun frameAt m =
          let
            val start = m * hop
            val buf = tab (frame, fn j =>
              let val idx = start + j
              in if idx < n then sub (x, idx) * sub (w, j) else 0.0 end)
          in
            Fft.rfft buf
          end
      in
        List.tabulate (nf, frameAt)
      end

    fun synthesize ({ frame, hop, window } : params) spectra =
      let
        val w = Window.make window frame
        val nf = List.length spectra
      in
        if nf = 0 then Array.fromList []
        else
          let
            val outLen = (nf - 1) * hop + frame
            val acc = Array.array (outLen, 0.0)
            val norm = Array.array (outLen, 0.0)
            fun place (m, spec) =
              let
                val time = Fft.irfft spec
                val start = m * hop
                fun loop j =
                  if j >= frame then ()
                  else
                    let
                      val pos = start + j
                      val wj = sub (w, j)
                    in
                      Array.update (acc, pos, sub (acc, pos) + wj * sub (time, j));
                      Array.update (norm, pos, sub (norm, pos) + wj * wj);
                      loop (j + 1)
                    end
              in
                loop 0
              end
            val _ = List.foldl (fn (s, m) => (place (m, s); m + 1)) 0 spectra
          in
            tab (outLen, fn i =>
              let val d = sub (norm, i)
              in if Real.abs d < 1e~12 then 0.0 else sub (acc, i) / d end)
          end
      end

    fun magnitude spec = tab (len spec, fn i => C.abs (sub (spec, i)))

    fun spectrogram p x = List.map magnitude (analyze p x)
  end
end
