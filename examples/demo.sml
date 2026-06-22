(* demo.sml - design an RBJ biquad lowpass and a windowed-sinc FIR lowpass,
   print their frequency-response magnitudes and tap lists. Deterministic:
   identical output on every run and under both MLton and Poly/ML. *)

structure W = Dsp.Window
structure BQ = Dsp.Biquad

fun fmtReal n r =
  let
    val s = Real.fmt (StringCvt.FIX (SOME n)) r
    val s = if String.isPrefix "~" s then "-" ^ String.extract (s, 1, NONE) else s
    (* collapse a negative zero ("-0.000") to positive: the only digits are 0,
       so the sign is rounding noise and must not differ across compilers *)
    val isZero =
      CharVector.all (fn c => c = #"0" orelse c = #"." orelse c = #"-") s
  in
    if isZero andalso String.isPrefix "-" s then String.extract (s, 1, NONE) else s
  end

val rate = 48000
val lp = BQ.lowpass { cutoff = 1000.0, q = 0.707, gainDb = 0.0, rate = rate }

val () = print "Biquad lowpass (cutoff 1000 Hz, Q 0.707, rate 48000):\n"
val () = print ("  coeffs b = [" ^ fmtReal 6 (#b0 lp) ^ ", " ^ fmtReal 6 (#b1 lp)
                ^ ", " ^ fmtReal 6 (#b2 lp) ^ "]\n")
val () = print ("  coeffs a = [" ^ fmtReal 6 (#a0 lp) ^ ", " ^ fmtReal 6 (#a1 lp)
                ^ ", " ^ fmtReal 6 (#a2 lp) ^ "]\n")
val () = print "  |H(f)| magnitude response:\n"
val () =
  List.app
    (fn f =>
       let val m = #magnitude (BQ.freqResponse lp { freq = f, rate = rate })
       in print ("    f = " ^ StringCvt.padLeft #" " 6 (Int.toString (Real.round f))
                 ^ " Hz  |H| = " ^ fmtReal 4 m ^ "\n") end)
    [0.0, 500.0, 1000.0, 2000.0, 8000.0, 20000.0]

val () = print "\nWindowed-sinc FIR lowpass (cutoff 800 Hz, 11 taps, Hamming, rate 8000):\n"
val fir = Dsp.Fir.lowpass { cutoff = 800.0, taps = 11, window = W.Hamming, rate = 8000 }
val () = print ("  taps = [\n")
val () =
  List.app
    (fn i => print ("    h[" ^ Int.toString i ^ "] = " ^ fmtReal 6 (Array.sub (fir, i)) ^ "\n"))
    (List.tabulate (Array.length fir, fn i => i))
val () = print "  ]\n"
val () = print ("  sum of taps (DC gain) = " ^ fmtReal 6 (Array.foldl op+ 0.0 fir) ^ "\n")

val () = print "\nDirect vs FFT convolution of [1,-2,3,0.5] * [0.5,1,-0.5]:\n"
val a = Array.fromList [1.0, ~2.0, 3.0, 0.5]
val b = Array.fromList [0.5, 1.0, ~0.5]
fun showArr arr =
  "[" ^ String.concatWith ", "
          (List.map (fn i => fmtReal 4 (Array.sub (arr, i)))
                    (List.tabulate (Array.length arr, fn i => i))) ^ "]"
val () = print ("  direct      = " ^ showArr (Dsp.convolve (a, b)) ^ "\n")
val () = print ("  fft-based   = " ^ showArr (Dsp.fftConvolve (a, b)) ^ "\n")
