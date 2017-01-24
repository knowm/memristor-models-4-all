## About 

| Model Class | Model Language | Simulator |
|---|---|---|
|Well-formed|Verilog-A|Xyce|

## Instructions

See: <http://knowm.org/well-posed-memristor-modeling-with-xyce-and-verilog-a/>


    cd .../memristor-models-4-all/Generic/Verilog-A
    buildxyceplugin -o hys *.va /usr/local/lib
    Xyce -plugin /usr/local/lib/hys.so test_hys_forward.cir
    Xyce -plugin /usr/local/lib/hys.so test_hys_reverse.cir
    Xyce -plugin /usr/local/lib/hys.so test_hys_transient.cir

    gnuplot
    plot '.../memristor-models-4-all/Generic/Verilog-A/test_hys_forward.cir.prn' using 2:3 with lines title "forward DC Sweep", '.../memristor-models-4-all/Generic/Verilog-A/test_hys_reverse.cir.prn' using 2:3 with lines title "reverse DC Sweep", '.../memristor-models-4-all/Generic/Verilog-A/test_hys_transient.cir.prn' using 3:4 with lines title "Transient"

    Xyce -plugin /usr/local/lib/hys.so test_hys_homotopy.cir
    plot '.../memristor-models-4-all/Generic/Verilog-A/test_hys_homotopy.cir.HOMOTOPY.prn' using 3:4 with lines title "I"
    plot '.../memristor-models-4-all/Generic/Verilog-A/test_hys_homotopy.cir.HOMOTOPY.prn' using 3:5 with lines title "S"

