## About 

| Model Class | Model Language | Simulator |
|---|---|---|
|Ill-formed|Verilog-A|Xyce|

## Instructions

    cd .../memristor-models-4-all/Joglekar_Linear_Resistance_Switch/Verilog-A
    buildxyceplugin -o memristor *.va /usr/local/lib
    Xyce -plugin /usr/local/lib/memristor.so memristor_sim.cir
    gnuplot
    plot '.../memristor-models-4-all/Joglekar_Linear_Resistance_Switch/Verilog-A/memristor_sim.cir.prn' using 3:4 with lines title "I-V"