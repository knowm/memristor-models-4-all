## About 

| Model Class | Model Language | Simulator |
|---|---|---|
|Ill-formed|Verilog-A|Xyce|

## Instructions

    cd ~/workspaces/workspace_knowm/memristor-models-4-all/Metastable_Switch/Verilog-A
    buildxyceplugin -o metastableswitch *.va /usr/local/lib
    Xyce -plugin /usr/local/lib/metastableswitch.so memristor_sim.cir
    gnuplot
    plot '~/workspaces/workspace_knowm/memristor-models-4-all/Metastable_Switch/Verilog-A/memristor_sim.cir.prn' using 3:4 with lines title "I-V"

## Note

This is a work in progress until Xyce/ADMS has implemented the `distNormal` function.