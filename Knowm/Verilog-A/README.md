## About 

| Model Class | Model Language | Simulator |
|---|---|---|
|Ill-formed|Verilog-A|Xyce|

## Resources

1. [Build Xyce from Source for ADMS Verilog-A Model Integration](http://knowm.org/build-xyce-from-source-for-adms-verilog-a-model-integration/)
1. [Well-posed Memristor Modeling with Xyce and Verilog-A](http://knowm.org/well-posed-memristor-modeling-with-xyce-and-verilog-a/)

## Instructions

    cd ~/workspaces/workspace_knowm/memristor-models-4-all/Knowm/Verilog-A
    buildxyceplugin -o metastableswitch *.va /usr/local/lib
    Xyce -plugin /usr/local/lib/metastableswitch.so memristor_sim.cir
    gnuplot
    plot '~/workspaces/workspace_knowm/memristor-models-4-all/Knowm/Verilog-A/memristor_sim.cir.prn' using 3:4 with lines title "I-V"


