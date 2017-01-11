## About 

| Model Class | Model Language | Simulator |
|---|---|---|
|Ill-formed|Native Xyce Device|Xyce|

## Instructions

See: <http://knowm.org/native-memristor-device-development-in-xyce/>

## Sources

1. Original Source of SPICE Model: <http://www.radioeng.cz/fulltexts/2009/09_02_210_214.pdf>

## How to Run Model

    Xyce /path/to/.../joglekar.cir
    gnuplot
    plot '/path/to/.../joglekar.cir.prn' using 3:4 with lines title "I-V"
