## About 

| Model Class | Model Language | Simulator |
|---|---|---|
|Ill-formed|Native Xyce Device|Xyce|

## Instructions

See: <http://knowm.org/native-memristor-device-development-in-xyce/>

## Sources

1. Original Source of Yakopcic Model: <https://cyakopcic1.files.wordpress.com/2014/02/a-memristor-device-model.pdf>

## How to Run Model

    Xyce /path/to/.../yakopcic.cir
    gnuplot
    plot '/path/to/.../yakopcic.cir.prn' using 3:4 with lines title "I-V"
