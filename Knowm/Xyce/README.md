## About 

| Model Class | Model Language | Simulator |
|---|---|---|
|Ill-formed|Native Xyce Device|Xyce|

## Instructions

See: <http://knowm.org/the-mean-metastable-switch-model-in-xyce/>

## Sources

1. Original Source of MSS Model: <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085175>

## How to Run Model

    Xyce /path/to/.../knowm.cir
    gnuplot
    plot '/path/to/.../knowm.cir.prn' using 3:4 with lines title "I-V"



