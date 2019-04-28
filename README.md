## Introduction

This project aims to be a collaborative effort to collect under one roof all the memristor models published to date. We welcome contribution of your own models or implementations of other people's models. We also welcome any suggestions and/or corrections. This is a collaborative effort. 

This project may be useful for you whether you are trying to simulate other people's existing models or as a way to start your own model development, as these working examples can help get you up to speed quickly. If it was useful for you please consider adding your own working example by submitting a pull request with the model files and a README file. We hope that eventually this repository can be used to quantitatively compare and benchmark memristor models in terms of model accuracy, model robustness and simulation speed to name a few possibilities.

*WE DO NOT PROVIDE FREE TECHNICAL SUPPORT* 

## Model Classes

In order to bring some sort of ordering to the different types of memristor models and their implementations, here are three model classes we can use for easier classification.

### SPICE Subcircuit

Around 2008, the first memristor models designed for SPICE started showing up. A [paper by Biolek et al.](http://www.radioeng.cz/fulltexts/2009/09_02_210_214.pdf) showed how to implement a memristor model in SPICE as a `.subckt`. I call this a "hacked" model because in order to implement it you are constrained in using only existing SPICE elements such as resistors, capacitors, and various specialized voltage and current sources, hacking together a model from these pieces. In a post on knowm.org: [The Joglekar Resistance Switch Memristor Model in LTSpice](http://knowm.org/the-joglekar-resistance-switch-memristor-model-in-ltspice/), the SPICE subcircuit method is explained in great detail, demonstrating a working model and simulation in LTSpice.

### Transient Limited

In another class of models, the memristor is usually defined in a more flexible and less rigid and less "hacked" manner in a language such as Verliog-A, MATLAB or ModSPEC. In such simulation environments, you are given more or less a blank slate to formulate your model using the language's full API and methods. Inside your model implementation, you are given access to the global voltage drop across the device and from there you define the current though the device. You are also given the current simulation time, and you can use this to calculate how long a certain voltage was applied since the last time-step, important in modeling memristors given their memory behavior. It only works for transient analysis, not DC sweeps.

### Complete

In this well-posed model class, the simulation will work for both transient and DC sweep. Problems are avoided by following a strict set of guidelines for avoiding common pitfalls.

In May, 2016, Wang and Roychowdhury published a paper titled [Well-Posed Models of Memristive Devices](https://arxiv.org/abs/1605.04897), in which they plainly list out all of the problems with the current approaches and implementations of memristor models:

> Common ill-posedness mechanisms in models include division-by-zero errors, often due to expressions like 1/(x−a), which become unbounded (and feature a “doubly infinite” jump) at x=a; the use of log(),  1/(x−a) or sqrt() without ensuring that their arguments are always positive, regardless of bias; the fundamental misconception that non-real (i.e., complex) numbers or infinity are legal values for device models (they are not!); and “sharp”/“pointy” functions like |x|, whose derivatives are not continuous.  Another key aspect of well-posedness is that the model’s equations must produce mathematically valid outputs for any mathematically valid input to the device. Possibly the most basic kind of input is one that is constant (“DC”) for all time. DC solutions (“operating points”) are fundamental in circuit design; they are typically the kind of solution a designer seeks first of all, and are used as starting points for other analyses like transient, small signal AC, etc. If a model’s equations do not produce valid DC outputs given DC inputs, it fails a very fundamental well-posedness requirement. For example, the equation d/dt o(t) = i(t) is ill posed, since no DC (constant) solution for the output o(t) is possible if the input i(t) is any non-zero constant. Such ill posedness is typically indicative of some fundamental physical principle being violated; for example, in the case of d/dt o(t) = i(t), the system is not strictly stable. Indeed, a well-posed model that is properly written and implemented should work consistently in every analysis (including DC, transient, AC, etc.).

## Copyright

Since this is a collection of people's published work and original work, we think the best keep track of copyrights by making sure to include in the model files the original copyright header or, if missing, adding one that makes a best effort to source the original copyright holder. If unknown, at least link to a source where the model was found. The README files are copyrighted by the original and contributor authors.

