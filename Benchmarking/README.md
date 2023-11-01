# Benchmarking

## Introduction

Whenever a change is made to _aLENS_ there is a possibility of increasing or decreasing the speed of simulations. To ensure that the changes made to _aLENS_ are not slowing down the simulations, we have a set of benchmarks that we run after every change. The benchmarks are typically run on a single node with varying MPI ranks and threads. Changes may cause slow down in one section of the _aLENS_ time step, but may speed up another section. Therefore, we have a series of benchmarks that run several computationally intensive tests based on physical systems that, when combined, cover the entire _aLENS_ time step. If a slow down is detected, there are specific tests designed to examine separate functionalities of the code to pinpoint the cause of the slow down. These two categories are called _general_ and _specific_ benchmarks, respectively.

Benchmarking takes place with systems large enough that you see speed up scaling with different numbers of threads and MPI ranks. Often changes will cause simulations to lose these scaling properties. Therefore, we run benchmarks with varying number of threads and MPI ranks to ensure that the scaling properties are maintained.

## Features to test

- [ ] Rod-Rod steric interactions\*
- [ ] Sphere-Sphere steric interactions\*
- [ ] Sphere-Rod steric interactions
- [ ] Rod-Rod binding kinetics\*
- [ ] Sphere-Sphere binding kinetics\*
- [ ] Rod-Sphere binding kinetics
- [ ] Rod-rod hinge constraint\*
- [ ] Rod-rod bending constraint\*
- [ ] rod-rod extending constraint
- [ ] Sphere-rod extending constraint
- [ ] Sphere-sphere extending constraint\*
- [ ] Tri-sphere bending constraint\*

\* = tested in general benchmarks

## Running Benchmarks

# General Benchmarks

## Benchmark 1 (Rod-Crosslinker Aster Gas)

Functionality tested:

- Periodic boundary conditions
- Rod steric interactions
- Motor walking (with end pausing)
- Binding kinetics to rods

Description: N rods are placed in a periodic box with M crosslinking motors Motors bind to, unbind from, and crosslink. Motors move to rod plus ends with a force-velocity while crosslinking. Because of the motor end pausing, the rods form dense asters with motors concentrated at there center.

## Benchmark 2 (Collapsed semi-flexible, sticky bead-spring filament)

Functionality tested:

- Open boundary conditions
- Sphere steric interactions
- Sphere-sphere extending constraint
- Tri-bead bending constraints
- Binding kinetics to beads with one head permanently bound

Description: A single semi-flexible bead-spring filament of N beads is held together by crosslinkers that have their first head permanently bound to a bead. Each bead has exactly one crosslinker head permanently bound but can have heads from other crosslinkers transiently bound to it. Semi-flexibility is achieved by a bending constraint between three consecutive beads.

## Benchmark 3 (Confined semi-flexible filament with crosslinkers)

Functionality tested:

- Spherically confined boundary conditions
- Rod steric interactions
- Rod binding kinetics
- Rod-rod hinge constraint
- Rod-rod bending constraint

Description: N filaments are confined in a spherical shell of size R. The filaments are semi-flexible and inextensible thanks to the pinning constraint. M crosslinkers are also confined to the interior of the shell. Crosslinkers bind to, unbind from, and transiently crosslink filaments creating a gelled network of filaments.

# Specific Benchmarks
