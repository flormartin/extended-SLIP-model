# Extended SLIP Model
This repository shows the work of the practical course *Dodo Alive* at the [Munich Institute of Robotics and Machine Intelligence](https://www.mirmi.tum.de/).

##Introduction

In a group project with three students, the goal was to find optimal control parameters to produce periodical walking gaits in an extended single-legged SLIP model. The work is based on a model presented by [Peuker et al. (2012)](https://ieeexplore.ieee.org/abstract/document/6290742) that is used for comparison to the simple SLIP model in humanoid walking.

## Extensions

In this work, the SLIP model is physically extended by a knee with an additional *leg mass*. Unlike the original model, the extended model includes *non-conservative forces* in the form of two mass-spring-damper systems about the hip and knee. Lastly, two-legged walking of the model is tested with an estimation for physical parameters of a Dodo.

## Actuation

The underlying paper suggests an alternating hip actuation policy between stance and flight phase. As an extension to this policy, hip actuation using a fourier type input is tested and compared to the more rigid actuation.

## Software
MatLab, CasADi


