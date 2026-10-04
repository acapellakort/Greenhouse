# YourPackage Documentation

## Introduction
greenhouse_inference is a Julia-based tool developed to perform Bayesian inference on model
parameters in greenhouse crop environments. 

Greenhouse crop growth is influenced by numerous factors, including temperature, humidity, CO₂ concentration, vapor pressure, and the operational controls of greenhouse systems like heating, ventilation, and shading. To analyze these parameters, greenhouse_inference uses data input from two essential CSV files:

1. Growth Environment Data:

Contains greenhouse conditions and control settings, with columns such as:

Time: Time stamp of data recording.
T1: Crop temperature (in Kelvin, often unavailable).
T2: Air temperature (in Kelvin).
C1: CO₂ concentration (in mg/m³).
V1: Vapor pressure (in Pa).
RH: Relative humidity (in %).
Tpipe: Heating pipe temperature (in Kelvin).

U1–U12: Control functions for various greenhouse systems, including:
U1: Thermal screen
U2: Fan-pad system
U3: Mechanical cooling
U4: Air heater
U5: Shading screen
U6: Side windows opening
U7: Forced ventilation
U8: Roof windows control
U9: Fog system
U10: CO₂ source
U11: Pipe temperature
U12: Light control

1. Weather Data:

Provides external environmental data relevant to greenhouse conditions, with columns such as:

Time: Time stamp of data recording.
Temperature: Ambient temperature 2 meters above ground (°C).
Relative Humidity: Humidity 2 meters above ground (%).
Cloud Cover: Total, high, medium, and low cloud cover (%) at different atmospheric layers.
Shortwave Radiation: Solar radiation at the surface (W/m²).
Wind Data: Wind speed and direction at multiple elevations, including:
10 m above ground
80 m above ground
900 mb pressure level

With this structured data, greenhouse_inference enables users to analyze observational data, incorporate control system impacts, and produce probabilistic estimates for model parameters in greenhouse crop growth models. Through Bayesian inference, it supports data-informed decision-making, aiding in the optimization of greenhouse operations and the advancement of crop management practices


## Prerequisites
- Julia version 1.10.5 or above
- Required packages: 
`CSV`, `DataFrames`,  `Unitful`, `Dates`, `LaTeXStrings`, `SparseArrays`
`DifferentialEquations`,  `JTwalk`, `Statistics`

## Installation
Provide a step-by-step guide on installing the software.
To install YourPackage, run:
_julia_

- using Pkg
- Pkg.add("YourPackage")

## Basic Usage
Describe essential commands or functions, with simple examples.

The dictionary `QoI_dict` inlcudes all constants names and prior supports
available for the inference. 

QoI_dict = Dict("constant name"    => [min, known_value max ],...)  

To SetUp the inference, the variables to be infeered are given in the vector:
`global QoI = ["alpha1", "alpha2", "alpha4",  "beta2", "gamma3", "gamma4", "nu4", "nu4CO2", "nu7", "eta1", "tau1", "psi2" ]`   

The MCMC chain length and burn-in are set up in the variables:
chain_length = 50000
burn_in      =  15000


## Advanced Usage (optional) 
Show more complex or specific functionalities.

## Troubleshooting 
Include common issues and solutions.

## References 
Link to any additional resources or documentation.
