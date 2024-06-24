---
bibliography:
- assets/references.bib
title: Surface Flux Transport Model User Manual
---

# Theoretical background

## Magnetic field evolution on the solar surface

The surface flux transport (SFT) model, which solves the radial
component of the magnetic field on the solar surface, has demonstrated
remarkable effectiveness in simulating the dynamics of the large-scale
magnetic field on the photosphere. The governing equation can be written
as,

```{math}
\frac{\partial B}{\partial t} = \frac{D}{R_\odot^2}\left[\frac{\partial}{\partial s}\left((1-s^2)\frac{\partial B}{\partial s}\right) + \frac{1}{1-s^2}\frac{\partial^2B}{\partial\phi^2}\right] - \frac{\partial}{\partial s}\left[\frac{v_s(s)}{R_\odot}\sqrt{1-s^2} B\right] - \Omega(s)\frac{\partial B}{\partial\phi}
```

where {math}`s` = sin{math}`\theta`, D is the magnetic diffusivity, {math}`\Omega (s)`
is the angular velocity in east-west direction on a sine-latitude grid,
{math}`v_s (s)` is the flow profile along north-south direction on a
sine-latitude grid and {math}`\phi` is the longitude. As the velocity profiles
involved in transporting the magnetic flux on the photosphere is a
function of latitude only, we can simplify this equation by taking
average in the longitudinal direction which will improve the
computational efficiency and provide us a lesser parameter space to
comprehensively explore the dynamics. After averaging the {math}`B_r`,
```{math}
\overline{B}(s,t)=(2\pi)^{-1}\int_0^{2\pi}B(s,\phi,t)\,\mathrm{d}\phi.
```
Using this reduced form of {math}`B_r`, we can re-write the 
equation as,
```{math}
\frac{\partial\overline{B}}{\partial t} = \frac{\partial}{\partial s}\left[\frac{D}{R_\odot^2}(1-s^2)\frac{\partial\overline{B}}{\partial s} - \frac{v_s(s)}{R_\odot}\sqrt{1-s^2}\overline{B}\right],
```
where the north-south velocity is, {math}`v_s(s) = D_us(1-s^2)^{p/2}.` In
this form the velocity, {math}`s` = sin{math}`\theta` and {math}`D_u` controls the
amplitude of the function.

```{figure} MC_flow.png
:scale: 80 %
:align: "center"
:alt: Meridional flow profile

Example meridional flow profile.
```

## BMR modeling algorithm

We introduce a source function to model the approximating bipolar
magnetic region (BMR) for an observed SHARP. The location of the center
of the BMR we use the locations of the positive and negative polarity
positions {math}`(s_+, \phi_+)` and {math}`(s_-,\phi_-)` on the computational grid.
Here $s$ denotes sine-latitude and $\phi$ denotes (Carrington)
longitude. We compute,

1. centroid of the BMR,\
\
```{math}
    s_0 = \frac12(s_+ + s_-),\qquad \phi_0 = \frac12(\phi_+ + \phi_-)
    \label{eqn:center}
``` 
2. polarity separation, which is the heliographic angle,\
\
```{math}
    \rho = \arccos\left[s_+s_- + \sqrt{1-s_+^2}\sqrt{1 - s_-^2}\cos(\phi_+-\phi_-) \right]
    \label{eqn:separation}
```
3. the tilt angle with respect to the equator, given by,\
\
```{math}
    \gamma = \arctan\left[\frac{\arcsin(s_+) - \arcsin(s_-)}{\sqrt{1-s_0^2}(\phi_- - \phi_+)}\right]
    \label{eqn:tilt}
```
