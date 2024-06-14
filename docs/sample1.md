---
bibliography:
- assets/references.bib
title: Surface Flux Transport Model User Manual
---

# The Basics []{#chapter:basics label="chapter:basics"}

## Magnetic field evolution on the solar surface

The surface flux transport (SFT) model, which solves the radial
component of the magnetic field on the solar surface, has demonstrated
remarkable effectiveness in simulating the dynamics of the large-scale
magnetic field on the photosphere. The governing equation can be written
as,

$$\frac{\partial B}{\partial t} = \frac{D}{R_\odot^2}\left[\frac{\partial}{\partial s}\left((1-s^2)\frac{\partial B}{\partial s}\right) + \frac{1}{1-s^2}\frac{\partial^2B}{\partial\phi^2}\right] - \frac{\partial}{\partial s}\left[\frac{v_s(s)}{R_\odot}\sqrt{1-s^2} B\right]- \Omega(s)\frac{\partial B}{\partial\phi},
\label{eqn:sft2d}$$ where $s =$ sin$\theta$, $D$ is the magnetic
diffusivity, $\Omega (s)$ is the angular velocity in east-west direction
on a sine-latitude grid, $v_s (s)$ is the flow profile along north-south
direction on a sine-latitude grid and $\phi$ is the longitude. As the
velocity profiles involved in transporting the magnetic flux on the
photosphere is a function of latitude only, we can simplify this
equation by taking average in the longitudinal direction which will
improve the computational efficiency and provide us a lesser parameter
space to comprehensively explore the dynamics. After averaging the
$B_r$,
$$\overline{B}(s,t)=(2\pi)^{-1}\int_0^{2\pi}B(s,\phi,t)\,\mathrm{d}\phi.
\label{eqn:sft1d}$$ Using this reduced form of $B_r$, we can re-write
equation-[\[eqn:sft2d\]](#eqn:sft2d){reference-type="ref"
reference="eqn:sft2d"} as,
$$\frac{\partial\overline{B}}{\partial t} = \frac{\partial}{\partial s}\left[\frac{D}{R_\odot^2}(1-s^2)\frac{\partial\overline{B}}{\partial s} - \frac{v_s(s)}{R_\odot}\sqrt{1-s^2}\overline{B}\right],
\label{eqn:evol}$$ where the north-south velocity is,
$$v_s(s) = D_us(1-s^2)^{p/2}.
\label{eqn:merid}$$ In this form the velocity, $s =$ sin$\theta$ and
$D_u$ controls the amplitude of the function.

![Example flow profile in the north-south direction using
equation-[\[eqn:merid\]](#eqn:merid){reference-type="ref"
reference="eqn:merid"}.](assets/MC_flow250.pdf){#fig:merid width="70%"}

Presently the code solves
equation-[\[eqn:evol\]](#eqn:evol){reference-type="ref"
reference="eqn:evol"} on a sine latitude grid,
figure-[1.2](#fig:grid){reference-type="ref" reference="fig:grid"} shows
an example of the grid point distribution on a ($s$,$\phi$) grid
containing (10,20) points.

![Example grid visualization with 10, 20 points in sine latitude ($s$)
and longitude ($\phi$ direction
respectively.](assets/grid_structure_sft1D.pdf){#fig:grid width="85%"}

### Bipolar approximation

We introduce a source function to model the approximating bipolar
magnetic region (BMR) for an observed SHARP. The location of the center
of the BMR we use the locations of the positive and negative polarity
positions $(s_+, \phi_+)$ and $(s_-,\phi_-)$ on the computational grid.
Here $s$ denotes sine-latitude and $\phi$ denotes (Carrington)
longitude. We compute,

::: itemize
centroid of the BMR,\
\
$$\begin{aligned}
    s_0 = \frac12(s_+ + s_-),\qquad \phi_0 = \frac12(\phi_+ + \phi_-)
    \label{eqn:center}
    
\end{aligned}$$ polarity separation, which is the heliographic angle,\
\
$$\begin{aligned}
    \rho = \arccos\left[s_+s_- + \sqrt{1-s_+^2}\sqrt{1 - s_-^2}\cos(\phi_+-\phi_-) \right]
    \label{eqn:separation}
    
\end{aligned}$$ the tilt angle with respect to the equator, given by,\
\
$$\begin{aligned}
    \gamma = \arctan\left[\frac{\arcsin(s_+) - \arcsin(s_-)}{\sqrt{1-s_0^2}(\phi_- - \phi_+)}\right]
    \label{eqn:tilt}
    
\end{aligned}$$
:::

Together with the unsigned flux, $|\Phi|$, these parameters define the
BMR for our chosen functional form. For an untilted BMR centered at
$s=\phi=0$, this functional form is defined as
$$B(s,\phi) = F(s,\phi) = -B_0\frac{\phi}{\rho}\exp\left[-\frac{\phi^2 + 2\arcsin^2(s)}{(a\rho)^2}\right],
\label{eqn:bmr}$$ where the amplitude $B_0$ is scaled to match the
corrected flux of the observed region on the computational grid. To
account for the location $(s_0,\phi_0)$ and tilt $\gamma$ of a general
region, we set $B(s,\phi) = F(s',\phi')$, where $(s',\phi')$ are
spherical coordinates in a frame where the region is centered at
$s'=\phi'=0$ and untilted. Figure
[1.3](#fig:example-bmr){reference-type="ref"
reference="fig:example-bmr"} shows an example BMR.

![A BMR centered at 25$^{\circ}$ latitude and 22.5$^{\circ}$ longitude
with a tilt of 30$^{\circ}$ with respect to the equator calculated using
equation-[\[eqn:bmr\]](#eqn:bmr){reference-type="ref"
reference="eqn:bmr"}](assets/example_bmr.png){#fig:example-bmr
width="50%"}

### Coordinate transformation

The coordinate transformation from the frame $(s,\phi)$ where the BMR is
centred at $s=\phi=0$ to the frame $(s',\phi')$ where it is centred at
$(s_0,\phi_0)$ with tilt $\gamma$. This amounts to a rotation, which is
easiest to express in Cartesian coordinates
$$x = \cos\phi\sqrt{1-s^2}, \quad y=\sin\phi\sqrt{1-s^2},\quad z=s.$$
Multiplying by the rotation matrices for the sequence of rotations
indicated in Figure [1.4](#fig:rotation-bmr){reference-type="ref"
reference="fig:rotation-bmr"} shows that Cartesian coordinates in the
rotated frame are $$\begin{aligned}
\nonumber
\left[
\begin{array}{c}
x'\\
y'\\
z'
\end{array}
\right] &=&
\left[
\begin{array}{ccc}
1\; & 0\; & 0\\
0\; & \cos\gamma\; & -\sin\gamma\\
0\; & \sin\gamma\; & \cos\gamma
\end{array}
\right]
\left[
\begin{array}{ccc}
\cos\lambda_0\; & 0\; & \sin\lambda_0\\
0\; & 1\; & 0\\
-\sin\lambda_0\; & 0\; & \cos\lambda_0
\end{array}
\right]
\cdot\\ 
&&\cdot
\left[
\begin{array}{ccc}
\cos\phi_0\; & \sin\phi_0\; & 0\\
-\sin\phi_0\; & \cos\phi_0\; & 0\\
0 & 0 & 1
\end{array}
\right]
\left[
\begin{array}{c}
x\\
y\\
z
\end{array}
\right],
\label{eqn:rotation-trans}
\end{aligned}$$ where $s_0=\sin\lambda_0$. From these we determine
$\phi'=\arctan(y'/x')$ and $s' = z'$.

The parameter $a$ in
equation-[\[eqn:bmr\]](#eqn:bmr){reference-type="ref"
reference="eqn:bmr"} controls the size of the BMR relative to the
separation, $\rho$, of the original polarity centroids. For given values
of $\lambda_0$, $\gamma$, and $\rho$, and $B_0$ chosen to give the
required magnetic flux, the parameter $a$ may be chosen to control the
axial dipole moment of the BMR. A good match to the axial dipole moment
of the original SHARP is obtained with $a=0.56$, and the same value
works for every region.

![A BMR centered at 45$^{\circ}$ latitude and 45$^{\circ}$ longitude
with a tilt of 45$^{\circ}$ with respect to the equator is shown in (a).
Coordinate transformation operations along $\phi$,$\lambda$ and $\gamma$
are shown in (b), (c) and (d). The corresponding transformation matris
is given by
equation-[\[eqn:rotation-trans\]](#eqn:rotation-trans){reference-type="ref"
reference="eqn:rotation-trans"}](assets/rotation_bmr.png){#fig:rotation-bmr
width="85%"}
