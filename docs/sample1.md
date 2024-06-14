---
bibliography:
- assets/references.bib
title: Surface Flux Transport Model User Manual
---

# The Basics 
## Magnetic field evolution on the solar surface

The surface flux transport (SFT) model, which solves the radial
component of the magnetic field on the solar surface, has demonstrated
remarkable effectiveness in simulating the dynamics of the large-scale
magnetic field on the photosphere. The governing equation can be written
as,

```{math}
a^2 + b^2 = c^2
```


```{math}
\frac{\partial B}{\partial t} = \frac{D}{R_\odot^2}\left[\frac{\partial}{\partial s}\left((1-s^2)\frac{\partial B}{\partial s}\right) + \frac{1}{1-s^2}\frac{\partial^2B}{\partial\phi^2}\right] - \frac{\partial}{\partial s}\left[\frac{v_s(s)}{R_\odot}\sqrt{1-s^2} B\right]
```


$\frac{\partial B}{\partial t} = \frac{D}{R_\odot^2}\left[\frac{\partial}{\partial s}\left((1-s^2)\frac{\partial B}{\partial s}\right) + \frac{1}{1-s^2}\frac{\partial^2B}{\partial\phi^2}\right] - \frac{\partial}{\partial s}\left[\frac{v_s(s)}{R_\odot}\sqrt{1-s^2} B\right]- \Omega(s)\frac{\partial B}{\partial\phi}$
