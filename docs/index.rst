.. Example documentation master file, created by
   sphinx-quickstart on Sat Sep 23 20:35:12 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Surface Flux Transport 1D documentation
=======================================

We can only observe half of the solar surface routinely. However, a global observation is
essential for understanding and forecasting solar activity. Sun is magnetically active and this 
activity cycle lasts for about 11 years during which the magnetic field distribution on the solar
surface -- the photosphere, evolves and generates a vast range of events e.g., emergece and 
decay of sunspots, solar flares, coronal mass ejections, flow of solar wind and energetic particles
into the heliosphere etc. In order to understand these we can use theoretical models to compute
the global magnetic field configuration on the photosphere and later use them to drive different 
models. Surface flux transport model solves the radial component of the magnetic induction equation 
on the solar surface by using prescribed flow profiles and source terms. This distribution provides
a numerical model of SFT developed in FORTRAN by folowing `Yeates (2020) <https://doi.org/10.1007/s11207-020-01688-y>`_.

Link to the GitHub repository: `https://sr-dash.github.io/SFT-1D/ <https://sr-dash.github.io/SFT-1D/>`_.

Link to the resources used in `Yeates (2020) <https://doi.org/10.1007/s11207-020-01688-y>`_: `sharps-bmrs <https://github.com/antyeates1983/sharps-bmrs>`_.

Here is a short guide on how to build the SFT 1D code, install the dependencies and check the output.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart.md
   notebooks/plotting-results.ipynb
   sft-theory.md

..
   some-feature.md
   another-feature.md
   sample0.md
   sample1.md
