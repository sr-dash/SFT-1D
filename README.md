# Surface Flux Transport (SFT) 1D

[![Build and Test](https://github.com/sr-dash/SFT-1D/actions/workflows/main.yml/badge.svg)](https://github.com/sr-dash/SFT-1D/actions/workflows/main.yml)
[![Sphinx Doc](https://github.com/sr-dash/SFT-1D/actions/workflows/documentation.yml/badge.svg)](https://github.com/sr-dash/SFT-1D/actions/workflows/documentation.yml)
[![HitCount](https://hits.dwyl.com/sr-dash/SFT-1D.svg?style=flat-square&show=unique)](http://hits.dwyl.com/sr-dash/SFT-1D)



Surface Flux Transport modeling in one dimension (latitude). 

## Documentation

A general introcuction to the model and preliminary documentation can be found [here](https://sr-dash.github.io/SFT-1D/). 

The input dataset can be downloaded from Zenodo at [https://zenodo.org/records/13830728](https://zenodo.org/records/13830728).

These pages are currently under development. 

For ease of installing the dependecies, supporting packages and building the executable please follow the steps mentioned below. 

* File reading/writing operatins in the SFT model uses NETCDF library. For a fresh installation of the libraries and it's dependencies, run the following:

```shell
$ bash install_netcdf.sh
```

* For building the executable for SFT model, run the following:

```shell
$ make
```

* For plotting the results a sample plotting routine is provided in the repository. To install the dependencies please run the following command before executing the python script.

```shell
$ pip install -r requirements.txt
```

Feel free to explore the documentation, and if you have any questions or feedback, please [reach out](dashs@hawaii.edu).
