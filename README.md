# lsdo_project_template

<!---
[![Python](https://img.shields.io/pypi/pyversions/lsdo_project_template)](https://img.shields.io/pypi/pyversions/lsdo_project_template)
[![Pypi](https://img.shields.io/pypi/v/lsdo_project_template)](https://pypi.org/project/lsdo_project_template/)
[![Coveralls Badge][13]][14]
[![PyPI version][10]][11]
[![PyPI Monthly Downloads][12]][11]
-->

[![GitHub Actions Test Badge](https://github.com/suspensemax/LSDO_Motor/actions/workflows/actions.yml/badge.svg)](https://github.com/LSDO_Motor/LSDO_Motor/actions)
[![Forks](https://img.shields.io/github/forks/suspensemax/LSDO_Motor.svg)](https://github.com/suspensemax/LSDO_Motor/network)
[![Issues](https://img.shields.io/github/issues/suspensemax/LSDO_Motor.svg)](https://github.com/suspensemax/LSDO_Motor/issues)

This package provides solvers for low-fidelity PMSMs(Permanent Magnet Synchronous Motors) in the context of NASA's Technical Challenge 1, under the NASA Undergraduate Leadership initiative. These solvers performs optimization using [CSDL](!https://lsdolab.github.io/csdl/). 

The Low-Fidelity Motor Model uses two sub-models to analyze the motors in the context of ULI: Sizing and Analysis. This takes into account the geometry of the motors as well as the Inductance and other properties of a motor, and using a Newtonion Solver will optimize the given control methods.  

This repository contains low-fidelity examples, documentation, and comprehensive tests for ease of use as well as scalable development. 

# Installation

## Installation instructions for users
For direct installation with all dependencies, run on the terminal or command line
```sh
pip install git+https://github.com/suspensemax/lsdo_project_template.git
```
Once PyPy is set up. 
For right now: 
```sh
git clone https://github.com/suspensemax/LSDO_Motor.git
```

If you want users to install a specific branch, run
```sh
pip install git+https://github.com/suspensemax/LSDO_Motor.git@branch
```
As of right now, to install a certain branch: 
```sh
git clone https://github.com/suspensemax/LSDO_Motor.git
git checkout {branch}
```
to utilize a certain branch

**Enabled by**: `packages=find_packages()` in the `setup.py` file.

## Installation instructions for developers
To install `LSDO_Motor`,clone the repository using *git*.
On the terminal or command line, run
```sh
git clone https://github.com/suspensemax/LSDO_Motor.git
# pip install -e ./LSDO_Motor
```

# For Developers
For details on documentation, refer to the README in `docs` directory.

For details on testing/pull requests, refer to the README in `tests` directory.

# License
This project is licensed under the terms of the **MIT License**

Copyright (c) [2023] [Luca Scotzniosky]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
