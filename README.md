# pyzeo

Python wrapper for the [Zeo++ library](http://zeoplusplus.org). Based on the latest released version 0.3.

## Installation

We provide pre-built wheels for several platforms (including Apple Silicon). So in most cases you will want to install the package with `pip`:

```sh
pip install pyzeo
```

### Installing from source

If you wish to install the package directly from the source code, you can do so with:

```sh
git clone https://github.com/nomad-coe/pyzeo
cd pyzeo
pip install .
```

Note that this type of installation will require a separate, platfrom-dependent
compilation step. These are some common problems you may encounter in the compilation step:

- **fatal error: Python.h: No such file or directory**: The package depends on
   C/C++ extensions that are compiled during the setup. For the compilation to
   work you will need to install the *pythonX.X-dev*-package, where X.X is the
   python version you use. E.g. for python 3.9 on Ubuntu this package could be
   installed with:

   .. code-block:: sh

       sudo apt install python3.9-dev

 - **Installation errors on MacOS**: The package depends on C++ extensions that
   are compiled during the setup. If experiencing problems with setup on MacOS,
   you may need to install the Xcode Command Line tools package. This can be
   done with:

   .. code-block:: sh

       xcode-select --install

### Cython

By default the package comes with pre-built Cython binding code. The cython
wrapper definitions live inside src/pyzeo. These bindings can be recreated
by first setting `USE_CYTHON=True` in setup.py, and then recreating the bindings
with:

```sh
python setup.py build_ext --inplace --force
```

Remember to disable cython afterwards by setting `USE_CYTHON=False` in setup.py.

> **Maintainers:** see [scripts/RELEASE.md](scripts/RELEASE.md) for the full release process.

### License

The python wrapper code is licensed under Apache 2.0. [Zeo++
library](http://zeoplusplus.org) and [Voro++
library](https://math.lbl.gov/voro++/) which are included in the source code
have the following license:

#### Zeo++

> Zeo++ Copyright (c) 2011, The Regents of the University
> of California, through Lawrence Berkeley National Laboratory (subject
> to receipt of any required approvals from the U.S. Dept. of Energy).
> All rights reserved.

> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are
> met:

> (1) Redistributions of source code must retain the above copyright
> notice, this list of conditions and the following disclaimer.

> (2) Redistributions in binary form must reproduce the above copyright
> notice, this list of conditions and the following disclaimer in the
> documentation and/or other materials provided with the distribution.

> (3) Neither the name of the University of California, Lawrence
> Berkeley National Laboratory, U.S. Dept. of Energy nor the names of
> its contributors may be used to endorse or promote products derived
> from this software without specific prior written permission.

> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
> "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
> LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
> A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
> OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
> SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
> LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
> DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
> THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
> (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

> You are under no obligation whatsoever to provide any bug fixes,
> patches, or upgrades to the features, functionality or performance of
> the source code ("Enhancements") to anyone; however, if you choose to
> make your Enhancements available either publicly, or directly to
> Lawrence Berkeley National Laboratory, without imposing a separate
> written license agreement for such Enhancements, then you hereby grant
> the following license: a  non-exclusive, royalty-free perpetual
> license to install, use, modify, prepare derivative works, incorporate
> into other computer software, distribute, and sublicense such
> enhancements or derivative works thereof, in binary and source code
> form.

#### Voro++

> Voro++ Copyright (c) 2008, The Regents of the University of California, through
> Lawrence Berkeley National Laboratory (subject to receipt of any required
> approvals from the U.S. Dept. of Energy). All rights reserved.

> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are met: 

> (1) Redistributions of source code must retain the above copyright notice, this
> list of conditions and the following disclaimer. 

> (2) Redistributions in binary form must reproduce the above copyright notice,
> this list of conditions and the following disclaimer in the documentation
> and/or other materials provided with the distribution. 

> (3) Neither the name of the University of California, Lawrence Berkeley
> National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
> be used to endorse or promote products derived from this software without
> specific prior written permission. 

> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
> ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
> WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
> DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
> ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
> (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
> LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
> ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
> (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
> SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

> You are under no obligation whatsoever to provide any bug fixes, patches, or
> upgrades to the features, functionality or performance of the source code
> ("Enhancements") to anyone; however, if you choose to make your Enhancements
> available either publicly, or directly to Lawrence Berkeley National
> Laboratory, without imposing a separate written license agreement for such
> Enhancements, then you hereby grant the following license: a  non-exclusive,
> royalty-free perpetual license to install, use, modify, prepare derivative
> works, incorporate into other computer software, distribute, and sublicense
> such enhancements or derivative works thereof, in binary and source code form.
