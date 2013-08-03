#### // photon


This is an experimental/learning project used to test various graphics-related techniques including but not limited to:
Rasterization, and Raytracing.




#### // building


* Pretty simple if you're familiar with CMake - I use it for included projects and even the dependencies use it.

1. clone and get all dependencies

	git clone git://github.com/apetrone/photon.git

	git submodule update --init

2. build soil

	cd deps/soil

	cmake -G <generator>

3. Build whichever project you want. Each are organized into their own subfolder. Here's an example for "raster" with CMake. After this you can use "make" or open the generated project in your favorite IDE.

	cd raster

	cmake -G <generator>

The resulting binaries should be in the root project folder.
These are command line applications -- invoke them from their project root.

#### // license (BSD 2-Clause)

Copyright (c) 2013, Adam Petrone
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.