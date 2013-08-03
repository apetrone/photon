#### // photon

license: (BSD 2-Clause)

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

The resulting binaries should be in the root project folder. Some generators don't do this - they place it in a Debug/Release subfolder.

These are command line applications -- invoke them from their project root. i.e. Make sure the current working directory is the project's own root folder. "./raster" or "./Debug/raster", if it's in a subfolder.
