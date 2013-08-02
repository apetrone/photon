#### // photon


This is an experimental/learning project used to test various graphics-related techniques including but not limited to:
Rasterization, and Raytracing.




#### // building

* clone and get all dependencies

	git clone git://github.com/apetrone/photon.git

	git submodule update --init

* build soil

	cd deps/soil

	cmake -G "Unix Makefiles"
	
	make

* build raster

	cd raster

	cmake -G "Unix Makefiles"

	make