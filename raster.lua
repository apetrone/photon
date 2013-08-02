solution "raster"
configurations {"Debug", "Release"}
project "raster"
language "C++"
location "build"
targetdir "bin"


kind "ConsoleApp"

files
{
	"src/raster.cpp"
}

includedirs
{
	"glm",
	"soil/include",
	"soil/src/",
}

libdirs
{
	"soil/lib/"
}

links
{
	"soil",
}