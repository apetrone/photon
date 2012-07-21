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
	"xwl/include"
}

libdirs
{
	"soil/lib/",
	"xwl/lib/x86/debug"
}

links
{
	"soil",
	"xwl"
}

configuration { "linux" }
	links
	{
		"GL"
	}

configuration { "macosx" }
	defines { "__MACH__" }
	linkoptions {
		"-framework Cocoa",
		"-framework OpenGL"
	}