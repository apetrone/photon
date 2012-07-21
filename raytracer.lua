solution "raytracer"
configurations {"Debug", "Release"}
project "raytracer"
language "C++"
location "build"
targetdir "bin"


kind "ConsoleApp"

files
{
	"src/raytracer.cpp",
	"src/common.hpp"
}

includedirs
{
	"glm",
	"soil/include",
	"soil/src/",
	"xwl/include",
	"src/"
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