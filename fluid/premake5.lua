workspace "Fluid Simulation"
    location "generated"
    language "C++"
    architecture "x86_64"
    toolset "clang"

    configurations { "debug", "release" }

    filter { "configurations:debug" }
        symbols "On"

    filter { "configurations:release" }
        optimize "On"

    filter "system:linux"
        links{ "dl", "pthread"}

        defines{ "_X11" }

    filter "system:windows"
        defines { "_WINDOWS" }

    filter "system:MAC"
        defines { "_MAC" }

    filter { }
    flags { "ProfileGC", "ShowGC" }

    targetdir ("build/bin/%{prj.name}/%{cfg.longname}")
    objdir ("build/obj/%{prj.name}/%{cfg.longname}")

include "deps/glad.lua"
include "deps/glfw.lua"
include "deps/glm.lua"
include "deps/stb.lua"


project "Fluid"
    kind "WindowedApp"

    includedirs
    {
        "src/",
        "deps/stb/"
    }

    files {
        "src/**.cpp",
        "src/**.hpp"
    }

    links { "GLAD", "GLFW", "STB" }

