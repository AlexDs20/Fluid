workspace "Fluid Simulation"
    location "generated"
    language "C++"
    architecture "x86_64"
    toolset "clang"

    configurations { "debug", "release", "dev" }

    filter { "configurations:debug" }
        symbols "On"

    filter { "configurations:release" }
        optimize "On"

    filter { "configurations:dev" }
        optimize "On"
        defines{ "_OMP" }
        links{ "omp" }

    filter "system:linux"
        defines{ "_X11" }

    filter { }

    targetdir ("build/bin/%{prj.name}/%{cfg.longname}")
    objdir ("build/obj/%{prj.name}/%{cfg.longname}")

include "deps/glad.lua"
include "deps/glfw.lua"
include "deps/glm.lua"
include "deps/stb.lua"


project "Fluid"
    kind "WindowedApp"
    openmp "On"

    includedirs
    {
        "src/",
        "deps/stb/"
    }

    files {
        "src/*.cpp",
        "src/*.hpp",
        "src/Math/**",
        "src/Physics/**",
        "src/Renderer/**",
    }

    links { "GLAD", "GLFW", "STB", "GLM" }
