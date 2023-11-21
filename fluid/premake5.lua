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
        symbols "On"

    filter { "configurations:dev" }
        optimize "On"
        defines{ "_OMP" }
        links{ "omp" }
        symbols "On"

    filter "system:linux"
        defines{ "_X11" }

    filter { }

    targetdir ("build/bin/%{prj.name}/%{cfg.longname}")
    objdir ("build/obj/%{prj.name}/%{cfg.longname}")

include "deps/glad.lua"
include "deps/glfw.lua"
include "deps/glm.lua"

project "Fluid"
    kind "WindowedApp"
    openmp "On"

    includedirs
    {
        "src/",
        "src/Console/",
        "src/Input/",
        "src/Math/",
        "src/Message/",
        "src/Physics/",
        "src/Renderer/",
        "deps/stb/",
        "deps/glad/include/",
        "deps/glfw/include/",
        "deps/glm/"
    }

    files {
        "src/*.cpp",
        "src/Math/**",
        "src/Physics/**",
        "src/Renderer/**",
        "src/Message/**",
    }

    buildoptions { "-mavx2" }

    links { "GLAD", "GLFW", "GLM" }
