workspace "Fluid Simulation"
    location "generated"
    language "C++"
    architecture "x86_64"

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

    targetdir ("build/bin/%{prj.name}/%{cfg.longname}")
    objdir ("build/obj/%{prj.name}/%{cfg.longname}")

project "Fluid"
    kind "ConsoleApp"

    includedirs
    {
        "src/"
    }

    files "src/**.cpp"
    files "src/**.hpp"
