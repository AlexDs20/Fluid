project "GLM"
    kind "SharedLib"
    language "C"
    architecture "x86_64"

    includedirs { "glm/" }

    files
    {
        "glm/glm/**"
    }

    filter "system:linux"
        pic "On"

        systemversion "latest"
        staticruntime "On"

        defines
        {
            "_GLM_X11"
        }

    filter "system:windows"
        systemversion "latest"
        staticruntime "On"

        defines
        {
            "_GLM_WIN32",
            "_CRT_SECURE_NO_WARNINGS"
        }

    filter "configurations:debug"
        runtime "debug"
        symbols "on"

    filter "configurations:release"
        runtime "release"
        optimize "on"
