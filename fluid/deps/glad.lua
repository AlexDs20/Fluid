project "GLAD"
    kind "SharedLib"
    language "C"
    architecture "x86_64"

    includedirs { "glad/include" }

    files { "glad/src/glad.c" }

    filter "system:linux"
      pic "On"

      systemversion "latest"
      staticruntime "On"

      defines { "_GLAD_X11" }

    filter "configurations:debug"
      runtime "debug"
      symbols "On"

    filter "configurations:release"
      runtime "release"
      optimize "On"
