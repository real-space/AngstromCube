if(HAS_RAPIDXML)
    include(FetchContent)

    FetchContent_Declare(
        rapidxml
  #     GIT_REPOSITORY https://github.com/discord/rapidxml.git
  #     GIT_TAG        2ae4b2888165a393dfb6382168825fddf00c27b9 # version 1.13
        GIT_REPOSITORY https://github.com/Fe-Bell/RapidXML.git
        GIT_TAG        v116  # version 1.16
    )
    FetchContent_MakeAvailable(rapidxml)
    # after configuring headers for rapidXML are found in _deps/rapidxml-src/RapidXML
endif(HAS_RAPIDXML)


if(HAS_RAPIDJSON)
    include(FetchContent)

    FetchContent_Declare(
        rapidjson
        GIT_REPOSITORY https://github.com/Tencent/rapidjson.git
        GIT_TAG        v1.1.0
    )
    FetchContent_MakeAvailable(rapidjson)
    # after configuring headers for rapidJSON are found in _deps/rapidjson-src
endif(HAS_RAPIDJSON)


if(HAS_TFQMRGPU)
    include(FetchContent)

    FetchContent_Declare(
        tfQMRgpu
        GIT_REPOSITORY https://github.com/real-space/tfQMRgpu.git
        # GIT_TAG        v0.9  # version 0.9
        GIT_TAG "8d56049929db70a8095c7f20bcdad0bfc534a02a"
    )
    FetchContent_MakeAvailable(tfQMRgpu)
    # after configuring headers for tfQMRgpu are found in _deps/tfQMRgpu-src
endif(HAS_TFQMRGPU)
