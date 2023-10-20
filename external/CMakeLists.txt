if(HAS_RAPIDXML)
    include(FetchContent)

    FetchContent_Declare(
        rapidxml
        GIT_REPOSITORY https://github.com/discord/rapidxml.git
        GIT_TAG        2ae4b2888165a393dfb6382168825fddf00c27b9 # version 1.13
    )
    FetchContent_MakeAvailable(rapidxml)
    # after configuring headers for rapidXML are found in _deps/rapidxml-src
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


