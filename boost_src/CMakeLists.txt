add_library(
    IMSRGBoostZip
    gzip.cpp
    zlib.cpp
)
target_link_libraries(
    IMSRGBoostZip
    PUBLIC
    Boost::boost
    ZLIB::ZLIB
)