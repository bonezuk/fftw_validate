cmake_minimum_required(VERSION 3.10)

project(FFTWValidate)

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	include_directories(AFTER "/home/bonez/Development/fftw/fftw_install/include" )
	add_library(fftw STATIC IMPORTED)
	set_property(TARGET fftw PROPERTY IMPORTED_LOCATION "/home/bonez/Development/fftw/fftw_install/lib/libfftw3.a" )

#	find_package(PkgConfig REQUIRED)
#	pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
#	add_library(fftw ALIAS PkgConfig::FFTW)
else (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	include_directories(AFTER "D:/Development/DSD/fftw-test/fftw-3.3.10/api" )
	add_library(fftw SHARED IMPORTED)
	set_property(TARGET fftw PROPERTY IMPORTED_LOCATION "D:/Development/DSD/fftw-test/fftw-3.3.10/build_2022/Release/fftw3.dll" )
	set_property(TARGET fftw PROPERTY IMPORTED_IMPLIB "D:/Development/DSD/fftw-test/fftw-3.3.10/build_2022/Release/fftw3.lib" )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")


# Add the executable
add_executable(fftw_validate fftw_validate.cxx)

# Include directories and link libraries
target_include_directories(fftw_validate PRIVATE fftw)
target_link_libraries(fftw_validate PRIVATE fftw)
