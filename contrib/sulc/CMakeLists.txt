SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ADD_CUSTOM_TARGET(sulc
	DEPENDS TorsionalDetailedPatchySwapInteraction 
)

# without this bit clang-based compilation fails with some "undefined symbols" errors during the linking stage
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	link_libraries(${ARGV0} "-undefined dynamic_lookup")
endif()

# this removes the prefix from each library. I didn't find any other way of doing it for all 
# libraries at once if the library type is MODULE
function(add_library_no_prefix target source)
	# ${ARGN} holds the list of arguments past the last expected argument, which is "source"
	add_library(${target} MODULE EXCLUDE_FROM_ALL ${source} ${ARGN})
	set_target_properties(${target} PROPERTIES PREFIX "")
endfunction()

# Observables

# Interactions
add_library_no_prefix(TorsionalDetailedPatchySwapInteraction src/Interactions/TorsionalDetailedPatchySwapInteraction.cpp)

# Particles
add_library_no_prefix(TorsionalPatchyParticle src/Particles/TorsionalPatchyParticle.cpp)

IF(CUDA)
	set(CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -rdynamic")
	find_package("CUDA")
	
	# same thing but for CUDA libs
	function(cuda_add_library_no_prefix target source)
		cuda_add_library(${target} MODULE EXCLUDE_FROM_ALL ${source} ${ARGN})
		target_link_libraries(${target} ${CUDA_LIBRARIES})
		set_target_properties(${target} PROPERTIES PREFIX "")
	endfunction()

	cuda_add_library_no_prefix(CUDATorsionalDetailedPatchySwapInteraction src/Interactions/CUDATorsionalDetailedPatchySwapInteraction.cu src/Interactions/TorsionalDetailedPatchySwapInteraction.cpp)
	
	ADD_DEPENDENCIES(sulc CUDATorsionalDetailedPatchySwapInteraction )
ENDIF(CUDA)
