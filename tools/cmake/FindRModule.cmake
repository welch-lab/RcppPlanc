# FindRModule.cmake: find a specific R module. from MLPACK
#[[ ----
Copyright (c) 2007-2023, mlpack contributors (see https://github.com/mlpack/mlpack/blob/master/COPYRIGHT.txt)
Modifications (c) 2023 Andrew Robbins, Welch Lab, University of Michigan
All rights reserved.

Redistribution and use of mlpack in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

------ ]]
function(find_r_module module)
    string(TOUPPER ${module} module_upper)

    if(NOT R_${module_upper})
        if(ARGC GREATER 1)
            # Not required but we have version constraints.
            set(VERSION_REQ ${ARGV1})
        endif()

        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        execute_process(COMMAND ${RSCRIPT_EXECUTABLE} "-e" "find.package('${module}')"
            RESULT_VARIABLE _${module}_status
            OUTPUT_VARIABLE _${module}_location
            ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

        # Some cleanup in location of R Module.
        string(REGEX MATCHALL "\".*\"" _${module}_location "${_${module}_location}")
        string(REGEX REPLACE "\"" "" _${module}_location "${_${module}_location}")

        if(NOT _${module}_status)
            # Now we have to check the version.
            if(VERSION_REQ)
                execute_process(COMMAND ${RSCRIPT_EXECUTABLE} "-e" "packageVersion('${module}')"
                    RESULT_VARIABLE _version_status
                    OUTPUT_VARIABLE _version_compare
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

                # Different versions of R may enclose the version number in different
                # delimiters.  Sometimes, semicolons show up too.
                string(REGEX MATCHALL "[‘'][0-9._]*[’']" _version_compare "${_version_compare}")
                string(REGEX REPLACE ";" "" _version_compare "${_version_compare}")
                string(REGEX REPLACE "[‘']" "" _version_compare "${_version_compare}")
                string(REGEX REPLACE "[’']" "" _version_compare "${_version_compare}")

                # Compare the version of the package using compareVersion().
                execute_process(COMMAND ${RSCRIPT_EXECUTABLE} "-e"
                    "compareVersion('${_version_compare}', '${VERSION_REQ}')"
                    RESULT_VARIABLE _compareVersion_status
                    OUTPUT_VARIABLE _compareVersion_result
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

                # Extract compareVersion() result i.e. 1 -> Newer, 0 -> Equal and -1 -> Later.
                string(REGEX REPLACE "\\[\\1\\]" "" _compareVersion_result "${_compareVersion_result}")

                if("${_compareVersion_result}" GREATER "-1")
                    set(R_${module_upper}
                        "${_${module}_location} (found suitable version \"${_version_compare}\", minimum required is \"${VERSION_REQ}\")"
                        CACHE STRING "Location of R module ${module}"
                    )
                else()
                    message(WARNING "Unsuitable version of R module ${module} (${VERSION_REQ} or greater required).")
                endif()
            else()
                # No version requirement so we are done.
                set(R_${module_upper} ${_${module}_location} CACHE STRING "Location of R module ${module}")
            endif()
        endif()
    endif()

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(R_${module} DEFAULT_MSG R_${module_upper})
endfunction()
