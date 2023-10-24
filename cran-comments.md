# Version 1.0.0

## Test Environments

* Local:
  * Windows 11, R 4.3.1 (x86_64-w64-mingw32)
  * r-hub/containers/ubuntu-release (R 4.3.1)
  * r-hub/containers/ubuntu-next (R 4.3.1 Patched r85316)
  * r-hub/containers/ubuntu-gcc12 (R-Devel r85316) -- additional packages: libdeflate-dev
  * r-hub/containers/ubuntu-clang (R-Devel r85316) -- additional packages: libdeflate-dev, libstdc++-12-dev
  * r-hub/containers/atlas (Fedora 38, R-Devel r85316)
  * r-hub/containers/gcc13 (Fedora 38, R-Devel r85316)
  * r-hub/containers/intel (Fedora 38, R-Devel r85316)
  * r-hub/containers/mkl (Fedora 38, R-Devel r85316)
  * r-hub/containers/nold  (Ubuntu 22.04, R-Devel r85316) -- additional packages: libdeflate-dev
  * r-hub/containers/clang-asan (Ubuntu 22.04, R-Devel r85316) -- additional packages: libdeflate-dev
* win-builder:
  * Windows Server 2022, Devel
* R-hub Builder (<https://builder.r-hub.io>)
  * Windows Server 2022 R-Devel
  * Ubuntu Linux 20.04.1 LTS, R-devel with rchk


## R CMD check results

``` <!-- language: lang-none -->

Found the following significant warnings:
  D:/temp/Rtmpk1MpYi/R.INSTALL2f5205985409f/RcppPlanc/src/_deps/hwloc-src/hwloc/topology-windows.c:1282:46: warning: format '%u' expects argument of type 'unsigned int', but argument 4 has type 'DWORD' {aka 'long unsigned int'} [-Wformat=]
  D:/temp/Rtmpk1MpYi/R.INSTALL2f5205985409f/RcppPlanc/src/_deps/hwloc-src/hwloc/topology-windows.c:1282:49: warning: format '%u' expects argument of type 'unsigned int', but argument 5 has type 'DWORD' {aka 'long unsigned int'} [-Wformat=]
  D:/temp/Rtmpk1MpYi/R.INSTALL2f5205985409f/RcppPlanc/src/_deps/hwloc-src/hwloc/topology-windows.c:1282:52: warning: format '%u' expects argument of type 'unsigned int', but argument 6 has type 'DWORD' {aka 'long unsigned int'} [-Wformat=]
  D:/temp/Rtmpk1MpYi/R.INSTALL2f5205985409f/RcppPlanc/src/_deps/hwloc-src/hwloc/topology-windows.c:1286:29: warning: passing argument 2 of 'GetComputerNameA' from incompatible pointer type [-Wincompatible-pointer-types]
```

This is from one of our linked dependencies. If need be, they can be patched out, but these code paths are never called.
This warning is not present on UNIX-alikes.

``` <!-- language: lang-none -->

/root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:88:10: warning: explicit comparison with infinity in fast floating point mode [-Wtautological-constant-compare]
  /root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:98:10: warning: explicit comparison with infinity in fast floating point mode [-Wtautological-constant-compare]
  /root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:134:10: warning: explicit comparison with NaN in fast floating point mode [-Wtautological-constant-compare]
  /root/R/x86_64-pc-linux-gnu-library/4.4/RcppArmadillo/include/armadillo_bits/arma_cmath.hpp:144:10: warning: explicit comparison with NaN in fast floating point mode [-Wtautological-constant-compare]
See ‘/check/RcppPlanc.Rcheck/00install.out’ for details.
```

This warning is exclusive to the intel compiler and is an issue with RcppArmadillo.

``` <!-- language: lang-none -->

❯ Found the following directories with names of version control directories:
    ./src/_deps/highfive-src/.git
```

  These are not included in the tarball, they are downloaded at configure time from version control tags.

``` <!-- language: lang-none -->
* checking R files for syntax errors ... WARNING
Warning in Sys.setlocale("LC_CTYPE", "en_US.UTF-8") :
  OS reports request to set locale to "en_US.UTF-8" cannot be honored
```

The R-Hub ubuntu containers appear to not have locales installed. At least, some of them.

``` <!-- language: lang-none -->
* checking data for non-ASCII characters ... [1s/0s] WARNING
  Error loading dataset 'ctrl.sparse':
   Error in .requirePackage(package) :
    unable to find required package 'Matrix'

  Error loading dataset 'stim.sparse':
   Error in .requirePackage(package) :
    unable to find required package 'Matrix'

  The dataset(s) may use package(s) not declared in Depends/Imports.
  ```

Matrix is declared in Suggests. This error does not appear to occur on R-devel.

``` <!-- language: lang-none -->
❯ checking line endings in Makefiles ... WARNING
  Found the following Makefile(s) with CR or CRLF line endings:
    src/Makefile
    src/_deps/highfive-src/src/benchmarks/Makefile
    src/_deps/highfive-subbuild/Makefile
    src/_deps/hwloc-build/Makefile
    src/_deps/hwloc-subbuild/Makefile
  Some Unix 'make' programs require LF line endings.
```

  These are not included in the tarball. They are generated at configure time.

``` <!-- language: lang-none -->
❯ checking for GNU extensions in Makefiles ... WARNING
  Found the following file(s) containing GNU extensions:
    src/Makefile
    src/_deps/highfive-src/src/benchmarks/Makefile
    src/_deps/highfive-subbuild/Makefile
    src/_deps/hwloc-build/Makefile
    src/_deps/hwloc-subbuild/Makefile
  Portable Makefiles do not use GNU extensions such as +=, :=, $(shell),
  $(wildcard), ifeq ... endif, .NOTPARALLEL See section 'Writing portable
  packages' in the 'Writing R Extensions' manual.
```

These are generated at configure time by CMake. I imagine that it does not use GNU extensions when run on a system that does not support them.

``` <!-- language: lang-none -->
❯ checking pragmas in C/C++ headers and code ... WARNING
  Files which contain non-portable pragma(s)
    'src/_deps/highfive-src/deps/catch2/extras/catch_amalgamated.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reusable_string_stream.hpp'
    'src/_deps/hwloc-src/hwloc/topology-windows.c'
  Files which contain pragma(s) suppressing diagnostics:
    'src/_deps/highfive-src/deps/catch2/examples/231-Cfg-OutputStreams.cpp'
    'src/_deps/highfive-src/deps/catch2/extras/catch_amalgamated.cpp'
    'src/_deps/highfive-src/deps/catch2/extras/catch_amalgamated.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/catch_template_test_macros.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/catch_test_case_info.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/catch_test_spec.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reusable_string_stream.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_template_test_registry.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_macro_impl.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_registry.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_spec_parser.hpp'
    'src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_floating_point.cpp'
    'src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X02-DisabledMacros.cpp'
    'src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Compilation.tests.cpp'
    'src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Condition.tests.cpp'
    'src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Exception.tests.cpp'
    'src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Generators.tests.cpp'
    'src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Message.tests.cpp'
    'src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Tricky.tests.cpp'
    'src/_deps/hwloc-src/hwloc/topology-windows.c'
    'src/_deps/hwloc-src/utils/lstopo/lstopo-draw.c'
    'src/_deps/hwloc-src/utils/lstopo/lstopo-windows.c'
```

These are vendored source files. If need be, we can patch them.


  ``` <!-- language: lang-none -->

❯ checking installed package size ... NOTE
    installed size is 33.0Mb
    sub-directories of 1Mb or more:
      libs  30.9Mb
```

Static linkage of quite a few libraries will do this. Reported size is from my Windows machine. Unix-alike installs
use dynamic linkage and embedded hwloc and will be smaller. Linux builds are typically ~5MB even with statically linked
hwloc.

``` <!-- language: lang-none -->

❯ checking top-level files ... NOTE
  Non-standard files/directories found at top level:
    'CMakeLists.txt' 'cmake'
```

This is CMake standard practice. If need be, they can be moved to src or another directory.

``` <!-- language: lang-none -->
❯ checking line endings in C/C++/Fortran sources/headers ... NOTE
  Found the following sources/headers with CR or CRLF line endings:
    src/CMakeFiles/3.24.3/CompilerIdC/CMakeCCompilerId.c
    src/CMakeFiles/3.24.3/CompilerIdCXX/CMakeCXXCompilerId.cpp
    src/CMakeFiles/FindOpenMP/OpenMPCheckVersion.c
    src/CMakeFiles/FindOpenMP/OpenMPCheckVersion.cpp
    src/CMakeFiles/FindOpenMP/OpenMPTryFlag.c
    src/CMakeFiles/FindOpenMP/OpenMPTryFlag.cpp
    src/_deps/highfive-src/deps/catch2/examples/010-TestCase.cpp
    src/_deps/highfive-src/deps/catch2/examples/020-TestCase-1.cpp
    src/_deps/highfive-src/deps/catch2/examples/020-TestCase-2.cpp
    src/_deps/highfive-src/deps/catch2/examples/030-Asn-Require-Check.cpp
    src/_deps/highfive-src/deps/catch2/examples/100-Fix-Section.cpp
    src/_deps/highfive-src/deps/catch2/examples/110-Fix-ClassFixture.cpp
    src/_deps/highfive-src/deps/catch2/examples/120-Bdd-ScenarioGivenWhenThen.cpp
    src/_deps/highfive-src/deps/catch2/examples/210-Evt-EventListeners.cpp
    src/_deps/highfive-src/deps/catch2/examples/231-Cfg-OutputStreams.cpp
    src/_deps/highfive-src/deps/catch2/examples/300-Gen-OwnGenerator.cpp
    src/_deps/highfive-src/deps/catch2/examples/301-Gen-MapTypeConversion.cpp
    src/_deps/highfive-src/deps/catch2/examples/302-Gen-Table.cpp
    src/_deps/highfive-src/deps/catch2/examples/310-Gen-VariablesInGenerators.cpp
    src/_deps/highfive-src/deps/catch2/examples/311-Gen-CustomCapture.cpp
    src/_deps/highfive-src/deps/catch2/extras/catch_amalgamated.cpp
    src/_deps/highfive-src/deps/catch2/extras/catch_amalgamated.hpp
    src/_deps/highfive-src/deps/catch2/fuzzing/NullOStream.cpp
    src/_deps/highfive-src/deps/catch2/fuzzing/NullOStream.h
    src/_deps/highfive-src/deps/catch2/fuzzing/fuzz_TestSpecParser.cpp
    src/_deps/highfive-src/deps/catch2/fuzzing/fuzz_XmlWriter.cpp
    src/_deps/highfive-src/deps/catch2/fuzzing/fuzz_textflow.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_benchmark.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_benchmark_all.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_chronometer.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_chronometer.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_clock.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_constructor.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_environment.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_estimate.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_execution_plan.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_optimizer.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_outlier_classification.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/catch_sample_analysis.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_analyse.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_benchmark_function.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_benchmark_function.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_complete_invoke.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_estimate_clock.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_measure.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_repeat.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_run_for_at_least.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_run_for_at_least.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_stats.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_stats.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/benchmark/detail/catch_timing.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_all.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_approx.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_approx.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_assertion_info.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_assertion_result.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_assertion_result.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_config.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_config.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_get_random_seed.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_get_random_seed.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_message.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_message.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_registry_hub.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_section_info.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_session.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_session.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_tag_alias.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_tag_alias_autoregistrar.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_tag_alias_autoregistrar.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_template_test_macros.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_test_case_info.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_test_case_info.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_test_macros.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_test_spec.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_test_spec.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_timer.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_timer.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_tostring.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_tostring.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_totals.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_totals.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_translate_exception.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_version.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_version.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/catch_version_macros.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generator_exception.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generator_exception.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generators.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generators.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generators_adapters.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generators_all.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generators_random.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generators_random.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/generators/catch_generators_range.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_all.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_capture.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_capture.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_config.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_config.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_enum_values_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_exception.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_exception.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_generatortracker.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_generatortracker.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_registry_hub.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_registry_hub.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_reporter.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_reporter.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_reporter_factory.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_reporter_factory.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_reporter_registry.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_reporter_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_tag_alias_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_testcase.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/interfaces/catch_interfaces_testcase.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_assertion_handler.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_assertion_handler.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_case_insensitive_comparisons.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_case_insensitive_comparisons.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_case_sensitive.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_clara.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_clara.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_commandline.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_commandline.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_compare_traits.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_compiler_capabilities.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_config_android_logwrite.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_config_counter.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_config_uncaught_exceptions.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_config_wchar.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_console_colour.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_console_colour.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_console_width.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_container_nonmembers.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_context.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_context.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_debug_console.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_debug_console.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_debugger.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_debugger.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_decomposer.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_decomposer.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_enforce.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_enforce.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_enum_values_registry.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_enum_values_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_errno_guard.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_errno_guard.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_exception_translator_registry.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_exception_translator_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_fatal_condition_handler.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_fatal_condition_handler.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_floating_point_helpers.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_floating_point_helpers.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_getenv.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_getenv.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_istream.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_istream.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_lazy_expr.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_lazy_expr.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_leak_detector.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_leak_detector.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_list.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_list.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_logical_traits.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_main.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_message_info.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_message_info.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_meta.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_move_and_forward.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_noncopyable.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_optional.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_output_redirect.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_output_redirect.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_parse_numbers.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_parse_numbers.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_platform.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_polyfills.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_polyfills.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_preprocessor.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_preprocessor_remove_parens.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_random_number_generator.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_random_number_generator.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_random_seed_generation.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_random_seed_generation.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reporter_registry.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reporter_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reporter_spec_parser.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reporter_spec_parser.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_result_type.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_result_type.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reusable_string_stream.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_reusable_string_stream.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_run_context.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_run_context.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_section.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_section.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_sharding.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_singletons.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_singletons.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_source_line_info.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_source_line_info.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_startup_exception_registry.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_startup_exception_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_stdstreams.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_stdstreams.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_stream_end_stop.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_string_manip.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_string_manip.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_stringref.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_stringref.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_tag_alias_registry.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_tag_alias_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_template_test_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_case_info_hasher.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_case_info_hasher.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_case_registry_impl.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_case_registry_impl.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_case_tracker.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_case_tracker.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_failure_exception.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_failure_exception.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_macro_impl.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_registry.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_registry.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_spec_parser.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_test_spec_parser.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_textflow.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_textflow.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_to_string.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_uncaught_exceptions.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_uncaught_exceptions.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_unique_name.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_unique_ptr.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_void_type.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_wildcard_pattern.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_wildcard_pattern.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_windows_h_proxy.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_xmlwriter.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/internal/catch_xmlwriter.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_all.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_container_properties.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_container_properties.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_contains.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_exception.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_exception.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_floating_point.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_floating_point.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_predicate.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_predicate.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_quantifiers.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_quantifiers.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_string.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_string.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_templated.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_templated.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/catch_matchers_vector.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/internal/catch_matchers_impl.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/matchers/internal/catch_matchers_impl.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_automake.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_automake.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_common_base.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_common_base.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_compact.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_compact.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_console.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_console.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_cumulative_base.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_cumulative_base.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_event_listener.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_event_listener.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_helpers.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_helpers.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_junit.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_junit.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_multi.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_multi.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_registrars.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_registrars.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_sonarqube.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_sonarqube.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_streaming_base.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_streaming_base.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_tap.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_tap.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_teamcity.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_teamcity.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_xml.cpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporter_xml.hpp
    src/_deps/highfive-src/deps/catch2/src/catch2/reporters/catch_reporters_all.hpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X01-PrefixedMacros.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X02-DisabledMacros.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X03-DisabledExceptions-DefaultHandler.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X04-DisabledExceptions-CustomHandler.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X05-DeferredStaticChecks.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X10-FallbackStringifier.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X11-DisableStringification.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X12-CustomDebugBreakMacro.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X21-PartialTestCaseEvents.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X22-BenchmarksInCumulativeReporter.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X23-CasingInReporterNames.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X24-ListenerStdoutCaptureInMultireporter.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X25-ListenerCanAskForCapturedStdout.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X26-ReporterPreferencesForPassingAssertionsIsRespected.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X27-CapturedStdoutInTestCaseEvents.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X28-ListenersGetEventsBeforeReporters.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X29-CustomArgumentsForReporters.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X30-BazelReporter.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X31-DuplicatedTestCases.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X32-DuplicatedTestCasesDifferentTags.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X33-DuplicatedTestCaseMethods.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X34-DuplicatedTestCaseMethodsDifferentFixtures.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X35-DuplicatedReporterNames.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X90-WindowsHeaderInclusion.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X91-AmalgamatedCatch.cpp
    src/_deps/highfive-src/deps/catch2/tests/ExtraTests/X92-NoTests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Clara.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/CmdLine.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/CmdLineHelpers.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/ColourImpl.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Details.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/FloatingPoint.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/GeneratorsImpl.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/InternalBenchmark.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Parse.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/PartTracker.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/RandomNumberGeneration.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Reporters.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Sharding.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Stream.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/String.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/StringManip.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Tag.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/TestCaseInfoHasher.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/TestSpec.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/TestSpecParser.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/TextFlow.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/ToString.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Traits.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/UniquePtr.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/IntrospectiveTests/Xml.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/TestRegistrations.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/TimingTests/Sleep.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Approx.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/BDD.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Benchmark.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Class.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Compilation.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Condition.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Decomposition.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/EnumToString.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Exception.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Generators.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Matchers.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/MatchersRanges.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Message.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Misc.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringByte.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringChrono.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringGeneral.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringOptional.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringPair.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringTuple.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringVariant.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringVector.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/ToStringWhich.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/Tricky.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/UsageTests/VariadicMacros.tests.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/helpers/parse_test_spec.cpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/helpers/parse_test_spec.hpp
    src/_deps/highfive-src/deps/catch2/tests/SelfTest/helpers/type_with_lit_0_comparisons.hpp
    src/_deps/highfive-src/deps/catch2/third_party/clara.hpp
    src/_deps/highfive-src/deps/catch2/tools/misc/coverage-helper.cpp
    src/_deps/highfive-src/doc/poster/example1_hdf5.cpp
    src/_deps/highfive-src/doc/poster/example1_highfive.cpp
    src/_deps/highfive-src/doc/poster/example3.cpp
    src/_deps/highfive-src/doc/poster/example6.cpp
    src/_deps/highfive-src/doc/poster/example_boost.cpp
    src/_deps/highfive-src/doc/poster/example_boost_ublas.cpp
    src/_deps/highfive-src/doc/poster/example_easy_highfive.cpp
    src/_deps/highfive-src/doc/poster/example_eigen.cpp
    src/_deps/highfive-src/doc/poster/example_props.cpp
    src/_deps/highfive-src/include/highfive/H5Attribute.hpp
    src/_deps/highfive-src/include/highfive/H5DataSet.hpp
    src/_deps/highfive-src/include/highfive/H5DataSpace.hpp
    src/_deps/highfive-src/include/highfive/H5DataType.hpp
    src/_deps/highfive-src/include/highfive/H5Easy.hpp
    src/_deps/highfive-src/include/highfive/H5Exception.hpp
    src/_deps/highfive-src/include/highfive/H5File.hpp
    src/_deps/highfive-src/include/highfive/H5FileDriver.hpp
    src/_deps/highfive-src/include/highfive/H5Group.hpp
    src/_deps/highfive-src/include/highfive/H5Object.hpp
    src/_deps/highfive-src/include/highfive/H5PropertyList.hpp
    src/_deps/highfive-src/include/highfive/H5Reference.hpp
    src/_deps/highfive-src/include/highfive/H5Selection.hpp
    src/_deps/highfive-src/include/highfive/H5Utility.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Annotate_traits.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Annotate_traits_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Attribute_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Converter_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5DataSet_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5DataType_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Dataspace_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Exception_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5FileDriver_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5File_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Friends.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Iterables_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Node_traits.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Node_traits_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Object_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Path_traits.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Path_traits_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5PropertyList_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5ReadWrite_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Reference_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Selection_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Slice_traits.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Slice_traits_misc.hpp
    src/_deps/highfive-src/include/highfive/bits/H5Utils.hpp
    src/_deps/highfive-src/include/highfive/bits/H5_definitions.hpp
    src/_deps/highfive-src/include/highfive/h5easy_bits/H5Easy_Eigen.hpp
    src/_deps/highfive-src/include/highfive/h5easy_bits/H5Easy_misc.hpp
    src/_deps/highfive-src/include/highfive/h5easy_bits/H5Easy_opencv.hpp
    src/_deps/highfive-src/include/highfive/h5easy_bits/H5Easy_public.hpp
    src/_deps/highfive-src/include/highfive/h5easy_bits/H5Easy_scalar.hpp
    src/_deps/highfive-src/include/highfive/h5easy_bits/H5Easy_vector.hpp
    src/_deps/highfive-src/include/highfive/h5easy_bits/H5Easy_xtensor.hpp
    src/_deps/highfive-src/src/benchmarks/hdf5_bench.cpp
    src/_deps/highfive-src/src/benchmarks/hdf5_bench_improved.cpp
    src/_deps/highfive-src/src/benchmarks/highfive_bench.cpp
    src/_deps/highfive-src/src/examples/boost_multi_array_2D.cpp
    src/_deps/highfive-src/src/examples/boost_multiarray_complex.cpp
    src/_deps/highfive-src/src/examples/boost_ublas_double.cpp
    src/_deps/highfive-src/src/examples/compound_types.cpp
    src/_deps/highfive-src/src/examples/create_attribute_string_integer.cpp
    src/_deps/highfive-src/src/examples/create_dataset_double.cpp
    src/_deps/highfive-src/src/examples/create_dataset_half_float.cpp
    src/_deps/highfive-src/src/examples/create_datatype.cpp
    src/_deps/highfive-src/src/examples/create_extensible_dataset.cpp
    src/_deps/highfive-src/src/examples/create_page_allocated_files.cpp
    src/_deps/highfive-src/src/examples/easy_attribute.cpp
    src/_deps/highfive-src/src/examples/easy_dumpoptions.cpp
    src/_deps/highfive-src/src/examples/easy_load_dump.cpp
    src/_deps/highfive-src/src/examples/eigen_matrix.cpp
    src/_deps/highfive-src/src/examples/hl_hdf5_inmemory_files.cpp
    src/_deps/highfive-src/src/examples/parallel_hdf5_collective_io.cpp
    src/_deps/highfive-src/src/examples/parallel_hdf5_independent_io.cpp
    src/_deps/highfive-src/src/examples/read_write_dataset_string.cpp
    src/_deps/highfive-src/src/examples/read_write_fixedlen_string.cpp
    src/_deps/highfive-src/src/examples/read_write_raw_ptr.cpp
    src/_deps/highfive-src/src/examples/read_write_single_scalar.cpp
    src/_deps/highfive-src/src/examples/read_write_vector_dataset.cpp
    src/_deps/highfive-src/src/examples/read_write_vector_dataset_references.cpp
    src/_deps/highfive-src/src/examples/readme_snippet.cpp
    src/_deps/highfive-src/src/examples/renaming_objects.cpp
    src/_deps/highfive-src/src/examples/select_by_id_dataset_cpp11.cpp
    src/_deps/highfive-src/src/examples/select_partial_dataset_cpp11.cpp
    src/_deps/highfive-src/tests/test_dependent_library/include/simpleton.hpp
    src/_deps/highfive-src/tests/test_dependent_library/src/otherton.cpp
    src/_deps/highfive-src/tests/test_dependent_library/src/simpleton.cpp
    src/_deps/highfive-src/tests/unit/test_all_types.cpp
    src/_deps/highfive-src/tests/unit/tests_high_five.hpp
    src/_deps/highfive-src/tests/unit/tests_high_five_base.cpp
    src/_deps/highfive-src/tests/unit/tests_high_five_easy.cpp
    src/_deps/highfive-src/tests/unit/tests_high_five_multi_dims.cpp
    src/_deps/highfive-src/tests/unit/tests_high_five_parallel.cpp
    src/_deps/highfive-src/tests/unit/tests_import_public_headers.cpp
  Some Unix compilers require LF line endings.
```

  These are not included in the tarball, they are downloaded at configure time by git. On UNIX-alikes this note is limited to
  `src/_deps/hwloc-src/contrib/android/AndroidApp/lstopo/src/main/cpp/lib.c`

``` <!-- language: lang-none -->
❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    '-freciprocal-math' '-funsafe-math-optimizations'
```

These compilation flags are checked at configure time and will not be used if they are not available.

``` <!-- language: lang-none -->
❯ checking compiled code ... NOTE
  Note: information on .o files for x64 is not available
  File 'C:/Users/andrew/AppData/Local/Temp/RtmpMLm9gJ/filec79641d0cf49/RcppPlanc.Rcheck/RcppPlanc/libs/x64/RcppPlanc.dll':
    Found '_assert', possibly from 'assert' (C)
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'rand', possibly from 'rand' (C)
    Found 'srand', possibly from 'srand' (C)

  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs nor [v]sprintf. The detected symbols are linked into
  the code but might come from libraries and not actually be called.

  See 'Writing portable packages' in the 'Writing R Extensions' manual.

  ```

This is likely from the exception handling routines of included libraries.
If need be, they can be patched out, but there are enough guards on our R
wrappers that they will never be called.

  ``` <!-- language: lang-none -->
❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

```

Per <https://github.com/r-hub/rhub/issues/560> this appears to be a bug.

  ``` <!-- language: lang-none -->
❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException''

```

Per <https://github.com/r-hub/rhub/issues/503> this appears to be a bug.

  ``` <!-- language: lang-none -->
❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable

```

This is an Rhub issue. It works locally.

``` <!-- language: lang-none -->
❯ checking for future file timestamps ... NOTE
  unable to verify current time

```

Per <https://github.com/r-hub/rhub/issues/440> this appears to be a bug.

``` <!-- language: lang-none -->
❯ checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
         user system elapsed
  inmf 15.049  0.115   1.532

```

Will vary by platform. This was recorded before the example arguments were lowered substantially.

``` <!-- language: lang-none -->
* checking tests ...
  Running 'testthat.R' [3344s/219s]
Running R code in 'testthat.R' had CPU time 15.3 times elapsed time
* checking re-building of vignette outputs ... [110s/12s] NOTE
Re-building vignettes had CPU time 9.4 times elapsed time
 [3346s/220s] NOTE

```

Only occurs in ubuntu-based docker containers. Likely the default openblas-pthread conflicting with openmp.

```c++
  Direct leak of 16 byte(s) in 1 object(s) allocated from:
#0 0x555ffa2c6318 in __interceptor_calloc (/opt/R/devel-asan/lib/R/bin/exec/R+0xb8318) (BuildId: f9ace1cd1e557acb45d82ef6d66e646229ab6542)
#1 0x7f81f807ea0f  (/usr/lib/x86_64-linux-gnu/libc++abi.so.1+0x25a0f) (BuildId: 47782717d7fe8d90800ba4d11684a17cfa249698)
#2 0x7f815f03f2a4 in planc::INMF<arma::Mat<double>>::INMF(std::__1::vector<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>, std::__1::allocator<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>>>&, unsigned int, double, bool) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/common/inmf.hpp:125:19
#3 0x7f815f035237 in planc::BPPINMF<arma::Mat<double>>::BPPINMF(std::__1::vector<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>, std::__1::allocator<std::__1::unique_ptr<arma::Mat<double>, std::__1::default_delete<arma::Mat<double>>>>>&, unsigned int, double) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/nmf/bppinmf.hpp:137:82
#4 0x7f815f035237 in Rcpp::Vector<19, Rcpp::PreserveStorage> runINMF<arma::Mat<double>>(std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>, unsigned int, double, unsigned int, bool) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/rcppplanc_nmf.cpp:572:23
#5 0x7f815ed901f5 in bppinmf_dense(std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>> const&, unsigned int, double, unsigned int, bool, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<arma::Mat<double>>) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/rcppplanc_nmf.cpp:632:16
#6 0x7f815eda1d36 in bppinmf(Rcpp::Vector<19, Rcpp::PreserveStorage>, unsigned int, double, unsigned int, bool, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double>>>>, Rcpp::Nullable<arma::Mat<double>>) /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/rcppplanc_nmf.cpp:666:16
#7 0x7f815f1cc172 in _RcppPlanc_bppinmf /check/RcppPlanc.Rcheck/00_pkg_src/RcppPlanc/src/RcppExports.cpp:134:34
#8 0x7f81f8283262 in R_doDotCall /tmp/R-devel/src/main/dotcode.c
```
This appears to be a libc++ issue with exception throwing per https://stackoverflow.com/questions/61015622/throw-with-clang-shows-possible-memory-leak.

```<!-- language: lang-none -->
Function Rcpp::Vector<19, Rcpp::PreserveStorage> onlineINMF_S1_mem<arma::Mat<double> >(std::__1::vector<arma::Mat<double>, std::__1::allocator<arma::Mat<double> > >, unsigned int, double, unsigned int, unsigned int, unsigned int, bool)
  [UP] ignoring variable <unnamed var:   %36 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %40 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %44 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %49 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %53 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  [UP] ignoring variable <unnamed var:   %57 = alloca %struct.SEXPREC*, align 8> as it has address taken, results will be incomplete 
  and similar rchk errors
```
These appear to be issues with rchk's handling of Armadillo's template metaprogramming.


0 errors ✔ | 5 warnings ✖ | 7 notes ✖

## revdepcheck results

None as this is a new package.
