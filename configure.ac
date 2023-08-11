AC_PREREQ([2.71])
AC_INIT([RcppPlanc],[0.8.0.9000],[robbiand@med.umich.edu])
AC_CONFIG_AUX_DIR([tools])
m4_include([tools/m4/ax_require_defined.m4])
m4_include([tools/m4/ax_openmp.m4])
m4_include([tools/m4/ax_blas.m4])
m4_include([tools/m4/ax_lapack.m4])
m4_include([tools/m4/ax_cxx_compile_stdcxx.m4])
m4_include([tools/m4/ax_append_flag.m4])
m4_include([tools/m4/ax_check_compile_flag.m4])
m4_include([tools/m4/ax_append_compile_flags.m4])
m4_include([tools/m4/ax_lib_hdf5.m4])

# Find the compiler and compiler flags used by R.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS="`"${R_HOME}/bin/R" CMD config CPPFLAGS`"
FC=`"${R_HOME}/bin/R" CMD config FC`
FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`

: ${LDFLAGS=`echo $LDFLAGS`}
: ${CXXFLAGSENV=`echo $CXXFLAGS`}
#get compilers

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX
AC_PROG_CC



#get libraries

AX_OPENMP([CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"])
AX_BLAS([], AC_MSG_ERROR([could not find your BLAS library], 1))
AX_LAPACK([], AC_MSG_ERROR([could not find your LAPACK library], 1))

AX_LIB_HDF5()

if test "$with_hdf5" = "yes"; then
AC_MSG_RESULT([Found hdf5 with version: $HDF5_VERSION])
else
  AC_MSG_ERROR([hdf5 does not seem to be installed on your platform.
  Please install the hdf5 library.
  The required HDF5 library files can be installed as follows:
      - Debian-based (e.g. Debian >= 8.0, Ubuntu >= 15.04): 'sudo apt-get install libhdf5-dev'
        - Old Debian-based (e.g Debian < 8.0, Ubuntu < 15.04): Install from source (see INSTALL)
        - OS X using Homebrew: 'brew install hdf5'
        - RPM-based (e.g Fedora): 'sudo yum install hdf5-devel'])
fi


#check for headers
AC_CHECK_HEADER([flexiblas/cblas.h], [CPPFLAGS="$CPPFLAGS -DHAVE_FLEXIBLAS_CBLAS_H"],[
  AC_CHECK_HEADER([openblas/cblas.h], [CPPFLAGS="$CPPFLAGS -DHAVE_OPENBLAS_CBLAS_H"], [
    AC_CHECK_HEADER([cblas.h])
  ])
])
#Checks for typedefs, structures, and compiler characteristics.
AX_CXX_COMPILE_STDCXX([17], [], mandatory)
AX_APPEND_COMPILE_FLAGS([-funsafe-math-optimizations -ffinite-math-only -fno-trapping-math -fno-math-errno -fno-signed-zeros -freciprocal-math -march=native -mtune=native])
#AC_C_INLINE
#AC_TYPE_INT32_T
#AC_TYPE_INT64_T
#AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
#AC_CHECK_FUNCS([floor pow sqrt])

PKG_CPPFLAGS="-Inmf -Innls -Icommon -DBUILD_SPARSE $HDF5_CPPFLAGS $CXXFLAGSENV"
PKG_CXXFLAGS="\$(SHLIB_OPENMP_CXXFLAGS) $HDF5_CFLAGS"
PKG_LIBS="\$(SHLIB_OPENMP_CXXFLAGS) $LAPACK_LIBS $BLAS_LIBS $FLIBS $HDF5_LDFLAGS $HDF5_LIBS -lhdf5_cpp -lhdf5_hl_cpp"

# Write the flags into the src/Makevars file.
AC_SUBST([PKG_CPPFLAGS], ["${CPPFLAGS} ${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${LDFLAGS} ${PKG_LIBS}"])
AC_SUBST([PKG_CXXFLAGS], ["${X86_FEATURE_CFLAGS} ${PKG_CXXFLAGS}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${PKG_CPPFLAGS}
    cxxflags: ${CXXFLAGS} ${PKG_CXXFLAGS}
    libs:     ${PKG_LIBS}

  --------------------------------------------------
"
