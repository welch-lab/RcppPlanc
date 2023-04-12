AC_PREREQ([2.71])
AC_INIT([RcppPlanc],[0.8.0.9000],[robbiand@med.umich.edu])
AC_CONFIG_AUX_DIR([tools])
m4_include([tools/m4/ax_openmp.m4])
m4_include([tools/m4/ax_blas.m4])
m4_include([tools/m4/ax_lapack.m4])
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

#get compilers

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX
AC_PROG_CC


#get libraries

AX_OPENMP([CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"])
AX_BLAS([], AC_MSG_ERROR([could not find your BLAS library], 1))
AX_LAPACK([], AC_MSG_ERROR([could not find your LAPACK library], 1))



#check for headers
#Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_INT32_T
#AC_TYPE_INT64_T
#AC_C_RESTRICT
#AC_TYPE_SIZE_T
#AC_TYPE_UINT64_T

# Checks for library functions.
#AC_CHECK_FUNCS([floor pow sqrt])

PKG_CPPFLAGS="-Inmf -Innls -Icommon -DBUILD_SPARSE"
PKG_CXXFLAGS="\$(SHLIB_OPENMP_CXXFLAGS)"
PKG_LIBS="\$(SHLIB_OPENMP_CXXFLAGS) $LAPACK_LIBS $BLAS_LIBS $FLIBS"

# Write the flags into the src/Makevars file.
AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS}"])
AC_SUBST([PKG_CXXFLAGS], ["${PKG_CXXFLAGS}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${PKG_CPPFLAGS}
    cxxflags: ${CXXFLAGS} ${PKG_CXXFLAGS}
    libs:     ${PKG_LIBS}

  --------------------------------------------------
"
