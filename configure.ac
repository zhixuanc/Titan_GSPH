#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([particler], [0.5.1], [dkumar@buffalo.edu])
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/constants.h])
AC_CONFIG_HEADER([config.h])
AC_PREFIX_DEFAULT([${PWD}])
AC_LANG_CPLUSPLUS

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_RANLIB

#-----------------------------------------------------------------------
#  Check for 2-D or 3-d version
#-----------------------------------------------------------------------
AC_MSG_CHECKING([build 3-D code])
AC_ARG_ENABLE(3d,
              AS_HELP_STRING([--disable-3d],
                             [build 2-d code, instead of default 3-d]))
              
if test "x$enable_3d" != "xno"; then
  AC_DEFINE([THREE_D],[],[build 3-d code])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

#-----------------------------------------------------------------------
#  Check for multi process build
#-----------------------------------------------------------------------
AC_MSG_CHECKING([build parallel code])
AC_ARG_ENABLE(parallel,
              AS_HELP_STRING([--enable-parallel],
                             [builds parallel code, requires some MPI installtion]))
if test "x$enable_parallel" != "xno"; then
  AC_DEFINE([MULTI_PROC],[],[build parallel code])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

# parallel-IO
AC_ARG_ENABLE([parallel-IO],
              AS_HELP_STRING([--enable-parallel-IO],
                       [enables single output file by all procs. Needs hardware support]))
AC_MSG_CHECKING([enable parallel-IO])
if test "x$enable_parallel_IO" = "xyes" && test "x$enable_parallel" != "xno"; then
  AC_DEFINE([PARALLEL_IO], [], [support parallel-IO])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

#-----------------------------------------------------------------------
#  Check for mpi installation
#-----------------------------------------------------------------------
AC_ARG_WITH([mpi],
            AS_HELP_STRING([--with-mpi],
                           [mpich2 is required for parallel build]))

# check for MPI installation
if test "x$enable_parallel" = "xyes"; then
    # check for MPICH PKG_CONFIG
    PKG_CHECK_MODULES([MPI],[mpich2-ch3],[mpich2='yes'],[mpich2='no'])

    # check for OPENMPI PKG_CONFIG
    if test "x$mpich2" = "xno"; then
        PKG_CHECK_MODULES([MPI],[openmpi],[openmpi='yes'],[openmpi='no'])
    fi

    # if build is debug, drop all flags except IGNORE_SEEK and include path
    if test "x$enable_debug" = "xyes" ; then
        for flag in $MPI_CFLAGS; do
            MPI_CFLAGS="$MPI_CFLAGS $(echo $flag | grep IGNORE_CXX_SEEK)"
            MPI_CFLAGS="$MPI_CFLAGS $(echo $flag | grep include)"
        done
    fi

    if test "x$mpich2" = "xyes" -o  "x$openmpi" = "xyes"; then
        CPPFLAGS="$CPPFLAGS $MPI_CFLAGS"
        CFLAGS="$CFLAGS $MPI_CFLAGS"
        CXXFLAGS="$CXXFLAGS $MPI_CFLAGS"
        LDFLAGS="$LDFLAGS $MPICH_LIBS"
    elif test "x$with_mpi" != "xno"; then
        CPPFLAGS="$CPPFLAGS -I/$with_mpi/include"
        CFLAGS="$CFLAGS -I/$with_mpi/include"
        CXXFLAGS="$CXXFLAGS -I/$with_mpi/include"
        LDFLAGS="$LDFLAGS -L/$with_mpi/lib"
    fi

    # check mpi header
    AC_CHECK_HEADERS([mpi.h])

    # check LIBATH to set lib flags
    AC_CHECK_LIB([mpich],[MPI_Init],,[mpich='no'])
    if test "x$mpich" = "xno"; then
        AC_CHECK_LIB([mpi],[MPI_Init],,[openmpi='no'])
    fi
    if test "x$mpich" = "xno" -a "x$openmpi" = "x$no"; then
        AC_MSG_FAILURE([can't find any mpi library])
    fi
fi
AM_CONDITIONAL(PARALLEL_BUILD, test "x$enable_parallel" != "xno")

# ---------------------------------------------------------------------------
# Check for --enable-debug
# ---------------------------------------------------------------------------
AC_MSG_CHECKING([to enable debug build])
AC_ARG_ENABLE(debug,
              AS_HELP_STRING([--enable-debug],
                             [enable debugging (disabled by default)]),
              [HAVE_DEBUG=$enableval], [HAVE_DEBUG=yes])

if test "x$enable_debug" = "xyes" ; then

# Try to strip out optimization.
  CXXFLAGS=`echo "$CXXFLAGS " | sed "s/-O2 //g"`
  CFLAGS=`echo "$CFLAGS " | sed "s/-O2 //g"`
  CXXFLAGS=`echo "$CXXFLAGS" | sed "s/-O //g"`
  CFLAGS=`echo "$CFLAGS" | sed "s/-O //g"`

# remove existing -g to avoid duplicates.
  CXXFLAGS=`echo "$CXXFLAGS" | sed "s/-g //"`
  CFLAGS=`echo "$CFLAGS" | sed "s/-g //"`

#  -g and -DDEBUG
  CXXFLAGS="$CXXFLAGS -O0 -g -DDEBUG"
  CFLAGS="$CFLAGS -O0 -g -DDEBUG -Wall"

  AC_MSG_RESULT([yes, modified CFLAGS="$CFLAGS"])
else
  AC_MSG_RESULT([no, CFLAGS="$CFLAGS"])
fi

# ---------------------------------------------------------------------------
# Check for lapack
# ---------------------------------------------------------------------------
AC_CHECK_LIB([lapack],[dgesv_])
AC_CHECK_LIB([z],[gzdopen])
#-----------------------------------------------------------------------
# Checks for HDF5 libraries.
#-----------------------------------------------------------------------
AC_ARG_WITH([hdf5],
            AS_HELP_STRING([--with-hdf5],
                            [hdf5 is required to build particler]))

if test -n "$with_hdf5"; then
   if test "x$with_hdf5" != "xno"; then
      CPPFLAGS="$CPPFLAGS -I$with_hdf5/include"
      CFLAGS="$CFLAGS -I$with_hdf5/include"
      CXXFLAGS="$CXXFLAGS -I$with_hdf5/include"
      LDFLAGS="$LDFLAGS -L$with_hdf5/lib"
   else
      AC_MSG_FAILURE([can't proceed without hdf5 installtion])
   fi
fi
AC_CHECK_LIB([hdf5],[H5open],,
             [AC_MSG_FAILURE([hdf5 is needed to build this code. Install hdf5 and try again.])])
#-----------------------------------------------------------------------
# Checks for GDAL libraries.
#-----------------------------------------------------------------------
AC_ARG_WITH([gdal],
            [AS_HELP_STRING([--with-gdal-config],
                            [provide full-path to gdal-config])
            ])

if test -n "$with_gdal" && test "x$with_gdal" != "xno"; then
  GDALINC=`$with_gdal --cflags`
  GDALLIB=`$with_gdal --libs`
  CPPFLAGS="$CPPFLAGS $GDALINC"
  LDFLAGS="$LDFLAGS $GDALLIB"
fi


AC_CHECK_HEADERS([gdal.h])
AC_CHECK_LIB([gdal],[GDALOpen],,)
AM_CONDITIONAL(HAVE_GDAL, test -n "$with_gdal" && test "x$with_gdal" != "xno")

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([hdf5.h malloc.h])
# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 gisapi/Makefile
                 datstr/Makefile
                 multiproc/Makefile
                 src/Makefile
                 preproc/Makefile
                ])
AC_OUTPUT
