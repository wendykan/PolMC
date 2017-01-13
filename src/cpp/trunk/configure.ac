# Process this file with autoconf to produce a configure script.
AC_INIT(polmc, 4.0, dccote@novajo.ca)
AC_CONFIG_SRCDIR([MuellerM.h])
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_LANG(C++)

# Checks for libraries (and add them to options)
AC_CHECK_LIB([m], [abs])

# Where is CoreFoundation ?
AC_ARG_WITH(corefoundation-includedir,AC_HELP_STRING([--with-corefoundation-includedir=value],
                      [indicate where CoreFoundation keeps its header files]) )
AC_ARG_WITH(corefoundation-libdir,AC_HELP_STRING([--with-corefoundation-libdir=value],
                      [indicate where CoreFoundation keeps its library file]) )

# Custom checks: we want to scan $libdir and $includedir
# but they are not passed by default when calling
# AC_CHECK_LIB.  Hence, we roll our own check by adding them
# to LDFLAGS and CPPFLAGS

# First, modify compiler and linker flags
SAVE_CXXFLAGS="$CXXFLAGS"
SAVE_LDFLAGS="$LDFLAGS"
eval test -d $with_corefoundation_includedir && CXXFLAGS="-I$with_corefoundation_includedir $CXXFLAGS"
eval test -d $with_corefoundation_libdir && LDFLAGS="-L$with_corefoundation_libdir $LDFLAGS"


# corefoundation support:does the user want CF support?
AC_ARG_ENABLE(corefoundation,AC_HELP_STRING([--enable-corefoundation],
                       [enable CoreFoundation for complex input files (default is YES)]),
	AC_DEFINE(__COREFOUNDATION,1, [ Set to 1 to use Core Foundation or Core Foundation Lite code but leave undefined for not using ] ),
	AC_DEFINE(__COREFOUNDATION,0, [ Set to 0 to use Core Foundation or Core Foundation Lite code but leave undefined for not using ] ) )
#If nothing was provided, set default to "yes"
if test -z "$enable_corefoundation"; then
                enable_corefoundation="yes"
fi

#Should check, but too complicated for now

have_corefoundation="yes"

if test "$enable_corefoundation" = yes; then
        if test "$have_corefoundation" = yes; then
                LIBS="-lCoreFoundation $LIBS"
                LDFLAGS="-L$with_corefoundation_libdir $LDFLAGS"
                LIBADD="-L$with_corefoundation_libdir $LIBADD"
                CXXFLAGS="-I$with_corefoundation_includedir $CXXFLAGS"
                AC_MSG_NOTICE([XML support via CoreFoundation Lite is enabled by user])
        else
                AC_MSG_ERROR([You must have CoreFoundation Lite installed. 
You must also provide prefix, libdir and includedir for the configure script, for instance:
./configure --enable-corefoundation--with-corefoundation-libdir=$HOME/lib --with-corefoundation-includedir=$HOME/include ])
        fi
else
        AC_MSG_NOTICE([XML support via CF-Lite is not required by user])
fi

# Return to previous state
#CXXFLAGS="$SAVE_CXXFLAGS"
#LDFLAGS="$SAVE_LDFLAGS"
#LIBS="-lpthread -lm $LIBS"

# Where is SCATMECH?
AC_ARG_WITH(scatmech-includedir,AC_HELP_STRING([--with-scatmech-includedir=value],
                      [indicate where SCATMECH keeps its header files]) )
AC_ARG_WITH(scatmech-libdir,AC_HELP_STRING([--with-scatmech-libdir=value],
                      [indicate where SCATMECH keeps its library file]) )

# Custom checks: we want to scan $libdir and $includedir
# but they are not passed by default when calling
# AC_CHECK_LIB.  Hence, we roll our own check by adding them
# to LDFLAGS and CPPFLAGS

# First, modify compiler and linker flags
SAVE_CXXFLAGS="$CXXFLAGS"
SAVE_LDFLAGS="$LDFLAGS"
eval test -d $with_scatmech_includedir && CXXFLAGS="-I$with_scatmech_includedir $CXXFLAGS"
eval test -d $with_scatmech_libdir && LDFLAGS="-L$with_scatmech_libdir $LDFLAGS"

#Then actually try to compile a bit of code and see if it works

# SCATMECH support:do we have SCATMECH installed?
AC_MSG_CHECKING([whether SCATMECH is installed])
AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[#include "askuser.h"]],[[SCATMECH::AskUser("",0);]]),
		[ have_scatmech="yes"] ,[ have_scatmech="no"] )
AC_MSG_RESULT($have_scatmech)

# SCATMECH support:does the user want SCATMECH support?
AC_ARG_ENABLE(scatmech,AC_HELP_STRING([--enable-scatmech],
                       [enable SCATMECH for polarization Monte Carlo (default is NO)]),
						AC_DEFINE(WITH_POLARIZATION,1, Set to 1 to use polarization Monte Carlo ),
						AC_DEFINE(WITH_POLARIZATION,0, Set to 0 to use intensity only Monte Carlo) )
#If nothing was provided, set default to "no"
if test -z "$enable_scatmech"; then
		enable_scatmech="no"
fi

# Return to previous state
CXXFLAGS="$SAVE_CXXFLAGS"
LDFLAGS="$SAVE_LDFLAGS"
LIBS="-lpthread -lm $LIBS"


# We are done and we know whether or not we have SCATMECH

if test "$enable_scatmech" = yes; then
	if test "$have_scatmech" = yes; then
		LIBS="-lscatmech $LIBS"
		LDFLAGS="-L$with_scatmech_libdir $LDFLAGS"
		LIBADD="-L$with_scatmech_libdir $LIBADD"
		CXXFLAGS="-I$with_scatmech_includedir $CXXFLAGS"
		AC_DEFINE(NO_SCATMECH_NAMESPACE, 1, Set to 1 so that no SCATMECH namespace used)
		AC_MSG_NOTICE([Polarization support via SCATMECH is enabled by user])
	else
		AC_MSG_ERROR([You must have the SCATMECH library installed for polarization support.
Get it at http://physics.nist.gov/scatmech, or do not enable scatmech support.
You must also provide prefix, libdir and includedir for the configure script, for instance:
./configure --enable-scatmech --with-scatmech-libdir=$HOME/lib --with-scatmech-includedir=$HOME/include ])
	fi
else
	AC_MSG_NOTICE([Polarization support via SCATMECH is not required by user])
fi

# Debugging support:does the user want debugging support?
AC_ARG_ENABLE(debug,AC_HELP_STRING([--enable-debug],
                       [enable debugging code]),
						AC_DEFINE(__MYDEBUG,1, [ Set to 1 to use debugging code but leave undefine for no debugging code ] ) )
#If nothing was provided, set default to "no"
if test -z "$enable_debug"; then
		enable_debug="no"
fi

if test "$enable_debug" = yes; then
		CXXFLAGS="-g $CXXFLAGS"
		AC_MSG_NOTICE([Debugging code enabled by user])
else
		AC_MSG_NOTICE([Debugging code is not activated by user])
fi

# Profiling support:does the user want profiling support?
AC_ARG_ENABLE(profile,AC_HELP_STRING([--enable-profile],
                       [enable profiling code]) )
#If nothing was provided, set default to "no"
if test -z "$enable_profile"; then
		enable_profile="no"
fi

if test "$enable_profile" = yes; then
		CXXFLAGS="-p $CXXFLAGS"
		AC_MSG_NOTICE([Profiling code enabled by user])
else
		AC_MSG_NOTICE([Profiling code is not activated by user])
fi


# What random number generator do we want?
AC_ARG_WITH(random-generator,AC_HELP_STRING([--with-random-generator=value],
                       [pick random number generator: mtrand or numrecipes]) )
#If nothing was provided, set default to "mtrand"
if test -z "$with_random_generator"; then
		with_random_generator="mtrand"
fi

if test "$with_random_generator" = mtrand; then
		AC_DEFINE(USE_MTRAND,1, [ Set to 1 to use mtrand random number generator ] )
		AC_MSG_NOTICE([Using mtrand for random number generator])
else
	 if test "$with_random_generator" = numrecipes; then
		AC_DEFINE(USE_NUMRECIPES_RAND,1, [ Set to 1 to use  random number generator from numerical recipes] )
		AC_MSG_NOTICE([Using random number generator from numerical recipes])
	fi
fi

# Do we have libxml-2.0 for MCGenericObject support?
AC_CHECK_PROG(with_libxml2, "xml2-config", "yes", "no")
if test "$with_libxml2" = yes; then
        CXXFLAGS="`xml2-config --cflags` $CXXFLAGS"
        LDFLAGS="`xml2-config --libs ` $LDFLAGS"
	AC_DEFINE(HAVE_LIBXML2,1,[ When libxml-2.0 is available, this is set to one and support for generic object is activated ])
	AC_MSG_NOTICE([Support for generic object defined via XML is activated])
else
	AC_MSG_NOTICE([Support for generic object defined via XML is not activated since libxml-2.0 is not installed on this system])
fi

AM_CONDITIONAL(XML_SUPPORT, test "$with_libxml2" = yes)


AC_CHECK_TYPES([long long])

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_CHECK_DECL(runtime_error,
                [AC_MSG_ERROR([We need `runtime_error'!])])

# Checks for library functions.
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([sqrt])
AC_CHECK_FUNCS([nan])
AC_CHECK_FUNCS([isnan])
AC_CHECK_FUNCS([isfinite])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
