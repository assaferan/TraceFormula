dnl
dnl Process this file with autoconf to produce a configure script.
dnl 

dnl
dnl Setup autoconf
dnl
AC_PREREQ([2.69])
AC_INIT([TraceFormula],[1.0.0-beta.1],[https://github.com/assaferan/TraceFormula],[TraceFormula],[https://github.com/assaferan/TraceFormula])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([cnf])

# AC_CANONICAL_HOST is needed to access the 'host_os' variable    
AC_CANONICAL_HOST

build_linux=no
build_windows=no
build_mac=no

# Detect the target system
case "${host_os}" in
    linux*)
        build_linux=yes
        ;;
    cygwin*|mingw*)
        build_windows=yes
        ;;
    darwin*)
        build_mac=yes
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

AM_CONDITIONAL([LINUX], [test "$build_linux" = "yes"])
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])
AM_CONDITIONAL([OSX], [test "$build_mac" = "yes"])

dnl
dnl C is the language
dnl
: ${CFLAGS=""}
AC_PROG_CC
AC_LANG([C])

dnl
dnl Setup automake
dnl
AM_INIT_AUTOMAKE([1.11 -Wall -Werror foreign subdir-objects])
AM_SILENT_RULES([yes])
AM_PROG_AR
LT_INIT
LT_INIT

dnl
dnl Setting the location of pari
dnl
if test "${PARI_PATH+set}" != set; then
     export PARI_PATH="/usr/local/"
fi

echo PARI_PATH=${PARI_PATH}

AC_ARG_VAR(PARI_PATH, [location of PARI/GP installation])

dnl
dnl Finally, generate the Makefiles and output everything
dnl
AC_CONFIG_FILES([Makefile src/Makefile])
AM_PROG_CC_C_O

AC_OUTPUT()
