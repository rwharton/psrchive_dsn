dnl @synopsis SWIN_PACKAGE_OPTIONS(
dnl 
AC_DEFUN([SWIN_PACKAGE_OPTIONS],
[
  AC_PROVIDE([SWIN_PACKAGE_OPTIONS])

  AC_ARG_WITH([[$1]-dir],
              AC_HELP_STRING([--with-[$1]-dir=DIR],
                             [[$1] is in DIR]))

  AC_ARG_WITH([[$1]-include-dir],
              AC_HELP_STRING([--with-[$1]-include-dir=DIR],
                             [[$1] header files are in DIR]))

  AC_ARG_WITH([[$1]-lib-dir],
              AC_HELP_STRING([--with-[$1]-lib-dir=DIR],
                             [[$1] library is in DIR]))

  if test x"$with_[$1]_dir" = x"no" ||
     test x"$with_[$1]_include-dir" = x"no" ||
     test x"$with_[$1]_lib_dir" = x"no"; then

    # user disabled [$1]. Leave cache alone.
    have_[$1]="user disabled"

  else

    # "yes" is not a specification
    if test x"$with_[$1]_dir" = xyes; then
      with_[$1]_dir=""
    fi

    if test x"$with_[$1]_include_dir" = xyes; then
      with_[$1]_include_dir=""
    fi

    if test x"$with_[$1]_lib_dir" = xyes; then
      with_[$1]_lib_dir=""
    fi

    # _include_dir and _lib_dir take precedence over _dir:
    if test x"$with_[$1]_include_dir" = x &&
       test x"$with_[$1]_dir" != x; then
      if test -d $with_[$1]_dir/include/[$1]; then
        with_[$1]_include_dir=$with_[$1]_dir/include/[$1]
      else
        with_[$1]_include_dir=$with_[$1]_dir/include
      fi
    fi
    if test x"$with_[$1]_lib_dir" = x &&
       test x"$with_[$1]_dir" != x; then
      with_[$1]_lib_dir=$with_[$1]_dir/lib
    fi

  fi

])

dnl @synopsis SWIN_PACKAGE_LIB_OPTIONS(
dnl 
AC_DEFUN([SWIN_PACKAGE_LIB_OPTIONS],
[
  AC_PROVIDE([SWIN_PACKAGE_LIB_OPTIONS])

  AC_ARG_WITH([[$1]-dir],
              AC_HELP_STRING([--with-[$1]-dir=DIR],
                             [[$1] is in DIR]))

  if test x"$with_[$1]_dir" = x"no"; then

    # user disabled [$1]. Leave cache alone.
    have_[$1]="user disabled"

  else

    # "yes" is not a specification
    if test x"$with_[$1]_dir" = xyes; then
      with_[$1]_dir=""
    fi

    # _dir overrides include_dir and lib_dir    
    if test x"$with_[$1]_dir" != x; then
      with_[$1]_lib_dir=$with_[$1]_dir
    fi

  fi

])

dnl @synopsis SWIN_PACKAGE_FIND(name,file)
dnl 
AC_DEFUN([SWIN_PACKAGE_FIND],
[
  AC_PROVIDE([SWIN_PACKAGE_FIND])

  TEST_INCS=`echo $CPPFLAGS | awk '{ for (i=1; i<=NF; i++) printf ("%s\n", $i);}' | grep '^-I' | sed -e 's/-I//g'`
  TEST_LIBS=`echo $LDFLAGS | awk '{ for (i=1; i<=NF; i++) printf ("%s\n", $i);}' | grep '^-L' | sed -e 's/-L//g'`
  TEST_PACKAGES=`echo $TEST_INCS $TEST_LIBS $LD_LIBRARY_PATH $PACKAGES | sed -e 's/:/ /g'`

  AC_CACHE_VAL([swin_cv_[$1]_found],
  [swin_cv_[$1]_found=""
  if test x"$TEST_PACKAGES" != x; then
    for cf_file in `find -L $TEST_PACKAGES -name "[$2]" 2> /dev/null`; do
      cf_path=`dirname $cf_file`
      swin_cv_[$1]_found="$swin_cv_[$1]_found $cf_path"
    done
  fi])
])

dnl @synopsis SWIN_PACKAGE_TRY_COMPILE(name,includes,function body)
dnl 
AC_DEFUN([SWIN_PACKAGE_TRY_COMPILE],
[
  AC_PROVIDE([SWIN_PACKAGE_TRY_COMPILE])

  cf_include_path_list="$with_[$1]_include_dir $swin_cv_[$1]_found [$4] ."

  ac_save_CPPFLAGS="$CPPFLAGS"

  AC_CACHE_VAL([swin_cv_[$1]_include_dir],
  [swin_cv_[$1]_include_dir=""
  for cf_dir in $cf_include_path_list; do
    CPPFLAGS="-I$cf_dir $ac_save_CPPFLAGS"
    AC_TRY_COMPILE([$2], [$3], have_[$1]=yes, have_[$1]=no)
    if test $have_[$1] = yes; then
      swin_cv_[$1]_include_dir="$cf_dir"
      break
    fi
  done])

  if test x"$swin_cv_[$1]_include_dir" != x; then
    have_[$1]=yes
    if test x"$swin_cv_[$1]_include_dir" != x.; then
      [$1]_CFLAGS="-I$swin_cv_[$1]_include_dir"
    else
      [$1]_CFLAGS=""
      swin_cv_[$1]_include_dir=""
    fi
  fi
  CFLAGS="$ac_save_CFLAGS"
])

dnl @synopsis SWIN_PACKAGE_TRY_LINK(name,includes,function body,ldflags)
dnl 
AC_DEFUN([SWIN_PACKAGE_TRY_LINK],
[
  AC_PROVIDE([SWIN_PACKAGE_TRY_LINK])

  swin_search_path="$swin_cv_[$1]_found"

  #
  # If a path/include/ directory is found, then add path/lib/ to the search
  #
  if test x"$swin_cv_[$1]_include_dir" != x; then
    swin_base=`dirname $swin_cv_[$1]_include_dir`
    if test x"$swin_base" != x. -a -d $swin_base/lib; then
      swin_search_path="$swin_search_path $swin_base/lib"
    fi
    if test x"$swin_base" != x. -a -d $swin_base/lib64; then
      swin_search_path="$swin_search_path $swin_base/lib64"
    fi
  fi

  cf_lib_path_list="$with_[$1]_lib_dir $swin_search_path [$5] ."

  ac_save_CFLAGS="$CFLAGS"
  CFLAGS="$[$1]_CFLAGS $ac_save_CFLAGS"

  ac_save_LIBS="$LIBS"

  for cf_dir in $cf_lib_path_list; do
    LIBS="-L$cf_dir [$4] $ac_save_LIBS"
    AC_TRY_LINK([$2], [$3], have_[$1]=yes, have_[$1]=no)
    if test $have_[$1] = yes; then
      if test x"$cf_dir" = x.; then
        [$1]_LIBS="[$4]"
      else
        [$1]_LIBS="-L$cf_dir [$4]"
      fi
      break
    fi
  done

  CFLAGS="$ac_save_CFLAGS"
  LIBS="$ac_save_LIBS"

])

