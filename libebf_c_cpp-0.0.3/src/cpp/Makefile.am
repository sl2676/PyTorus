lib_LTLIBRARIES = libebf_cpp.la
libebf_cpp_la_SOURCES = ebftable.hpp ebftable.cpp ebf.cpp ebf.hpp ebfvector.hpp
libebf_cpp_la_LDFLAGS = -version-info 0:0:0
include_HEADERS = ebftable.hpp ebf.hpp ebfvector.hpp

bin_PROGRAMS=ebf_libcpp_test
# ebftk_SOURCES= timer.cpp  timer.hpp  toolkit.cpp
# ebftk_LDADD=libebf_cpp.la
# TESTS = ebf_libcpp_test
# check_PROGRAMS=ebf_libcpp_test
ebf_libcpp_test_SOURCES= ebf_libcpp_test.cpp ebf_demo.cpp ebf_demo.hpp timer.hpp timer.cpp
ebf_libcpp_test_LDADD=libebf_cpp.la

# ebfdatadir=$(datadir)"ebf/"
# AM_CPPFLAGS = -DEBFDATADIR=\"$(ebfdatadir)\"


