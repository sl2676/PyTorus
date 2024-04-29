lib_LTLIBRARIES = libebf_c.la
libebf_c_la_SOURCES = ebftable.h ebftable.c ebf.c ebf.h
libebf_c_la_LDFLAGS = -version-info 0:0:0
include_HEADERS= ebftable.h ebf.h

# TESTS = ebf_libc_test
# check_PROGRAMS=ebf_libc_test
bin_PROGRAMS=ebf_libc_test
ebf_libc_test_SOURCES= ebf_libc_test.c ebf_demo.c ebf_demo.h
ebf_libc_test_LDADD=libebf_c.la
# ebfdatadir=$(prefix)"/../data"
# AM_CFLAGS = -DEBFDATADIR=\"$(ebfdatadir)\"

