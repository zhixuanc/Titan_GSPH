#
# automake template for map directory
#
# $Id:$
#

bin_PROGRAMS = preprocess

preprocess_SOURCES = buckstr.h \
                     createfunky.cc \
                     preprocess.h preprocess.cc  \
                     hilbert.h hilbert.c \
                     gisapi.h inc_plane.cc

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/gisapi
           

preprocess_LDADD = $(top_builddir)/src/libgmfghdf5.a \
                   $(top_builddir)/gisapi/libgisapi.a
