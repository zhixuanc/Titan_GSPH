#
# automake template for map directory
#
# $Id:$
#

noinst_LIBRARIES = libbgmesh.a

libbgmesh_a_SOURCES = bgmesh.h \
                      bnd_image.h bnd_image.cc \
                      bucket.h bucket.cc \
                      hashtab.h hashtab.cc \
                      hilbert.h hilbert.c \
                      mesh_update.cc \
                      particle.h particle.cc \
                      properties.h

libbgmesh_a_CPPFLAGS = -I$(top_srcdir)/src \
                       -I$(top_srcdir)/multiproc

