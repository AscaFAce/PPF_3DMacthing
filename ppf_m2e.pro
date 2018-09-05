QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        main.cpp \
    ppf3d_public.cpp \
    t_hash_int.cpp


INCLUDEPATH += /usr/include/pcl-1.8 \
               /usr/include/eigen3  \
               /usr/include/vtk-5.10 \
               /usr/include/boost

LIBS += /usr/lib/libvtk*.so \
        /usr/lib/x86_64-linux-gnu/libboost_*.so \
        /usr/lib/libpcl_*.so

HEADERS += \
    ppf3d_public.h \
    hash_murmur86.hpp \
    hash_murmur64.hpp \
    hash_murmur.hpp \
    t_hash_int.hpp \
    stdafx.h \
    icp_copy.hpp

