TEMPLATE = lib

CONFIG -= app_bundle
CONFIG -= qt

# Directories
DESTDIR = dist/
OBJECTS_DIR = build/


QMAKE_CXXFLAGS_RELEASE += -fpermissive
QMAKE_CXXFLAGS_DEBUG += -fpermissive
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3


#QMAKE_CXXFLAGS+= -fopenmp
#QMAKE_LFLAGS +=  -fopenmp

CONFIG += c++11
QMAKE_CXXFLAGS += -std=c++11

mac : {
QMAKE_CXXFLAGS += -stdlib=libc++
}

QMAKE_LFLAGS += -lstdc++
QMAKE_LFLAGS += -dynamiclib


INCLUDEPATH += /usr/local/include
#for boost
INCLUDEPATH += /opt/local/include
LIBS += -L/opt/local/lib
LIBS += -L/usr/local/lib

mac : {
LIBS += -lboost_system-mt -lboost_filesystem-mt -lboost_program_options
}

unix : {
LIBS += -lboost_system -lboost_filesystem -lboost_program_options

}



SOURCES += source/main.cpp \
    source/simplicialcomplex.cpp \
    source/topsimplex.cpp \
    source/vertex.cpp \
    source/implicit_representation.cpp \
    source/io_functions.cpp \
    source/Timer.cpp

HEADERS += \
    source/simplicialcomplex.h \
    source/topsimplex.h \
    source/vertex.h \
    source/Usage.h \
    source/Timer.h



