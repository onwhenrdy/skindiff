TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

unix*{
    QMAKE_CXXFLAGS += -pedantic -Werror
    QMAKE_CXXFLAGS += -Wno-sign-compare
    QMAKE_CXXFLAGS_WARN_ON -= -W3

    QMAKE_CXXFLAGS_RELEASE -= -O2 -O1
    QMAKE_CXXFLAGS_RELEASE += -DNDEBUG -O3 -mtune=native -march=native

    QMAKE_LFLAGS_RELEASE -= -01
}


Release::DIR = release
Release:DESTDIR = $$DIR/build
Release:OBJECTS_DIR = $$DIR/obj
Release:MOC_DIR = $$DIR/moc
Release:RCC_DIR = $$DIR/rcc
Release:UI_DIR = $$DIR/ui

Debug::DIR = debug
Debug:DESTDIR = $$DIR/build
Debug:OBJECTS_DIR = $$DIR/obj
Debug:MOC_DIR = $$DIR/moc
Debug:RCC_DIR = $$DIR/rcc
Debug:UI_DIR = $$DIR/ui

SOURCES += main.cpp \
    ../../tdmatrix.cpp \
    tdmatrixtests.cpp \
    thomasalgorithmtests.cpp

HEADERS += \
    ../../tdmatrix.h
