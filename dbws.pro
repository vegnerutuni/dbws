TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    parameter.cpp \
    phase.cpp \
    spacegroup.cpp \
    diffractogram.cpp \
    dbwsexception.cpp

HEADERS += \
    param_inc.h \
    tabelas.h \
    parameter.h \
    phase.h \
    spacegroup.h \
    diffractogram.h \
    dbwsexception.h

