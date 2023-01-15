DEFINES+=DATDIR=\\\"$$PWD/dat\\\"
#DEFINES+=DATDIR=\\\"$$PWD\\\\"

QT += core
QT -= gui

CONFIG += c++11


TARGET = SpikingNetLarvaMBON
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

LIBS += `pkg-config gsl --libs`

#src/STPD_mod.cpp \

SOURCES += \
    src/CFSNeuron.cpp \
    src/CNeuron.cpp \
    src/CRSNeuron.cpp \
    src/CRSNeuron_test.cpp \
    src/SpikingNetwork.cpp \
    src/IFNeuron.cpp \
    src/INeuron.cpp \
    src/PoissonNeuron.cpp \
    src/PoissonSource.cpp \
    src/synapseEnsemble.cpp \
    src/synapticTransmission.cpp \
    src/main.cpp \
    src/synapseSW.cpp \
    src/isynapse.cpp \
    src/isynapseensemble.cpp

HEADERS += \
    src/CFSNeuron.h \
    src/CNeuron.h \
    src/CRSNeuron.h \
    src/IFNeuron.h \
    src/INeuron.h \
    src/PoissonNeuron.h \
    src/PoissonSource.h \
    src/stdafx.h \
    src/synapseEnsemble.h \
    src/synapticTransmission.h \
    src/synapseSW.h \
    src/isynapse.h \
    src/isynapseensemble.h

DISTFILES += \
    dat/plot.py
