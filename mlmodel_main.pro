

TARGET = CLASS
TEMPLATE = app
SOURCES += main.cpp \
    epmodel.cpp \
    mainwindow.cpp \
    model.cpp \
    modelrun.cpp \
    plotwindow.cpp \
    subplot.cpp \
    modeloutput.cpp \
    switches.cpp \
    defaultinput.cpp \
    landsoil.cpp \
    modelchem.cpp \
    modelinput.cpp
HEADERS += mainwindow.h \
    epmodel.h \
    mlm_class.h \
    modeloutput.h \
    modelinput.h \
    model.h \
    modelrun.h \
    plotwindow.h \
    subplot.h \
    landsoil.h \
    modelchemtypes.h \
    modelchem.h
FORMS += mainwindow.ui \
    plotwindow.ui \
    subplot.ui

greaterThan(QT_MAJOR_VERSION, 4)
{
  QT += widgets 
  QT += printsupport
}

#For adding git hash to about window:
REVISION = $$system(git --git-dir $$PWD/.git --work-tree $$PWD describe --tags)
#REVISION = $$system(git describe --tags)
DEFINES += GITHASH=\\\"$$REVISION\\\"
