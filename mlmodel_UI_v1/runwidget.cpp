#include <QtGui>
#include <QString>
#include "runwidget.h"
#include "mlm_main.h"
#include <iostream>

const int IdRole = Qt::UserRole;

runWidget::runWidget(rundata *thisrun, QWidget *parent) : QWidget(parent)
{
    runButton = new QPushButton(tr("run"));
    tminInput = new QLineEdit(tr(""));
    tminLabel = new QLabel(tr("t-min:"));
    tminLabel->setBuddy(tminInput);
    tmaxInput = new QLineEdit(tr(""));
    tmaxLabel = new QLabel(tr("t-max:"));
    tmaxLabel->setBuddy(tmaxInput);
    ampInput = new QLineEdit(tr("1"));
    ampLabel = new QLabel(tr("Amplitude:"));
    ampLabel->setBuddy(ampInput);

    connect(runButton, SIGNAL(clicked()), this, SLOT(createData(thisrun)));

    QGridLayout *mainLayout = new QGridLayout;

    mainLayout->addWidget(tminLabel, 1,0, Qt::AlignRight);
    mainLayout->addWidget(tminInput, 1,1);
    mainLayout->addWidget(tmaxLabel, 2,0, Qt::AlignRight);
    mainLayout->addWidget(tmaxInput, 2,1);

    mainLayout->addWidget(ampLabel, 3,0, Qt::AlignRight);
    mainLayout->addWidget(ampInput, 3,1);

    mainLayout->addWidget(runButton, 4,1);
    setLayout(mainLayout);

    setWindowTitle(tr("Initiate model-run"));
}

void runWidget::createData(rundata *thisrun)
{
    double tmin = tminInput->text().toDouble();
    double tmax = tmaxInput->text().toDouble();
    double amp = ampInput->text().toDouble();
    double tsteps = (tmax-tmin)*60;
    //modelruns[widgetnum] = new rundata;

    *thisrun->h = 101.;

    std::cout << "Here we can read h: " << thisrun->h << std::endl;
}


