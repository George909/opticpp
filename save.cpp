#include "save.h"
#include "ui_save.h"

save::save(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::save)
{
    ui->setupUi(this);
}

save::~save()
{
    delete ui;
}
