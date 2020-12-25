#ifndef SAVE_H
#define SAVE_H

#include <QDialog>

namespace Ui {
class save;
}

class save : public QDialog
{
    Q_OBJECT

public:
    explicit save(QWidget *parent = nullptr);
    ~save();

private:
    Ui::save *ui;
};

#endif // SAVE_H
