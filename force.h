#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <array>
#include <complex>
#include "qcustomplot/qcustomplot.h"
#include <QDialog>
#include "optics/optics.h"

namespace Ui {
class force;
}

class force : public QDialog
{
    Q_OBJECT

public:
    explicit force(QWidget *parent = nullptr);
    ~force();

    void set_image(std::vector<std::array<std::complex<double>, 3>>& field, optics::environment env, double step );

private:
    Ui::force *ui;
};

#endif // FORCE_H
