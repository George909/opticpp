#ifndef POLARIZATION_H
#define POLARIZATION_H

#include <complex>
#include <array>
#include <vector>

#include <QDialog>
#include <QImage>
#include <QPainter>

namespace Ui {
class Polarization;
}

class Polarization : public QDialog
{
    Q_OBJECT

public:
    explicit Polarization(QWidget *parent = nullptr);

    QPixmap set_image(QImage& img,
                   std::vector<std::array<std::complex<double>, 3>>& field);

    ~Polarization();

private:
    Ui::Polarization *ui;
    QPixmap image;
};

#endif // POLARIZATION_H
