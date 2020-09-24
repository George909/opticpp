#include "polarization.h"
#include "ui_polarization.h"

Polarization::Polarization(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Polarization)
{
    ui->setupUi(this);
}


void Polarization::set_image(QImage& img,
                             std::vector<std::array<std::complex<double>, 3>>& field)
{
    double s0, s1, s2, s3;

    double delta{};
    double hi, ratio, angle;
    int k{};

    this->image = QPixmap::fromImage(img);
    QPainter painter(&image);
    painter.setPen(QPen(Qt::green));

    for (int i = 16; i < image.height(); i += 16) {
        for (int j = 16; j < image.width(); j += 16) {
            k = i * image.width() + j;
            delta = std::arg(field[k][0]) - std::arg(field[k][1]);

            s0 = std::pow(std::abs(field[k][0]), 2.) + std::pow(std::abs(field[k][1]), 2.);
            s1 = std::pow(std::abs(field[k][0]), 2.) - std::pow(std::abs(field[k][1]), 2.);
            s2 = 2 * std::abs(field[k][0]) * std::abs(field[k][1]) * std::cos(delta);
            s3 = 2 * std::abs(field[k][0]) * std::abs(field[k][1]) * std::sin(delta);

            hi = std::asin(s3 / s0) / 2.;
            ratio = std::tan(hi);
            angle = (std::atan2(s2, s1)) / 2. / M_PI * 180;

            double arrow_angle = 0.45;
            painter.save();
            painter.translate(i, j);
            painter.rotate(-angle);
            painter.translate(-i, -j);
            painter.drawEllipse(QPoint(i, j), 8, static_cast<int>(8 * ratio));
            painter.restore();

            if (s3 < 0) {
                painter.drawLine(i, j, i + 4, j);
                painter.drawLine(i - 4, j,
                             i + 4 - 4 * std::cos(arrow_angle),
                             j - 4 * std::sin(arrow_angle));

                painter.drawLine(i - 4, j,
                             i + 4 - 4 * std::cos(arrow_angle),
                             j + 4 * std::sin(arrow_angle));
            }
            else {
                painter.drawLine(i, j, i - 4, j);
                painter.drawLine(i + 4, j,
                                 i - 4 + 4 * std::cos(arrow_angle),
                                 j - 4 * std::sin(arrow_angle));

                painter.drawLine(i + 4,
                                 j, i - 4 + 4 * std::cos(arrow_angle),
                                 j + 4 * std::sin(arrow_angle));
            }
        }
    }
    ui->image->setPixmap(this->image);
}

Polarization::~Polarization()
{
    delete ui;
}
