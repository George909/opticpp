#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>

#include <vector>
#include <complex>
#include <array>
#include <algorithm>
#include <execution>
#include "polarization.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_actionSettings_triggered();

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_polarization_clicked();

    void on_frenel_clicked();

private:
    Ui::MainWindow *ui;

    enum class dimension {X,Y,Z};

    void field_to_images(const std::vector<std::array<std::complex<double>, 3>>& field,
                         QImage& a, QImage& p, dimension d);
    void field_to_images(const std::vector<std::complex<double>>& field,
                         QImage& a, QImage& p);

    size_t image_w;
    size_t image_h;
    std::vector<std::complex<double>> input_field;
    QImage input_a, input_p, output_a, output_p,
           output_ax, output_px, output_ay,
           output_py, output_az, output_pz;

    std::vector<std::array<std::complex<double>, 3>> output_field;
    Polarization p;
};
#endif // MAINWINDOW_H
