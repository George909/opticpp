#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "optics/optics.h"
#include <chrono>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->image_h = 128;
    this->image_w = 128;   
    input_a = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    input_p = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_a = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_p = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_ax = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_px = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_ay = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_py = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_az = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
    output_pz = QImage(image_w, image_h, QImage::Format::Format_Grayscale8);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionSettings_triggered()
{

}

void MainWindow::field_to_images(const std::vector<std::array<std::complex<double>, 3>> &field,
                                 QImage &a, QImage &p, MainWindow::dimension d)
{
    std::complex<double> c = (*std::max_element(std::execution::par,
                                               field.begin(),
                                               field.end(),
                                               [=](const auto& c1, const auto& c2) -> bool {
                                                     return std::abs(c1[(size_t)d]) < std::abs(c2[(size_t)d]);}))[(size_t)d];
    double max_a = std::abs(c);

    std::transform(std::execution::par,
                   field.begin(),
                   field.end(),
                   a.bits(),
                   [=](const auto& c) {
                       return std::abs(c[(size_t)d]) / max_a * 255.;});

    std::transform(std::execution::par,
                   field.begin(),
                   field.end(),
                   p.bits(),
                   [=](const auto& c){
        return (std::arg(c[(size_t)d]) + M_PI) / 2. / M_PI * 255.;});
}

void MainWindow::field_to_images(const std::vector<std::complex<double>> &field,
                                 QImage &a, QImage &p)
{
    double max_a = std::abs(*std::max_element(std::execution::par,
                                              field.begin(),
                                              field.end(),
                                              [](const auto& c1, const auto& c2) -> bool {
                                                    return std::abs(c1) < std::abs(c2);}));

    std::transform(std::execution::par,
                   field.begin(),
                   field.end(),
                   a.bits(),
                   [max_a](const auto& c) {
                       return std::abs(c) / max_a * 255.;});

    std::transform(std::execution::par,
                   field.begin(),
                   field.end(),
                   p.bits(),
                   [](const std::complex<double>& c){
                       return (std::arg(c) + M_PI) / 2. / M_PI * 255.;});
}

//CALCULTE INPUT FIELD
void MainWindow::on_pushButton_clicked() {

    this->input_field = optics::beam(image_w, image_h, 
                                     optics::some_ampl(ui->sigma->text().toDouble(), 0.25), 
                                     optics::vortex_phase(ui->t_charge->text().toDouble(), 
                                                          ui->phi->text().toDouble()));

    field_to_images(input_field, input_a, input_p);

    ui->input_a->setPixmap(QPixmap::fromImage(input_a));
    ui->input_p->setPixmap(QPixmap::fromImage(input_p));
}

//CALCULATE OUTPUT FIELD
void MainWindow::on_pushButton_2_clicked()
{
    optics::rw_info inf;
    inf.na = this->ui->na->text().toDouble();
    inf.lambda = 532e-9;
    inf.in_h = this->ui->size_input_y->text().toDouble();
    inf.in_w = this->ui->size_input_x->text().toDouble();
    inf.out_h = this->ui->size_output_y->text().toDouble();
    inf.out_w = this->ui->size_output_x->text().toDouble();

    this->output_field.resize(image_h * image_w);
    auto start = std::chrono::steady_clock::now();
    switch (this->ui->pol_type->currentIndex()) {
        case 0 :
            optics::richards_wolf<optics::linear>(input_field, inf, image_h, image_w, output_field);
            break;
        case 1 :
            optics::richards_wolf<optics::rcircular>(input_field, inf, image_h, image_w, output_field);
            break;
        case 2 :
            optics::richards_wolf<optics::lcircular>(input_field, inf, image_h, image_w, output_field);
            break;
        case 3 :
            optics::richards_wolf<optics::radial>(input_field, inf, image_h, image_w, output_field);
            break;
        case 4 :
            optics::richards_wolf<optics::azimutal>(input_field, inf, image_h, image_w, output_field);
            break;
    }
    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    qDebug() << time.count() << Qt::endl;

    field_to_images(output_field, output_ax, output_px, dimension::X);
    field_to_images(output_field, output_ay, output_py, dimension::Y);
    field_to_images(output_field, output_az, output_pz, dimension::Z);

    auto compare_vec3 ([=](const auto& c1, const auto& c2) -> bool {
        return std::abs(c1[0]) + std::abs(c1[1]) + std::abs(c1[2]) <
               std::abs(c2[0]) + std::abs(c2[1]) + std::abs(c2[2]);});

    std::array<std::complex<double>,3> c = (*std::max_element(std::execution::par,
                                               output_field.begin(),
                                               output_field.end(),
                                               compare_vec3));

    double max_a = std::abs(c[0]) + std::abs(c[1]) + std::abs(c[2]);

    std::transform(std::execution::par,
                   output_field.begin(),
                   output_field.end(),
                   output_a.bits(),
                   [=](const auto& c) {
                       return (std::abs(c[0]) + std::abs(c[1]) + std::abs(c[2])) / max_a * 255.;});

    ui->output_a->setPixmap(QPixmap::fromImage(output_a));
    ui->output_ax->setPixmap(QPixmap::fromImage(output_ax));
    ui->output_px->setPixmap(QPixmap::fromImage(output_px));
    ui->output_ay->setPixmap(QPixmap::fromImage(output_ay));
    ui->output_py->setPixmap(QPixmap::fromImage(output_py));
    ui->output_az->setPixmap(QPixmap::fromImage(output_az));
    ui->output_pz->setPixmap(QPixmap::fromImage(output_pz));
}

void MainWindow::on_polarization_clicked(){
    p.set_image(this->output_a, this->output_field);
    p.show();
}

void MainWindow::on_frenel_clicked()
{
    optics::frenel_inf inf;

    inf.d = this->ui->na->text().toDouble();
    inf.lambda = 532e-9;
    inf.in_h = this->ui->size_input_y->text().toDouble();
    inf.in_w = this->ui->size_input_x->text().toDouble();
    inf.out_h = this->ui->size_output_y->text().toDouble();
    inf.out_w = this->ui->size_output_x->text().toDouble();

    std::vector<std::complex<double>> field = optics::frenel(input_field, inf, image_w, image_h);

    field_to_images(field, output_a, output_p);

    ui->output_a->setPixmap(QPixmap::fromImage(output_a));
    ui->output_i->setPixmap(QPixmap::fromImage(output_p));
}
