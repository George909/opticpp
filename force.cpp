#include "force.h"
#include "ui_force.h"

force::force(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::force)
{
    ui->setupUi(this);
}

force::~force()
{
    delete ui;
}

void force::set_image(std::vector<std::array<std::complex<double>, 3>>& field, optics::environment env, double step) {
// configure axis rect: 
    ui->plot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming

    ui->plot->axisRect()->setupFullAxesBox(true);
    ui->plot->xAxis->setLabel("x");
    ui->plot->yAxis->setLabel("y");
 
// set up the QCPColorMap:
    QCPColorMap *colorMap = new QCPColorMap(ui->plot->xAxis, ui->plot->yAxis);
    int nx = 128;
    int ny = 128;
    std::vector<double> force(nx*ny);
    optics::frc(field, force, nx, ny, step, env);
    colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
    colorMap->data()->setRange(QCPRange(-step * nx / 2, step * nx / 2), QCPRange(-step * ny / 2, step * ny / 2)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
// now we assign some data, by accessing the QCPColorMapData instance of the color map:
    double z;
    for (int xIndex=0; xIndex<nx; ++xIndex) {
        for (int yIndex=0; yIndex<ny; ++yIndex) {
            z = force[yIndex * nx + xIndex];
            colorMap->data()->setCell(xIndex, yIndex, z);
  }
}
 
    // add a color scale:
    QCPColorScale *colorScale = new QCPColorScale(ui->plot);
    ui->plot->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
    colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    colorMap->setColorScale(colorScale); // associate the color map with the color scale
    colorScale->axis()->setLabel("Force");
 
// set the color gradient of the color map to one of the presets:
    colorMap->setGradient(QCPColorGradient::gpPolar);
// we could have also created a QCPColorGradient instance and added own colors to
// the gradient, see the documentation of QCPColorGradient for what's possible.
 
// rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:
    colorMap->rescaleDataRange();
 
// make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
    QCPMarginGroup *marginGroup = new QCPMarginGroup(ui->plot);
    ui->plot->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
 
// rescale the key (x) and value (y) axes so the whole color map is visible:
    ui->plot->rescaleAxes();
};
