#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Timer.h>

// Qt headers
#include <QtGui>
#include <QString>

#include <fstream>
#include <limits>
#include <vector>
#include <fstream>
#include <list>

#include "include/internal/BeltramiKleinTraits.h"
#include "include/internal/PoincareDiskTraits.h"
#include "internal/RandomDomainGenerator.h"
#include "internal/Qt/Input/WASDInput.h"
#include "internal/Qt/Input/PointInput.h"
#include "include/internal/Qt/Input/PolylineInput.h"
#include "internal/Qt/GraphicItems/RoutingScenarioGraphicsItem.h"
#include "internal/RoutingScenario.h"
#include "internal/Qt/Input/TriangulationSelectTriangle.h"
#include "ui_Main.h"

typedef CGAL::Poincare_disk_traits<> Poincare_disk_traits;
typedef CGAL::Beltrami_klein_traits<> Beltrami_klein_traits;

//either choose Poincare_disk_traits or Beltrami_klein_traits
typedef Beltrami_klein_traits K;

typedef K::FT FT;
typedef K::Point_2 Point_2;
typedef K::Circle_2 Circle_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;

typedef CGAL::Delaunay_mesh_face_base_2<K> Face_base;
typedef CGAL::Triangulation_vertex_base_2<K> Vertex_base;
typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
typedef CGAL::No_constraint_intersection_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> Triangulation;

typedef CGAL::Qt::Routing_scenario<Triangulation> RoutingScenario;
typedef CGAL::Qt::Routing_scenario_graphics_item<Triangulation> RoutingScenarioGraphicsItem;
typedef Triangulation::Vertex_handle Vertex_handle;

class MainWindow :
        public CGAL::Qt::DemosMainWindow,
        public Ui::Routing_on_the_hyperbolic_plane {
    Q_OBJECT

    RoutingScenario routingScenario;
    RoutingScenarioGraphicsItem *routingGraphicsItem;
    Circle_2 p_disk = Circle_2(Point_2(0, 0), 1);
    QGraphicsEllipseItem *disk;
    QGraphicsScene scene;
    CGAL::Qt::PointInput<Triangulation> *pi;
    CGAL::Qt::TriangulationSelectTriangle<Triangulation> *select_triangle;
    CGAL::Qt::Polyline_input<Triangulation> *polyInput;
    CGAL::Qt::WASDInput<Triangulation> *wasd_input;

public:
    MainWindow();

    int path_position = -1;

public Q_SLOTS:
    void processInput(const CGAL::Object &o);
    void transformation();
    void changedActiveNode(const CGAL::Object &o);
    void changedDestinationNode(const CGAL::Object &o);

    void clear();
    void load_routing_scenario(const QString &);

    //all gui actions
    //options:
    void on_actionShowTriangulation_toggled(bool checked);
    void on_actionShowTriangulationBetween_toggled(bool checked);
    void on_actionShowFaces_toggled(bool checked);
    void on_actionShowConstraints_toggled(bool checked);
    void on_actionShowPoints_toggled(bool checked);
    void on_actionShowPointToPointVisibility_toggled(bool checked);
    void on_actionShowPointToAllVisibility_toggled(bool checked);
    void on_actionShowVisibilityGraph_toggled(bool checked);
    void on_actionShowDijkstraTree_toggled(bool checked);
    void on_actionShowDecomposition_toggled(bool checked);
    void on_actionShowOrigin_toggled(bool checked);
    void on_actionColorUnitDisk_triggered(bool checked);
    void on_actionShowUnitCircle_triggered(bool checked);

    //edit:
    void on_actionInsertConstraints_toggled(bool checked);
    void on_actionInsertPoint_toggled(bool checked);
    void on_actionInsertRandomPoints_triggered();
    void on_actionLoadRoutingScenario_triggered();
    void on_actionSetPointOnObstacle_toggled(bool checked);

    //tools:
    void on_actionExportPNG_triggered();
    void on_actionExportSVG_triggered();
    void on_actionSaveRoutingScenario_triggered();
    void on_actionClear_triggered();
    void on_actionRecenter_triggered();
    void on_actionDragMode_toggled(bool checked);
    void on_actionRedraw_triggered();
    void on_actionWalkPathForward_triggered();
    void on_actionWalkPathBackward_triggered();
    void on_actionTestSuite_triggered();
    void on_actionTranslateToPoint_triggered();
    void on_actionWalk_toggled(bool checked);

    //compute:
    void on_actionComputeVisibilityGraph_triggered();
    void on_actionDijkstra_triggered();
    void on_actionFindPath_toggled(bool checked);
    void on_actionGenerateRandomDomain_triggered();
    void on_pathOptimize_released();
};

void MainWindow::on_actionShowUnitCircle_triggered(const bool checked) {
    disk->setVisible(checked);
}

void MainWindow::on_actionColorUnitDisk_triggered(const bool checked) {
    if (checked) {
        disk->setBrush(QBrush(QColor(0, 150, 0)));
        routingGraphicsItem->epen.setBrush(::Qt::yellow);
        routingGraphicsItem->set_obstacle_brush(::Qt::blue);
    } else {
        disk->setBrush(::Qt::NoBrush);
        routingGraphicsItem->epen.setBrush(::Qt::darkGreen);
        routingGraphicsItem->set_obstacle_brush(QColor(140, 140, 140, 100));
    }
    routingGraphicsItem->repaint();
}

MainWindow::MainWindow()
    : DemosMainWindow() {
    setupUi(this);

    this->graphicsView->setAcceptDrops(false);

    // Add Poincaré disk
    qreal origin_x = CGAL::to_double(p_disk.center().x());
    qreal origin_y = CGAL::to_double(p_disk.center().y());
    qreal radius = std::sqrt(CGAL::to_double(p_disk.squared_radius()));
    qreal diameter = std::sqrt(CGAL::to_double(4 * p_disk.squared_radius()));
    qreal left_top_corner_x = origin_x - radius;
    qreal left_top_corner_y = origin_y - radius;
    qreal width = diameter, height = diameter;

    disk = new QGraphicsEllipseItem(left_top_corner_x, left_top_corner_y, width, height);
    QPen pen(::Qt::black, 0.000);
    disk->setPen(pen);
    scene.addItem(disk);

    routingScenario = RoutingScenario();
    routingGraphicsItem = new RoutingScenarioGraphicsItem(&routingScenario);

    //add routing_graphics_item to scene
    scene.addItem(routingGraphicsItem);

    //connecting input signals with routing_scenario and routing_graphics_item
    pi = new CGAL::Qt::PointInput<Triangulation>(&scene, routingGraphicsItem, this);
    QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
                     this, SLOT(processInput(CGAL::Object)));
    QObject::connect(pi, SIGNAL(changedActiveNode(CGAL::Object)),
                     this, SLOT(changedActiveNode(CGAL::Object)));
    QObject::connect(pi, SIGNAL(changedDestinationNode(CGAL::Object)),
                     this, SLOT(changedDestinationNode(CGAL::Object)));

    select_triangle = new CGAL::Qt::TriangulationSelectTriangle<Triangulation>(&scene, &routingScenario, routingGraphicsItem, this);
    QObject::connect(select_triangle, SIGNAL(changedActiveNode(CGAL::Object)),
                     this, SLOT(changedActiveNode(CGAL::Object)));
    QObject::connect(select_triangle, SIGNAL(changedDestinationNode(CGAL::Object)),
                     this, SLOT(changedDestinationNode(CGAL::Object)));

    polyInput = new CGAL::Qt::Polyline_input<Triangulation>(this, &scene, routingGraphicsItem);
    QObject::connect(polyInput, SIGNAL(generate(CGAL::Object)),
                     this, SLOT(processInput(CGAL::Object)));

    QObject::connect(this->actionQuit, SIGNAL(triggered()),
                     this, SLOT(close()));

    wasd_input = new CGAL::Qt::WASDInput<Triangulation>(routingGraphicsItem, this);
    QObject::connect(wasd_input, SIGNAL(transformed()), this, SLOT(transformation()));

    auto *ag1 = new QActionGroup(this);
    ag1->addAction(this->actionInsertPoint);
    ag1->addAction(this->actionSetPointOnObstacle);
    ag1->addAction(this->actionInsertConstraints);
    ag1->addAction(this->actionDragMode);

    this->actionInsertPoint->setChecked(true);
    this->actionSetPointOnObstacle->setChecked(false);
    this->actionDragMode->setChecked(false);
    this->actionInsertConstraints->setChecked(false);
    this->actionWalk->setChecked(false);
    ag1->setExclusionPolicy(QActionGroup::ExclusionPolicy::ExclusiveOptional);

    this->actionShowTriangulation->setChecked(false);
    this->actionShowFaces->setChecked(false);
    this->actionShowConstraints->setChecked(false);
    this->actionShowPoints->setChecked(false);
    this->actionShowOrigin->setChecked(false);

    auto *ag2 = new QActionGroup(this);
    ag2->addAction(this->actionShowVisibilityGraph);
    ag2->addAction(this->actionShowPointToPointVisibility);
    ag2->addAction(this->actionShowPointToAllVisibility);
    ag2->addAction(this->actionShowDijkstraTree);
    ag2->addAction(this->actionFindPath);
    ag2->setExclusionPolicy(QActionGroup::ExclusionPolicy::ExclusiveOptional);

    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowPointToPointVisibility->setChecked(false);
    this->actionShowPointToAllVisibility->setChecked(false);
    this->actionShowDijkstraTree->setChecked(false);

    auto *ag3 = new QActionGroup(this);
    ag3->addAction(this->actionShowTriangulation);
    ag3->addAction(this->actionShowTriangulationBetween);
    ag3->setExclusionPolicy(QActionGroup::ExclusionPolicy::ExclusiveOptional);

    // Setup the scene and the view
    scene.setItemIndexMethod(QGraphicsScene::NoIndex);
    scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
    this->graphicsView->setScene(&scene);
    this->graphicsView->setMouseTracking(true);

    // Turn the vertical axis upside down
    this->graphicsView->scale(1, -1);

    //Point_2 p = Point_2(0.9999, 0.999999999);
    //std::cout<< sizeof(Triangulation::Face) << std::endl;
    //std::cout << std::setprecision(100) << "tanh(14)" << std::tanh(18.715) << std::endl;
    //std::cout << nextafter(1.0, 0) << std::endl;

    // The navigation adds zooming and translation functionality to the
    // QGraphicsView
    this->addNavigation(this->graphicsView);

    this->setupStatusBar();
    this->setupOptionsMenu();
    this->addAboutCGAL();
}

void MainWindow::processInput(const CGAL::Object &o) {
    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowDijkstraTree->setChecked(false);
    this->actionShowPointToAllVisibility->setChecked(false);
    this->actionShowPointToPointVisibility->setChecked(false);
    this->actionFindPath->setChecked(false);
    std::vector<Point_2> points;
    if (CGAL::assign(points, o)) {
        routingScenario.insert_obstacle(points.begin(), points.end(), true);
        routingGraphicsItem->changed();
    }
    else {
        Point_2 p;
        if (CGAL::assign(p, o)) {
            routingScenario.insert_point(p);
            routingGraphicsItem->changed();
        }
    }
}

void MainWindow::transformation() {
    routingGraphicsItem->transformation();
}

void MainWindow::changedActiveNode(const CGAL::Object &o) {
    Vertex_handle vh;
    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowDijkstraTree->setChecked(false);
    this->actionShowPointToAllVisibility->setChecked(false);
    this->actionShowPointToPointVisibility->setChecked(false);
    this->actionFindPath->setChecked(false);
    if (CGAL::assign(vh, o)) {
        if (vh != routingScenario.start_node_handle && vh != routingScenario.destination_node_handle) {
            routingScenario.set_point_to_start(vh);
            routingGraphicsItem->repaint();
        }
    } else {
        Point_2 p;
        if (CGAL::assign(p, o)) {
            routingScenario.set_start_point(p);
            routingGraphicsItem->changed();
        }
    }
}

void MainWindow::changedDestinationNode(const CGAL::Object &o) {
    Vertex_handle vh;
    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowPointToAllVisibility->setChecked(false);
    this->actionShowPointToPointVisibility->setChecked(false);
    this->actionFindPath->setChecked(false);
    if (CGAL::assign(vh, o)) {
        if (vh != routingScenario.start_node_handle && vh != routingScenario.destination_node_handle) {
            routingScenario.set_point_to_destination(vh);
            routingGraphicsItem->repaint();
        }
    } else {
        Point_2 p;
        if (CGAL::assign(p, o)) {
            routingScenario.set_destination_point(p);
            routingGraphicsItem->changed();
        }
    }
}

void MainWindow::on_actionExportPNG_triggered() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    "Export to PNG", ".", "PNG (*.png)\n");
    if (!fileName.isNull()) {
        QPixmap pixMap = this->graphicsView->grab();
        pixMap.save(fileName);
    }
}

void MainWindow::on_actionExportSVG_triggered() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    "Export to SVG", ".", "SVG (*.svg)\n");
    //for svg
    QSvgGenerator generator;
    generator.setFileName(fileName);
    generator.setSize(QSize(static_cast<int>(this->graphicsView->viewport()->width()),
        static_cast<int>(this->graphicsView->viewport()->width())));
    generator.setViewBox(QRectF(-10,-10, this->graphicsView->viewport()->width() + 20,
        this->graphicsView->viewport()->width() + 20));
    QPainter painter;
    painter.begin(&generator);
    scene.clearSelection();
    scene.render(&painter);
    painter.end();
}

void MainWindow::on_actionInsertPoint_toggled(const bool checked) {
    if (checked) {
        scene.installEventFilter(pi);
    } else {
        scene.removeEventFilter(pi);
    }
}

void MainWindow::on_actionSetPointOnObstacle_toggled(const bool checked) {
    if (checked) {
        scene.addItem(select_triangle->active_triangle);
        scene.installEventFilter(select_triangle);
    } else {
        scene.removeEventFilter(select_triangle);
        scene.removeItem(select_triangle->active_triangle);
        select_triangle->active_triangle->setVisible(false);
    }
}

void MainWindow::on_actionInsertConstraints_toggled(const bool checked) {
    if (checked) {
        scene.installEventFilter(polyInput);
    } else {
        scene.removeEventFilter(polyInput);
    }
}

void MainWindow::on_actionWalk_toggled(const bool checked) {
    if (checked) {
        this->speedSpinBox->setDisabled(true);
        this->actionShowOrigin->setChecked(true);
        this->actionInsertPoint->setDisabled(true);
        this->actionInsertConstraints->setDisabled(true);
        this->actionSetPointOnObstacle->setChecked(false);
        this->actionSetPointOnObstacle->setDisabled(true);
        wasd_input->set_walking_speed(this->speedSpinBox->value());
        scene.installEventFilter(wasd_input);
    } else {
        this->actionShowOrigin->setChecked(false);
        this->speedSpinBox->setDisabled(false);
        this->actionInsertPoint->setDisabled(false);
        this->actionInsertConstraints->setDisabled(false);
        this->actionSetPointOnObstacle->setDisabled(false);
        scene.removeEventFilter(wasd_input);
        this->actionRecenter->trigger();
    }
}

void MainWindow::on_actionClear_triggered() {
    clear();
}

void MainWindow::clear() {
    routingScenario.clear();
    routingGraphicsItem->clear();
    this->actionShowPointToPointVisibility->setChecked(false);
    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowPointToAllVisibility->setChecked(false);
    this->actionShowDijkstraTree->setChecked(false);
    this->actionSetPointOnObstacle->setChecked(false);
    this->actionFindPath->setChecked(false);
    on_actionRecenter_triggered();
}

void MainWindow::on_actionShowTriangulation_toggled(const bool checked) {
    routingGraphicsItem->triangulation_graphics_item->set_visible_edges(checked);
    if (checked) {
        statusBar()->showMessage(QString("Triangulation has %1 nodes and %2 faces.").arg(
                                     routingScenario.number_of_vertices()).arg(routingScenario.t->number_of_faces()),
                                 4000);
    }
}

void MainWindow::on_actionDragMode_toggled(const bool checked) {
    if (checked) {
        this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
    } else {
        this->graphicsView->setDragMode(QGraphicsView::NoDrag);
    }
}

void MainWindow::on_actionShowFaces_toggled(const bool checked) {
    routingGraphicsItem->set_show_obstacle_interior(checked);
}

void MainWindow::on_actionShowConstraints_toggled(const bool checked) {
    routingGraphicsItem->triangulation_graphics_item->set_visible_constraints(checked);
}

void MainWindow::on_actionShowPoints_toggled(const bool checked) {
    statusBar()->showMessage(QString("There are %1 nodes.").arg(routingScenario.number_of_vertices()), 5000);
    routingGraphicsItem->triangulation_graphics_item->set_visible_nodes(checked);
}

void MainWindow::on_actionInsertRandomPoints_triggered() {
    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowDijkstraTree->setChecked(false);
    this->actionShowPointToAllVisibility->setChecked(false);
    this->actionShowPointToPointVisibility->setChecked(false);
    this->actionFindPath->setChecked(false);
    on_actionRecenter_triggered();
    bool ok = false;
    const int number_of_points = QInputDialog::getInt(this, "Number of random points",
                                                      "Enter number of random points", 5000,
                                                      0, (std::numeric_limits<int>::max)(), 1, &ok);
    if (!ok) {
        return;
    }
    CGAL::Timer timer;
    timer.start();
    auto generator = CGAL::Qt::Random_domain_generator<Triangulation>(&routingScenario);
    if(this->blueNoise->checkState()) {
        generator.blue_noise(this->candidatesSpinBox->value(), number_of_points, this->radiusSpinBox->value());
        std::cout << "blue noise took: " << timer.time() << " seconds." << std::endl;
        statusBar()->showMessage(QString("Blue noise: %1 seconds.").arg(timer.time()), 6000);
        timer.reset();
    } else {
        const std::vector<Point_2> points = generator.inverse_sampling(number_of_points, this->radiusSpinBox->value());
        timer.stop();
        std::cout << "sampling points took: " << timer.time() << " seconds." << std::endl;
        statusBar()->showMessage(QString("Sampling points: %1 seconds.").arg(timer.time()), 6000);
        timer.reset();
        timer.start();

        routingScenario.insert_points(points.begin(), points.end());

        timer.stop();
        statusBar()->showMessage(QString("Inserting points: %1 seconds.").arg(timer.time()), 6000);
        std::cout << "inserting points took: " << timer.time() << " seconds." << std::endl;
        timer.reset();
    }

    timer.start();
    routingScenario.remove_unconstrained_points_in_obstacle_interior();
    timer.stop();
    std::cout << "remove_unconstrained_points_in_obstacle_interior took: " << timer.time() << " seconds." << std::endl;
    statusBar()->showMessage(QString("remove_unconstrained_points_in_obstacle_interior: %1 seconds.").arg(timer.time()), 6000);
    timer.reset();

    routingGraphicsItem->changed();
}

void MainWindow::on_actionLoadRoutingScenario_triggered() {
    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowPointToPointVisibility->setChecked(false);
    this->actionShowPointToAllVisibility->setChecked(false);

    QString fileName = QFileDialog::getOpenFileName(this, "Open file", ".");
    if (!fileName.isEmpty()) {
        if (!fileName.isEmpty()) {
            if (routingScenario.number_of_vertices() > 2) {
                QMessageBox msgBox(QMessageBox::Warning,
                                   "Open new polygonal domain",
                                   "Do you really want to clear the current scene?",
                                   (QMessageBox::Yes | QMessageBox::No),
                                   this);
                int ret = msgBox.exec();
                if (ret == QMessageBox::Yes) {
                    clear();
                } else
                    return;
            }
            load_routing_scenario(fileName);
        }

        actionRecenter->trigger();
        routingGraphicsItem->changed();
    }
}

void MainWindow::load_routing_scenario(const QString &filename) {
    std::ifstream ifs(qPrintable(filename));
    std::vector<Point_2> obstacle;
    int length;
    std::string model;
    ifs >> model;
    bool b;
    int projection = -1;
    if(typeid(K) == typeid(Beltrami_klein_traits) && model == "P") {
        std::cout << "Translating points from Poincare disk model to Beltrami-Klein model" << std::endl;
        projection = 1;
    }
    if(typeid(K) == typeid(Poincare_disk_traits) && model == "B") {
        std::cout << "Translating points from Beltrami-Klein model to Poincare disk model" << std::endl;
        projection = 0;
    }

    CGAL::Timer timer;
    timer.start();

    while (ifs >> length) {
        for(int i = 0; i < length; ++i) {
            double x;
            ifs >> x;
            double y;
            ifs >> y;

            if(projection == -1){
                obstacle.push_back(Point_2(x, y));
            } else if (projection == 0) {
                obstacle.push_back(routingScenario.beltrami_klein_to_poincare(Point_2(x, y)));
            } else {
                obstacle.push_back(routingScenario.poincare_to_beltrami_klein(Point_2(x, y)));
            }
        }
        routingScenario.insert_obstacle(obstacle.begin(), obstacle.end(), true);
        obstacle.clear();
    }
    routingGraphicsItem->changed();

    timer.stop();
    statusBar()->showMessage(QString("Loading routing scenario took: %1 seconds").arg(timer.time()), 4000);
}

void MainWindow::on_actionSaveRoutingScenario_triggered() {
    QString fileName = QFileDialog::getSaveFileName(this, "Save file", ".");
    if (!fileName.isEmpty()) {
        std::ofstream ofs(qPrintable(fileName));
        if (typeid(K) == typeid(Beltrami_klein_traits)) {
            ofs << "B" << std::endl;
        } else {
            ofs << "P" << std::endl;
        }

        CGAL::Timer timer;
        timer.start();

        std::list<std::vector<Vertex_handle>> obstacles = routingScenario.get_obstacles();
        for (std::vector<Vertex_handle> obstacle: obstacles) {
            const int length = obstacle.size();
            ofs << length << std::endl;
            for(const Vertex_handle vh : obstacle) {
                ofs << std::fixed << std::setprecision(14) << vh->point().x() << " " << vh->point().y() << std::endl;
            }
        }

        timer.stop();
        statusBar()->showMessage(QString("Saving routing scenario took: %1 seconds").arg(timer.time()), 4000);
    }
}

void MainWindow::on_actionRecenter_triggered() {
    qreal origin_x = CGAL::to_double(p_disk.center().x());
    qreal origin_y = CGAL::to_double(p_disk.center().y());
    qreal radius = std::sqrt(CGAL::to_double(p_disk.squared_radius()));
    qreal diameter = std::sqrt(CGAL::to_double(4 * p_disk.squared_radius()));
    qreal scale = 1.1;

    this->graphicsView->setSceneRect(origin_x - radius, origin_y - radius, diameter, diameter);
    this->graphicsView->fitInView(origin_x - scale * radius, origin_y - scale * radius,
                                  scale * diameter, scale * diameter,
                                  Qt::KeepAspectRatio);
    routingGraphicsItem->reset_transformation();
}

void MainWindow::on_actionShowPointToPointVisibility_toggled(const bool checked) {
    if (checked) {
        if (routingScenario.start_node_handle != nullptr && routingScenario.destination_node_handle != nullptr) {
            bool b = routingScenario.compute_intersected_faces();
            if (b) {
                statusBar()->showMessage(QString("Nodes see each other."), 4000);
            } else {
                statusBar()->showMessage(QString("View is blocked by an obstacle."), 4000);
            }
            routingGraphicsItem->set_show_point_to_point_visibility(true);
        } else {
            statusBar()->showMessage(QString("Set start/destination node first."), 4000);
            this->actionShowPointToPointVisibility->setChecked(false);
        }
    } else {
        routingGraphicsItem->set_show_point_to_point_visibility(false);
    }
}

void MainWindow::on_actionShowVisibilityGraph_toggled(const bool checked) {
    if (checked) {
        if (routingScenario.defined_visibility_graph) {
            routingGraphicsItem->set_show_visibility_graph(checked);
        } else {
            statusBar()->showMessage(QString("First compute visibility graph."), 4000);
            routingGraphicsItem->set_show_visibility_graph(false);
            this->actionShowVisibilityGraph->setChecked(false);
        }
    } else {
        routingGraphicsItem->set_show_visibility_graph(false);
    }
}

void MainWindow::on_actionComputeVisibilityGraph_triggered() {
    this->actionFindPath->setChecked(false);
    this->actionShowVisibilityGraph->setChecked(false);
    this->actionShowDijkstraTree->setChecked(false);
    if (routingScenario.t->dimension() >= 2) {
        CGAL::Timer timer;
        timer.start();

        int number_of_orientation_tests = 0;
        if(this->comboBox->currentIndex() == 0) {
            number_of_orientation_tests = routingScenario.build_visibility_graph();
        }
        if(this->comboBox->currentIndex() == 1) {
            routingScenario.use_triangulation_as_visibility_graph();
        }
        if(this->comboBox->currentIndex() == 2) {
            routingScenario.build_visibility_graph_naive();
        }
        timer.stop();
        statusBar()->showMessage(
            QString("Building: %1 seconds, #Nodes: %2, #Edges: %3").arg(timer.time()).arg(
                routingScenario.number_of_vertices()).arg(routingScenario.adjacencies.size() / 2), 8000);
        std::cout << "Visibility-graph build took: " << timer.time() << " seconds. #Edges: " << routingScenario.adjacencies.size()
                / 2 << "." << std::endl;
        std::cout << "Vertex_index_map size: " << routingScenario.vertex_index_map.size() << std::endl;
        if(number_of_orientation_tests != 0) {
            std::cout << "Average time orientation test: " << timer.time()/number_of_orientation_tests << std::endl;
        }
        timer.reset();
    } else {
        statusBar()->showMessage(QString("Already computed."), 4000);
    }
}

void MainWindow::on_actionShowPointToAllVisibility_toggled(const bool checked) {
    if (checked) {
        if (routingScenario.start_node_handle != nullptr && routingScenario.t->dimension() >= 2) {
            routingScenario.compute_visible_vertices_from_start_node();
            statusBar()->showMessage(QString("Active_node can see %1 nodes. ").arg(
                                         routingScenario.visibles_start_node.size()), 4000);
            routingGraphicsItem->set_show_point_to_all_visibility(true);
        } else {
            statusBar()->showMessage(QString("Set start node first."), 4000);
            this->actionShowPointToAllVisibility->setChecked(false);
        }
    } else {
        routingGraphicsItem->set_show_point_to_all_visibility(false);
    }
}

void MainWindow::on_actionDijkstra_triggered() {
    if (routingScenario.defined_visibility_graph) {
        if (!routingScenario.defined_dijkstra) {
            //choose random node if start node not defined
            if (routingScenario.start_node_handle == nullptr) {
                routingScenario.set_point_to_start(routingScenario.index_vertex_map[0]);
            }
            CGAL::Timer timer;
            timer.start();
            routingScenario.dijkstra();
            timer.stop();
            //routingScenario.print_out_pred();
            statusBar()->showMessage(QString("Dijkstra: %1 seconds.").arg(timer.time()), 6000);
        } else {
            statusBar()->showMessage(QString("Already computed."), 4000);
        }
    } else {
        statusBar()->showMessage(QString("Compute visibility graph first."), 4000);
    }
}

void MainWindow::on_actionShowDijkstraTree_toggled(const bool checked) {
    if (checked) {
        if (routingScenario.defined_dijkstra) {
            routingGraphicsItem->set_show_dijkstra_tree(true);
        } else {
            statusBar()->showMessage(QString("Compute Dijkstra first."), 4000);
            this->actionShowDijkstraTree->setChecked(false);
        }
    } else {
        routingGraphicsItem->set_show_dijkstra_tree(false);
    }
}

void MainWindow::on_actionFindPath_toggled(const bool checked) {
    path_position = -1;
    if (checked) {
        if (routingScenario.start_node_handle != nullptr && routingScenario.destination_node_handle != nullptr) {
            if (routingScenario.defined_dijkstra) {
                const bool reachable = routingScenario.get_path_from_dijkstra();
                if (reachable) {
                    const double path_length = routingScenario.get_path_length();
                    routingGraphicsItem->set_show_path(true);
                    statusBar()->showMessage(QString("The path is %1 long.").arg(path_length), 8000);
                } else {
                    statusBar()->showMessage(QString("Destination node unreachable from start node."), 4000);
                    this->actionFindPath->setChecked(false);
                }
            } else {
                if (routingScenario.defined_visibility_graph) {
                    CGAL::Timer timer;
                    timer.start();
                    const bool reachable = routingScenario.a_star();
                    timer.stop();
                    if (reachable) {
                        routingGraphicsItem->set_show_path(true);
                        double length = routingScenario.get_path_length();
                        statusBar()->showMessage(
                            QString("A*: %1 seconds, path length: %2.").arg(timer.time()).arg(length), 8000);
                    } else {
                        statusBar()->showMessage(
                            QString("A*: %1 seconds. Destination node unreachable from start node.").arg(timer.time()),
                            8000);
                        this->actionFindPath->setChecked(false);
                    }
                } else {
                    statusBar()->showMessage(QString("Compute visibility graph or Dijkstra first."), 4000);
                    this->actionFindPath->setChecked(false);
                }
            }
        } else {
            statusBar()->showMessage(QString("Set start/destination node first."), 4000);
            this->actionFindPath->setChecked(false);
        }
    } else {
        routingGraphicsItem->set_show_path(false);
    }
}

void MainWindow::on_actionGenerateRandomDomain_triggered() {
    if (routingScenario.number_of_vertices() > 2) {
        QMessageBox msgBox(QMessageBox::Warning,
                           "Generate new polygonal domain",
                           "Do you really want to clear the current scene?",
                           (QMessageBox::Yes | QMessageBox::No),
                           this);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes) {
            clear();
        } else
            return;
    }
    //this->actionWalk->setChecked(false);

    auto generator = CGAL::Qt::Random_domain_generator<Triangulation>(&routingScenario);
    CGAL::Timer timer;
    timer.start();
    generator.generate_random_domain(this->numberPointsSpinBox->value() * 1000, this->radiusSpinBox->value(),
                       this->thresholdSpinBox->value(), this->erosionSpinBox->value(),
        this->dilationSpinBox->value(), this->minSpinBox->value(),this->checkBox->checkState(),
        this->blueNoise->checkState(), this->candidatesSpinBox->value());
    timer.stop();
    this->actionSetPointOnObstacle->setChecked(true);
    statusBar()->showMessage(
        QString("Generating domain took: %1 seconds. Number of vertices: %2.").arg(timer.time()).arg(
            routingScenario.number_of_vertices()), 8000);
    actionRecenter->trigger();
    routingGraphicsItem->changed();
}

void MainWindow::on_actionTranslateToPoint_triggered() {
    if (routingScenario.start_node_handle != nullptr) {
        routingGraphicsItem->set_focus(routingScenario.start_point);
        statusBar()->showMessage(QString("Computed translation to active_point."), 4000);
    } else {
        statusBar()->showMessage(QString("Set active_point first."), 4000);
    }
    Q_EMIT(transformation());
}

void MainWindow::on_actionShowOrigin_toggled(const bool checked) {
    routingGraphicsItem->set_show_origin(checked);
}

void MainWindow::on_actionRedraw_triggered() {
    routingGraphicsItem->set_approximation_radius(this->approxSpinBox->value());
    routingGraphicsItem->changed();
}

void MainWindow::on_actionWalkPathForward_triggered() {
    if (routingScenario.defined_path) {
        path_position = path_position + this->pathStepSpinBox->value();
        if (path_position >= routingScenario.get_indices_path().size()) {
            path_position = 0;
        }
        routingGraphicsItem->set_focus(routingScenario.get_point_on_path(path_position));
        Q_EMIT(transformation());
    } else {
        statusBar()->showMessage(QString("Compute path first."), 4000);
    }
}

void MainWindow::on_actionWalkPathBackward_triggered() {
    if (routingScenario.defined_path) {
        path_position = path_position - this->pathStepSpinBox->value();
        if (path_position < 0) {
            path_position = routingScenario.get_indices_path().size() - 1;
        }
        routingGraphicsItem->set_focus(routingScenario.get_point_on_path(path_position));
        Q_EMIT(transformation());
    } else {
        statusBar()->showMessage(QString("Compute path first."), 4000);
    }
}

void MainWindow::on_actionTestSuite_triggered() {
    routingScenario.test_suite();
}

void MainWindow::on_actionShowDecomposition_toggled(const bool checked) {
    routingGraphicsItem->set_show_decomposition(checked);
}

void MainWindow::on_actionShowTriangulationBetween_toggled(bool checked) {
    routingGraphicsItem->triangulation_graphics_item->set_show_triangulation_between_obstacles(checked);
}

void MainWindow::on_pathOptimize_released() {
    if(routingScenario.defined_path) {
        routingScenario.path_optimization();
        routingGraphicsItem->repaint();
        double length = routingScenario.get_path_length();
        statusBar()->showMessage(QString("Path length: %1.").arg(length), 4000);
    }
}

#include "main.moc"

int main(int argc, char **argv) {
    QApplication app(argc, argv);

    app.setOrganizationDomain("geometryfactory.com");
    app.setOrganizationName("GeometryFactory");
    app.setApplicationName("Routing on the hyerbolic-plane demo");

    // Import resources from libCGAL (QT5).
    CGAL_QT_INIT_RESOURCES;

    //example independent use
    /*typedef Beltrami_klein_traits K;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Face_base;
    typedef CGAL::Triangulation_vertex_base_2<K> Vertex_base;
    typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
    typedef CGAL::No_constraint_intersection_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> T;
    typedef CGAL::Qt::Routing_scenario<T> RoutingScenario;
    RoutingScenario routing_scenario = RoutingScenario();
    auto rg = CGAL::Qt::Random_domain_generator<T>(&routing_scenario);
    rg.generate_random_domain(20000000, 0.999, 0.01, 5, 5, 5, true, false);
    std::cout << "Number of vertices: " << routing_scenario.number_of_vertices() << std::endl;

    routing_scenario.set_start_point(CGAL::ORIGIN);
    routing_scenario.set_destination_point(Point_2(0, 0.99));

    CGAL::Timer timer;
    timer.start();
    routing_scenario.use_triangulation_as_visibility_graph();
    timer.stop();
    std::cout << "Building visibility graph took: " << timer.time() << " seconds." << std::endl;
    timer.reset();
    //routing_scenario.use_triangulation_as_visibility_graph();
    //routing_scenario.build_visibility_graph_naive();
    timer.start();
    routing_scenario.a_star();
    timer.stop();
    std::cout << "A* took: " << timer.time() << " seconds." << std::endl;
    timer.reset();

    std::cout << "FINISH" << std::endl;*/

    MainWindow mainWindow;
    mainWindow.show();
    mainWindow.on_actionRecenter_triggered();

    return app.exec();
}
