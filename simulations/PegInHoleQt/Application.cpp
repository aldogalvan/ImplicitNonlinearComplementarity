
#include "Application.h"

//------------------------------------------------------------------------------
#include <QFile>
#include <QString>
#include <QMessageBox>
//------------------------------------------------------------------------------
using namespace std;
using namespace chai3d;
//------------------------------------------------------------------------------

ApplicationWidget::ApplicationWidget (QWidget *parent)
{
    //--------------------------------------------------------------------------
    // INITIALIZATION
    //--------------------------------------------------------------------------

    // initialize variables
    m_parent  = (Interface*)(void*)parent;
    m_running = false;
    m_timer = new QBasicTimer;
    m_mouseMoveCamera = false;


    // reset frequency counters
    m_graphicRate.reset();
    m_hapticRate.reset();
    m_simulationRate.reset();


    //--------------------------------------------------------------------------
    // WORLD - CAMERA - LIGHTING
    //--------------------------------------------------------------------------

    // create a new world.
    m_world = new cWorld();

    // set the background color of the environment
    m_world->m_backgroundColor.setBlack();

    // create a camera and insert it into the virtual world
    m_camera = new cCamera(m_world);
    m_world->addChild(m_camera);

    // define a basis in spherical coordinates for the camera
    m_camera->setSphericalReferences(cVector3d(0,0,0),    // origin
                                     cVector3d(0,0,1),    // zenith direction
                                     cVector3d(1,0,0));   // azimuth direction

    m_camera->setSphericalDeg(4.0,    // spherical coordinate radius
                              0,      // spherical coordinate azimuth angle
                              0);     // spherical coordinate polar angle

    // set the near and far clipping planes of the camera
    m_camera->setClippingPlanes (0.01, 20.0);

    // create a light source
    m_light = new cSpotLight (m_world);

    // add light to camera
    m_camera->addChild(m_light);

    // enable light source
    m_light->setEnabled(true);

    // define the direction of the light beam
    m_light->setDir(-1.0,-1.0,-0.5);


    //--------------------------------------------------------------------------
    // HAPTIC DEVICES / TOOLS
    //--------------------------------------------------------------------------

    // create a haptic device handler
    cHapticDeviceHandler handler;

    // get access to the first available haptic device found
    handler.getDevice (m_device, 0);

    // open a connection with the haptic device
    m_device->open();

    // retrieve information about the current haptic device
    cHapticDeviceInfo info = m_device->getSpecifications();

    // retrieve information about the current haptic device
    cHapticDeviceInfo hapticDeviceInfo = m_device->getSpecifications();

    simulation = new PegInHole(m_world);

    //--------------------------------------------------------------------------
    // WIDGETS
    //--------------------------------------------------------------------------

    // create a background
    cBackground* background = new cBackground();
    m_camera->m_backLayer->addChild(background);

    // set background properties
    background->setCornerColors(cColorf(1.00f, 1.00f, 1.00f),
                                cColorf(1.00f, 1.00f, 1.00f),
                                cColorf(0.85f, 0.85f, 0.85f),
                                cColorf(0.85f, 0.85f, 0.85f));
};

//------------------------------------------------------------------------------

ApplicationWidget::~ApplicationWidget ()
{
    stop();
    delete m_world;
    delete m_timer;
}

//------------------------------------------------------------------------------

bool
ApplicationWidget::loadModel (string filename)
{
    // remove object from the world, so we can now safely modify it
    m_worldLock.acquire();
    m_modelLock.acquire();
    m_modelLock.release();
    m_worldLock.release();

    // reset camera and tool
    m_camera->setSphericalDeg(4.0,    // spherical coordinate radius
                              0,      // spherical coordinate azimuth angle
                              0);     // spherical coordinate polar angle


    // failure
    QMessageBox::warning(this, "ModelViewer", QString("Failed to load model %1").arg(filename.c_str()), QMessageBox::Ok);
    return false;
}

//------------------------------------------------------------------------------

void* ApplicationWidget::hapticThread ()
{
    // acquire run lock
    m_runLock.acquire();

    // update state
    m_running = true;

    while (m_running)
    {
        simulation->update_haptics();

        // update frequency counter
        m_hapticRate.signal(1);
    }

    // disable forces
    m_device->setForceAndTorqueAndGripperForce (cVector3d(0.0, 0.0, 0.0),
                                                cVector3d(0.0, 0.0, 0.0),
                                                0.0);

    // update state
    m_running = false;

    // release run lock
    m_runLock.release();

    // exit thread
    return (NULL);
}

//------------------------------------------------------------------------------

void* ApplicationWidget::simulationThread()
{
    // acquire run lock
    m_runLock.acquire();

    // update state
    m_running = true;

    cPrecisionClock clock;

    while (m_running)
    {

        clock.stop();
        double dt = clock.getCurrentTimeSeconds();
        clock.start(true);

        simulation->step(dt);

        // update frequency counter
        m_simulationRate.signal(1);
    }


    // update state
    m_running = false;

    // release run lock
    m_runLock.release();

    // exit thread
    return (NULL);
}

//------------------------------------------------------------------------------

void ApplicationWidget::initializeGL ()
{
#ifdef GLEW_VERSION
    glewInit ();
#endif

    // enable anti-aliasing
    QGLWidget::setFormat(QGLFormat(QGL::SampleBuffers));
}

//------------------------------------------------------------------------------

void ApplicationWidget::paintGL ()
{
    if (!m_running) return;

    m_worldLock.acquire();

    // render world
    m_camera->renderView(m_width, m_height);

    // wait until all GL commands are completed
    glFinish();

    m_graphicRate.signal(1);

    m_worldLock.release();
}

//------------------------------------------------------------------------------

void ApplicationWidget::resizeGL (int a_width,  int a_height)
{
    m_worldLock.acquire ();

    m_width = a_width;
    m_height = a_height;

    m_worldLock.release ();
}

//------------------------------------------------------------------------------

int ApplicationWidget::start ()
{
    // start graphic rendering
    m_timer->start(25, this);

    // start haptic thread
    m_hapticThread.start (_hapticThread, CTHREAD_PRIORITY_HAPTICS, this);

    // start simulation thread
    m_simulationThread.start(_simulationThread, CTHREAD_PRIORITY_SIMULATION, this);

    return(0);
}

//------------------------------------------------------------------------------

int ApplicationWidget::stop ()
{
    // stop the simulation thread and wait it to join
    m_running = false;
    m_runLock.acquire();
    m_runLock.release();


    // stop graphics
    m_timer->stop();

    return 0;
}

//------------------------------------------------------------------------------

void ApplicationWidget::wheelEvent (QWheelEvent *event)
{
    double radius = m_camera->getSphericalRadius() + (double)(event->delta())*5e-4;
    m_camera->setSphericalRadius(radius);
}

//------------------------------------------------------------------------------

void ApplicationWidget::mousePressEvent(QMouseEvent *event)
{
    m_mouseX = event->pos().x();
    m_mouseY = event->pos().y();
    m_mouseMoveCamera = true;
}

//------------------------------------------------------------------------------

void ApplicationWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (m_mouseMoveCamera)
    {
        int x = event->pos().x();
        int y = event->pos().y();

        // compute mouse motion
        int dx = x - m_mouseX;
        int dy = y - m_mouseY;
        m_mouseX = x;
        m_mouseY = y;

        // compute new camera angles
        double polarDeg = m_camera->getSphericalPolarDeg() -0.5 * dy;
        double azimuthDeg = m_camera->getSphericalAzimuthDeg() - 0.5 * dx;

        // assign new angles
        m_camera->setSphericalPolarDeg(cClamp(polarDeg, 1.0, 179.0));
        m_camera->setSphericalAzimuthDeg(azimuthDeg);

        // line up tool with camera
        // m_tool->setLocalRot(m_camera->getLocalRot());
    }
}

//------------------------------------------------------------------------------

void ApplicationWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_mouseMoveCamera = false;
}

//------------------------------------------------------------------------------

void _hapticThread (void *arg)
{
    ((ApplicationWidget*)arg)->hapticThread();
}

//------------------------------------------------------------------------------

void _simulationThread  (void *arg)
{
    ((ApplicationWidget*)arg)->simulationThread();
}