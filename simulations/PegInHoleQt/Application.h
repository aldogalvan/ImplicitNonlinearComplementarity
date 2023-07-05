
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_APPLICATION_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_APPLICATION_H

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "PegInHole.hpp"
//------------------------------------------------------------------------------
#include "Interface.h"
//------------------------------------------------------------------------------
#include <QBasicTimer>
#include <QWheelEvent>
//------------------------------------------------------------------------------
#include <string>
//------------------------------------------------------------------------------

void _hapticThread (void *arg);
void _simulationThread(void *arg);

//------------------------------------------------------------------------------

class ApplicationWidget : public QGLWidget
{

    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    ApplicationWidget (QWidget *parent);
    virtual ~ApplicationWidget ();


    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:

    bool loadModel(std::string filename);
    int start();
    int stop();
    void waitForStop();
    void* hapticThread();
    void* simulationThread();
    bool isRunning() { return m_running; }

    double getGraphicRate() { return (m_graphicRate.getFrequency()); }
    double getHapticRate() { return  (m_hapticRate.getFrequency()); }
    double getSimulationRate() {return (m_simulationRate.getFrequency()); }


    //--------------------------------------------------------------------------
    // PROTECTED METHODS:
    //--------------------------------------------------------------------------

protected:

    void initializeGL();
    void resizeGL(int a_width, int a_height);
    void paintGL();
    void wheelEvent(QWheelEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void timerEvent(QTimerEvent *event) { updateGL(); }


    //--------------------------------------------------------------------------
    // PUBLIC MEMBERS:
    //--------------------------------------------------------------------------

public:

    // the main simulation
    PegInHole* simulation;

    // application control
    Interface* m_parent;
    chai3d::cMutex m_worldLock;
    chai3d::cMutex m_modelLock;
    chai3d::cMutex m_runLock;
    chai3d::cThread m_hapticThread;
    chai3d::cThread m_simulationThread;
    chai3d::cFrequencyCounter m_graphicRate;
    chai3d::cFrequencyCounter m_hapticRate;
    chai3d::cFrequencyCounter m_simulationRate;

    QBasicTimer *m_timer;
    bool m_running;
    int m_width;
    int m_height;
    int m_mouseX;
    int m_mouseY;
    bool m_mouseMoveCamera;

    // CHAI3D world
    chai3d::cGenericHapticDevicePtr m_device;
    chai3d::cWorld* m_world;
    chai3d::cCamera* m_camera;
    chai3d::cSpotLight* m_light;
};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_APPLICATION_H
