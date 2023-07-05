

#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_INTERFACE_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_INTERFACE_H


#include "ui_Interface.h"
//------------------------------------------------------------------------------
#include "chai3d.h"
//------------------------------------------------------------------------------
#include <QFileSystemModel>
#include <QGLWidget>
#include <QLabel>
#include <QMainWindow>
#include <QMessageBox>
#include <QTimer>
#include <QShortcut>
//------------------------------------------------------------------------------
class ApplicationWidget;
//------------------------------------------------------------------------------

namespace Ui
{
    class InterfaceClass;
}

//------------------------------------------------------------------------------

class Interface : public QMainWindow
{
    Q_OBJECT

    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    Interface (QWidget *parent = 0, Qt::WindowFlags flags = 0);
    ~Interface ();


    //--------------------------------------------------------------------------
    // PRIVATE MEMBERS - UI:
    //--------------------------------------------------------------------------

private:

    Ui::InterfaceClass ui;
    QShortcut *EscKey;
    QShortcut *FKey;
    QShortcut *SKey;
    QShortcut *QKey;
    QTimer *StatusTimer;
    ApplicationWidget *Application;
    QLabel GraphicRate;
    QLabel SimulationRate;
    QLabel HapticRate;
    QFileSystemModel *dirModel;


    //--------------------------------------------------------------------------
    // PRIVATE MEMBERS:
    //--------------------------------------------------------------------------

private:

    int AbortRequest;


    //--------------------------------------------------------------------------
    // PRIVATE SLOTS:
    //--------------------------------------------------------------------------

private slots:

    void EnterFullScreen();
    void ExitFullScreen();
    void ToggleFullScreen();
    void SetFullScreen(bool fs);
    void ToggleSettings();
    void ShowSettings(bool show);
    void UpdateStatus();


    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:

    int  Start();
    void Stop();
};


#endif //IMPLICITNONLINEARCOMPLEMENTARITY_INTERFACE_H
