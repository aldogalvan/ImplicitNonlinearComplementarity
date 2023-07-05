#include <iostream>
#include "chai3d.h"
#include <Eigen/Dense>
#include "src/collision.hpp"
#include "src/implicit_lcp.hpp"
#include "imgui.h"
//------------------------------------------------------------------------------
#include "extras/GLFW/include/GLFW/glfw3.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
//------------------------------------------------------------------------------
using namespace chai3d;
using namespace std;
using namespace Eigen;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// GENERAL SETTINGS
//------------------------------------------------------------------------------

// stereo Mode
/*
    C_STEREO_DISABLED:            Stereo is disabled
    C_STEREO_ACTIVE:              Active stereo for OpenGL NVDIA QUADRO cards
    C_STEREO_PASSIVE_LEFT_RIGHT:  Passive stereo where L/R images are rendered next to each other
    C_STEREO_PASSIVE_TOP_BOTTOM:  Passive stereo where L/R images are rendered above each other
*/
cStereoMode stereoMode = C_STEREO_DISABLED;

enum MouseStates
{
    MOUSE_IDLE,
    MOUSE_MOVE_CAMERA
};

enum HapticStates
{
    HAPTIC_IDLE,
    HAPTIC_SELECTION
};


// fullscreen mode
bool fullscreen = false;

// mirrored display
bool mirroredDisplay = false;


//------------------------------------------------------------------------------
// DECLARED VARIABLES
//------------------------------------------------------------------------------

// the device position
cVector3d device_pos;

// the god object position
cVector3d god_object_pos;

// the haptic device velocity
cVector3d device_velocity;

// the sphere
cMultiMesh* object1;

// the plane
cMultiMesh* object2;

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cSpotLight *light;

// a haptic device handler
cHapticDeviceHandler* handler;

// a pointer to the current haptic device
cGenericHapticDevicePtr hapticDevice;

// a label to display the haptic device model
cLabel* labelHapticDeviceModel;

// a label to display the position [m] of the haptic device
cLabel* labelHapticDevicePosition;

// a global variable to store the position [m] of the haptic device
cVector3d hapticDevicePosition;

// a font for rendering text
cFontPtr font;

// a label to display the rate [Hz] at which the simulation is running
cLabel* labelRates;


// a flag to indicate if the haptic simulation currently running
bool simulationRunning = false;

// a flag to indicate if the haptic simulation has terminated
bool simulationFinished = true;

// a frequency counter to measure the simulation graphic rate
cFrequencyCounter freqCounterGraphics;

// a frequency counter to measure the simulation haptic rate
cFrequencyCounter freqCounterHaptics;

// a frequency counter to measure the simulation haptic rate
cFrequencyCounter freqCounterSimulation;

// mouse state
MouseStates mouseState = MOUSE_IDLE;

// last mouse position
double mouseX, mouseY;

// haptic thread
cThread* hapticsThread;

// simulation thread
cThread* simulationThread;

// a handle to window display context
GLFWwindow* window = NULL;

// current width of window
int width  = 0;

// current height of window
int height = 0;

// swap interval for the display context (vertical synchronization)
int swapInterval = 1;


//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

// callback when the window display is resized
void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height);

// callback when an error GLFW occurs
void errorCallback(int error, const char* a_description);

// callback when a key is pressed
void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);

// callback to handle mouse click
void mouseButtonCallback(GLFWwindow* a_window, int a_button, int a_action, int a_mods);

// callback to handle mouse motion
void mouseMotionCallback(GLFWwindow* a_window, double a_posX, double a_posY);

// callback to handle mouse scroll
void mouseScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY);


// this function renders the scene
void updateGraphics(void);

// this function contains the main haptics simulation loop
void updateHaptics(void);

void updateSimulation(void);

// this function closes the application
void close(void);

void fillEigenMatrices( cMultiMesh* mesh, MatrixXd& vertices, MatrixXi& triangles);

int main(int argc, char* argv[])
{
    //--------------------------------------------------------------------------
    // INITIALIZATION
    //--------------------------------------------------------------------------

    cout << endl;
    cout << "-----------------------------------" << endl;
    cout << "CHAI3D" << endl;
    cout << "Demo: 03-analytics" << endl;
    cout << "Copyright 2003-2016" << endl;
    cout << "-----------------------------------" << endl << endl << endl;
    cout << "Keyboard Options:" << endl << endl;
    cout << "[1] - Enable/Disable potential field" << endl;
    cout << "[2] - Enable/Disable damping" << endl;
    cout << "[f] - Enable/Disable full screen mode" << endl;
    cout << "[m] - Enable/Disable vertical mirroring" << endl;
    cout << "[q] - Exit application" << endl;
    cout << endl << endl;


    //--------------------------------------------------------------------------
    // OPEN GL - WINDOW DISPLAY
    //--------------------------------------------------------------------------

    // initialize GLFW library
    if (!glfwInit())
    {
        cout << "failed initialization" << endl;
        cSleepMs(1000);
        return 1;
    }

    // set error callback
    glfwSetErrorCallback(errorCallback);

    // compute desired size of window
    const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    int w = 0.8 * mode->height;
    int h = 0.5 * mode->height;
    int x = 0.5 * (mode->width - w);
    int y = 0.5 * (mode->height - h);

    // set OpenGL version
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    // set active stereo mode
    if (stereoMode == C_STEREO_ACTIVE)
    {
        glfwWindowHint(GLFW_STEREO, GL_TRUE);
    }
    else
    {
        glfwWindowHint(GLFW_STEREO, GL_FALSE);
    }

    // create display context
    window = glfwCreateWindow(w, h, "CHAI3D", NULL, NULL);
    if (!window)
    {
        cout << "failed to create window" << endl;
        cSleepMs(1000);
        glfwTerminate();
        return 1;
    }

    // get width and height of window
    glfwGetWindowSize(window, &width, &height);

    // set position of window
    glfwSetWindowPos(window, x, y);

    // set key callback
    glfwSetKeyCallback(window, keyCallback);

    // set resize callback
    glfwSetWindowSizeCallback(window, windowSizeCallback);

    // set mouse position callback
    glfwSetCursorPosCallback(window, mouseMotionCallback);

    // set mouse button callback
    glfwSetMouseButtonCallback(window, mouseButtonCallback);

    // set mouse scroll callback
    glfwSetScrollCallback(window, mouseScrollCallback);

    // set current display context
    glfwMakeContextCurrent(window);

    // sets the swap interval for the current display context
    glfwSwapInterval(swapInterval);

#ifdef GLEW_VERSION
    // initialize GLEW library
    if (glewInit() != GLEW_OK)
    {
        cout << "failed to initialize GLEW library" << endl;
        glfwTerminate();
        return 1;
    }
#endif


    // Initialize ImGui
//    const char* glsl_version = "#version 140";
//    IMGUI_CHECKVERSION();
//    ImGui::CreateContext();
//    ImGuiIO& io = ImGui::GetIO(); (void)io;
//    ImGui::StyleColorsDark();
//    // Setup Platform/Renderer bindings
//    ImGui_ImplGlfw_InitForOpenGL(window, true);
//    ImGui_ImplOpenGL3_Init(glsl_version);

    //--------------------------------------------------------------------------
    // WORLD - CAMERA - LIGHTING
    //--------------------------------------------------------------------------

    // create a new world.
    world = new cWorld();

    // set the background color of the environment
    world->m_backgroundColor.setWhite();

    // create a camera and insert it into the virtual world
    camera = new cCamera(world);
    world->addChild(camera);

    // position and orient the camera
    camera->set(cVector3d(0.5, 0.0, 0.0),    // camera position (eye)
                cVector3d(0.0, 0.0, 0.0),    // look at position (target)
                cVector3d(0.0, 0.0, 1.0));   // direction of the (up) vector

    // set the near and far clipping planes of the camera
    camera->setClippingPlanes(0.01, 10.0);

    // set stereo mode
    camera->setStereoMode(stereoMode);

    // set stereo eye separation and focal length (applies only if stereo is enabled)
    camera->setStereoEyeSeparation(0.005);
    camera->setStereoFocalLength(0.5);

    // set vertical mirrored display mode
    camera->setMirrorVertical(mirroredDisplay);

    // create a directional light source
    light = new cSpotLight(world);

    // insert light source inside world
    world->addChild(light);

    // enable light source
    light->setEnabled(true);

    // define direction of light beam
    light->setLocalPos(cVector3d(1,0,0));
    light->setDir(-1, 0, 0.0);

    object1  = new cMultiMesh();
    object1->loadFromFile("/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/model.obj");
    object1->getMesh(0)->m_material->setRed();
    object1->m_material->setShininess(0.75);
    object1->setLocalPos(cVector3d(0.2,0,0.7));
    world->addChild(object1);

    object2 = new cMultiMesh();
    object2->loadFromFile("/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/model.obj");
    object2->setLocalPos(cVector3d(0.0,0.0,0.0));
    object2->m_material->setShininess(0.75);
    object2->getMesh(0)->m_material->setBlueCyan();

    world->addChild(object2);

    world->computeGlobalPositions();

    //--------------------------------------------------------------------------
    // HAPTIC DEVICE
    //--------------------------------------------------------------------------

    // create a haptic device handler
    handler = new cHapticDeviceHandler();

    // get a handle to the first haptic device
    handler->getDevice(hapticDevice, 0);

    // open a connection with the haptic device
    hapticDevice->open();

    // retrieve information about the current haptic device
    cHapticDeviceInfo info = hapticDevice->getSpecifications();

    //--------------------------------------------------------------------------
    // WIDGETS
    //--------------------------------------------------------------------------

    // create a font
    font = NEW_CFONTCALIBRI20();

    // create a label to display the haptic device model
    labelHapticDeviceModel = new cLabel(font);
    camera->m_frontLayer->addChild(labelHapticDeviceModel);
    labelHapticDeviceModel->setText(info.m_modelName);

    // create a label to display the position of haptic device
    labelHapticDevicePosition = new cLabel(font);
    camera->m_frontLayer->addChild(labelHapticDevicePosition);

    // create a label to display the haptic and graphic rate of the simulation
    labelRates = new cLabel(font);
    camera->m_frontLayer->addChild(labelRates);
    labelRates->m_fontColor.setBlack();

    //--------------------------------------------------------------------------
    // START SIMULATION
    //--------------------------------------------------------------------------

    // create a thread which starts the main haptics rendering loop
    hapticsThread = new cThread();
    hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

    // create a thread which starts the main simulation loop
    simulationThread = new cThread();
    simulationThread->start(updateSimulation, CTHREAD_PRIORITY_SIMULATION);

    // setup callback when application exits
    atexit(close);


    //--------------------------------------------------------------------------
    // MAIN GRAPHIC LOOP
    //--------------------------------------------------------------------------

    // call window size callback at initialization
    windowSizeCallback(window, width, height);

    while (!glfwWindowShouldClose(window))
    {
        // get width and height of window
        glfwGetWindowSize(window, &width, &height);

        // render graphics
        updateGraphics();

        // swap buffers
        glfwSwapBuffers(window);

        // process events
        glfwPollEvents();

        // signal frequency counter
        freqCounterGraphics.signal(1);
    }

    // close window
    glfwDestroyWindow(window);

    // terminate GLFW library
    glfwTerminate();

    // exit
    return (0);
}

//------------------------------------------------------------------------------

void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height)
{
    // update window size
    width  = a_width;
    height = a_height;

    // update position of label
    labelHapticDevicePosition->setLocalPos(20, width - 60, 0);

    // update position of label
    labelHapticDeviceModel->setLocalPos(20, height - 40, 0);

}

//------------------------------------------------------------------------------

void errorCallback(int a_error, const char* a_description)
{
    cout << "Error: " << a_description << endl;
}

//------------------------------------------------------------------------------

void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
    // filter calls that only include a key press
    if ((a_action != GLFW_PRESS) && (a_action != GLFW_REPEAT))
    {
        return;
    }

    // option - exit
    if ((a_key == GLFW_KEY_ESCAPE) || (a_key == GLFW_KEY_Q))
    {
        glfwSetWindowShouldClose(a_window, GLFW_TRUE);
    }

    // option - enable/disable force field
    if (a_key == GLFW_KEY_1)
    {

    }

    // option - enable/disable damping
    if (a_key == GLFW_KEY_2)
    {

    }

    // option - toggle fullscreen
    if (a_key == GLFW_KEY_F)
    {
        // toggle state variable
        fullscreen = !fullscreen;

        // get handle to monitor
        GLFWmonitor* monitor = glfwGetPrimaryMonitor();

        // get information about monitor
        const GLFWvidmode* mode = glfwGetVideoMode(monitor);

        // set fullscreen or window mode
        if (fullscreen)
        {
            glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
            glfwSwapInterval(swapInterval);
        }
        else
        {
            int w = 0.8 * mode->height;
            int h = 0.5 * mode->height;
            int x = 0.5 * (mode->width - w);
            int y = 0.5 * (mode->height - h);
            glfwSetWindowMonitor(window, NULL, x, y, w, h, mode->refreshRate);
            glfwSwapInterval(swapInterval);
        }
    }

    // option - toggle vertical mirroring
    if (a_key == GLFW_KEY_M)
    {
        mirroredDisplay = !mirroredDisplay;
        camera->setMirrorVertical(mirroredDisplay);
    }
}


//------------------------------------------------------------------------------

void mouseButtonCallback(GLFWwindow* a_window, int a_button, int a_action, int a_mods)
{
    if (a_button == GLFW_MOUSE_BUTTON_RIGHT && a_action == GLFW_PRESS)
    {
        // store mouse position
        glfwGetCursorPos(window, &mouseX, &mouseY);

        // update mouse state
        mouseState = MOUSE_MOVE_CAMERA;
    }

    else
    {
        // update mouse state
        mouseState = MOUSE_IDLE;
    }
}

//------------------------------------------------------------------------------

void mouseMotionCallback(GLFWwindow* a_window, double a_posX, double a_posY)
{
    if (mouseState == MOUSE_MOVE_CAMERA)
    {
        // compute mouse motion
        int dx = a_posX - mouseX;
        int dy = a_posY - mouseY;
        mouseX = a_posX;
        mouseY = a_posY;

        // compute new camera angles
        double azimuthDeg = camera->getSphericalAzimuthDeg() - 0.5 * dx;
        double polarDeg = camera->getSphericalPolarDeg() - 0.5 * dy;

        // assign new angles
        camera->setSphericalAzimuthDeg(azimuthDeg);
        camera->setSphericalPolarDeg(polarDeg);

    }
}

//------------------------------------------------------------------------------

void mouseScrollCallback(GLFWwindow* a_window, double a_offsetX, double a_offsetY)
{
    double r = camera->getSphericalRadius();
    r = cClamp(r + 0.1 * a_offsetY, 0.5, 3.0);
    camera->setSphericalRadius(r);
}

//------------------------------------------------------------------------------

void close(void)
{
    // stop the simulation
    simulationRunning = false;

    // wait for graphics and haptics loops to terminate
    while (!simulationFinished) { cSleepMs(100); }

    // close haptic device
    hapticDevice->close();

    // delete resources
    delete hapticsThread;
    delete world;
    delete handler;
}

//------------------------------------------------------------------------------

void updateGraphics(void)
{

    /////////////////////////////////////////////////////////////////////
    // UPDATE WIDGETS
    /////////////////////////////////////////////////////////////////////

    // update position data
    labelHapticDevicePosition->setText( hapticDevicePosition.str(3));

    // update haptic and graphic rate data
    labelRates->setText(cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " +
                        cStr(freqCounterSimulation.getFrequency(), 0) + " Hz / " +
                        cStr(freqCounterHaptics.getFrequency(), 0) + " Hz");

    // update position of label
    labelRates->setLocalPos((int)(0.5 * (width - labelRates->getWidth())), 15);


    /////////////////////////////////////////////////////////////////////
    // RENDER SCENE
    /////////////////////////////////////////////////////////////////////

    // update shadow maps (if any)
    world->updateShadowMaps(false, mirroredDisplay);

    // render world
    camera->renderView(width, height);

    // wait until all OpenGL commands are completed
    glFinish();

    // check for any OpenGL errors
    GLenum err;
    err = glGetError();
    if (err != GL_NO_ERROR) cout << "Error:  %s\n" << gluErrorString(err);
}

//------------------------------------------------------------------------------

void updateSimulation(void)
{
    // simulation in now running
    simulationRunning = true;
    simulationFinished = false;

    cPrecisionClock clock;

    MatrixXd object1_vertices;
    MatrixXi object1_triangles;
    fillEigenMatrices(object1, object1_vertices, object1_triangles);

    VectorXd q_constrained_1(7); q_constrained_1 << 0, 0, 0, 1, 0, 0, 0;
    q_constrained_1.block<3,1>(0,0) = device_pos.eigen();
    VectorXd q_unconstrained_1(7); q_unconstrained_1 << 0, 0, 0, 1, 0, 0, 0;
    q_unconstrained_1.block<3,1>(0,0) = device_pos.eigen();
    VectorXd u_constrained_1(6); u_constrained_1.setZero();
    VectorXd u_unconstrained_1(6); u_unconstrained_1.setZero();

    MatrixXd object2_vertices;
    MatrixXi object2_triangles;
    fillEigenMatrices(object2, object2_vertices, object2_triangles);

    VectorXd q_constrained_2(7); q_constrained_2 << 0, 0, 0, 1, 0, 0, 0;
    VectorXd q_unconstrained_2(7); q_unconstrained_2 << 0, 0, 0, 1, 0, 0, 0;
    VectorXd u_constrained_2(6); u_constrained_2.setZero();
    VectorXd u_unconstrained_2(6); u_unconstrained_2.setZero();

    // the matrix
    VectorXd M = VectorXd::Ones(6);

    // main haptic simulation loop
    while (simulationRunning) {

        clock.stop();
        double dt = clock.getCurrentTimeSeconds();
        clock.start(true);

        // compute collisions
        vector<Contact*> collisions;
        if (CollisionDetector::findCollisionsRigidRigid(q_constrained_1, q_unconstrained_1, object1_vertices, object1_triangles,
                                              q_constrained_2, q_unconstrained_2, object2_vertices, object2_triangles,
                                              collisions))
        {
//            object1_velocity_constrained.setZero();
//            object1_position_constrained = object1_position_unconstrained;
//            object1_velocity_unconstrained.setZero();

            ImplicitLCP::setup_implicit_lcp_rigid(q_constrained_1, q_unconstrained_1,
                                                  u_constrained_1, u_unconstrained_1, M, collisions, dt);
        }
        else
        {
            u_constrained_1 = u_unconstrained_1;
            q_constrained_1 = q_unconstrained_1;
        }

        // the the global variable for the god object position
        god_object_pos = cVector3d(q_constrained_1.block<3,1>(0,0));

        // update the device velocity
        q_unconstrained_1.block<3,1>(0,0) = device_pos.eigen();
        u_unconstrained_1.block<3,1>(0,0) = (q_unconstrained_1.block<3,1>(0,0) -
                q_constrained_1.block<3,1>(0,0));


        for (auto ptr : collisions)
        {
            delete [] ptr;
        }

        object1->setLocalPos(q_constrained_1.block<3,1>(0,0));
        object1->setLocalRot(Quaterniond(q_constrained_1(3),q_constrained_1(4),
                                         q_constrained_1(5),q_constrained_1(6)).toRotationMatrix());

        // update frequency counter
        freqCounterSimulation.signal(1);
    }

    // exit haptics thread
    simulationFinished = true;
}
void updateHaptics(void) {

    // simulation in now running
    simulationRunning = true;
    simulationFinished = false;

    cPrecisionClock clock;

    int scale_factor = 5; // the scaling
    cVector3d start_pos = cVector3d(0,0,0.5); // the starting position

    double k = 2000;// the stiffness
    double b = 0.5; //  the damping

    // main haptic simulation loop
    while (simulationRunning) {

        clock.stop();
        double dt = clock.getCurrentTimeSeconds();
        clock.start(true);

        hapticDevice->getPosition(device_pos);
        device_pos *= scale_factor;
        device_pos += start_pos;

        hapticDevice->getLinearVelocity(device_velocity);
        device_velocity *= scale_factor;

        Vector3d force(0,0,0);

//        force = (device_pos.eigen() - god_object_pos.eigen())*k;
        hapticDevice->setForce(force);

        // update frequency counter
        freqCounterHaptics.signal(1);
    }

    // exit haptics thread
    simulationFinished = true;
}

void fillEigenMatrices( cMultiMesh* mesh, MatrixXd& vertices, MatrixXi& triangles) {

    // Get the number of vertices and triangles in the mesh
    int numVertices = mesh->getNumVertices();
    int numTriangles = mesh->getNumTriangles();

    // Resize the Eigen matrices to accommodate the vertices and triangles
    vertices.resize(numVertices, 3);
    triangles.resize(numTriangles, 3);

    // Fill the vertices matrix
    for (int i = 0; i < numVertices; ++i) {

        const cVector3d& vertex = mesh->getVertexPos(i).eigen();
        vertices.row(i) << vertex.x(), vertex.y(), vertex.z();
    }

    // Fill the triangles matrix
    for (int i = 0; i < numTriangles; ++i) {
        cVector3d triangle (mesh->getMesh(0)->m_triangles->getVertexIndex0(i),
                            mesh->getMesh(0)->m_triangles->getVertexIndex1(i),
                            mesh->getMesh(0)->m_triangles->getVertexIndex2(i));
        triangles.row(i) << triangle.x(), triangle.y(), triangle.z();
    }

}