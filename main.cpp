#include <iostream>
#include "chai3d.h"
#include <Eigen/Dense>
//------------------------------------------------------------------------------
#include "extras/GLFW/include/GLFW/glfw3.h"
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

struct collisionInfo
{
    Vector3d collisionPt;
    Vector3d collisionNormal;
    double collisionDepth;
};

struct LCP
{
    MatrixXd A;
    MatrixXd b;
    MatrixXd x;
};

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

// the sphere
cShapeSphere* god_object_sphere;
cVector3d god_object;
cVector3d haptic_device;

// the triangle
cShapeLine* triangle[3];

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cDirectionalLight *light;

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

// a global variable to store the velocity [m/s] of the haptic device
cVector3d hapticDeviceVelocity;

// a font for rendering text
cFontPtr font;

// a label to display the rate [Hz] at which the simulation is running
cLabel* labelRates;

// a line representing the velocity of the haptic device
cShapeLine* velocity;


// a flag to indicate if the haptic simulation currently running
bool simulationRunning = false;

// a flag to indicate if the haptic simulation has terminated
bool simulationFinished = true;

// a frequency counter to measure the simulation graphic rate
cFrequencyCounter freqCounterGraphics;

// a frequency counter to measure the simulation haptic rate
cFrequencyCounter freqCounterHaptics;

// mouse state
MouseStates mouseState = MOUSE_IDLE;

// last mouse position
double mouseX, mouseY;

// haptic thread
cThread* hapticsThread;

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

// this function closes the application
void close(void);

bool sphereTriangleCollision(double sphereRadius, const cVector3d& spherePosition,
                             const cVector3d& triangleVertex1, const cVector3d& triangleVertex2,
                             const cVector3d& triangleVertex3, collisionInfo& colInfo);
bool isPointInsideTriangle(const cVector3d& point, const cVector3d& triangleVertex1,
                           const cVector3d& triangleVertex2, const cVector3d& triangleVertex3);
VectorXd lemkes_algorithm(const MatrixXd& M, const VectorXd& q);
Matrix<double, 3, 8> computeFrictionConeJacobian(const Vector3d& contactNormal);

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


    //--------------------------------------------------------------------------
    // WORLD - CAMERA - LIGHTING
    //--------------------------------------------------------------------------

    // create a new world.
    world = new cWorld();

    // set the background color of the environment
    world->m_backgroundColor.setBlack();

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
    light = new cDirectionalLight(world);

    // insert light source inside world
    world->addChild(light);

    // enable light source
    light->setEnabled(true);

    // define direction of light beam
    light->setDir(-1.0, 0.0, 0.0);

    // draw a triangle
    triangle[0] = new cShapeLine(cVector3d(0,0,0),cVector3d(0,0,0.2));
    triangle[1] = new cShapeLine(cVector3d(0,0,0.2),cVector3d(0,0.2,0.2));
    triangle[2] = new cShapeLine(cVector3d(0,0.2,0.2),cVector3d(0,0,0.0));

    world->addChild(triangle[0]);
    world->addChild(triangle[1]);
    world->addChild(triangle[2]);

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

    //


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

    //--------------------------------------------------------------------------
    // START SIMULATION
    //--------------------------------------------------------------------------

    // create a thread which starts the main haptics rendering loop
    hapticsThread = new cThread();
    hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

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

void updateHaptics(void)
{
    // simulation in now running
    simulationRunning  = true;
    simulationFinished = false;

    // 6 DOF God Object
    god_object_sphere = new cShapeSphere(0.01);
    world->addChild(god_object_sphere);

    cVector3d haptic_device;
    hapticDevice->getPosition(haptic_device);
    cVector3d god_object = haptic_device;

    Vector3d qdot; qdot.setZero();

    double frictionCoeff = 0.5;

    cPrecisionClock clock;

    // main haptic simulation loop
    while(simulationRunning)
    {
        clock.stop();
        double dt = clock.getCurrentTimeSeconds();
        clock.start(true);

        hapticDevice->getPosition(haptic_device);

        cVector3d vertex1(0,0,0); cVector3d vertex2(0,0,0.1); cVector3d vertex3(0,0.1,0.1);

        collisionInfo colInfo;

        // LCP formulation
        if (sphereTriangleCollision(0.01,haptic_device,
                                vertex1,vertex2,vertex3,colInfo))
        {
            LCP lcp;

            int numCollisions = 1;
            int numObjects = 1;
            int numDOF = 3*numObjects;

            // mass matrix
            MatrixXd M(3,3); M.setIdentity();

            // friction cone basis matrix
            MatrixXd D(8,8); D.setIdentity();

            // N matrix
            MatrixXd N(numDOF,numCollisions);
            MatrixXd J_n(3,3); J_n.setIdentity();
            J_n(0,0) = colInfo.collisionNormal(0);
            J_n(1,0) = colInfo.collisionNormal(1);
            J_n(2,0) = colInfo.collisionNormal(2);

            N.block<3,1>(0,0) = J_n.transpose()*colInfo.collisionNormal;

            // B matrix
            MatrixXd B(numDOF,numCollisions*8);
            auto Jb = computeFrictionConeJacobian(colInfo.collisionNormal);
            B.block<3,8>(0,0) = Jb*D;

            // E matrix
            MatrixXd E(8,1); E.setOnes();

            // A Matrix
            MatrixXd A(10,10);
            A.block<1,1>(0,0) = dt *N.transpose()*M.inverse()*N;
            A.block<1,8>(0,1) = dt * N.transpose() * M.inverse() * B;
            A(0,2) = 0 ;

            A.block<8,1>(1,0) = dt * B.transpose() * M.inverse() * N;
            A.block<8,8>(1,1) = dt * B.transpose() * M.inverse() * B;
            A.block<8,1>(1,2) = E;

            A(9,0) = frictionCoeff;
            A.block<1,8>(9,1) = E.transpose();
            A(9,9) = 0;

            // q vector
            VectorXd q(10);
            q.head(1) = N.transpose()*M.inverse()*M*qdot;
            auto test = B.transpose()*M.inverse()*M*qdot;
            q.head<8>(1) = VectorXd(8);//B.transpose()*M.inverse()*M*qdot;
            q(9) = 0;
        }
        else
        {
            god_object = haptic_device;
        }

        god_object_sphere->setLocalPos(god_object);
        hapticDevice->setForce(cVector3d(0,0,0));

        // update frequency counter
        freqCounterHaptics.signal(1);
    }

    // exit haptics thread
    simulationFinished = true;
}

//------------------------------------------------------------------------------

bool sphereTriangleCollision(double sphereRadius, const cVector3d& spherePosition,
                             const cVector3d& triangleVertex1, const cVector3d& triangleVertex2,
                             const cVector3d& triangleVertex3, collisionInfo& contactInfo)
{
    // Compute the normal of the triangle
    cVector3d triangleNormal = cNormalize(cCross(triangleVertex2 - triangleVertex1,
                                                 triangleVertex3 - triangleVertex1));

    // Compute the signed distance from the sphere center to the plane of the triangle
    double signedDistance = cDot(triangleNormal, spherePosition - triangleVertex1);

    // Check if the sphere is above or below the triangle plane
    if (signedDistance > sphereRadius || signedDistance < -sphereRadius)
    {
        // No collision
        return false;
    }

    // Compute the projection of the sphere center onto the triangle plane
    cVector3d projection = spherePosition - signedDistance * triangleNormal;

    // Check if the projection is inside the triangle
    if (!isPointInsideTriangle(projection, triangleVertex1, triangleVertex2, triangleVertex3))
    {
        // No collision
        return false;
    }

    // Compute the contact point and penetration depth
    contactInfo.collisionPt = projection.eigen();
    contactInfo.collisionDepth = sphereRadius - signedDistance;
    contactInfo.collisionNormal = (projection - spherePosition).eigen();
    contactInfo.collisionNormal.normalize();

    // Collision detected
    return true;
}

bool isPointInsideTriangle(const cVector3d& point, const cVector3d& triangleVertex1,
                           const cVector3d& triangleVertex2, const cVector3d& triangleVertex3)
{
    // Compute the barycentric coordinates of the point with respect to the triangle
    cVector3d edge0 = triangleVertex2 - triangleVertex1;
    cVector3d edge1 = triangleVertex3 - triangleVertex1;
    cVector3d edge2 = point - triangleVertex1;
    double dot00 = cDot(edge0, edge0);
    double dot01 = cDot(edge0, edge1);
    double dot02 = cDot(edge0, edge2);
    double dot11 = cDot(edge1, edge1);
    double dot12 = cDot(edge1, edge2);
    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if the point is inside the triangle
    return (u >= 0.0) && (v >= 0.0) && (u + v <= 1.0);
}

VectorXd lemkes_algorithm(const MatrixXd& M, const VectorXd& q) {
    int n = q.size();
    VectorXd z = VectorXd::Zero(n);
    VectorXd w = VectorXd::Zero(n);
    VectorXi basis = VectorXi::LinSpaced(n, n, 2 * n - 1);

    int i = 0;

    while (true) {
        std::cout << i << std::endl;
        int entering_index;
        w.maxCoeff(&entering_index);
        int entering_var = basis(entering_index);

        if (w(entering_index) <= 0)
            break;

        VectorXd d = M.colPivHouseholderQr().solve(-M.col(entering_var));
        int leaving_index = -1;
        double min_ratio = std::numeric_limits<double>::max();
        for (int j = 0; j < n; j++) {
            if (d(j) > 0 && z(j) / d(j) < min_ratio) {
                leaving_index = j;
                min_ratio = z(j) / d(j);
            }
        }
        int leaving_var = basis(leaving_index);

        double theta = std::min(z(leaving_index) / d(leaving_index), 1.0);
        z += theta * d;
        z(leaving_index) = 0;

        double phi = std::min(w(entering_index) / M(leaving_index, entering_var), 1.0);
        w += phi * M.col(entering_var);
        w(entering_index) = 0;

        basis(entering_index) = leaving_var;
        basis(leaving_index) = entering_var;
        i++;
    }

    return z.head(n);
}


Matrix<double, 3, 8> computeFrictionConeJacobian(const Vector3d& contactNormal)
{
    // Compute two orthogonal vectors in the tangential plane
    Vector3d tangent1, tangent2;
    if (std::abs(contactNormal.x()) < std::abs(contactNormal.y()))
    {
        tangent1 << 0, -contactNormal.z(), contactNormal.y();
    }
    else
    {
        tangent1 << -contactNormal.z(), 0, contactNormal.x();
    }
    tangent1.normalize();
    tangent2 = contactNormal.cross(tangent1);

    // Compute the angles for the friction cone vertices
    const int numVertices = 8;
    const double coneAngle = M_PI / 4.0;
    VectorXd angles = VectorXd::LinSpaced(numVertices, 0, 2 * M_PI - coneAngle).array() + coneAngle / 2.0;
    // Construct the Jacobian matrix for the friction cone
    Matrix<double, 3, 8> J_friction;
    for (int i = 0; i < numVertices; ++i)
    {
        double angle = angles(i);
        Vector3d vertexDirection = std::cos(angle) * tangent1 + std::sin(angle) * tangent2;
        J_friction.col(i) = vertexDirection;
    }

    return J_friction;
}
