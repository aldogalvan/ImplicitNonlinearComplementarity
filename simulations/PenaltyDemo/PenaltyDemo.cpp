
#include "PenaltyDemo.hpp"

struct TrajectoryResult {
    double position;
    double velocity;
};

CollisionResult checkSensorCollision(Sensor* sensor) {
    CollisionResult result;

    // Check if the start and end positions of the sensor cross the z = 0 plane.
    if ((sensor->globalStartPos.z() >= 0 && sensor->globalNormal.z() <= 0) ||
        (sensor->globalStartPos.z() <= 0 && sensor->globalNormal.z() >= 0)) {
        result.collision = true;
        result.distance = fabs(sensor->globalStartPos.z());
        result.p0 = Eigen::Vector3d(sensor->globalStartPos.x(), sensor->globalStartPos.y(), 0.0);
        result.normal = Vector3d(0,0,1);
        result.sensor = sensor;
    } else {
        result.collision = false;
        result.distance = 0.0;
        result.p0 = Eigen::Vector3d::Zero();
        result.normal = Vector3d{0.0, 0.0, 0.0};
        result.sensor = nullptr;
    }

    return result;
}

TrajectoryResult sinusoidalTrajectory(double amplitude, double frequency, double phase, double timeOffset, double time) {
    TrajectoryResult result;

    // Calculate position as a function of time
    result.position = amplitude * sin(2 * M_PI * frequency * (time - timeOffset) + phase);

    // Calculate velocity as the derivative of position
    // Using the derivative of sin(u) = cos(u) * du/dt
    result.velocity = 2 * M_PI * amplitude * frequency * cos(2 * M_PI * frequency * (time - timeOffset) + phase);

    return result;
}

MatrixXd Q(Quaterniond q)
{
    MatrixXd Q = MatrixXd::Zero(4,3);
    Q(0,0) = 0.5*q.w();
    Q(0,1) = 0.5*q.z();
    Q(0,2) = -0.5*q.y();
    Q(1,0) = -0.5*q.z();
    Q(1,1) = 0.5*q.w();
    Q(1,2) = 0.5*q.x();
    Q(2,0) = 0.5*q.y();
    Q(2,1) = -0.5*q.x();
    Q(2,2) = 0.5*q.w();
    Q(3,0) = -0.5*q.x();
    Q(3,1) = -0.5*q.y();
    Q(3,2) = -0.5*q.z();
    return Q;
}
MatrixXd pQpx()
{
    MatrixXd ret(4,3);
    ret(1,2) = 1;
    ret(2,1) = -1;
    ret(3,0) = -1;
    return 0.5*ret;
}

MatrixXd pQpy()
{
    MatrixXd ret(4,3);
    ret(0,2) = -1;
    ret(2,0) = 1;
    ret(3,1) = -1;
    return 0.5*ret;
}

MatrixXd pQpz()
{
    MatrixXd ret(4,3);
    ret(0,1) = 1;
    ret(1,0) = -1;
    ret(3,2) = -1;
    return 0.5*ret;
}

MatrixXd pQpw()
{
    MatrixXd ret(4,3);
    ret(0,0) = 1;
    ret(1,1) = 1;
    ret(2,2) = 1;
    return 0.5*ret;
}

MatrixXd pRpx(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = q.x(); ret(0,1) = q.y(); ret(0,2) = q.z();
    ret(1,0) = q.y(); ret(1,1) = -q.x(); ret(1,2) = -q.w();
    ret(2,0) = q.z(); ret(2,1) = q.w(); ret(2,2) = -q.x();
    return 2*ret;
}

MatrixXd pRpy(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = -q.y(); ret(0,1) = q.x(); ret(0,2) = q.w();
    ret(1,0) = q.x(); ret(1,1) = q.y(); ret(1,2) = q.z();
    ret(2,0) = -q.w(); ret(2,1) = q.z(); ret(2,2) = -q.y();
    return 2*ret;
}

MatrixXd pRpz(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = -q.z(); ret(0,1) = -q.w(); ret(0,2) = q.x();
    ret(1,0) = q.w(); ret(1,1) = -q.z(); ret(1,2) = q.y();
    ret(2,0) = q.x(); ret(2,1) = q.y(); ret(2,2) = q.z();
    return 2*ret;
}

MatrixXd pRpw(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = q.w(); ret(0,1) = -q.z(); ret(0,2) = q.y();
    ret(1,0) = q.z(); ret(1,1) = q.w(); ret(1,2) = -q.x();
    ret(2,0) = -q.y(); ret(2,1) = q.x(); ret(2,2) = q.w();
    return 2*ret;
}

// Function to check collision between a ray and a triangle defined by its vertices.
CollisionResult checkRayTriangleCollision(const Sensor& ray, const Eigen::Vector3d& vertexA,
                                          const Eigen::Vector3d& vertexB, const Eigen::Vector3d& vertexC) {
    CollisionResult result;

    // Calculate the vectors for the two edges of the triangle.
    Eigen::Vector3d edge1 = vertexB - vertexA;
    Eigen::Vector3d edge2 = vertexC - vertexA;

    // Calculate the normal to the triangle.
    Eigen::Vector3d normal = edge1.cross(edge2);
    normal.normalize();

    // Calculate the direction vector of the ray.
    Eigen::Vector3d dir = ray.normal;
    dir.normalize();

    // Calculate the vector from the ray's origin to one of the triangle's vertices.
    Eigen::Vector3d h = ray.startPos - vertexA;

    // Calculate the determinant to determine if the ray and triangle are parallel.
    double a = edge1.dot(dir.cross(edge2));
    if (a > -1e-6 && a < 1e-6) {
        result.collision = false;  // Ray and triangle are parallel.
        return result;
    }

    // Calculate parameters for barycentric coordinates.
    double f = 1.0 / a;
    double u = f * h.dot(dir.cross(edge2));
    double v = f * dir.dot(h.cross(edge1));

    // Check if the intersection point is within the triangle.
    if (u >= 0.0 && v >= 0.0 && u + v <= 1.0) {
        // Calculate the distance to the intersection point along the ray.
        double t = f * edge2.dot(h);

        // Check if the intersection point is within the ray's length.
        if (t >= 0.0 && t <= ray.len) {
            result.collision = true;
            result.distance = t;
            result.p0 = vertexA + u * edge1 + v * edge2;
        }
    }
}

void PenaltyDemo::initialize()
{
    // create the meshes
    spoon = new PenaltyGodObject(m_devicePtr,true,0,"/home/aldo/ImplicitNonlinearComplementarity/resources/RigidBodyDemo/cube.obj");
    spoon->scaleMesh(1);
    spoon->getMesh(0)->m_material->setBlue(); spoon->getMesh(0)->m_material->setShininess(0.5);
    m_world->addChild(spoon);
    spoon->set_local_pos(Vector3d(0,0,0.2));
    spoon->compute_inertia_matrix();
    spoon->wrapMesh();

//    cMesh* ground =  new cMesh();
//    cCreatePlane(ground,5,5);
//    m_world->addChild(ground);

    // create the  mug
//     mug = new PenaltyObject(false,0,"/home/aldo/ImplicitNonlinearComplementarity/resources/RigidBodyDemo/cube.obj");
//     mug->scaleMesh(2);
//     mug->getMesh(0)->m_material->setRed();  mug->getMesh(0)->m_material->setShininess(0.5);
//     m_world->addChild( mug);
//     mug->compute_inertia_matrix();
//     mug->wrapMesh();
//     mug->computeAllNormals();
//     mug->wrapMesh();
//     mug->setWireMode(false);

    if (m_readTrajectory)
    {
        inputFile = new std::ifstream("recording.txt");
        if (!inputFile->is_open()) {
            std::cerr << "Failed to open the file for reading." << std::endl;
        }
    }

    if (m_recordTrajectory)
    {
        outputFile = new std::ofstream("recording.txt");
        if (!outputFile->is_open()) {
            std::cerr << "Failed to open the file for writing." << std::endl;
            return;
        }
    }

    if (m_recordData)
    {
        dataFile = new std::ofstream("bdf2.txt");
        if (!dataFile->is_open()) {
            std::cerr << "Failed to open the file for writing." << std::endl;
            return;
        }
    }
}

void PenaltyDemo::updateHaptics(Vector3d& f)
{
    spoon->updateFromDevice();
    f = Vector3d(0,0,0);
}

void PenaltyDemo::step(double h)
{

    // update the  mug

//    VectorXd spoon_start = spoon->get_configuration();
//    VectorXd spoon_end = spoon->get_configuration();
//    const MatrixXd& spoon_vertices = spoon->vertices();
//    const MatrixXi& spoon_triangles = spoon->triangles();
//    VectorXd  mug_start =  mug->get_configuration();
//    VectorXd  mug_end =  mug->get_configuration();
//    const MatrixXd&  mug_vertices =  mug->vertices();
//    const MatrixXi&  mug_triangles =  mug->triangles();


    vector<Collision> potentialCollisions;
//    if (1)
//    {
//        CollisionDetector::findCollisionsBroadPhase(spoon_start, spoon_end, spoon_vertices, spoon_triangles,
//                                                    mug_start, mug_end, mug_vertices, mug_triangles,
//                                                    potentialCollisions);
//    }

    linearizedBDF1Solver(h,potentialCollisions);
    h_minus_1 = h;
}

void PenaltyDemo::updateGraphics()
{
    spoon->update_mesh_position();
//    mug->update_mesh_position();
}

void PenaltyDemo::linearizedBDF1Solver(double h, vector<Collision> potentialCollisions)
{
    VectorXd x_d = spoon->x_d; Quaterniond q_d = spoon->q_d; q_d.normalize();
    VectorXd v_d = spoon->xdot_d; VectorXd omega_d = spoon->omega_d;

    cPrecisionClock clock;
    clock.start(true);

    if (m_readTrajectory)
    {
        auto state = sinusoidalTrajectory(0.1,0.7,0,0,simTime);
        x_d(2) = state.position;
        x_d(1) = state.position;
        x_d(2) = state.velocity;
        x_d(1) = state.velocity;
    }

    MatrixXd I = MatrixXd::Identity(3,3);
    Vector3d x_c = spoon->x; Quaterniond q_c = spoon->q; q_c.normalize();
    Vector3d v_c = spoon->xdot; Vector3d omega_c = spoon->omega;
    double kc = spoon->Kc; double bc = spoon->Bc;
    double kt = spoon->Kt; double bt = spoon->Bt;
    double m = spoon->mass;
    MatrixXd J_s = spoon->Ibody; MatrixXd Jinv_s = spoon->IbodyInv;
    Vector3d c = VectorXd::Zero(3); Quaterniond q = Quaterniond::Identity();
    MatrixXd L = q_c*J_s*q_c.inverse()*omega_c; Matrix3d R = q_c.toRotationMatrix();
    Quaterniond dq = q_d*q.inverse()*q_c.inverse();
    Vector3d u_c = 2*acos(q.w())*dq.vec();
    Quaterniond q_inv = q.inverse();
    MatrixXd C = Vector4d(dq.x(),dq.y(),dq.z(),dq.w())*RowVector4d(q_inv.x(),q_inv.y(),q_inv.z(),q_inv.w());

    ////////// Coupling Jacobian /////////////////////////////////
    MatrixXd J_c = MatrixXd::Zero(13,13);

    // the virtual coupling force
    Vector3d F_c = kc*(x_d-x_c) + bc*(v_d - v_c); //! CHECK!
    VectorXd T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c); //! CHECK

    // Derivative of angular velocity wrt quaternion
    VectorXd pwpqx = (pRpx(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpx(q_c).transpose())*L; //! CHECK
    VectorXd pwpqy = (pRpy(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpy(q_c).transpose())*L; //! CHECK
    VectorXd pwpqz = (pRpz(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpz(q_c).transpose())*L; //! CHECK
    VectorXd pwpqw = (pRpw(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpw(q_c).transpose())*L; //! CHECK

    // derivatives of force with respect to quaternion
    VectorXd pFcpqx= -kc*pRpx(q_c)*c + bc * vec2skew(c)*pwpqx; //! CHECK!
    VectorXd pFcpqy= -kc*pRpy(q_c)*c + bc * vec2skew(c)*pwpqy; //! CHECK!
    VectorXd pFcpqz= -kc*pRpz(q_c)*c + bc * vec2skew(c)*pwpqz; //! CHECK!
    VectorXd pFcpqw= -kc*pRpw(q_c)*c + bc * vec2skew(c)*pwpqw; //! CHECK!

    // Derivative of linear velocity with respect to quaternion
    MatrixXd pupq = 2*acos(dq.w())*C.block<3,4>(0,0) - 2/(sqrt(1-dq.w()*dq.w()))*dq.vec()*C.block<1,4>(3,0); //! CHECK!
    if (pupq.hasNaN())
        pupq.setZero();

    // Elements of Jacobian
    MatrixXd pFcpx = -kc*MatrixXd::Identity(3,3); //! CHECK!
    MatrixXd pTcpx = vec2skew(R*c)*pFcpx; //! CHECK!
    MatrixXd pFcpq(3,4);
    pFcpq.col(0) = pFcpqx; //! CHECK!
    pFcpq.col(1) = pFcpqy; //! CHECK!
    pFcpq.col(2) = pFcpqz; //! CHECK!
    pFcpq.col(3) = pFcpqw; //! CHECK!
    MatrixXd pTcpq(3,4);

    pTcpq.col(0) = vec2skew(R*c)*pFcpqx - vec2skew(F_c)*pRpx(q_c)*c + kt*pupq.col(0) - bt*pwpqx; //! Check!
    pTcpq.col(1) = vec2skew(R*c)*pFcpqy - vec2skew(F_c)*pRpy(q_c)*c + kt*pupq.col(1) - bt*pwpqy; //! Check!
    pTcpq.col(2) = vec2skew(R*c)*pFcpqz - vec2skew(F_c)*pRpz(q_c)*c + kt*pupq.col(2) - bt*pwpqz; //! Check!
    pTcpq.col(3) = vec2skew(R*c)*pFcpqw - vec2skew(F_c)*pRpw(q_c)*c + kt*pupq.col(3) - bt*pwpqw; //! Check!
    MatrixXd pFcpP = spoon->Bc*I/m; //! Check!
    MatrixXd pTcpP = vec2skew(R*c)*pFcpP; //! Check!
    MatrixXd pFcpL = bc*vec2skew(c)*R*(MatrixXd::Identity(3,3)/m)*R.transpose();
    MatrixXd pTcpL = (bc*vec2skew(R*c)*vec2skew(c) - bt*I)*R*(I/m)*R.transpose();

    J_c.block<3,3>(7,0) = pFcpx;
    J_c.block<3,3>(10,0) = pTcpx;
    J_c.block<3,4>(7,3) = pFcpq;
    J_c.block<3,4>(10,3) = pTcpq;
    J_c.block<3,3>(7,7) = pFcpP;
    J_c.block<3,3>(10,7) = pTcpP;
    J_c.block<3,3>(7,10) = pFcpL;
    J_c.block<3,3>(10,10) = pTcpL;

    /////////////////////  END OF COUPLING JACOBIAN ////////////////////

    ///////////////////// PENALTY JACOBIAN ////////////////////////////

    // create the jacobians and force vectors
    MatrixXd J_p = MatrixXd::Zero(13,13);
    VectorXd F_p = VectorXd::Zero(3); VectorXd T_p = VectorXd::Zero(3);

    vector<CollisionResult> collisions = computeCollisions(potentialCollisions);

    if (0)
    {

        for (int i = 0; i < collisions.size(); i++)
        {
            double k = collisions[i].sensor->k;
            double b = collisions[i].sensor->b;
            Vector3d p0 = collisions[i].p0;
            Vector3d r = collisions[i].sensor->startPos;
            Vector3d n = collisions[i].normal;
            Matrix3d N = n * n.transpose();

            Vector3d F_p_ = -k*N*(x_c + R*r - p0) - k* - b*N*(v_c + omega_c.cross(r));
            Vector3d T_p_ = Vector3d::Zero();
            F_p += F_p_; T_p += T_p_;

            MatrixXd pFppx = -k*N;
            MatrixXd pTppx = vec2skew(R*r)*pFppx;
            MatrixXd pFppq = MatrixXd::Zero(3,4);
            pFppq.col(0) = -k*N*pRpx(q_c)*r - b*N*vec2skew(omega_c)*pRpx(q_c)*r + b*N*vec2skew(R*r)*pwpqx;
            pFppq.col(1) = -k*N*pRpy(q_c)*r - b*N*vec2skew(omega_c)*pRpy(q_c)*r + b*N*vec2skew(R*r)*pwpqy;
            pFppq.col(2) = -k*N*pRpz(q_c)*r - b*N*vec2skew(omega_c)*pRpz(q_c)*r + b*N*vec2skew(R*r)*pwpqz;
            pFppq.col(3) = -k*N*pRpw(q_c)*r - b*N*vec2skew(omega_c)*pRpw(q_c)*r + b*N*vec2skew(R*r)*pwpqw;
            MatrixXd pTppq = MatrixXd::Zero(3,4);
            pTppq.col(0) = vec2skew(R*r)*pFppq.col(0) - vec2skew(F_p_) * pRpx(q_c)*r;
            pTppq.col(1) = vec2skew(R*r)*pFppq.col(1) - vec2skew(F_p_) * pRpy(q_c)*r;
            pTppq.col(2) = vec2skew(R*r)*pFppq.col(2) - vec2skew(F_p_) * pRpz(q_c)*r;
            pTppq.col(3) = vec2skew(R*r)*pFppq.col(3) - vec2skew(F_p_) * pRpw(q_c)*r;
            MatrixXd pFppP = -(b/m)*N;
            MatrixXd pTppP = vec2skew(R*r)*pFppP;
            MatrixXd pFppL = b*N* vec2skew(R*r)*R*(I/m)*R.transpose();
            MatrixXd pTppL = vec2skew(R*r)*pFppL;

            J_p.block<3,3>(7,0) += pFppx;
            J_p.block<3,3>(10,0) += pTppx;
            J_p.block<3,4>(7,3) += pFppq;
            J_p.block<3,4>(10,3) += pTppq;
            J_p.block<3,3>(7,7) += pFppP;
            J_p.block<3,3>(10,7) += pTppP;
            J_p.block<3,3>(7,10) += pFppL;
            J_p.block<3,3>(10,10) += pTppL;

            cout <<" HEY " << flush;
        }
    }


    ///////////////////// END OF PENALTY JACOBIAN ////////////////////////////

    ////////// RIGID BODY JACOBIAN /////////////////////////////////

    MatrixXd Q_ = Q(q_c);
    MatrixXd pQpx_ = pQpx();
    MatrixXd pQpy_ = pQpy();
    MatrixXd pQpz_ = pQpz();
    MatrixXd pQpw_ = pQpw();

    MatrixXd pqdotpL = Q_*R*(I/m)*R.transpose();
    MatrixXd pqdotpq(4,4);
    pqdotpq.col(0) = pQpx_*omega_c + Q_*pwpqx;
    pqdotpq.col(1) = pQpy_*omega_c + Q_*pwpqy;
    pqdotpq.col(2) = pQpz_*omega_c + Q_*pwpqz;
    pqdotpq.col(3) = pQpw_*omega_c + Q_*pwpqw;

    MatrixXd J_r = MatrixXd::Zero(13,13);
    J_r.block<3,3>(0,10) = I / m;
    J_r.block<4,4>(3,3) = pqdotpq;
    J_r.block<4,3>(3,10) = pqdotpL;

    ////////// END RIGID BODY JACOBIAN /////////////////////////////////


    ///////////////////// LINEARIZED IMPLICIT SOLVER ////////////////////////////

    MatrixXd A = MatrixXd::Zero(13,13);
    A = MatrixXd::Identity(13,13) - h*(J_r + J_c + J_p);

    VectorXd b = VectorXd::Zero(13);
    b.block<3,1>(0,0) = h*v_c;
    Quaterniond temp = Quaterniond(omega_c(0),omega_c(1),omega_c(2),0)*q_c;
    b.block<4,1>(3,0) = h*0.5*Vector4d(temp.x(),temp.y(),temp.z(),temp.w());
    b.block<3,1>(7,0) = h*(F_c + F_p);
    b.block<3,1>(10,0) = h*(T_c + T_p);

    VectorXd y = A.colPivHouseholderQr().solve(b);

    // Peform integration
    x_c = y.block<3,1>(0,0) + x_c;
    q_c = sumQuaternions(Quaterniond(y(3),y(4),y(5),y(6)), q_c);
    q_c.normalize();
    v_c = y.block<3,1>(7,0)/m + v_c;
    omega_c = spoon->IbodyInv*y.block<3,1>(10,0) + omega_c;

    // update new states
    spoon->x = x_c;
    spoon->q = q_c;
    spoon->q.normalize();
    spoon->update_inertia_matrix();
    spoon->xdot = v_c;
    spoon->omega = omega_c;

    F_c = kc*(x_d-x_c) + bc*(v_d - v_c); //! CHECK!
    T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c); //! CHECK

    if (m_recordData)
    {
        if (!dataFile->is_open())
        {
            std::cerr << "Failed to open the file for writing." << std::endl;
        }

        *dataFile << simTime << " " << x_d(0) << " " << x_d(1) << " " << x_d(2) << " " <<
                  q_d.x() << " " << q_d.y() << " " << q_d.z() << " " << q_d.w() << " " <<
                  x_c(0) << " " << x_c(1) << " " << x_c(2) << " " <<
                  q_c.x() << " " << q_c.y() << " " << q_c.z() << " " << q_c.w() << " " <<
                  F_c(0) << " " << F_c(1) << " " << F_c(2) << " " <<
                  T_c(0) << " " << T_c(1) << " " << T_c(2) << endl;

    }


    simTime += h;
}

void PenaltyDemo::BDF1Solver(double h, vector<Collision> potentialCollisions)
{
    VectorXd x_d = spoon->x_d; Quaterniond q_d = spoon->q_d; q_d.normalize();
    VectorXd v_d = spoon->xdot_d; VectorXd omega_d = spoon->omega_d;

    if (m_readTrajectory)
    {
        auto state = sinusoidalTrajectory(0.1,0.7,0,0,simTime);
        x_d(2) = state.position;
        x_d(1) = state.position;
        x_d(2) = state.velocity;
        x_d(1) = state.velocity;
    }

    MatrixXd I = MatrixXd::Identity(3,3);
    Vector3d x_c = spoon->x; Quaterniond q_c = spoon->q; q_c.normalize();
    Vector3d v_c = spoon->xdot; Vector3d omega_c = spoon->omega;
    double kc = spoon->Kc; double bc = spoon->Bc;
    double kt = spoon->Kt; double bt = spoon->Bt;
    double m = spoon->mass;
    MatrixXd J_s = spoon->Ibody; MatrixXd Jinv_s = spoon->IbodyInv;
    Vector3d c = VectorXd::Zero(3); Quaterniond q = Quaterniond::Identity();
    MatrixXd L = q_c*J_s*q_c.inverse()*omega_c; Matrix3d R = q_c.toRotationMatrix();
    Quaterniond dq = q_d*q_c.inverse()*q_c.inverse()*q;
    Vector3d u_c = 2*acos(q.w())*dq.vec();
    Quaterniond q_inv = q.inverse();
    MatrixXd C = Vector4d(dq.x(),dq.y(),dq.z(),dq.w())*RowVector4d(q_inv.x(),q_inv.y(),q_inv.z(),q_inv.w());


    // the virtual coupling force
    Vector3d F_c = kc*(x_d-x_c) + bc*(v_d - v_c); //! CHECK!
    VectorXd T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c); //! CHECK

    int it = 0;
    while(it < maxIt)
    {
        ////////// Coupling Jacobian /////////////////////////////////
        MatrixXd J_c = MatrixXd::Zero(13,13);

        // Derivative of angular velocity wrt quaternion
        VectorXd pwpqx = (pRpx(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpx(q_c).transpose())*L; //! CHECK
        VectorXd pwpqy = (pRpy(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpy(q_c).transpose())*L; //! CHECK
        VectorXd pwpqz = (pRpz(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpz(q_c).transpose())*L; //! CHECK
        VectorXd pwpqw = (pRpw(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpw(q_c).transpose())*L; //! CHECK

        // derivatives of force with respect to quaternion
        VectorXd pFcpqx= -kc*pRpx(q_c)*c + bc * vec2skew(c)*pwpqx; //! CHECK!
        VectorXd pFcpqy= -kc*pRpy(q_c)*c + bc * vec2skew(c)*pwpqy; //! CHECK!
        VectorXd pFcpqz= -kc*pRpz(q_c)*c + bc * vec2skew(c)*pwpqz; //! CHECK!
        VectorXd pFcpqw= -kc*pRpw(q_c)*c + bc * vec2skew(c)*pwpqw; //! CHECK!

        // Derivative of linear velocity with respect to quaternion
        MatrixXd pupq = 2*acos(dq.w())*C.block<3,4>(0,0) - 2/(sqrt(1-dq.w()*dq.w()))*dq.vec()*C.block<1,4>(3,0); //! CHECK!
        if (pupq.hasNaN())
            pupq.setZero();

        // Elements of Jacobian
        MatrixXd pFcpx = -kc*MatrixXd::Identity(3,3); //! CHECK!
        MatrixXd pTcpx = vec2skew(R*c)*pFcpx; //! CHECK!
        MatrixXd pFcpq(3,4);
        pFcpq.col(0) = pFcpqx; //! CHECK!
        pFcpq.col(1) = pFcpqy; //! CHECK!
        pFcpq.col(2) = pFcpqz; //! CHECK!
        pFcpq.col(3) = pFcpqw; //! CHECK!
        MatrixXd pTcpq(3,4);

        pTcpq.col(0) = vec2skew(R*c)*pFcpqx - vec2skew(F_c)*pRpx(q_c)*c + kt*pupq.col(0) - bt*pwpqx; //! Check!
        pTcpq.col(1) = vec2skew(R*c)*pFcpqy - vec2skew(F_c)*pRpy(q_c)*c + kt*pupq.col(1) - bt*pwpqy; //! Check!
        pTcpq.col(2) = vec2skew(R*c)*pFcpqz - vec2skew(F_c)*pRpz(q_c)*c + kt*pupq.col(2) - bt*pwpqz; //! Check!
        pTcpq.col(3) = vec2skew(R*c)*pFcpqw - vec2skew(F_c)*pRpw(q_c)*c + kt*pupq.col(3) - bt*pwpqw; //! Check!
        MatrixXd pFcpP = spoon->Bc*I/m; //! Check!
        MatrixXd pTcpP = vec2skew(R*c)*pFcpP; //! Check!
        MatrixXd pFcpL = bc*vec2skew(c)*R*(MatrixXd::Identity(3,3)/m)*R.transpose();
        MatrixXd pTcpL = (bc*vec2skew(R*c)*vec2skew(c) - bt*I)*R*(I/m)*R.transpose();

        J_c.block<3,3>(7,0) = pFcpx;
        J_c.block<3,3>(10,0) = pTcpx;
        J_c.block<3,4>(7,3) = pFcpq;
        J_c.block<3,4>(10,3) = pTcpq;
        J_c.block<3,3>(7,7) = pFcpP;
        J_c.block<3,3>(10,7) = pTcpP;
        J_c.block<3,3>(7,10) = pFcpL;
        J_c.block<3,3>(10,10) = pTcpL;

        /////////////////////  END OF COUPLING JACOBIAN ////////////////////

        ///////////////////// PENALTY JACOBIAN ////////////////////////////

        // create the jacobians and force vectors
        //    VectorXd F_p = VectorXd::Zero(3); VectorXd T_p = VectorXd::Zero(4);
        //
        //    // Penalty force
        //    MatrixXd F_p  =
        //
        //    MatrixXd J_p = ;

        ///////////////////// END OF PENALTY JACOBIAN ////////////////////////////

        ////////// RIGID BODY JACOBIAN /////////////////////////////////

        MatrixXd Q_ = Q(q_c);
        MatrixXd pQpx_ = pQpx();
        MatrixXd pQpy_ = pQpy();
        MatrixXd pQpz_ = pQpz();
        MatrixXd pQpw_ = pQpw();

        MatrixXd pqdotpL = Q_*R*(I/m)*R.transpose();
        MatrixXd pqdotpq(4,4);
        pqdotpq.col(0) = pQpx_*omega_c + Q_*pwpqx;
        pqdotpq.col(1) = pQpy_*omega_c + Q_*pwpqy;
        pqdotpq.col(2) = pQpz_*omega_c + Q_*pwpqz;
        pqdotpq.col(3) = pQpw_*omega_c + Q_*pwpqw;

        MatrixXd J_r = MatrixXd::Zero(13,13);
        J_r.block<3,3>(0,10) = I / m;
        J_r.block<4,4>(3,3) = pqdotpq;
        J_r.block<4,3>(3,10) = pqdotpL;

        ////////// END RIGID BODY JACOBIAN /////////////////////////////////


        ///////////////////// LINEARIZED IMPLICIT SOLVER ////////////////////////////

        MatrixXd A = MatrixXd::Zero(13,13);
        A = MatrixXd::Identity(13,13) - h*(J_r + J_c);

        VectorXd b = VectorXd::Zero(13);
        b.block<3,1>(0,0) = h*v_c;
        Quaterniond temp = Quaterniond(omega_c(0),omega_c(1),omega_c(2),0)*q_c;
        b.block<4,1>(3,0) = h*0.5*Vector4d(temp.x(),temp.y(),temp.z(),temp.w());
        b.block<3,1>(7,0) = h*F_c;
        b.block<3,1>(10,0) = h*T_c;

        VectorXd y = A.colPivHouseholderQr().solve(b);

        // Peform integration
        x_c = 0.75*y.block<3,1>(0,0) + x_c;
        q_c = sumQuaternions(Quaterniond(y(3),y(4),y(5),y(6)), q_c);
        q_c.normalize();
        v_c = 0.75*y.block<3,1>(7,0)/m + v_c;
        omega_c = q_c*J_s*q_c.inverse()*y.block<3,1>(10,0) + omega_c;


        // Update the constants
        L = q_c*J_s*q_c.inverse()*omega_c; R = q_c.toRotationMatrix();
        dq = q_d*q_c.inverse()*q.inverse();
        u_c = 2*acos(q.w())*dq.vec();
        C = Vector4d(dq.x(),dq.y(),dq.z(),dq.w())*RowVector4d(q_inv.x(),q_inv.y(),q_inv.z(),q_inv.w());
        F_c = kc*(x_d - x_c) + bc*(v_d - v_c);
        T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c);

        if (y.norm() < 1e-3)
            break;

        it++;

    }

    spoon->x = x_c;
    spoon->q = q_c;
    spoon->q.normalize();
    spoon->update_inertia_matrix();
    spoon->xdot = v_c;
    spoon->omega = omega_c;


    simTime += h;

}
void PenaltyDemo::linearizedBDF2Solver(double h, vector<Collision> potentialCollisions)
{

    MatrixXd J_s = spoon->Ibody; MatrixXd Jinv_s = spoon->IbodyInv;
    VectorXd x_d = spoon->x_d; Quaterniond q_d = spoon->q_d; q_d.normalize();
    VectorXd v_d = spoon->xdot_d; VectorXd omega_d = spoon->omega_d;

    if (m_readTrajectory)
    {
        auto state = sinusoidalTrajectory(0.1,0.7,0,0,simTime);
        x_d(2) = state.position;
        x_d(1) = state.position;
        x_d(2) = state.velocity;
        x_d(1) = state.velocity;
    }

    double phi = h / h_minus_1;
    double k1 = (1+phi)*(1+phi)/(1+2*phi);
    double k2 = phi*phi/(1+2*phi);
    double k3 = h*(1+phi)/(1+2*phi);

    MatrixXd I = MatrixXd::Identity(3,3);
    VectorXd x_c_minus_1 = spoon->x_minus_1; Quaterniond q_c_minus_1 = spoon->q_minus_1;
    VectorXd x_c = spoon->x; Quaterniond q_c = spoon->q; q_c.normalize();
    VectorXd v_c_minus_1 = spoon->xdot_minus_1; VectorXd omega_c_minus_1 = spoon->omega_minus_1;
    VectorXd v_c = spoon->xdot; VectorXd omega_c = spoon->omega;
    double kc = spoon->Kc; double bc = spoon->Bc;
    double kt = spoon->Kt; double bt = spoon->Bt;
    double m = spoon->mass;
    Vector3d c = VectorXd::Zero(3); Quaterniond q = Quaterniond::Identity();
    MatrixXd L = q_c*J_s*q_c.inverse()*omega_c; Matrix3d R = q_c.toRotationMatrix();
    Quaterniond dq = q_d*q_c.inverse()*q.inverse();
    Vector3d u_c = 2*acos(q.w())*dq.vec();
    Quaterniond q_inv = q.inverse();
    MatrixXd C = Vector4d(dq.x(),dq.y(),dq.z(),dq.w())*RowVector4d(q_inv.x(),q_inv.y(),q_inv.z(),q_inv.w());

    // the virtual coupling force
    Vector3d F_c = kc * ( x_d - x_c ) + bc * ( v_d - v_c ); //! CHECK!
    VectorXd T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c); //! CHECK

    ////////// Coupling Jacobian /////////////////////////////////
    MatrixXd J_c = MatrixXd::Zero(13,13);

    // Derivative of angular velocity wrt quaternion
    VectorXd pwpqx = (pRpx(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpx(q_c).transpose())*L; //! CHECK
    VectorXd pwpqy = (pRpy(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpy(q_c).transpose())*L; //! CHECK
    VectorXd pwpqz = (pRpz(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpz(q_c).transpose())*L; //! CHECK
    VectorXd pwpqw = (pRpw(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpw(q_c).transpose())*L; //! CHECK

    // derivatives of force with respect to quaternion
    VectorXd pFcpqx= -kc*pRpx(q_c)*c + bc * vec2skew(c)*pwpqx; //! CHECK!
    VectorXd pFcpqy= -kc*pRpy(q_c)*c + bc * vec2skew(c)*pwpqy; //! CHECK!
    VectorXd pFcpqz= -kc*pRpz(q_c)*c + bc * vec2skew(c)*pwpqz; //! CHECK!
    VectorXd pFcpqw= -kc*pRpw(q_c)*c + bc * vec2skew(c)*pwpqw; //! CHECK!

    // Derivative of linear velocity with respect to quaternion
    MatrixXd pupq = 2*acos(dq.w())*C.block<3,4>(0,0) - 2/(sqrt(1-dq.w()*dq.w()))*dq.vec()*C.block<1,4>(3,0); //! CHECK!
    if (pupq.hasNaN())
        pupq.setZero();

    // Elements of Jacobian
    MatrixXd pFcpx = -kc*MatrixXd::Identity(3,3); //! CHECK!
    MatrixXd pTcpx = vec2skew(R*c)*pFcpx; //! CHECK!
    MatrixXd pFcpq(3,4);
    pFcpq.col(0) = pFcpqx; //! CHECK!
    pFcpq.col(1) = pFcpqy; //! CHECK!
    pFcpq.col(2) = pFcpqz; //! CHECK!
    pFcpq.col(3) = pFcpqw; //! CHECK!
    MatrixXd pTcpq(3,4);

    pTcpq.col(0) = vec2skew(R*c)*pFcpqx - vec2skew(F_c)*pRpx(q_c)*c + kt*pupq.col(0) - bt*pwpqx; //! Check!
    pTcpq.col(1) = vec2skew(R*c)*pFcpqy - vec2skew(F_c)*pRpy(q_c)*c + kt*pupq.col(1) - bt*pwpqy; //! Check!
    pTcpq.col(2) = vec2skew(R*c)*pFcpqz - vec2skew(F_c)*pRpz(q_c)*c + kt*pupq.col(2) - bt*pwpqz; //! Check!
    pTcpq.col(3) = vec2skew(R*c)*pFcpqw - vec2skew(F_c)*pRpw(q_c)*c + kt*pupq.col(3) - bt*pwpqw; //! Check!
    MatrixXd pFcpP = spoon->Bc*I/m; //! Check!
    MatrixXd pTcpP = vec2skew(R*c)*pFcpP; //! Check!
    MatrixXd pFcpL = bc*vec2skew(c)*R*(MatrixXd::Identity(3,3)/m)*R.transpose();
    MatrixXd pTcpL = (bc*vec2skew(R*c)*vec2skew(c) - bt*I)*R*(I/m)*R.transpose();

    J_c.block<3,3>(7,0)  = pFcpx;
    J_c.block<3,3>(10,0) = pTcpx;
    J_c.block<3,4>(7,3)  = pFcpq;
    J_c.block<3,4>(10,3) = pTcpq;
    J_c.block<3,3>(7,7)  = pFcpP;
    J_c.block<3,3>(10,7) = pTcpP;
    J_c.block<3,3>(7,10) = pFcpL;
    J_c.block<3,3>(10,10) = pTcpL;

    /////////////////////  END OF COUPLING JACOBIAN ////////////////////

    ///////////////////// PENALTY JACOBIAN ////////////////////////////

    // create the jacobians and force vectors
//    VectorXd F_p = VectorXd::Zero(3); VectorXd T_p = VectorXd::Zero(4);
//
//    // Penalty force
//    MatrixXd F_p  =
//
//    MatrixXd J_p = ;

    ///////////////////// END OF PENALTY JACOBIAN ////////////////////////////

    ////////// RIGID BODY JACOBIAN /////////////////////////////////

    MatrixXd Q_ = Q(q_c);
    MatrixXd pQpx_ = pQpx();
    MatrixXd pQpy_ = pQpy();
    MatrixXd pQpz_ = pQpz();
    MatrixXd pQpw_ = pQpw();

    MatrixXd pqdotpL = Q_*R*(I/m)*R.transpose();
    MatrixXd pqdotpq(4,4);
    pqdotpq.col(0) = pQpx_*omega_c + Q_*pwpqx;
    pqdotpq.col(1) = pQpy_*omega_c + Q_*pwpqy;
    pqdotpq.col(2) = pQpz_*omega_c + Q_*pwpqz;
    pqdotpq.col(3) = pQpw_*omega_c + Q_*pwpqw;

    MatrixXd J_r = MatrixXd::Zero(13,13);
    J_r.block<3,3>(0,10) = I / m;
    J_r.block<4,4>(3,3) = pqdotpq;
    J_r.block<4,3>(3,10) = pqdotpL;

    ////////// END RIGID BODY JACOBIAN /////////////////////////////////

    ///////////////////// LINEARIZED IMPLICIT SOLVER ////////////////////////////

    MatrixXd A = MatrixXd::Zero(13,13);
    A = MatrixXd::Identity(13,13) - h*(J_r + J_c);

    VectorXd b = VectorXd::Zero(13);
    b.block<3,1>(0,0) = k3*v_c;
    Quaterniond temp = Quaterniond(omega_c(0),omega_c(1),omega_c(2),0)*q_c;
    b.block<4,1>(3,0) = k3*0.5*Vector4d(temp.x(),temp.y(),temp.z(),temp.w());
    b.block<3,1>(7,0) = k3*F_c;
    b.block<3,1>(10,0) = k3*T_c;

    VectorXd y = A.colPivHouseholderQr().solve(b);

    // save last states
    spoon->x_minus_1 = spoon->x; spoon->q_minus_1 = spoon->q;
    spoon->xdot_minus_1 = spoon->xdot; spoon->omega_minus_1 = spoon->omega;

    // Peform integration
    x_c = y.block<3,1>(0,0) + k1*x_c - k2*x_c_minus_1;
    q_c = sumQuaternions(Quaterniond(y(3),y(4),y(5),y(6)), q_c);
    q_c.normalize();
    v_c = y.block<3,1>(7,0)/m + k1*v_c - k2*v_c_minus_1;
    omega_c = spoon->IbodyInv*y.block<3,1>(10,0) + k1*omega_c - k2*omega_c_minus_1;

    F_c = kc * ( x_d - x_c ) + bc * ( v_d - v_c ); //! CHECK!
    T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c); //! CHECK

    // update new states
    spoon->x = x_c;
    spoon->q = q_c;
    spoon->q.normalize();
    spoon->update_inertia_matrix();
    spoon->xdot = v_c;
    spoon->omega = omega_c;

    if (m_recordData)
    {
        if (!dataFile->is_open())
        {
            std::cerr << "Failed to open the file for writing." << std::endl;
        }

        *dataFile << simTime << " " << x_d(0) << " " << x_d(1) << " " << x_d(2) << " " <<
                  q_d.x() << " " << q_d.y() << " " << q_d.z() << " " << q_d.w() << " " <<
                  x_c(0) << " " << x_c(1) << " " << x_c(2) << " " <<
                  q_c.x() << " " << q_c.y() << " " << q_c.z() << " " << q_c.w() << " " <<
                    F_c(0) << " " << F_c(1) << " " << F_c(2) << " " <<
               T_c(0) << " " << T_c(1) << " " << T_c(2) << endl;
    }

    simTime += h;
}

void PenaltyDemo::BDF2Solver(double h, vector<Collision> potentialCollisions)
{

    VectorXd x_d = spoon->x_d; Quaterniond q_d = spoon->q_d; q_d.normalize();
    VectorXd v_d = spoon->xdot_d; VectorXd omega_d = spoon->omega_d;

    if (m_readTrajectory)
    {
        auto state = sinusoidalTrajectory(0.1,0.7,0,0,simTime);
        x_d(2) = state.position;
        x_d(1) = state.position;
        x_d(2) = state.velocity;
        x_d(1) = state.velocity;
    }

    double phi = h / h_minus_1;
    double k1 = (1+phi)*(1+phi)/(1+2*phi);
    double k2 = phi*phi/(1+2*phi);
    double k3 = h*(1+phi)/(1+2*phi);

    VectorXd x_c_minus_1 = spoon->x_minus_1; Quaterniond q_c_minus_1 = spoon->q_minus_1;
    VectorXd x_c = spoon->x; Quaterniond q_c = spoon->q; q_c.normalize();
    VectorXd v_c_minus_1 = spoon->xdot_minus_1; VectorXd omega_c_minus_1 = spoon->omega_minus_1;
    VectorXd v_c = spoon->xdot; VectorXd omega_c = spoon->omega;
    MatrixXd I = MatrixXd::Identity(3,3);
    double kc = spoon->Kc; double bc = spoon->Bc;
    double kt = spoon->Kt; double bt = spoon->Bt;
    double m = spoon->mass;
    Vector3d c = VectorXd::Zero(3); Quaterniond q = Quaterniond::Identity();
    MatrixXd J_s = spoon->Ibody; MatrixXd Jinv_s = spoon->IbodyInv;
    MatrixXd L = q_c*J_s*q_c.inverse()*omega_c; Matrix3d R = q_c.toRotationMatrix();
    Quaterniond dq = q_d*q_c.inverse()*q.inverse();
    Vector3d u_c = 2*acos(q.w())*dq.vec();
    Quaterniond q_inv = q.inverse();
    MatrixXd C = Vector4d(dq.x(),dq.y(),dq.z(),dq.w())*RowVector4d(q_inv.x(),q_inv.y(),q_inv.z(),q_inv.w());

    // the virtual coupling force
    Vector3d F_c = kc*(x_d-x_c) + bc*(v_d - v_c); //! CHECK!
    VectorXd T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c); //! CHECK

    int it = 0;
    while(it < maxIt)
    {
        ////////// Coupling Jacobian /////////////////////////////////
        MatrixXd J_c = MatrixXd::Zero(13,13);

        // Derivative of angular velocity wrt quaternion
        VectorXd pwpqx = (pRpx(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpx(q_c).transpose())*L; //! CHECK
        VectorXd pwpqy = (pRpy(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpy(q_c).transpose())*L; //! CHECK
        VectorXd pwpqz = (pRpz(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpz(q_c).transpose())*L; //! CHECK
        VectorXd pwpqw = (pRpw(q_c)*(I/m)*R.transpose() + R*(I/m)*pRpw(q_c).transpose())*L; //! CHECK

        // derivatives of force with respect to quaternion
        VectorXd pFcpqx= -kc*pRpx(q_c)*c + bc * vec2skew(c)*pwpqx; //! CHECK!
        VectorXd pFcpqy= -kc*pRpy(q_c)*c + bc * vec2skew(c)*pwpqy; //! CHECK!
        VectorXd pFcpqz= -kc*pRpz(q_c)*c + bc * vec2skew(c)*pwpqz; //! CHECK!
        VectorXd pFcpqw= -kc*pRpw(q_c)*c + bc * vec2skew(c)*pwpqw; //! CHECK!

        // Derivative of linear velocity with respect to quaternion
        MatrixXd pupq = 2*acos(dq.w())*C.block<3,4>(0,0) - 2/(sqrt(1-dq.w()*dq.w()))*dq.vec()*C.block<1,4>(3,0); //! CHECK!

        if (pupq.hasNaN())
            pupq.setZero();

        // Elements of Jacobian
        MatrixXd pFcpx = -kc*MatrixXd::Identity(3,3); //! CHECK!
        MatrixXd pTcpx = vec2skew(R*c)*pFcpx; //! CHECK!
        MatrixXd pFcpq(3,4);
        pFcpq.col(0) = pFcpqx; //! CHECK!
        pFcpq.col(1) = pFcpqy; //! CHECK!
        pFcpq.col(2) = pFcpqz; //! CHECK!
        pFcpq.col(3) = pFcpqw; //! CHECK!
        MatrixXd pTcpq(3,4);

        pTcpq.col(0) = vec2skew(R*c)*pFcpqx - vec2skew(F_c)*pRpx(q_c)*c + kt*pupq.col(0) - bt*pwpqx; //! Check!
        pTcpq.col(1) = vec2skew(R*c)*pFcpqy - vec2skew(F_c)*pRpy(q_c)*c + kt*pupq.col(1) - bt*pwpqy; //! Check!
        pTcpq.col(2) = vec2skew(R*c)*pFcpqz - vec2skew(F_c)*pRpz(q_c)*c + kt*pupq.col(2) - bt*pwpqz; //! Check!
        pTcpq.col(3) = vec2skew(R*c)*pFcpqw - vec2skew(F_c)*pRpw(q_c)*c + kt*pupq.col(3) - bt*pwpqw; //! Check!
        MatrixXd pFcpP = bc*I/m; //! Check!
        MatrixXd pTcpP = vec2skew(R*c)*pFcpP; //! Check!
        MatrixXd pFcpL = bc*vec2skew(c)*R*(MatrixXd::Identity(3,3)/m)*R.transpose();
        MatrixXd pTcpL = (bc*vec2skew(R*c)*vec2skew(c) - bt*I)*R*(I/m)*R.transpose();

        J_c.block<3,3>(7,0) = pFcpx;
        J_c.block<3,3>(10,0) = pTcpx;
        J_c.block<3,4>(7,3) = pFcpq;
        J_c.block<3,4>(10,3) = pTcpq;
        J_c.block<3,3>(7,7) = pFcpP;
        J_c.block<3,3>(10,7) = pTcpP;
        J_c.block<3,3>(7,10) = pFcpL;
        J_c.block<3,3>(10,10) = pTcpL;

        /////////////////////  END OF COUPLING JACOBIAN ////////////////////

        ///////////////////// PENALTY JACOBIAN ////////////////////////////

        // create the jacobians and force vectors
        //    VectorXd F_p = VectorXd::Zero(3); VectorXd T_p = VectorXd::Zero(4);
        //
        //    // Penalty force
        //    MatrixXd F_p  =
        //
        //    MatrixXd J_p = ;

        ///////////////////// END OF PENALTY JACOBIAN ////////////////////////////

        ////////// RIGID BODY JACOBIAN /////////////////////////////////

        MatrixXd Q_ = Q(q_c);
        MatrixXd pQpx_ = pQpx();
        MatrixXd pQpy_ = pQpy();
        MatrixXd pQpz_ = pQpz();
        MatrixXd pQpw_ = pQpw();

        MatrixXd pqdotpL = Q_*R*(I/m)*R.transpose();
        MatrixXd pqdotpq(4,4);
        pqdotpq.col(0) = pQpx_*omega_c + Q_*pwpqx;
        pqdotpq.col(1) = pQpy_*omega_c + Q_*pwpqy;
        pqdotpq.col(2) = pQpz_*omega_c + Q_*pwpqz;
        pqdotpq.col(3) = pQpw_*omega_c + Q_*pwpqw;

        MatrixXd J_r = MatrixXd::Zero(13,13);
        J_r.block<3,3>(0,10) = I / m;
        J_r.block<4,4>(3,3) = pqdotpq;
        J_r.block<4,3>(3,10) = pqdotpL;

        ////////// END RIGID BODY JACOBIAN /////////////////////////////////


        ///////////////////// LINEARIZED IMPLICIT SOLVER ////////////////////////////

        MatrixXd A = MatrixXd::Zero(13,13);
        A = MatrixXd::Identity(13,13) - h*(J_r + J_c);
        // preconditioner
        A += 1e6*MatrixXd::Identity(13,13);

        VectorXd b = VectorXd::Zero(13);
        b.block<3,1>(0,0) = k3*v_c;
        Quaterniond temp = Quaterniond(omega_c(0),omega_c(1),omega_c(2),0)*q_c;
        b.block<4,1>(3,0) = k3*0.5*Vector4d(temp.x(),temp.y(),temp.z(),temp.w());
        b.block<3,1>(7,0) = k3*F_c;
        b.block<3,1>(10,0) = k3*T_c;

        VectorXd y = A.colPivHouseholderQr().solve(b);

        // Perform integration
        x_c = 0.75*y.block<3,1>(0,0) + k1*x_c - k2*x_c_minus_1;
        q_c = sumQuaternions(Quaterniond(y(3),y(4),y(5),y(6)), q_c);
        q_c.normalize();
        v_c = 0.75*y.block<3,1>(7,0)/m + k1*v_c - k2*v_c_minus_1;
        omega_c = q_c*J_s*q_c.inverse()*y.block<3,1>(10,0)+ k1*omega_c - k2*omega_c_minus_1;

        // Update the constants
        L = q_c*J_s*q_c.inverse()*omega_c; R = q_c.toRotationMatrix();
        dq = q_d*q_c.inverse()*q.inverse();
        u_c = 2*acos(q.w())*dq.vec();
        C = Vector4d(dq.x(),dq.y(),dq.z(),dq.w())*RowVector4d(q_inv.x(),q_inv.y(),q_inv.z(),q_inv.w());
        F_c = kc*(x_d-x_c) + bc*(v_d - v_c); //! CHECK!
        T_c = (R*c).cross(F_c) + kt*u_c + bt*(omega_d - omega_c); //! CHECK

        if (y.norm() < 1e-3)
            break;

        it++;

    }

    // save last states
    spoon->x_minus_1 = spoon->x; spoon->q_minus_1 = spoon->q;
    spoon->xdot_minus_1 = spoon->xdot; spoon->omega_minus_1 = spoon->omega;

    // update new states
    spoon->x = x_c;
    spoon->q = q_c;
    spoon->q.normalize();
    spoon->update_inertia_matrix();
    spoon->xdot = v_c;
    spoon->omega = omega_c;

    if (m_recordData)
    {
        if (!dataFile->is_open())
        {
            std::cerr << "Failed to open the file for writing." << std::endl;
            return;
        }

        *dataFile << simTime << " " << x_d(0) << " " << x_d(1) << " " << x_d(2) << " " <<
                    q_d.x() << " " << q_d.y() << " " << q_d.z() << " " << q_d.w() << " " <<
                    x_c(0) << " " << x_c(1) << " " << x_c(2) << " " <<
                    q_c.x() << " " << q_c.y() << " " << q_c.z() << " " << q_c.w() << " " <<
                    F_c(0) << " " << F_c(1) << " " << F_c(2) << " " <<
                    T_c(0) << " " << T_c(1) << " " << T_c(2) << endl;

    }

    simTime += h;
}

vector<CollisionResult> PenaltyDemo::computeCollisions(vector<Collision> potentialCollisions)
{
//    int numPotentialCollisions = potentialCollisions.size(); // Assuming all three vectors have the same size.

//    spoon->updateSensors();
    vector<vector<Sensor*>> sensors = spoon->m_sensors;
    vector<vector<CollisionResult>> temp_results;
    temp_results.resize(sensors.size());

    for (int i = 0; i < sensors.size(); i++) {
        temp_results[i].resize(sensors[i].size());
        for (int j = 0 ; j < sensors[i].size(); j++)
        {
            temp_results[i][j] = checkSensorCollision(sensors[i][j]);
        }
    }

    vector<CollisionResult> results;

    for (int i = 0; i < sensors.size(); i++) {
        for (int j = 0 ; j < sensors[i].size(); j++)
        {
            if (temp_results[i][j].collision == true) {
                results.emplace_back(temp_results[i][j]);
            }
            }
    }



//    auto& lvertices = spoon->vertices();
//    auto& ltriangles = spoon->triangles();
//    Vector3d lpos = spoon->x;
//    Quaterniond lq = spoon->q;
//    auto& rvertices = mug->vertices();
//    auto& rtriangles = mug->triangles();
//    Vector3d rpos = mug->x;
//    Quaterniond rq = mug->q;
//    vector<vector<Sensor*>> leftSensors = spoon->m_sensors;
//    vector<vector<Sensor*>> rightSensors = mug->m_sensors;
//    vector<CollisionResult> results;
//
//    // TODO: Parallelize
//    for (int i = 0 ; i < numPotentialCollisions; i++)
//    {
//        int leftIdx = potentialCollisions[i].collidingTriangle1;
//        int rightIdx = potentialCollisions[i].collidingTriangle2;
//
//        // check the left rays on the right triangle
//        vector<CollisionResult> leftCollisionResults;
//        for (int l = 0; l < leftSensors[leftIdx].size(); l++)
//        {
//            //
////            Vector3d vidx = rtriangles.row(rightIdx);
////            Vector3d a =
////            checkRayTriangleCollision(ray, )
//        }
//
//        // check the right rays on the left triangle
//        vector<CollisionResult> rightCollisionResults;
//        for (int r = 0; r < rightSensors[rightIdx].size(); r++)
//        {
//
//        }
//
//    }
//

    return results;
}