// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Linear Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr bool INIT_ON_FIRST_PREDICTION = false;
constexpr double INIT_POS_STD = 0;
constexpr double INIT_VEL_STD = 0; //15;
constexpr double ACCEL_STD = 0.1; //0.1;
constexpr double GPS_POS_STD = 3.0;
// -------------------------------------------------- //

void KalmanFilter::predictionStep(double dt)
{
    if (!isInitialised() && INIT_ON_FIRST_PREDICTION)
    {
        // Implement the State Vector and Covariance Matrix Initialisation in the
        // section below if you want to initialise the filter WITHOUT waiting for
        // the first measurement to occur. Make sure you call the setState() /
        // setCovariance() functions once you have generated the initial conditions.
        // Hint: Assume the state vector has the form [X,Y,VX,VY].
        // Hint: You can use the constants: INIT_POS_STD, INIT_VEL_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
            VectorXd state = Vector4d::Zero();
            MatrixXd cov = Matrix4d::Zero();

            // Assume the initial position is (X,Y) = (0,0) m
            // Assume the initial velocity is 5 m/s at 45 degrees (VX,VY) = (5*cos(45deg),5*sin(45deg)) m/s
            //state << 0, 0, 5.0*cos(M_PI/4), 5.0*sin(M_PI/4);
            state << 0, 0, 0, 0;

            cov(0,0) = INIT_POS_STD*INIT_POS_STD;
            cov(1,1) = INIT_POS_STD*INIT_POS_STD;
            cov(2,2) = INIT_VEL_STD*INIT_VEL_STD;
            cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;

            setState(state);
            setCovariance(cov);
        // ----------------------------------------------------------------------- //
    }

    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Prediction Step for the system in the  
        // section below.
        // Hint: You can use the constants: ACCEL_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        // F is the process model defined by the system dynamics.
        MatrixXd F = Matrix4d::Zero();
        F << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;   

        // Q is the covariance of the process model noise.
        MatrixXd Q = Matrix2d::Zero();    
        Q << pow(ACCEL_STD, 2), 0,
             0, pow(ACCEL_STD, 2);

        // L is the process model noise sensitivity. Often this is assumed to be the Identity matrix to simplify the equation.
        MatrixXd L(4,2);
        L << 0.5*pow(dt, 2), 0, 
             0, 0.5*pow(dt, 2),
             dt, 0,
             0, dt;

        // Perform the prediction step, calculating the priori state and covariance estimates.
        state = F*state; // no inputs so as simple as this.
        cov = F*cov*(F.transpose()) + L*Q*(L.transpose());


        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Update Step for the GPS Measurements in the 
        // section below.
        // Hint: Assume that the GPS sensor has a 3m (1 sigma) position uncertainty.
        // Hint: You can use the constants: GPS_POS_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE 
        
        // H is the measurement matrix
        MatrixXd H(2,4);
        H << 1, 0, 0, 0,
             0, 1, 0, 0;

        // R is the measurement noise covariance matrix. It tells us the uncertainty in the GPS measurement.
        MatrixXd R = Matrix2d::Zero();
        R << GPS_POS_STD*GPS_POS_STD, 0,
             0, GPS_POS_STD*GPS_POS_STD;

        // z is the 2x1 measurement vector.
        VectorXd z(2,1);
        z << meas.x, meas.y;

        // innovation is the error between the predicted measurements and the actual measurements we just received.
        VectorXd innovation = z - H*state;

        // S is the innovation covariance - the uncertainty in our innovation.
        MatrixXd S = H*cov*(H.transpose()) + R;

        // K is the Kalman Gain. This is the ratio of the innovation uncertainty and the current state uncertainty.
        // It allows us to provide a weighted estimate based on how much we trust the predicted state vs the current measurement.
        MatrixXd K = cov*(H.transpose())*S.inverse();

        state = state + K*innovation;
        cov = (MatrixXd::Identity(4,4) - K*H)*cov;

        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
    else
    {
        // Implement the State Vector and Covariance Matrix Initialisation in the
        // section below. Make sure you call the setState/setCovariance functions
        // once you have generated the initial conditions.
        // Hint: Assume the state vector has the form [X,Y,VX,VY].
        // Hint: You can use the constants: GPS_POS_STD, INIT_VEL_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
            VectorXd state = Vector4d::Zero();
            MatrixXd cov = Matrix4d::Zero();

            state << meas.x, meas.y, 0, 0;

            // R is the measurement noise covariance matrix. It tells us the uncertainty in the GPS measurement.
            MatrixXd R = Matrix2d::Zero();
            R << GPS_POS_STD*GPS_POS_STD, 0,
             0, GPS_POS_STD*GPS_POS_STD;

            cov << R(0,0), R(0,1), 0, 0,
                  R(1,0), R(1,1), 0, 0,
                  0, 0, pow(INIT_VEL_STD, 2), 0,
                  0, 0, 0, pow(INIT_VEL_STD, 2);


            setState(state);
            setCovariance(cov);
        // ----------------------------------------------------------------------- //
    }        
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance()
{
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0){pos_cov << cov(0,0), cov(0,1), cov(1,0), cov(1,1);}
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState()
{
    if (isInitialised())
    {
        VectorXd state = getState(); // STATE VECTOR [X,Y,VX,VY]
        double psi = std::atan2(state[3],state[2]);
        double V = std::sqrt(state[2]*state[2] + state[3]*state[3]);
        return VehicleState(state[0],state[1],psi,V);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt){predictionStep(dt);}
void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map){}
void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map){}

