// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Unscented Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 0.05;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 2;
constexpr double INIT_PSI_STD = 5.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;
// -------------------------------------------------- //

// ----------------------------------------------------------------------- //
// USEFUL HELPER FUNCTIONS


std::vector<VectorXd> generateSigmaPoints(VectorXd state, MatrixXd cov) {
    std::vector<VectorXd> sigmaPoints;

    unsigned int numStates = state.rows();
    double k = 3.0 - numStates;

    VectorXd x0 = state;

    sigmaPoints.push_back(x0);

    MatrixXd newCov = (numStates+k)*cov;
    MatrixXd sqrtCov = newCov.llt().matrixL();

    for(int i = 0; i < numStates; ++i)
    {
        VectorXd deltaX = sqrtCov.col(i);

        sigmaPoints.push_back(state + deltaX);
        sigmaPoints.push_back(state - deltaX);
    }

    return sigmaPoints;
}

std::vector<double> generateSigmaWeights(unsigned int numStates) {
    std::vector<double> weights;


    double k = 3.0 - numStates;

    double W0 = k/(numStates + k);
    weights.push_back(W0);

    double wi = 1/(2*(numStates+k));

    for (int i = 1; i <= 2*numStates; i++) {
        weights.push_back(wi);
    }

    return weights;
}

VectorXd normaliseState(VectorXd state)
{
    state(2) = wrapAngle(state(2));
    return state;
}

// ----------------------------------------------------------------------- //

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Update Step for the Lidar Measurements in the 
        // section below.
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // Hint: You can use the constants: LIDAR_RANGE_STD, LIDAR_THETA_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        BeaconData map_beacon = map.getBeaconWithId(meas.id); // Match Beacon with built in Data Association Id
        if (meas.id != -1 && map_beacon.id != -1) // Check that we have a valid beacon match
        {
           

        }
        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Prediction Step for the system in the  
        // section below.
        // Hint: You can use the constants: ACCEL_STD, GYRO_STD
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
        MatrixXd Q = Matrix2d::Zero();
        Q(0,0) = dt*GYRO_STD*GYRO_STD;
        Q(1,1) = dt*ACCEL_STD*ACCEL_STD;

        // Augment the state vector.
        int n_x = state.size();
        int aug_size = 2;
        int x_aug_size = n_x+aug_size;

        VectorXd x_aug = VectorXd::Zero(x_aug_size);
        x_aug.head(n_x) = state;

        // Augment the covariance matrix.
        MatrixXd cov_aug = MatrixXd::Zero(x_aug_size, x_aug_size);
        cov_aug.topLeftCorner(n_x, n_x) = cov;
        cov_aug.bottomRightCorner(aug_size, aug_size) = Q;

        
        // Generate Sigma Points
        std::vector<VectorXd> sigmaPoints = generateSigmaPoints(x_aug, cov_aug);
        std::vector<double> sigmaWeights = generateSigmaWeights(x_aug_size);

        // Transform Sigma Points with Process Model
        std::vector<VectorXd> sigma_points_predict;
        for (const auto& sigma_point : sigmaPoints)
        {
            double x = sigma_point(0);
            double y = sigma_point(1);
            double psi = sigma_point(2);
            double V = sigma_point(3);
            double psi_dot = gyro.psi_dot;
            double psi_dot_noise = sigma_point(4);
            double accel_noise = sigma_point(5);

            // Predict using process model
            double newX = x + dt*V*cos(psi);
            double newY = y + dt*V*sin(psi);
            double newPsi = psi + dt*psi_dot+psi_dot_noise;
            double newV = V + dt*accel_noise;

            VectorXd sigmaPointPredicted = Vector4d::Zero();
            sigmaPointPredicted(0,0) = newX;
            sigmaPointPredicted(1,0) = newY;
            sigmaPointPredicted(2,0) = newPsi;
            sigmaPointPredicted(3,0) = newV;

            sigma_points_predict.push_back(normaliseState(sigmaPointPredicted));

        }

        // Calculate mean of transformed sigma points.
        state = VectorXd::Zero(n_x);
        for (int i = 0; i < sigma_points_predict.size(); i++) {
            state += sigmaWeights[i]*sigma_points_predict[i];
        }

        // Calculate covariance of transformed sigma points.
        cov = MatrixXd::Zero(n_x, n_x);
        for (int i = 0; i < sigma_points_predict.size(); i++) {
            VectorXd diff = normaliseState(sigma_points_predict[i] - state);
            cov += sigmaWeights[i]* diff * diff.transpose(); 
        }


    

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

        VectorXd z = Vector2d::Zero();
        MatrixXd H = MatrixXd(2,4);
        MatrixXd R = Matrix2d::Zero();

        z << meas.x,meas.y;
        H << 1,0,0,0,0,1,0,0;
        R(0,0) = GPS_POS_STD*GPS_POS_STD;
        R(1,1) = GPS_POS_STD*GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        MatrixXd K = cov*H.transpose()*S.inverse();

        state = state + K*y;
        cov = (Matrix4d::Identity() - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        // You may modify this initialisation routine if you can think of a more
        // robust and accuracy way of initialising the filter.
        // ----------------------------------------------------------------------- //
        // YOU ARE FREE TO MODIFY THE FOLLOWING CODE HERE

        VectorXd state = Vector4d::Zero();
        MatrixXd cov = Matrix4d::Zero();

        state(0) = meas.x;
        state(1) = meas.y;
        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = GPS_POS_STD*GPS_POS_STD;
        cov(2,2) = INIT_PSI_STD*INIT_PSI_STD;
        cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;

        setState(state);
        setCovariance(cov);

        // ----------------------------------------------------------------------- //
    }             
}

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
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
        VectorXd state = getState(); // STATE VECTOR [X,Y,PSI,V,...]
        return VehicleState(state[0],state[1],state[2],state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt){}
