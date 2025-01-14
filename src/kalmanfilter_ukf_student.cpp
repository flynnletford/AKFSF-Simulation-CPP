// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Unscented Kalman Filter
//

// Notes on using the Unscented Kalman Filter (UKF):
// 1. Just like the EKF, the UKF is not guaranteed to converge.
// 2. UKF can diverge if the state is too far away from the truth as it still uses an approximation of the system.
// 3. To avoid convergence issues, try to initialise the full state and covariance from measurement data. 
//    Assuming zero with large uncertainty can lead to issues. The handleGPSMeasurement function provides an example
//    of initialising the position state values as well as the covariance matrix with default noise values. However,
//    the initial velocity and heading values are not set. If the UKF were to be started with the real-world system
//    having non-zero values for these two parameters, errors are to be expected.

#include "kalmanfilter.h"
#include "utils.h"

constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;

// The sigma points allow us to capture the distributions shape with a minimum number of points.
// Sigma points are based on the square root of the covariance, hence the name "sigma".
// For a state vector of dimensions n x 1, we generate 2n + 1 sigma points.
std::vector<VectorXd> generateSigmaPoints(VectorXd state, MatrixXd cov) {
    std::vector<VectorXd> sigmaPoints;

    unsigned int numStates = state.rows();
    double k = 3.0 - numStates;

    VectorXd x0 = state;

    sigmaPoints.push_back(x0);

    MatrixXd newCov = (numStates+k)*cov;

    // We can calculate the square root of a real positive definite matrix using the Cholesky Decomposition.
    // llt() represents the equation L*L.transpose(). The matrix L is the square root of the original matrix.
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

    // For a Gaussian distribution, the following weight has been shown to improve the accuracy.
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

VectorXd normaliseLidarMeasurement(VectorXd meas)
{
    meas(1) = wrapAngle(meas(1));
    return meas;
}


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
            int n_x = state.size();
            int n_z = 2;

            // Measurement Vector
            VectorXd z = VectorXd::Zero(n_z);
            z << meas.range, meas.theta;

            // Noise Covariance Matrix
            MatrixXd R = MatrixXd::Zero(n_z,n_z);
            R(0,0) = LIDAR_RANGE_STD*LIDAR_RANGE_STD;
            R(1,1) = LIDAR_THETA_STD*LIDAR_THETA_STD;

            // Augment the State Vector with Noise States
            int n_aug = n_x + n_z;

            VectorXd x_aug = VectorXd::Zero(n_aug);
            x_aug.head(n_x) = state;

            // Augment the Covariance Matrix with Noise States
            MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);
            P_aug.topLeftCorner(n_x,n_x) = cov;
            P_aug.bottomRightCorner(n_z,n_z) = R;

            // Generate sigma points
            std::vector<VectorXd> x_sig = generateSigmaPoints(x_aug, P_aug);
            std::vector<double> weights_sig = generateSigmaWeights(n_aug);

            // Transform sigma points with LiDAR measurement model.
            std::vector<VectorXd> z_sig;
            for (const auto& x : x_sig) {
                double Lx = map_beacon.x;
                double Ly = map_beacon.y;
                double rangeNoise = x[4];
                double headingNoise = x[5];

                double deltaX = Lx - x[0];
                double deltaY = Ly - x[1];

                double rangeHat = sqrt(deltaX*deltaX + deltaY*deltaY) + rangeNoise;
                double headingHat = atan2(deltaY, deltaX) - x[2] + headingNoise;

                VectorXd z_hat = VectorXd::Zero(n_z);
                z_hat << rangeHat, headingHat;
                z_sig.push_back(normaliseLidarMeasurement(z_hat));
            }
           
            // Calculate Measurement Mean
            VectorXd z_mean = VectorXd::Zero(n_z);
            for(unsigned int i = 0; i < z_sig.size(); i++){z_mean += weights_sig[i] * z_sig[i];}

            // Calculate Innovation Covariance
            // This represents how much uncertainty there is inside the measurement.
            MatrixXd S = MatrixXd::Zero(n_z,n_z);
            for(unsigned int i = 0; i < z_sig.size(); i++)
            {
                VectorXd diff = normaliseLidarMeasurement(z_sig[i] - z_mean);
                S += weights_sig[i] * diff * diff.transpose();
            }

            // Calculate the Cross Covariance
            // This represents how much uncertainty there is between the state and the measurement.
            MatrixXd Pxz = MatrixXd::Zero(n_aug, n_z);
            for (int i = 0; i < z_sig.size(); i++) {
                VectorXd x_diff = normaliseState(x_sig[i]  - x_aug);
                VectorXd z_diff = normaliseLidarMeasurement(z_sig[i]  - z_mean);
                Pxz += weights_sig[i] * x_diff * z_diff.transpose();
            }

            // Implement the UKF Update step equations
            VectorXd measurementInnovation = normaliseLidarMeasurement(z - z_mean);
            MatrixXd K = Pxz*S.inverse();
            VectorXd new_x_aug = x_aug + K*measurementInnovation;

            MatrixXd new_P_aug = P_aug - K*S*K.transpose();

            state = new_x_aug.head(n_x);
            cov = new_P_aug.topLeftCorner(n_x, n_x);
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

        // STEP 1: Augment the state vector and covariance matrix with noise.

        // Augment the state vector.
        MatrixXd Q = Matrix2d::Zero();
        Q(0,0) = GYRO_STD*GYRO_STD;
        Q(1,1) = ACCEL_STD*ACCEL_STD;

        int n_x = state.size();
        int n_w = 2;
        int n_aug = n_x+n_w;

        // Our augmented state will just be the state vector with [0; 0] appended as we assume a noise distribution
        // with zero mean for both gyroscope noise, and process model acceleration noise.
        VectorXd x_aug = VectorXd::Zero(n_aug);
        x_aug.head(n_x) = state;

        // Augment the covariance matrix.
        MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);
        P_aug.topLeftCorner(n_x, n_x) = cov;
        P_aug.bottomRightCorner(n_w, n_w) = Q;
        
        // STEP 2: Generate Sigma Points
        std::vector<VectorXd> sigmaPoints = generateSigmaPoints(x_aug, P_aug);
        std::vector<double> sigmaWeights = generateSigmaWeights(n_aug);

        // STEP 3: Transform Sigma Points with Process Model

        // This has the effect of predicting the states forward in time using the process model.
        // The transformed sigma points provide a prediction of the Gaussian distribution for the next time step.
        std::vector<VectorXd> sigma_points_predict;
        for (const auto& sigma_point : sigmaPoints)
        {
            double x = sigma_point(0);
            double y = sigma_point(1);
            double psi = sigma_point(2);
            double V = sigma_point(3);

            // Our model assume the input is the turn rate of the vehicle as measured by a gryoscope.
            double psi_dot = gyro.psi_dot;
            double psi_dot_noise = sigma_point(4);

            // We assume that the acceleration is unknown and so is modelled by process noise as a random variable.
            double accel_noise = sigma_point(5);

            // Predict using process model
            double newX = x + dt*V*cos(psi);
            double newY = y + dt*V*sin(psi);
            double newPsi = psi + dt*(psi_dot+psi_dot_noise);
            double newV = V + dt*accel_noise;

            VectorXd sigmaPointPredicted = Vector4d::Zero();
            sigmaPointPredicted << newX, newY, newPsi, newV;

            sigma_points_predict.push_back(sigmaPointPredicted);
        }

        // STEP 4: Calculate the mean and covariance of transformed sigma points.

        // 4a) Calculate mean of transformed sigma points.
        state = VectorXd::Zero(n_x);
        for (int i = 0; i < sigma_points_predict.size(); ++i) {
            state += sigmaWeights[i]*sigma_points_predict[i];
        }

        // 4b) Calculate covariance of transformed sigma points.
        cov = MatrixXd::Zero(n_x, n_x);
        for (int i = 0; i < sigma_points_predict.size(); ++i) {
            VectorXd diff = normaliseState(sigma_points_predict[i] - state);
            cov += sigmaWeights[i]* diff * diff.transpose(); 
        }

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the UKF update state would just produce the same result.
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
