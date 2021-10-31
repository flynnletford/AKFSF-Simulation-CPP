// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Extended Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double GYRO_BIAS_STD = 10.0;

constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;

constexpr double MAX_LIDAR_RANGE = 90.0;
// -------------------------------------------------- //

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
}

int stopCount = 0;

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
        // HINT: You can use the constants: LIDAR_RANGE_STD, LIDAR_THETA_STD
        // HINT: The mapped-matched beacon position can be accessed by the variables
        // map_beacon.x and map_beacon.y
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        double Px = state(0);
        double Py = state(1);

        std::vector<BeaconData> map_beacons_within_range = map.getBeaconsWithinRange(Px, Py, MAX_LIDAR_RANGE); // Match Beacon with built in Data Association Id
        
        if (map_beacons_within_range.size() == 0) {
            std::cout << "No map beacons found within range; skipping" << std::endl;
            return;
        }

        // if (stopCount > 30) {
        //     return;
        // }

        BeaconData map_beacon;
        int count = 0;
        double lowestCost;
        std::cout << "*******************************************" << std::endl;
        for(const auto& beacon : map_beacons_within_range) {
            double Lx = beacon.x;
            double Ly = beacon.y;

            double predictedRange = sqrt( pow(Lx - Px,2) + pow(Ly - Py,2) );
            double predictedBearing = wrapAngle(atan2( Ly - Py, Lx - Px ) - state(2));

            VectorXd v(2,1);
            v << meas.range - predictedRange, wrapAngle(meas.theta - predictedBearing);

            VectorXd costValues(2,1);
            costValues << meas.range - predictedRange, cos(wrapAngle(meas.theta - predictedBearing))*MAX_LIDAR_RANGE;
            double cost = sqrt(costValues.transpose()*costValues);
            if ( (count == 0) || (cost < lowestCost) ) {
                map_beacon = beacon;
                lowestCost = cost;

                if (map_beacons_within_range.size() > 3) {
                    std::cout << "count: " << count << std::endl;
                    std::cout << "predictedRange: " << predictedRange << ", predictedBearing: " << predictedBearing << std::endl;
                    std::cout << "meas.range: " << meas.range << ", meas.theta: " << meas.theta << std::endl;
                    std::cout << "v(0,0): " << v(0,0) << ", v(1,0): " << v(1,0) << std::endl;
                }
            }
            count++;
        }
        
        std::cout << "*******************************************" << std::endl;
        stopCount++;
        
        
        if (meas.id != -1 && map_beacon.id != -1)
        {           
            // The map matched beacon positions can be accessed using: map_beacon.x AND map_beacon.y
            double Px = state(0);
            double Py = state(1);
            double Lx = map_beacon.x;
            double Ly = map_beacon.y;


            double predictedRange = sqrt( pow(Lx - Px,2) + pow(Ly - Py,2) );
            double predictedBearing = wrapAngle(atan2( Ly - Py, Lx - Px ) - state(2));

            MatrixXd H(2,5);
            H <<    (Px - Lx)/predictedRange, (Py - Ly)/predictedRange, 0, 0, 0,
                    (Py - Ly)/-pow(predictedRange, 2), (Px - Lx)/pow(predictedRange, 2), -1, 0, 0;

            VectorXd v(2,1);
            v << meas.range - predictedRange, wrapAngle(meas.theta - predictedBearing);

            MatrixXd R(2,2);
            R << pow(LIDAR_RANGE_STD, 2), 0,
                 0, pow(LIDAR_THETA_STD, 2);

            MatrixXd S = H*cov*H.transpose() + R;

            MatrixXd K = cov*H.transpose()*S.inverse();

            MatrixXd I = MatrixXd::Identity(5,5);

            state = state + K*v;
            state(2) = wrapAngle(state(2));
            cov = (I - K*H)*cov;
        } else {
            std::cout << "No LiDAR association; skipping" << std::endl;
        }

        // ----------------------------------------------------------------------- //
        if (state(3) > 0) {
            std::cout << "Setting velocity here!!!" << std::endl;
        } else {

        setState(state);
        setCovariance(cov);
        }


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
        // HINT: Assume the state vector has the form [PX, PY, PSI, V].
        // HINT: Use the Gyroscope measurement as an input into the prediction step.
        // HINT: You can use the constants: ACCEL_STD, GYRO_STD
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        double gyro_noise = 0.0;

        VectorXd deltaState(5,1);
        deltaState << state(3)*cos(state(2)),
                      state(3)*sin(state(2)),
                      gyro.psi_dot - state(4) + gyro_noise,
                      0,
                      0;

        state = state + dt*deltaState;

        // Ensure angle is between -pi and pi
        state(2) = wrapAngle(state(2));


        MatrixXd stateJacobian(5,5);
        stateJacobian << 1, 0, -dt*state(3)*sin(state(2)), dt*cos(state(2)), 0,
                         0, 1, dt*state(3)*cos(state(2)), dt*sin(state(2)), 0,
                         0, 0, 1, 0, 0,
                         0, 0, 0, 1, 0,
                         0, 0, 0, 0, 1;

        MatrixXd Q(5,5);
        Q << 0,0,0,0,0,
             0,0,0,0,0,
             0,0,pow(dt,2)*pow(GYRO_STD, 2), 0,0,
             0, 0, 0, pow(dt,2)*pow(ACCEL_STD, 2), 0,
             0,0,0,0, pow(dt,2)*pow(GYRO_BIAS_STD, 2);


        cov = stateJacobian*cov*stateJacobian.transpose() + Q;

        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the EKF update state would just produce the same result.
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        VectorXd z = Vector2d::Zero();
        MatrixXd H = MatrixXd(2,5);
        MatrixXd R = Matrix2d::Zero();

        z << meas.x,meas.y;
        H << 1,0,0,0,0,0,1,0,0,0;
        R(0,0) = GPS_POS_STD*GPS_POS_STD;
        R(1,1) = GPS_POS_STD*GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;

        MatrixXd S = H * cov * H.transpose() + R;

        // Innovation checking: Only fuse gps data if it's below the threshold for a valid value.
        // This is the Chi Squared Distribution value.

        // Normalised Innovation Squared
        double NIS = (y.transpose()) * (S.inverse()) * y;
        //std::cout << "NIS: " << NIS << std::endl;

        if (NIS >= 11.07) {
            //std::cout << "NIS for GPS measurement is above threshold; ignoring" << std::endl;
            return;
        }

        MatrixXd K = cov*H.transpose()*S.inverse();

        state = state + K*y;
        cov = (MatrixXd::Identity(5,5) - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        int state_size = 4;
        int gyro_aug_size = 1;
        int n_aug = state_size + gyro_aug_size;

        VectorXd state = VectorXd::Zero(n_aug);
        MatrixXd cov = MatrixXd::Zero(n_aug, n_aug);

        state(0) = meas.x;  // x
        state(1) = meas.y;  // y
        state(2) = -1.5708; // heading          // For testing
        state(3) = -2.0;    // V                // For testing
        state(4) = 0;       // Gyro bias - initialise as zero.

        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = GPS_POS_STD*GPS_POS_STD;
        cov(2,2) = INIT_PSI_STD*INIT_PSI_STD;
        cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;
        cov(4,4) = GYRO_BIAS_STD*GYRO_BIAS_STD; // Initial uncertainty in gyroscope bias.

        setState(state);
        setCovariance(cov);
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
        VectorXd state = getState(); // STATE VECTOR [X,Y,PSI,V,...]
        return VehicleState(state[0],state[1],state[2],state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt){}
