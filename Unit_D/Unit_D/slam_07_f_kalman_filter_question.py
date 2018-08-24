# The full Kalman filter, consisting of prediction and correction step.
#
# slam_07_f_kalman_filter
# Claus Brenner, 12.12.2012
from lego_robot import *
from math import sin, cos, pi, atan2, sqrt
from numpy import *
from slam_d_library import get_observations, write_cylinders

def distance(p0, p1):   # For the purpose of putting out deviation.
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)
class ExtendedKalmanFilter:
    def __init__(self, state, covariance,
                 robot_width, scanner_displacement,
                 control_motion_factor, control_turn_factor,
                 measurement_distance_stddev, measurement_angle_stddev):
        # The state. This is the core data of the Kalman filter.
        self.state = state
        self.covariance = covariance

        # Some constants.
        self.robot_width = robot_width
        self.scanner_displacement = scanner_displacement
        self.control_motion_factor = control_motion_factor
        self.control_turn_factor = control_turn_factor
        self.measurement_distance_stddev = measurement_distance_stddev
        self.measurement_angle_stddev = measurement_angle_stddev

    @staticmethod
    def g(state, control, w):
        x, y, theta = state
        l, r = control
        if r != l:
            alpha = (r - l) / w
            rad = l/alpha
            g1 = x + (rad + w/2.)*(sin(theta+alpha) - sin(theta))
            g2 = y + (rad + w/2.)*(-cos(theta+alpha) + cos(theta))
            g3 = (theta + alpha + pi) % (2*pi) - pi
        else:
            g1 = x + l * cos(theta)
            g2 = y + l * sin(theta)
            g3 = theta

        return array([g1, g2, g3])

    @staticmethod
    def dg_dstate(state, control, w):
        theta = state[2]
        l, r = control
        if r != l:
            alpha = (r - l)/w
            R = l / alpha
            g1 = (R + (w / 2)) * (cos(theta + alpha) - cos(theta))
            g2 = (R + (w / 2)) * (sin(theta + alpha) - sin(theta))
            m = array([[1.0, 0.0, g1], [0.0, 1.0, g2], [0.0, 0.0, 1.0]])
            # Making a list of all rows, then give the list to array constructor.
        else:
            # This is for the special case r == l.
            g1 = -l * sin(theta)
            g2 = l * cos(theta)
            m = array([[1.0, 0.0, g1], [0.0, 1.0, g2], [0.0, 0.0, 1.0]])  # Replace this.
        return m

    @staticmethod
    def dg_dcontrol(state, control, w):

        theta = state[2]
        l, r = tuple(control)
        if r != l:
            alpha = (r - l) / w
            wr = (w * r) / ((r - l) ** 2)
            wl = (w * l) / ((r - l) ** 2)
            r2l = (r + l) / (2 * (r - l))

            g1_l = wr * (sin(theta + alpha) - sin(theta)) - r2l * cos(theta + alpha)
            g2_l = wr * (-cos(theta + alpha) + cos(theta)) - r2l * sin(theta + alpha)
            g3_l = - (1 / w)

            g1_r = -wl * (sin(theta + alpha) - sin(theta)) + r2l * cos(theta + alpha)
            g2_r = -wl * (-cos(theta + alpha) + cos(theta)) + r2l * sin(theta + alpha)
            g3_r = 1 / w

            m = array([[g1_l, g1_r], [g2_l, g2_r], [g3_l, g3_r]])
        else:

            # This is for the special case l == r.
            g1_l = .5 * (cos(theta) + (l / w) * sin(theta))
            g2_l = .5 * (sin(theta) - (l / w) * cos(theta))
            g3_l = - 1 / w

            g1_r = .5 * ((-l / w) * sin(theta) + cos(theta))
            g2_r = .5 * ((l / w) * cos(theta) + sin(theta))
            g3_r = 1 / w

            m = array([[g1_l, g1_r], [g2_l, g2_r], [g3_l, g3_r]])
        return m

    @staticmethod
    def get_error_ellipse(covariance):
        """Return the position covariance (which is the upper 2x2 submatrix)
           as a triple: (main_axis_angle, stddev_1, stddev_2), where
           main_axis_angle is the angle (pointing direction) of the main axis,
           along which the standard deviation is stddev_1, and stddev_2 is the
           standard deviation along the other (orthogonal) axis."""
        eigenvals, eigenvects = linalg.eig(covariance[0:2,0:2])
        angle = atan2(eigenvects[1,0], eigenvects[0,0])
        return (angle, sqrt(eigenvals[0]), sqrt(eigenvals[1]))        

    def predict(self, control):

        left, right = control
        alpha_1 = self.control_motion_factor
        alpha_2 = self.control_turn_factor

        g2l = (alpha_1 * left) ** 2 + (alpha_2 * (left - right)) ** 2
        g2r = (alpha_1 * right) ** 2 + (alpha_2 * (left - right)) ** 2

        Sigma_control = diag([g2l, g2r])
        Vt = self.dg_dcontrol(self.state, control, self.robot_width)
        VtT = Vt.T

        Sigma_covariance = self.covariance
        Gt = self.dg_dstate(self.state, control, self.robot_width)
        GtT = Gt.T

        self.covariance = dot(dot(Gt, Sigma_covariance), GtT) + dot(dot(Vt, Sigma_control), VtT)
        self.state = self.g(self.state, control, self.robot_width)  # It's also the expectation!
    @staticmethod
    def h(state, landmark, scanner_displacement):
        """Takes a (x, y, theta) state and a (x, y) landmark, and returns the
           measurement (range, bearing)."""
        dx = landmark[0] - (state[0] + scanner_displacement * cos(state[2]))
        dy = landmark[1] - (state[1] + scanner_displacement * sin(state[2]))
        r = sqrt(dx * dx + dy * dy)
        alpha = (atan2(dy, dx) - state[2] + pi) % (2*pi) - pi
        return array([r, alpha])
        # Put in the coordinates of the robot, return the coordinates of fixed landmarks in the local coordinate.

    @staticmethod
    def dh_dstate(state, landmark, scanner_displacement):

        x, y, theta = state
        x_m, y_m = landmark
        x_e = x + scanner_displacement * cos(theta)     # car to the scanner.
        y_e = y + scanner_displacement * sin(theta)
        delta_x = x_m - x_e
        delta_y = y_m - y_e
        q = (delta_x) ** 2 + (delta_y) ** 2

        dr_dx = -delta_x / sqrt(q)
        dr_dy = -delta_y / sqrt(q)
        dr_dtheta = (scanner_displacement / sqrt(q)) * \
                    (delta_x * sin(theta) - delta_y * cos(theta))
        dalpha_dx = delta_y / q
        dalpha_dy = -delta_x / q
        dalpha_dtheta = -(scanner_displacement / q) * \
                        (delta_x * cos(theta) + delta_y * sin(theta)) - 1
        return array([[dr_dx, dr_dy, dr_dtheta], [dalpha_dx, dalpha_dy, dalpha_dtheta]])

    def correct(self, measurement, landmark):
        """The correction step of the Kalman filter."""

        # --->>> Put your new code here.
        #
        # You will have to compute:
        # H, using dh_dstate(...).
        Ht = self.dh_dstate(self.state, landmark, self.scanner_displacement)
        HtT = Ht.T
        # Q, a diagonal matrix, from self.measurement_distance_stddev and
        #  self.measurement_angle_stddev (remember: Q contains variances).
        g1 = (self.measurement_distance_stddev) ** 2
        g2 = (self.measurement_angle_stddev) ** 2
        Q = diag([g1, g2])
        # K, from self.covariance, H, and Q.
        #  Use linalg.inv(...) to compute the inverse of a matrix.
        K = dot(self.covariance, dot(HtT, linalg.inv(dot(Ht, dot(self.covariance, HtT)) + Q)))
        # The innovation: it is easy to make an error here, because the
        #  predicted measurement and the actual measurement of theta may have
        #  an offset of +/- 2 pi. So here is a suggestion:
        #   innovation = array(measurement) -\
        #                self.h(self.state, landmark, self.scanner_displacement)
        #   innovation[1] = (innovation[1] + pi) % (2*pi) - pi
        innovation = array(measurement) - self.h(self.state, landmark, self.scanner_displacement)
        # The detected (r, alpha) minus the expected measurement based on the assumed position.
        innovation[1] = (innovation[1] + pi) % (2 * pi) - pi    # Edit one member of  an array.
        # Then, you'll have to compute the new self.state.
        mu_t = self.state + dot(K, innovation)
        self.state = mu_t
        # And finally, compute the new self.covariance. Use eye(3) to get a 3x3
        #  identity matrix.
        self.covariance = dot((eye(3) - dot(K, Ht)), self.covariance)     # Be sure to use dot function
        #
        # Hints:
        # dot(A, B) is the 'normal' matrix product (do not use: A*B).
        # A.T is the transposed of a matrix A (A itself is not modified).
        # linalg.inv(A) returns the inverse of A (A itself is not modified).
        # eye(3) returns a 3x3 identity matrix.

        # pass # Remove this.


    # @staticmethod
    # def compute_center(point_list):
    #     # Safeguard against empty list.
    #     if not point_list:
    #         return (0.0, 0.0)
    #     # If not empty, sum up and divide.
    #     sx = sum([p[0] for p in point_list])
    #     sy = sum([p[1] for p in point_list])
    #     return (sx / len(point_list), sy / len(point_list))
    #
    # @staticmethod
    # def estimate_transform(left_list, right_list, fix_scale=False):
    #     if len(left_list) < 4 or len(right_list) < 4:  # at least two points
    #         return None
    #     # Compute left and right center.
    #     lc = ExtendedKalmanFilter.compute_center(left_list)
    #     rc = ExtendedKalmanFilter.compute_center(left_list)
    #     l_i = [(subtract(l, lc)) for l in
    #            left_list]  # tuple subtract tuple, not only x but also y, l_i is a list
    #     r_i = [(subtract(r, rc)) for r in right_list]
    #     cs, ss, rr, ll = 0.0, 0.0, 0.0, 0.0
    #     for i in range(len(left_list)):
    #         cs += r_i[i][0] * r_i[i][0] + r_i[i][1] * r_i[i][1]
    #         ss += -(r_i[i][0] * l_i[i][1]) + r_i[i][1] * l_i[i][0]
    #         rr += (r_i[i][0] * r_i[i][0]) + (r_i[i][1] * r_i[i][1])
    #         ll += (l_i[i][0] * l_i[i][0]) + (l_i[i][1] * l_i[i][1])
    #
    #     if rr == 0.0 or ll == 0.0:
    #         return None
    #     if fix_scale:
    #         la = 1.0
    #     else:
    #         la = sqrt(rr / ll)
    #
    #     if cs == 0.0 or ss == 0.0:
    #         return None
    #     else:
    #         c = cs / sqrt(cs * cs + ss * ss)
    #         s = ss / sqrt(cs * cs + ss * ss)
    #
    #     tx = rc[0] - la * (c * lc[0] - s * lc[1])
    #     ty = rc[1] - (la * ((s * lc[0]) + (c * lc[1])))
    #
    #     return la, c, s, tx, ty
    #
    #
    # @staticmethod
    # def apply_transform(trafo, p):
    #     la, c, s, tx, ty = trafo
    #     lac = la * c
    #     las = la * s
    #     # print p
    #     x = lac * p[0] - las * p[1] + tx
    #     y = las * p[0] + lac * p[1] + ty
    #     # print (x, y)
    #     return (x, y)
    #
    # def correct_pose(self, trafo):
    #     la, c, s, tx, ty = trafo
    #     x, y = ExtendedKalmanFilter.apply_transform(trafo, (self.state[0], self.state[1]))
    #     theta = self.state[2] + atan2(s, c)
    #     self.state[0:3] = [x, y, theta]

if __name__ == '__main__':
    # Robot constants.
    scanner_displacement = 30.0
    ticks_to_mm = 0.349
    robot_width = 155.0

    # Cylinder extraction and matching constants.
    minimum_valid_distance = 20.0
    depth_jump = 100.0
    cylinder_offset = 90.0
    max_cylinder_distance = 300.0

    # Filter constants.
    control_motion_factor = 0.35  # Error in motor control.
    control_turn_factor = 0.6  # Additional error due to slip when turning.
    measurement_distance_stddev = 200.0  # Distance measurement error of cylinders.
    measurement_angle_stddev = 15.0 / 180.0 * pi  # Angle measurement error.

    # Measured start position.
    initial_state = array([1850.0, 1897.0, 213.0 / 180.0 * pi])
    # Covariance at start position.
    initial_covariance = diag([100.0**2, 100.0**2, (10.0 / 180.0 * pi) ** 2])

    max_deviation = 0   # For the purpose of putting out deviation.

    # Setup filter.
    kf = ExtendedKalmanFilter(initial_state, initial_covariance,
                              robot_width, scanner_displacement,
                              control_motion_factor, control_turn_factor,
                              measurement_distance_stddev,
                              measurement_angle_stddev)
    # Instance with initial state.

    # Read data.
    logfile = LegoLogfile()     # instance
    logfile.read("robot4_motors.txt")
    logfile.read("robot4_scan.txt")
    logfile.read("robot_arena_landmarks.txt")
    logfile.read("robot4_reference.txt")
    reference_cylinders = [l[1:3] for l in logfile.landmarks]      # X, Y and radius in float type
    # Get a new modified array from an existed array.

    # Loop over all motor tick records and all measurements and generate
    # filtered positions and covariances.
    # This is the Kalman filter loop, with prediction and correction.
    states = []
    covariances = []
    matched_ref_cylinders = []      # Magenta points.
    for i in xrange(len(logfile.motor_ticks)):      # For all steps.
        # Prediction.
        control = array(logfile.motor_ticks[i]) * ticks_to_mm
        kf.predict(control)     # Modify the instance kf to expected state.

        # Correction.
        observations = get_observations(
            logfile.scan_data[i],
            depth_jump, minimum_valid_distance, cylinder_offset,
            kf.state, scanner_displacement,
            reference_cylinders, max_cylinder_distance)     # Detected(also the measurement) and real landmark pairs.



        for j in xrange(len(observations)):
            kf.correct(*observations[j])        # Separate elements into several arguments.
        # measurement = observation[j][0], landmark = observation[j][1]
        # Log state, covariance, and matched cylinders for later output.
        states.append(kf.state)
        covariances.append(kf.covariance)
        matched_ref_cylinders.append([m[1] for m in observations])

        # det_car = []
        # real_car = []
        # trafo = []
        # for k in xrange(len(observations)):
        #     distance, angle = observations[k][0]
        #     x, y = distance * cos(angle), distance * sin(angle)
        #     scanner_pose = kf.state + [scanner_displacement * cos(kf.state[2]),
        #                                scanner_displacement * sin(kf.state[2]),
        #                                0.0]
        #     x, y = LegoLogfile.scanner_to_world(scanner_pose, (x, y))
        #     det_car.append(array([x, y]))
        #     real_car.append(array(observations[k][1]))
        #     trafo = ExtendedKalmanFilter.estimate_transform(det_car, real_car, fix_scale=True)
        # if trafo:
        #     # print kf.state[1:4]
        #     kf.correct_pose(trafo)

    # Write all states, all state covariances, and matched cylinders to file.
    f = open("kalman_prediction_and_correction.txt", "w")
    for i in xrange(len(states)):
        # Output the center of the scanner, not the center of the robot.
        print >> f, "F %f %f %f" % \
            tuple(states[i] + [scanner_displacement * cos(states[i][2]),
                               scanner_displacement * sin(states[i][2]),
                               0.0])

        # # For the purpose of putting out deviation.
        # p0 = tuple(states[i] + [scanner_displacement * cos(states[i][2]),
        #                        scanner_displacement * sin(states[i][2]),
        #                        0.0])
        # deviation = distance(p0, logfile.reference_positions[i])
        # if deviation > max_deviation:
        #     max_deviation = deviation

        # Convert covariance matrix to angle stddev1 stddev2 stddev-heading form
        e = ExtendedKalmanFilter.get_error_ellipse(covariances[i])
        print >> f, "E %f %f %f %f" % (e + (sqrt(covariances[i][2,2]),))
        # Also, write matched cylinders.
        write_cylinders(f, "W C", matched_ref_cylinders[i])        

    f.close()

    print max_deviation     # For the purpose of putting out deviation.
