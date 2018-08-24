# Compute the derivative of a scan.
# 03_b_scan_derivative
# Claus Brenner, 09 NOV 2012
from pylab import *
from lego_robot import *
import matplotlib.pyplot as plt


# Find the derivative in scan data, ignoring invalid measurements.
def compute_derivative(scan, min_dist):
    jumps = [ 0 ]
    for i in xrange(1, len(scan) - 1):
        # --->>> Insert your code here.
        # Compute derivative using formula "(f(i+1) - f(i-1)) / 2".
        # Do not use erroneous scan values, which are below min_dist.
        #jumps.append(i%20 * 10) # Replace this line, append derivative instead.
        l = scan[i-1]
        r = scan[i+1]
        if l > min_dist and r > min_dist:
            deriv = (r - l)/2.0
            jumps.append(deriv)
        else:
            jumps.append(0)
    jumps.append(0)
    return jumps


if __name__ == '__main__':

    minimum_valid_distance = 20.0

    # Read the logfile which contains all scans.
    logfile = LegoLogfile()
    logfile.read("robot4_scan.txt")

    # Pick one scan.
    scan_no = 8
    scan = logfile.scan_data[scan_no]

    # Compute derivative, (-1, 0, 1) mask.
    der = compute_derivative(scan, minimum_valid_distance)

    # Plot scan and derivative.
    title("Derivative of One Scan")
    #title("Plot of scan %d" % scan_no)
    # plot(scan)
    # plot(der)
    plt.plot(scan)
    plt.plot(der)


    show()
