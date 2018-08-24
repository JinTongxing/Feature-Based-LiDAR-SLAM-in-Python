# Histogram implementation of a bayes filter - combines
# convolution and multiplication of distributions, for the
# movement and measurement steps.
# 06_d_histogram_filter
# Claus Brenner, 28 NOV 2012
from pylab import plot, show, ylim
from distribution import *

def move(distribution, delta):
    """Returns a Distribution that has been moved (x-axis) by the amount of
       delta."""
    return Distribution(distribution.offset + delta, distribution.values)



# --->>> Copy your convolve(a, b) and multiply(a, b) functions here.
def convolve(a, b):
    """Convolve distribution a and b and return the resulting new distribution."""
    dist_lst = []
    offs = a.offset + b.offset
    for a_val in a.values:
        res = []
        for b_val in b.values:
            res.append(a_val * b_val)
        dist_lst.append(Distribution(offs, res))
        offs += 1
    c = Distribution.sum(dist_lst)
    return c

def multiply(a, b):
    """Multiply two distributions and return the resulting distribution."""
    start = min([a.start(), b.start()]) # return the min in a list
    stop = max([a.stop(), b.stop()])
    mul_dist = []
    for ix in range(start, stop):
        mul_dist.append(a.value(ix) * b.value(ix))
    d = Distribution(start, mul_dist)
    Distribution.normalize(d)
    return d


if __name__ == '__main__':
    arena = (0,2200)

    # Start position. Exactly known - a unit pulse.
    start_position = 10
    position = Distribution.unit_pulse(start_position)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         linestyle='steps')

    # Movement data.
    controls  =    [ 20 ] * 100

    # Measurement data. Assume (for now) that the measurement data
    # is correct. - This code just builds a cumulative list of the controls,
    # plus the start position.
    p = start_position
    measurements = []
    for c in controls:
        p += c
        measurements.append(p)  # [30, 50, 70, 90...]

    # This is the filter loop.
    for i in xrange(len(controls)):
        # Move, by convolution. Also termed "prediction".
        control = Distribution.triangle(controls[i], 10)
        position = convolve(position, control)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='b', linestyle='steps')

        # Measure, by multiplication. Also termed "correction".
        measurement = Distribution.triangle(measurements[i], 50)
        position = multiply(position, measurement)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='r', linestyle='steps')

    show()
