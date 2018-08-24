import numpy as np
from scipy import linalg
from lego_robot import *
from matplotlib import pyplot as plt
logfile = LegoLogfile()
logfile.read("ekf_slam_correction.txt")
world_ellipse = logfile.world_ellipses
ss = []
for we in world_ellipse:
    s = []
    for cylinder in we:
        s.append(np.pi*cylinder[1]*cylinder[2])
    ss.append(s)
ss_a = np.array(ss)
#print ss_a[:,0]
plt.figure()
for i in range(6):
    plt.plot(ss_a[:,i],label="cylinder_%d"%i,lw = 2)
plt.xlim(0,40)
plt.ylim(0,)
plt.xlabel('Steps')
plt.ylabel('Error Ellipses Area/mm2')
plt.title('Uncertainty Chage Over Time')
plt.legend()
plt.grid(linestyle="--")
plt.show()
