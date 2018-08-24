from lego_robot import *
import math
import numpy as np
import matplotlib.pyplot as plt
def distance(p0, p1):
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)
logfile = LegoLogfile()
logfile.read("test.txt")
landmarks = logfile.landmarks
#for landmark in landmarks:
#    print landmark[1:4]
landmark_l = [ landmark[1:3] for landmark in landmarks]
print landmark_l

logfile.read("ekf_slam_correction.txt")
wcs = logfile.world_cylinders
diss = []
for wc in wcs:
    dis = []    
    dis.append(distance(wc[0],landmark_l[0]))
    dis.append(distance(wc[1],landmark_l[4]))
    dis.append(distance(wc[2],landmark_l[1]))
    dis.append(distance(wc[3],landmark_l[2]))
    dis.append(distance(wc[4],landmark_l[3]))
    dis.append(distance(wc[5],landmark_l[5]))
    diss.append(dis)
# print diss
diss_array = np.array(diss)
#print diss_array[:,0]
plt.figure()
for i in range(0,6):
    plt.plot(diss_array[:,i],label="cylinder_%d"%i,lw=2)
plt.xlabel("Steps")
plt.ylabel("Error(mm)")
plt.xlim(0,300)
plt.ylim(0,100)
plt.title("Euclidean Error Distance")
plt.grid(linestyle="--")
plt.legend()
plt.show()
