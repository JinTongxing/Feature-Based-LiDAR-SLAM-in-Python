import numpy as np
c = []
a = (1, 2)
b = (3, 4)
print tuple(np.subtract(a,b))


def compute_center(point_list):
    # Safeguard against empty list.
    if not point_list:
        return (0.0, 0.0)
    # If not empty, sum up and divide.
    sx = sum([p[0] for p in point_list])
    sy = sum([p[1] for p in point_list])
    return (sx / len(point_list), sy / len(point_list))

