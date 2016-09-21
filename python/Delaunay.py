import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.spatial import Delaunay

def func(x1, x2, x3):
	return np.square(x1+x2-1.0) + np.square(x1-x2-2.0) + np.sqaure(float(x1-x3))

vertices = [[1,1,1], [1,-1,1], [-1,-1,1], [-1,1,1], [1,1,-1], [1,-1,-1], [-1,-1,-1], [-1,1,-1]]
points = np.array(vertices)
tri = Delaunay(points)

print points[tri.simplices]
#plt.show()
