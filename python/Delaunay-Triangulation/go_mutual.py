import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.spatial import Delaunay
from PIL import Image
from collections import defaultdict
from sklearn.metrics import mean_squared_error
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import MI

class Delaunay_d(object):

	def __init__(self, d, iteration):
		self.datafile = 'go_data_trail'
		self.file1 = "trail1.gif"
		self.file2 = "trail2.gif"
		self.img0 = Image.open(self.file1).convert('L')
		self.img1 = Image.open(self.file2).convert('L')
		self.width, self.height = self.img0.size
		self.globalbestPoint = [] # 0.3922227 , 0.19259966 , 0.23557716 , 0.43344293

		self.d, self.iteration = d, iteration
		self.Mn, self.gn, self.vn, self.bestRho = np.inf, np.inf, np.inf, np.inf
		self.bestPoint = []
		self.vertex_values = defaultdict(list)
		self.volumes = defaultdict(list)
		self.vertices = [[0, 0, 0, 0],
						[0, 0, 0, 1],
						[0, 0, 1, 0],
						[0, 1, 0, 0],
						[1, 0, 0, 0],
						[0, 0, 1, 1],
						[0, 1, 1, 0],
						[1, 1, 0, 0],
						[0, 1, 0, 1],
						[1, 0, 1, 0],
						[1, 0, 0, 1],
						[0, 1, 1, 1],
						[1, 1, 1, 0],
						[1, 0, 1, 1],
						[1, 1, 0, 1],
						[1, 1, 1, 1],
						[0.3, 0.5, 0.6, 0.9]]
		self.T = Delaunay(self.vertices, incremental=True, qhull_options='Qt Fv')
		self.T.coplanar
		self.simplices = np.array(self.vertices)[self.T.simplices]
		for vertex in self.vertices:
			# print self.func(vertex), vertex
			if self.Mn > self.func(vertex):
				self.globalbestPoint = vertex
			self.Mn = min(self.Mn, self.func(vertex))
		self.error = [0 for _ in xrange(self.iteration)]

	def ScaleRotateTranslate(self, image, angle, center=None, new_center=None, scale=None):
		if center is None:
			return image.rotate(angle)
		# angle = angle / 180.0 * math.pi
		nx, ny = x, y = center
		sx = sy = 1.0
		if new_center:
			nx, ny = new_center
		if scale:
			sx, sy = scale
		cosine = np.cos(np.deg2rad(angle))
		sine = np.sin(np.deg2rad(angle))
		a = cosine / sx
		b = sine / sx
		c = x-nx*a-ny*b
		d = -sine / sy
		e = cosine / sy
		f = y - nx*d - ny*e
		return image.transform(image.size, Image.AFFINE, (a,b,c,d,e,f), resample=Image.BICUBIC)

	def MSE(self, imageA, imageB):
		return MI.mutual_information(imageA, imageB)

	def func(self, p):
		# print tuple(p), self.vertex_values
		idx = tuple(p)
		if not self.vertex_values[idx]:
			x = self.width / 2.0 + (40 * p[0] - 20.0)
			y = self.height / 2.0 + (40 * p[1] - 20.0)
			# r = 0.24273422 * 40.0 - 20.0
			# s = 0.42331118 * 0.4 + 0.8
			r = 0.5 * 40.0 - 20.0
			s = 0.5 * 0.4 + 0.8
			x0, y0, x1, y1 = math.floor(x), math.floor(y), math.ceil(x), math.ceil(y)

			# clockwise
			if x0 == x1 and y0 != y1:
				v0 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x0, y0), (s, s)))
				v1 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x0, y1), (s, s)))
				rs = (y1-y)/(y1-y0)*v0 + (y-y0)/(y1-y0)*v1
			elif x0 != x1 and y0 == y1:
				v0 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x0, y0), (s, s)))
				v1 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x1, y0), (s, s)))
				rs = (x1-x)/(x1-x0)*v0 + (x-x0)/(x1-x0)*v1
			elif x0 == x1 and y0 == y1:
				rs = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x, y), (s, s)))
			else:
				v0 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x0, y0), (s, s)))
				v1 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x1, y0), (s, s)))
				v2 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x0, y1), (s, s)))
				v3 = self.MSE(self.img0, self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x1, y1), (s, s)))
				l1 = (x1-x)/(x1-x0)*v0 + (x-x0)/(x1-x0)*v1
				l2 = (x1-x)/(x1-x0)*v2 + (x-x0)/(x1-x0)*v3
				rs = (y1-y)/(y1-y0)*l1 + (y-y0)/(y1-y0)*l2
			# trans_B = self.ScaleRotateTranslate(self.img1, r, (self.width / 2.0, self.height/2.0), (x, y), (s, s))
			# rs = self.MSE(self.img0, trans_B)
			self.vertex_values[idx] = 1.0 / rs
		else:
			return self.vertex_values[idx]

		return 1.0 / rs

	def volume(self, simplex):
		idx = tuple(map(tuple, simplex))
		temp = self.volumes[idx]
		if temp:
			return temp
		matrix = np.copy(simplex)
		matrix = np.insert(matrix, 0, 1, axis=1)
		rs = 1.0/2/3/4 * abs(np.linalg.det(matrix))
		self.volumes[idx] = rs
		return rs

	def calRho(self):
		# print "In this iteration has: ", len(np.array(self.vertices)[self.T.simplices])
		for simplex in np.array(self.vertices)[self.T.simplices]:
			# if simplex visited before, then skip
			# idx = tuple(map(tuple, simplex))
			# if self.volumes[idx]:
			# 	continue
			v = self.volume(simplex)
			if v == 0.0 or np.isinf(v):
				continue
			if v < self.vn:
				self.vn = v
				# q = 34.346004
				q = 27.152900397563425
				# print 'hahahahahahahahahaahahahaha'
				self.gn = np.sqrt(q * self.vn * np.log(1.0 / self.vn))
			ff = np.sum( self.func(s) for s in simplex ) / (self.d+1)
			rho = v / (ff - self.Mn + self.gn)**2 
			if rho > self.bestRho:
				self.bestPoint = np.sum(simplex, axis=0) / (self.d + 1)
				# print 'Best Point this time is: ', self.bestPoint
				self.bestRho = rho
		return ff

	def iter(self):
		for n in range(self.iteration):
			self.bestRho = 0.0
			self.calRho()
			l = len(self.simplices)
			print "###################", n, "###################", "simplices: ", l
			if self.Mn > self.func(self.bestPoint):
				self.globalbestPoint = self.bestPoint
				# print "error: ", self.func(self.bestPoint), "Best Point: ", self.bestPoint, " Value: ", self.func(self.bestPoint)
			self.Mn = min(self.Mn, self.func(self.bestPoint))
			self.error[n] += self.Mn
			self.vertices.append(self.bestPoint)
			self.T.add_points([self.bestPoint])
			self.simplices = np.array(self.vertices)[self.T.simplices]

			# p = self.globalbestPoint
			# x = self.width / 2.0 + (100 * p[0] - 50.0)
			# y = self.height / 2.0 + (100 * p[1] - 50.0)
			# r = p[2] * 40.0 - 20.0
			# s = p[3] * 0.4 + 0.8
			# print "error: ", self.error[n],"Best Global: ", self.globalbestPoint,"\nBest Point: ", (x,y,r,s), " Value: ", self.func(self.bestPoint)
			print "error: ", self.error[n],"Best Global: ", self.globalbestPoint,"\nBest Point: ", self.bestPoint, " Value: ", self.func(self.bestPoint)
		np.savetxt(self.datafile, self.vertices)

	def plot(self):
		self.vertices = np.array(self.vertices)
		plt.triplot(self.vertices[:,0], self.vertices[:,1])
		# plt.plot(self.vertices[:,0], self.vertices[:,1], '.')
		plt.plot(self.globalbestPoint[0], self.globalbestPoint[1], 'o')
		plt.show()

	def show(self):
		img0 = Image.open(self.file1)
		img1 = Image.open(self.file2)
		p = self.globalbestPoint
		# p = [  0.42237259 , 0.22260757 , 0.24273422,  0.42331118]
		x = self.width / 2.0 + (40 * p[0] - 20.0)
		y = self.height / 2.0 + (40 * p[1] - 20.0)
		# r = 0.24273422 * 40.0 - 20.0
		# s = 0.42331118 * 0.4 + 0.8
		r = 0.5 * 40.0 - 20.0
		s = 0.5 * 0.4 + 0.8
		print (x,y,r,s)
		img1 = self.ScaleRotateTranslate(img1, r, (self.width / 2.0, self.height/2.0), (x, y), (s, s))
		# img1 = self.ScaleRotateTranslate(self.img1, r, None, (x, self.height-y), (s, s))
		img0.show()
		img1.show()
		img1.save('rs-trail1.png')
		rs = (Image.blend(img1, img0, 0.5))
		rs.save('rs-trail2.png')
		rs.show()

	def test2(self):
		file = open('go_data2')
		vertices = []
		for line in file.readlines():
			temp = []
			for i in line.split(' '):
				temp.append(float(i))
			vertices.append(temp)
		m = np.inf
		for i, arr in enumerate(vertices):
			arr = np.array(arr)
			val = self.func(arr)
			if m > val:
				m = val
				print i, arr, val

	def test3(self):
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		X = np.arange(-1, 1, 0.1)
		Y = np.arange(-1, 1, 0.1)
		X, Y = np.meshgrid(X, Y)
		zs = np.array( [self.func([x,y]) for x,y in zip(np.ravel(X), np.ravel(Y))])
		Z = zs.reshape(X.shape)
		surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
		                       linewidth=0, antialiased=False)
		ax.zaxis.set_major_locator(LinearLocator(10))
		ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		fig.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()

obj = Delaunay_d(4, 5000)
obj.iter()
obj.show()
obj.plot()
# obj.test3()


#best point for knee: 0.59379093  0.64182245  0.39307772  0.78009655	mim: 0.639451565899
#best point for 1.jpg: 0.42237259 , 0.22260757 , 0.24273422,  0.42331118	
# #best point for trail: 0.5250707   0.47523014  0.64535374  0.47940955	
################### 4999 ################### simplices:  124452go 
# error:  1.70939344944 Best Global:  [ 0.5250707   0.47523014  0.64535374  0.47940955]
# Best Point:  [ 0.92731563  0.40221884  0.57274906  0.82734122]  Value:  6.00865842989
# (151.00282783910586, 129.00920569681983, 0.0, 1.0)