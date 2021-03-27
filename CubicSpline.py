from config import math, makeGrid


class cubicSpline:

	def __init__(self, r, points: dict):
		"""
		Initialization of the cubic interpolation spline interface class
		
		points - list of values of the variable x
		values - corresponding x table values of the function
		a, b, c, d - dictionaries of partition coefficient values
		f - vector of the right parts of linear system
		"""
		self.values = []
		self.points = []
		for key in points:
			self.points.append(key)
			self.values.append(points[key])
		self.segmentsCount = len(points)  - 1
		self.a = dict()
		self.b = dict()
		self.c = dict()
		self.d = dict()
		self.f = dict()
		self.calculateCfcs(r)

	def calculateCfcs(self, r):
		"""Calculating coefficients
		"""

		hList, _ = makeGrid(self.points[0], self.points[-1], self.segmentsCount, r)

		# Formation of the primary LS view:
		for i in range(self.segmentsCount - 1):

			self.a[i+1] = hList[i] # replaces Alpha into this act
			self.b[i]   = 2*(hList[i] + hList[i+1]) # Beta
			self.d[i]   = hList[i+1] # Gamma
			self.f[i]   = 3*((self.values[i+2] - self.values[i+1]) / hList[i+1] 
				          - (self.values[i+1] - self.values[i]) / hList[i]) # Фи

		# run-through method: lite-version of Gaussian elimination (1st step)
		for i in range(1, self.segmentsCount - 1):
			self.b[i] -= self.a[i] / self.b[i - 1] * self.d[i - 1] # Сигма
			self.f[i] -= self.a[i] / self.b[i - 1] * self.f[i - 1] # Пси
		
		# (2nd step)
		self.c[self.segmentsCount-1] = self.f[self.segmentsCount-2] / self.b[self.segmentsCount-2]
		for i in range(self.segmentsCount - 2, 0, -1):
			self.c[i] = (self.f[i-1] - self.c[i+1]*self.d[i-1]) / self.b[i-1]
		
		self.c[0] = 0.0 # Краевое условие нулевой кривизны
		
		# Coefficients recovery
		for i in range(0, self.segmentsCount - 1):
			self.a[i] = self.values[i]
			self.b[i] = (self.values[i+1] - self.values[i])/hList[i] - (2*self.c[i] + self.c[i+1]) / 3
			self.d[i] = (self.c[i+1] - self.c[i]) / (3*hList[i])

		self.a[self.segmentsCount-1] = self.values[self.segmentsCount-1]
		self.b[self.segmentsCount-1] = (self.values[self.segmentsCount] - self.values[self.segmentsCount-1]) /\
			hList[self.segmentsCount - 1] - (2/3)*self.c[self.segmentsCount-1]*hList[self.segmentsCount - 1] 
		self.d[self.segmentsCount-1] = -self.c[self.segmentsCount-1]/(3*hList[i])

	def __call__(self, value):
		"""	Evaluating the spline value at the current point
		"""
		for i in range(self.segmentsCount):
			
			if (value > self.points[i] and value < self.points[i+1]) or \
			   (abs(value - self.points[i]) < math["eps"]) or \
			   (abs(value - self.points[i+1]) < math["eps"]):

				delta = value - self.points[i]
				result = {}
				result[0] = self.a[i] + self.b[i] * delta + self.c[i] * delta**2 + self.d[i] * delta**3
				result[1] = self.b[i] + 2 * self.c[i] * delta + 3 * self.d[i] * delta**2
				result[2] = 2 * self.c[i] + 6 * self.d[i] * delta

				return result

		raise Exception("Chosen point is out of any segment...")
