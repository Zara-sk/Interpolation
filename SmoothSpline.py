from config import math


class smoothSpline:

	def __init__(self, p, w, points: dict):
		""" Initializing the Smoothing Spline interface class
		
		points - list of values of the variable x
		values - corresponding x table values of the function
		segmentCount - number of segments
		alpha - expansion coefficients for the functional
		s1, s2, s3 - coefficients of the matrix diagonals
		p, w - smoothness coefficients of the approximation
		"""
		self.values = []
		self.points = []
		for key in points:
			self.points.append(key)
			self.values.append(points[key])

		self.segmentCount = len(points) - 1
		self.alpha = [0]*(self.segmentCount+1)
		self.s1    = [0]*(self.segmentCount+1)
		self.s2    = [0]*(self.segmentCount+1)
		self.s3    = [0]*(self.segmentCount+1)
		self.p = p
		self.w = w

		self.calculate()


	def calculate(self):

		# Calculating LS coefficients
		for i in range(self.segmentCount):

			ksi1 = self.findKsi(self.points[i], i)
			ksi2 = self.findKsi(self.points[i+1], i)

			# finding basic functions
			fi11 = self.findBasicFunction(1, ksi1)
			fi12 = self.findBasicFunction(2, ksi1)
			fi21 = self.findBasicFunction(1, ksi2)
			fi22 = self.findBasicFunction(2, ksi2)

			h = self.points[i+1] - self.points[i]

			# contributing to diagonal elements
			self.s1[i+1] += (1 - self.p) * self.w * (fi11*fi12 + fi21*fi22) - (1/h*self.p)
			self.s2[i]   += (1 - self.p) * self.w * (fi11*fi11 + fi21*fi21) + (1/h*self.p)
			self.s2[i+1] += (1 - self.p) * self.w * (fi12*fi12 + fi22*fi22) + (1/h*self.p)
			self.s3[i]   += (1 - self.p) * self.w * (fi12*fi11 + fi22*fi21) - (1/h*self.p)

			# сontribution to the right-hand side vector:
			self.alpha[i]   += (1 - self.p) * self.w * (fi11*self.values[i] + fi21*self.values[i+1])
			self.alpha[i+1] += (1 - self.p) * self.w * (fi12*self.values[i] + fi22*self.values[i+1])

		# Run-through (1st step)
		for i in range(1, self.segmentCount + 1):
			self.s2[i]    -= self.s1[i] / self.s2[i-1] * self.s3[i - 1]
			self.alpha[i] -= self.s1[i] / self.s2[i-1] * self.alpha[i-1]

		# Run-through (2nd step)
		self.alpha[self.segmentCount] /= self.s2[self.segmentCount]
		for i in range(self.segmentCount - 1, -1, -1):
			self.alpha[i] = (self.alpha[i] - self.alpha[i+1] * self.s3[i] / self.s2[i])


	def findKsi(self, x: float, segmentNumber: int):
		return 2 * (x - self.points[segmentNumber]) /\
				   (self.points[segmentNumber+1] - self.points[segmentNumber]) - 1


	def findBasicFunction(self, index: int, ksi: float):
		return 0.5 * (1 - ((-1)**index)*ksi)


	def __call__(self, value):
		"""	Evaluating the spline value at the current point
		"""
		for i in range(self.segmentCount):

			if (value > self.points[i] and value < self.points[i+1]) or \
			   (abs(value - self.points[i]) < math["eps"]) or \
			   (abs(value - self.points[i+1]) < math["eps"]):

				h = self.points[i+1] - self.points[i]
				ksi = self.findKsi(value, i)
				result = {}
				result[0] = self.alpha[i] * self.findBasicFunction(1, ksi) + self.alpha[i+1] * self.findBasicFunction(2, ksi)
				result[1] = self.alpha[i] * (-0.5) + self.alpha[i+1] / h 
				result[2] = 0
				return result

		raise("Chosen point is out of any segment...")