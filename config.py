""" File with auxiliary functions and parameters
"""

math = {
		"eps": 1e-7,
	   }

def makeGrid(a: float, b: float, n: int, r:float):
	""" Generating regular and adaptive grid partitions

	a, b - boundaries of an arbitrary segment
	n - number of split segments
	r - "discharge" ratio
	"""
	h1 = (b - a)/sum(r**i for i in range(n))
	hList = [h1 * r**i for i in range(n)]
	xList = [sum(hList[i] for i in range(j)) for j in range(n)]

	return hList, xList