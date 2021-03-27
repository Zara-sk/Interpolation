
if __name__ == '__main__':
	
	import matplotlib.pyplot as plt
	from numpy import arange, sin, cos

	from CubicSpline import cubicSpline 
	from SmoothSpline import smoothSpline

	"""Evaluating the inaccuracy"""

	def Splining(func, h: float, w: float, p: float, dif: int):
		table1 = {i: func(i) for i in arange(0, 10, h  )}
		table2 = {i: func(i) for i in arange(0, 10, h/2)}
		table3 = {i: func(i) for i in arange(0, 10, h/4)}

		spline1 = cubicSpline(1, table1)
		spline2 = cubicSpline(1, table2)
		spline3 = cubicSpline(1, table3)
		
		y1, y2, y3 = [], [], []
		
		for i in range(10):
			y1.append(spline1(w*i + p)[dif])
			y2.append(spline2(w*i + p)[dif])
			y3.append(spline3(w*i + p)[dif])

		return y1, y2, y3


	def test_derivate_p0(func, h: float, w: float, p: float, label: str):
		
		y1, y2, y3 = Splining(func, h, w, p, 0)

		print("Function:   ", label)
		print(f"for h={ h }:  ", max(abs(y1[i] - func(w*i + p)) for i in range(10)))
		print(f"for h={h/2}: ",  max(abs(y2[i] - func(w*i + p)) for i in range(10)))
		print(f"for h={h/4}:",   max(abs(y3[i] - func(w*i + p)) for i in range(10)))
		print("="*40)


	def test_derivate_p1(func, func_dif1, h: float, w: float, p: float, label: str):
		
		y1, y2, y3 = Splining(func, h, w, p, 1)

		print("Function:   ", label, "dy/dx")
		print(f"for h={ h }:  ", max(abs(y1[i] - func_dif1(i, w, p)) for i in range(10)))
		print(f"for h={h/2}: ",  max(abs(y2[i] - func_dif1(i, w, p)) for i in range(10)))
		print(f"for h={h/4}:",   max(abs(y3[i] - func_dif1(i, w, p)) for i in range(10)))
		print("="*40)
		

	def test_derivate_p2(func, func_dif2, h: float, w: float, p: float, label: str):
		
		y1, y2, y3 = Splining(func, h, w, p, 2)

		print("Function:   ", label, "d2y/dx2")
		print(f"for h={ h }:  ", max(abs(y1[i] - func_dif2(i, w, p)) for i in range(10)))
		print(f"for h={h/2}: ",  max(abs(y2[i] - func_dif2(i, w, p)) for i in range(10)))
		print(f"for h={h/4}:",   max(abs(y3[i] - func_dif2(i, w, p)) for i in range(10)))
		print("="*40)


	"""Test #1: f(x) = x^2"""
	x2_0 = lambda x: x**2 
	x2_1 = lambda x, a, b: 2*x*a**2+2*a*b
	x2_3 = lambda x, a, b: 2*a**2

	test_derivate_p0(x2_0, 0.1, 1, 0.14, "x^2")
	# for h=0.1:   0.04183178196587348
	# for h=0.05:  0.038000000000209866
	# for h=0.025: 0.014636729734741562
	
	test_derivate_p1(x2_0, x2_1, 0.1, 0.78, 0.49, "x^2")
	# for h=0.1:   2.4043999999997503
	# for h=0.05:  2.3543999999979217
	# for h=0.025: 2.329399999987027

	test_derivate_p2(x2_0, x2_3, 0.1, 0.88, 0.14, "x^2")
	# for h=0.1:   0.7153016151377547
	# for h=0.05:  0.4532619104571509
	# for h=0.025: 0.4518608589528599
	

	"""Test #2: f(x) = x^3"""
	x3_0 = lambda x: x**3
	x3_1 = lambda x, a, b: 3*a*(a*x + b)**2
	x3_2 = lambda x, a, b: 6*a**2*(a*x + b)

	test_derivate_p0(x3_0, 0.1, 1, 0.2, "x^3")
	# for h=0.1:   2.466041006996761
	# for h=0.05:  1.3062500005738684
	# for h=0.025: 0.6715312499968604

	test_derivate_p1(x3_0, x3_1, 0.005, 0.74, 0.13, "x^3")
	# for h=0.005:   15.702997999019644
	# for h=0.0025:  15.647110507876079
	# for h=0.00125: 15.619157370000877

	test_derivate_p2(x3_0, x3_2, 0.1, 0.96, 0.38, "x^3")
	# for h=0.1:   4.243030664220711
	# for h=0.05:  4.243007999347881
	# for h=0.025: 4.243007999845489


	"""Test #3: f(x) = x^4"""
	x4_0 = lambda x: x**4
	x4_1 = lambda x, a, b: 4*a*(a*x + b)**3
	x4_2 = lambda x, a, b: 12*a**2*(a*x + b)**2

	test_derivate_p0(x4_0, 0.005, 1, 0.28, "x^4")
	# for h=0.005:   2.5687878925500627
	# for h=0.0025:  1.2880837669526954
	# for h=0.00125: 0.6449647885410741

	test_derivate_p1(x4_0, x4_1, 0.005, 0.1, 0.91, "x^4")
	# for h=0.1:   3.193923599999476
	# for h=0.05:  2.534668600000097
	# for h=0.025: 2.2628892749994987

	test_derivate_p2(x4_0, x4_2, 0.1, 1, 0.22, "x^4")
	# for h=0.1:   0.005450198395806183
	# for h=0.05:  0.0022127018922201147
	# for h=0.025: 5.000468240723421e-05


	"""Test #4: f(x) = sinx"""
	sinx_0 = lambda x: sin(x)
	sinx_1 = lambda x, a, b: a*cos(a*x+b)
	sinx_2 = lambda x, a, b: -a**2*sin(a*x+b)

	test_derivate_p0(sinx_0, 0.1, 0.9, 0.1, "sin(x)")
	# for h=0.1:   0.04426814024894654
	# for h=0.05:  0.02349664272701435
	# for h=0.025: 0.012085937041755979

	test_derivate_p1(sinx_0, sinx_1, 0.005, 1, 0.41, "sinx")
	# for h=0.005:   0.49081372522355826
	# for h=0.0025:  0.4921823194986039
	# for h=0.00125: 0.4928662936090835

	test_derivate_p2(sinx_0, sinx_2, 0.01, 0.9, 0.23, "sinx")
	# for h=0.1:   0.00020914585513298611
	# for h=0.05:  8.60865765999197e-05
	# for h=0.025: 2.029435198114271e-06

	if 0:
		table = {i: sinx_0(i) for i in arange(0, 10, 0.005)}
		y1, y2, y3 = Splining(sinx_0, 0.005, 0.74, 0.22, 0)
		fig = plt.figure(figsize=[8,6])
		# plt.plot([(0.74*x + 0.22) for x in range(10)], [sinx_0(x, 0.74, 0.22) for x in range(10)], label="Точная")
		plt.plot([(0.74*x + 0.22) for x in range(10)], [sinx_0(0.74*x + 0.22) for x in range(10)], label="Точная")
		plt.plot([(0.74*x + 0.22) for x in range(10)], [y1[x] for x in range(10)], label="Апроксимирующая h = 0.1")
		plt.plot([(0.74*x + 0.22) for x in range(10)], [y2[x] for x in range(10)], label="Апроксимирующая h = 0.05")
		plt.plot([(0.74*x + 0.22) for x in range(10)], [y3[x] for x in range(10)], label="Апроксимирующая h= 0.0025")
		plt.legend()
		# plt.title("Первая производная функции x^4")
		plt.show()