rh = RecombinationHistory(Yp=0.0)
x_saha = range(rh.x_start, -5, 1000)

@test isapprox([Xe_saha_equation(rh, i) for i in x_saha], [Xe_saha_equation_with_He(rh, i) for i in x_saha], rtol = 1e-6)

rh = RecombinationHistory()
x = (rh.x_start, rh.x_end, 2000)
g = visibility_function_of_x(rh)

