import math
import numpy as np

def fun(x):
    f = 3*math.exp(x) + x**2 + 1.0/(x+2.0)
    return f

def ExactIntegration(x):
    f = 3*math.exp(x) + x**3/3 + np.log(x+2)
    return f

def gauss(ngp):
    """
    Get Gauss points in the parent element domain [-1, 1] and 
    the corresponding weights.
    
    Args: 
        ngp : (int) number of Gauss points.
        
    Returns: w,gp
        w  : weights.
        gp : Gauss points in the parent element domain.
    """
    if ngp == 1:
        gp = [0]
        w = [2]
    elif ngp == 2:
        gp = [-0.57735027, 0.57735027]
        w = [1, 1]
    elif ngp == 3:
        gp = [-0.7745966692, 0.7745966692, 0.0]
        w = [0.5555555556, 0.5555555556, 0.8888888889]
    elif ngp == 4:
        gp = [-0.8611363116, 0.8611363116, -0.3399810436, 0.3399810436]
        w = [0.3478548451, 0.3478548451, 0.6521451549, 0.6521451549]
    return w, gp

def lobatto(ngp):
    """
    Get Gauss Lobato points in the parent element domain [-1, 1] and 
    the corresponding weights.
    
    Args: 
        ngp : (int) number of Gauss points.
        
    Returns: w,gp
        w  : weights.
        gp : Gauss points in the parent element domain.
    """
    if ngp == 3:
        gp = [0, -1.0, 1.0]
        w = [1.3333333333, 0.3333333333, 0.3333333333]
    elif ngp == 4:
        gp = [-0.4472135955, 0.4472135955, -1.0, 1.0]
        w = [0.8333333333, 0.8333333333, 0.1666666667, 0.1666666667]
    return w, gp

def integrate(f, ngp, gp, w):
    return sum(f(gp[i])*w[i] for i in range(ngp))

if __name__ == "__main__":
    Ie = ExactIntegration(1)-ExactIntegration(-1)

    [w , gp] = gauss(1)   # extract Gauss points and weights
    I1 = integrate(fun, 1, gp, w)
    Error_1 = (I1-Ie)/Ie  # Relative error

    [w , gp] = gauss(2)   # extract Gauss points and weights
    I2 = integrate(fun, 2, gp, w)
    Error_2 = (I2-Ie)/Ie  # Relative error

    [w , gp] = gauss(3)   # extract Gauss points and weights
    I3 = integrate(fun, 3, gp, w)
    Error_3 = (I3-Ie)/Ie  # Relative error

    [w , gp] = lobatto(3)   # extract Gauss-Lobatto points and weights
    I3L = integrate(fun, 3, gp, w)
    Error_3L = (I3L-Ie)/Ie  # Relative error

    [w , gp] = lobatto(4)   # extract Gauss-Lobatto points and weights
    I4L = integrate(fun, 4, gp, w)
    Error_4L = (I4L-Ie)/Ie  # Relative error

    print("Exact integration: ", Ie)
    print("1-point Gauss Quadrature: ", I1, "  Relative error = ", Error_1)
    print("2-point Gauss Quadrature: ", I2, "  Relative error = ", Error_2)
    print("3-point Gauss Quadrature: ", I3, "  Relative error = ", Error_3)
    print("3-point Gauss Lobatto Quadrature: ", I3L, "  Relative error = ", Error_3L)
    print("4-point Gauss Lobatto Quadrature: ", I4L, "  Relative error = ", Error_4L)
