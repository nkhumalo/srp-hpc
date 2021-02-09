import numpy
import math

ll = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]
pi   = math.acos(-1.0)
pi12 = pi/2.0
pi14 = pi12/2.0
pi18 = pi14/2.0

# See Jacobi rotation at: https://mathworld.wolfram.com/JacobiRotationMatrix.html
#
#  For a matrix A: cot(2 Phi) = (A_qq - A_pp)/(2 A_qp)                         (1)
#
#  The rotation matrix is R  = [[cos(phi), sin(phi)] [-sin(phi), cos(phi)]]    (2)
#
#  From Eq.(1) we have 
#
#     tan(2 phi) = (2 A_qp)/(A_qq - A_pp)
#
#     phi = 1/2 * arctan(2 A_qp/(A_qq - A_pp))
#
#     c = cos(phi) = cos(1/2 * arctan(2 A_qp/(A_qq - A_pp)))
#     s = sin(phi) = sin(1/2 * arctan(2 A_qp/(A_qq - A_pp)))
#
# So, we have E(c2) where we assume that E is a quadratic function.
# For a given matrix A and damping factor dd (which is applied only to the off-
# diagonal elements we have that
#
#     c2 = cos^2(1/2 * arctan(2 * dd * A_qp/(A_qq - A_pp)))
#
# Then when we minimize E(c2) finding cmin we get
#
#     sqrt(cmin) = cos(1/2 * arctan(2 * dd * A_qp/(A_qq - A_pp)))
#     2*acos(sqrt(cmin)) = arctan(2 * dd * A_qp/(A_qq - A_pp))
#     tan(2*acos(sqrt(cmin))) = 2 * dd * A_qp/(A_qq - A_pp))
#     tan(2*acos(sqrt(cmin)))(A_qq - A_pp)/(2 A_qp) = dd

for dd in ll:
    #dd = -dd
    A    = numpy.array([[1.0, dd],[dd, 2.0]])
    e,v  = numpy.linalg.eig(A)
    dc2  = v[0][0]*v[0][0]
    c2   = math.cos(0.5*math.atan(2.0*dd))
    c2   = c2*c2
    #print("%f, %f, %f, %f" % (dd,math.log(dd),math.atan(math.log(dd))/pi12,(ang*ang)/pi18-1.0))
    #print("%f, %f, %f, %f" % (dd,math.log(dd),math.atan(math.log(dd))/pi12,(ang*ang)/(2.0*pi18*pi18)-1.0))
    print("%f, %f, %f" % (dd,c2,dc2))
    #print(dd)
    #print(e)
    #print(v)
    #print("")
