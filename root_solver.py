from mpmath import *
from numpy import complex
import numpy as np
from numpy import *

def complex_root_solver(x,para,step,converge,D):
	#x        : the root you want to solve
	#para     : other parameter u will use in the function D
	#step     : how big the step u will step
	#converge : converge condition when abs(D)<converge
	done = 0
	while(abs(D(x,para))>converge and done<4000):
		J1 = (D(x+1e-6,para)-D(x,para))/1e-6
		J2 = (D(x+1e-6j,para)-D(x,para))/1e-6
		J11 = J1.real
		J12 = J2.real
		J21 = J1.imag
		J22 = J2.imag
		J = [J11,J12],[J21,J22]
		Jacobian = mat(J)
		inver_J = Jacobian.I
		aa = [D(x,para).real,D(x,para).imag]
		a1 = mat(aa)
		a2 = a1.T
		bb = [x.real,x.imag]
		b1 = mat(bb)
		b2 = b1.T
		b2 = b2 - inver_J*a2*step
		x = complex(b2[0,0],b2[1,0])
		done = done + 1
	if( abs(D(x,para))>converge or D(x,para).real>1e20 or D(x,para).imag>1e20 or math.isnan(D(x,para).real)or math.isnan(D(x,para).imag)):
		return complex(0,0)
	else:
		print('done !, x is ',x,' and parameter is ',para,',and done is',done)
		return x


