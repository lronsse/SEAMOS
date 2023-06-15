import scipy
import scipy.integrate as spy
from scipy.integrate import quad

b = 3.04 # [m] wing span
RatioChords = 0.32 #aileron chord/wing chord
ailchord = 0.32*0.2688 #chord aileron
AilSpan = 0.32*3.04
AilSpan = round(AilSpan)
print('Aileron span: ',AilSpan,'m')
print('Aileron chord: ',ailchord)

c = 0.2688 #[m] source: Mathis
S_ref = 0.65# [m^2] wing area
tau = 0.4 #effectiveness

#initial b1 & b2: 0.868 & 1.368
b1 = 0.868
b2 = 1.3 #end of aileron (10% of the wing tip is free of aileron)
print('Ailerons are placed at a distance b1 =', b1,'m from the symmetrical axis and extend to', b2,'m')
# y = b1
d_alpha_upper = 20
d_alpha_lower = 15
d_alpha = 1/2*(d_alpha_upper+d_alpha_lower)

cl_alpha = 0.0984
cd_0 = 0.018

def f(y):
    return c*y

# Compute the integral
result_int1, error = spy.quad(f, b1, b2)
#print(result_int1)

#find the aileron control derivative
Cl_d_alpha = 2*cl_alpha*tau/(S_ref*b)*result_int1

#find the roll damping Coefficient CL_P
def f2(y):
    return c*y**2

result_int2, error = spy.quad(f2, 0, b/2)
Cl_P = -(4*(cl_alpha+cd_0)/(S_ref*b))*result_int2

#Find the roll rate dphi/dt
Pi = scipy.pi
Phi = 35 #
V = 20 #[m/s]
P = - Cl_d_alpha/Cl_P*d_alpha *(2*V/b)

t = (1/P)*Phi #[s]

print('The roll rate is',P,'degrees per second')
print('The time to reach a bank angle of',Phi,'degrees, is', t,'s')


