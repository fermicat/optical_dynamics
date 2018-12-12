import math
hbar = 1.0546e-34
epsilon0 = 8.8542e-12
echarge = 1.602176487e-19
m0 = 9.10938215e-31
kB = 1.3806504e-23
c0 = 3.0e8
eGauss = echarge/math.sqrt(4.0 * math.pi * epsilon0)

L = 3e-3 # 3 for 70 ps, 1 for 23 ps round-trip time
rl = 1
rr = 0.6

T1 = 1e-10
T2 = 50e-15
Diff = 0.0046
overlap = 0.1
dipole = 0.8e-9 * echarge
nr = 3.5
loss = 500
lambda0 = 3.2e-6

#rho_e = 2*(0.0309*m0)/(2*math.pi * hbar**2) * (hbar / T2)
#rho_e = 0.6e22

Rabbi2 = T2 * dipole**2 / hbar**2
print("Rabbi2 =", Rabbi2)

Etransition = hbar * c0 * 2.0 * math.pi / lambda0
kwave = (Etransition/hbar) / c0
totloss = loss - 1/(2*L) * math.log(rl*rr)
print ("totloss =", totloss)

gain_coef_worho = overlap * dipole**2 * (Etransition/hbar) * T2 / (2 * hbar * epsilon0 * nr * c0)
print("gain_coef_w/orho =", gain_coef_worho)

J = 0.8
DeltaN = 2 * J - 1
rho = totloss/(gain_coef_worho*DeltaN)

print("rho = ", rho)

print(5.62e22 * gain_coef_worho * DeltaN * 100 / totloss, "%")
# 10%


G = overlap * dipole * (Etransition/hbar) *rho /(nr * epsilon0 *c0)
print("G =", G)