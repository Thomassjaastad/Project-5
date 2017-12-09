from __future__ import division
from numpy import *
from matplotlib.pyplot import *

file_d = open("Diffusion_temperature.txt", "r")
#file = open("statistics.txt", "r")

Temp = []
E_p = []
E_k = []
E_tot = []
t = []
Dif = []

T_dif = []
D_final = []
D_ave = []
for line in file_d:
	numbers = line.split()
	T_dif.append(float(numbers[5]))
	D_ave.append(float(numbers[3]))
	D_final.append(float(numbers[1]))

"""
for line in file:
	numbers = line.split()
	E_k.append(float(numbers[0]))
	E_p.append(float(numbers[1]))
	E_tot.append(float(numbers[2]))
	Temp.append(float(numbers[3]))
	t.append(float(numbers[4]))
	Dif.append(float(numbers[5]))

file.close()
"""


#Converting to SI unit
"""
TempSI = array(Temp)*119.7
E_pSI = array(E_p)*1.651E-21
E_kSI = array(E_k)*1.651E-21
E_totSI = array(E_tot)*1.651E-21
tSI = array(t)*1.002E-13
DifSI = array(Dif)*1E-13
"""
  #For intial temp 200, 500, 700 ,1000 K

"""
subplot(2, 2, 1)
xlabel('$t[s]$', fontsize = 15)
ylabel('$E_k[J]$', fontsize = 16)
title('$Kinetic$ $Energy$ $SI$ $units$', fontsize = 14)
plot(tSI, E_kSI)


subplot(2,2,2)
xlabel('$t[s]$', fontsize = 15)
ylabel('$E_p[J]$', fontsize = 16)
title('$Potential$ $Energy$ $SI$ $units$', fontsize = 14)
plot(tSI, E_pSI, 'g')

subplot(2,2,3)
xlabel('$t[s]$', fontsize = 15)
ylabel('$E_{tot}[J]$', fontsize = 16)
title('$Total$ $Energy$ $SI$ $units$', fontsize = 14)
plot(tSI, E_totSI, 'r')

subplot(2,2,4)
xlabel('$t[s]$', fontsize = 15)
ylabel('$Temperature[K]$', fontsize = 16)
title('$Temperature$ $SI$ $units$', fontsize = 14)
plot(tSI, TempSI, 'y')
"""

plot(T_dif, D_final, 'o')
xlabel('$T[K]$', fontsize = 18)
ylabel('$D$', fontsize= 18)
show()

plot(T_dif,D_ave, 'o')
xlabel('$T[K]$', fontsize = 18)
ylabel('$D[m^2/s]$', fontsize= 18)
show()


"""
plot(tSI, Ratio_1000, label = 'T = 1000K')
hold('on')
plot(tSI, Ratio_700, label = 'T = 700K')
hold('on')
plot(tSI, Ratio_500, label = 'T = 500K')
hold('on')
plot(tSI, Ratio_200, label = 'T = 200K')

xlabel('$t[s]$', fontsize = 16)
ylabel('$T/T_i$', fontsize = 18)
title('$Temperature$ $Ratios$', fontsize = 20)
legend()
show()
"""