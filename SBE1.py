#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SBE1.py
#  
#  Copyright 2014 ainhi <ainhi@nhi>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


import math
import cmath
import numpy



def main():
	
	return 0


exp=math.exp
ln=math.log
sqrt=math.sqrt
pi=math.pi


T2 = 200 # Thoi gian khu pha
X0 = 0.5 # Chi0
Er = 4.2 # Nang luong lien ket
detuning = 100.0 # Nang luong troi cua photon
emax = 300.0 # emax
h = 658.5 # Hang so Planck
aB=125.0e-8

delta_t = 30.0 # delta(t) trong cong thuc tinh dE0
N = 80 # So doan chia cua nang luong
de = emax/N #delta (epsilon) 
dt = 2.0 # Buoc thoi gian
t0 = -3*delta_t # Thoi gian ban dau
tf = 17*delta_t # Thoi gian cuoi
M = int((tf-t0)/dt) # So buoc lap thoi gian


g_pola=de*sqrt(de)/(Er*sqrt(Er))/pi/pi
g_density=0.5*(de*sqrt(de))/(Er*sqrt(Er))/((aB*aB*aB)*(pi*pi))/10e+16
z0 = 0 + 0*1j
# Xay dung ham g(n,l), xay dung mang (0,0) cho g(n,l)

g = numpy.zeros ((N+1, N+1)) # N+1 hang, N+1 cot
for n in range(1,N+1):
	for l in range(1,N+1):
		if n==l:
			g[n][l]=0.0
		else:
			g[n][l]=ln((sqrt(n) + sqrt(l))/abs(sqrt(n)-sqrt(l)))*sqrt(Er*de/n)/pi
			
	
f = numpy.zeros((N+1),float)
p = numpy.zeros((N+1),complex)
l = numpy.zeros((N+1,N+1),complex)
for n in range(1,N+1):
	l[n][0] = l[n][1] = z0
def sump(n): #sump
	sump = 0.0 + 0*1j
	for l in range(1,N+1):
		sump = sump + g[n][l]*p[l]
	return sump	
	
def En(n): #sumf
	En = 0.0
	for l in range(1,N+1):
		En = En + 2.0*g[n][l]*f[l]
	return En
	
	

def ft(t,f,p,n):
	ORn = (1/h)*(0.5*h*sqrt(pi)*X0*(exp(-t**2/(delta_t**2)))/delta_t+sump(n))
	w = ORn*p.conjugate()
	ft=-2*(w.imag)
	return ft
	
		

def pt(t,f,p,n):
	ORn= (1/h)*(0.5*h*sqrt(pi)*X0*(exp(-t**2/(delta_t**2)))/delta_t+sump(n))
	pt = -1j*(n*de-detuning-En(n))*p/h+1j*(1.0-2.0*f)*ORn-p/T2
	return pt
	
k1 = k2 = k3 = k4 = 0.0 +0.0*1j
l1 = l2 = l3 =l4 = 0.0 + 0.0*1j

z = open('SBE100sum', 'w')
k = open('SBE100', 'w')

for n in range(1,N+1):
	f[n]=0.0 
	p[n]=0+0*1j
	k.write(str(t0))
	k.write('			')
	k.write(str(n*de))
	k.write('			')
	k.write(str(f[n]))
	k.write('			')
	k.write('			')
	k.write(str(abs(p[n])))
	k.write('\n')
k.write('\n')
sumpp = z0
sumf = 0.0
Sumf = 0.0	
z.write(str(t0))
z.write('			')
z.write(str(abs(sumpp)))
z.write('			')
z.write(str(Sumf))
z.write('\n')

for i in range(0,M+1):
	t = t0 + i*dt
	sumpp += sqrt(n)*p[n]*g_pola
	sumf += sqrt(de) *(f[n])
	Sumf += g_density*sumf
	for n in range(1,N+1):
		k1 = dt*ft(t,f[n],p[n],n)
		l1 = dt*pt(t,f[n],p[n],n)
		k2 = dt*ft(t+0.5*dt,f[n]+0.5*k1,p[n]+0.5*l1,n)
		l2 = dt*pt(t+0.5*dt,f[n]+0.5*k1,p[n]+0.5*l1,n)
		k3 = dt*ft(t+0.5*dt,f[n]+0.5*k2,p[n]+0.5*l2,n)
		l3 = dt*pt(t+0.5*dt,f[n]+0.5*k2,p[n]+0.5*l2,n)
		k4 = dt*ft(t+dt,f[n]+k3,p[n]+l3,n)
		l4 = dt*pt(t+dt,f[n]+k3,p[n]+l3,n)
		
		f[n] = f[n] + (k1 + 2*k2 + 2*k3 +k4 )/6
		p[n] = p[n] + (l1 + 2*l2 + 2*l3 + l4)/6
		k.write(str(t+dt))
		k.write('			')
		k.write(str(n*de))
		k.write('			')
		k.write(str(f[n]))
		k.write('			')
		k.write(str(abs(p[n])))
		k.write(str('\n'))
	k.write('\n')
		
	z.write(str(t+dt))
	z.write('			')
	z.write(str(abs(sumpp)))
	z.write('			')
	z.write(str(Sumf))
	z.write('\n')
	
				
		
k.close()
z.close()	


	
		



if __name__ == '__main__':
	main()

