# SSP 2021 OD
# Aaron Zhang
# 7/17/21

from math import sin, cos, sqrt, radians, degrees
import numpy as np
from ScalEqLagrange import roots, rho2
from fandg import *
from babyod import *

k = 0.0172020989484 #Gaussian gravitational constant
cAU = 173.144643267 #speed of light in au/(mean solar)day
eps = radians(23.4374) #Earth's obliquity
tepoch = 2459419.7916667

#relative path of input text file
with open("MoG /ZhangInput.txt" , "r") as f:
    lines = f.readlines()
    inputs = [[[],[],[],[],[],[]],
         [[],[],[],[],[],[]],
         [[],[],[],[],[],[]]]
    #user input if more than 3 observations are available
    newLines = []
    if len(lines) > 3:
       newLines.append(lines[int(input("Choose first data set: ")) - 1])
       newLines.append(lines[int(input("Choose second data set: ")) - 1])
       newLines.append(lines[int(input("Choose third data set: ")) - 1])
    else:
        newLines.append(lines[0])
        newLines.append(lines[1])
        newLines.append(lines[2])
    for i in range(len(newLines)):
        #nine columns after splitting
        nineValues = newLines[i].split(" ")
        times = nineValues[3].split(":")
        mins = float(times[1])/60
        sec = float(times[2])/3600
        #n = final time
        n = (float(times[0]) + mins + sec)
        Y = float(nineValues[0])
        M = float(nineValues[1])
        D = float(nineValues[2])
        J0 = (367 * Y) - int(7 * (Y+ int((M+9) / 12)) / 4) + int(275 * M / 9 ) + D + 1721013.5
        JD = J0 + n / 24
        inputs[i][0].append(JD)
        ra = nineValues[4].split(":")
        #conversion to decimal degrees then radians
        raValues = float(ra[0]) + float(ra[1])/60 + float(ra[2])/3600 
        raValuesDeg = raValues * 15
        raRad = radians(raValuesDeg)
        inputs[i][1].append(raRad)
        dec = nineValues[5].split(":")
        if float(dec[0]) < 0:
            decValues = float(dec[0]) - float(dec[1])/60 - float(dec[2])/3600
        else:
            decValues = float(dec[0]) + float(dec[1])/60 + float(dec[2])/3600
        decRad = radians(decValues)
        inputs[i][2].append(decRad)
        sunx = nineValues[6]
        suny = nineValues[7]
        sunz = nineValues[8]
        inputs[i][3].append(sunx)
        inputs[i][4].append(suny)
        inputs[i][5].append(sunz)
t0 = inputs[0][0][0]
t1 = inputs[1][0][0]
t2 = inputs[2][0][0]
p = [[],[],[]]
R = [[],[],[]]
for i in range(len(newLines)):
    RA = inputs[i][1]
    dec = inputs[i][2]
    #Eq 1
    p[i] = [cos(inputs[i][1][0]) * cos(inputs[i][2][0]), sin(inputs[i][1][0]) * cos(inputs[i][2][0]), sin(inputs[i][2][0])]
    R[i] = [float(inputs[i][3][0]), float(inputs[i][4][0]), float(inputs[i][5][0])]

#Mog 1 determing the roots and rho2
def mog1():
    #Step 1
    D0 = np.dot(p[0], np.cross(p[1], p[2]))
    D11 = np.dot(np.cross(R[0], p[1]), p[2])
    D12 = np.dot(np.cross(R[1], p[1]), p[2])
    D13 = np.dot(np.cross(R[2], p[1]), p[2])
    D21 = np.dot(np.cross(p[0], R[0]), p[2])
    D22 = np.dot(np.cross(p[0], R[1]), p[2])
    D23 = np.dot(np.cross(p[0], R[2]), p[2])
    D31 = np.dot(p[0], np.cross(p[1], R[0]))
    D32 = np.dot(p[0], np.cross(p[1], R[1]))
    D33 = np.dot(p[0], np.cross(p[1], R[2]))
    tau3 = k * (t2 - t1)
    tau1 = k * (t0 - t1)
    tau = tau3 - tau1
    A1 = tau3 / tau
    B1 = A1 / 6 * (tau**2 - tau3**2)
    A3 = -tau1 / tau
    B3 = A3 / 6 * (tau**2 - tau1**2)
    A = (A1 * D21 - D22 + A3 * D23) / -D0
    B = (B1 * D21 + B3 * D23) / -D0
    E = -2 * np.dot(p[1], R[1])
    F = np.dot(R[1], R[1])
    a = -1 * (A**2 + (A * E) + F)
    b = -1 * (2 * (A * B) + (B * E))
    c = -1 * (B**2)
    r2 = np.roots([1,0,a,0,0,b,0,0,c])
    return D0,D11,D12,D13,D21,D22,D23,D31,D32,D33,tau3,tau1,tau,A1,B1,A3,B3,A,B,E,F,a,b,c,r2
D0,D11,D12,D13,D21,D22,D23,D31,D32,D33,tau3,tau1,tau,A1,B1,A3,B3,A,B,E,F,a,b,c,r2 = mog1()
roots = []
rho2 = []
for i in range(len(r2)):
    if r2[i] > 0:
        roots.append(r2[i].real)
        rho2.append(A + (B / (r2[i].real **3)))
rTwo = roots[0]
if len(roots) > 0:
    print("Roots:", roots)
    print("Rho2:", rho2)
    print("Note: only index 0 works for the input values given in the OC input file")
    rTwo = roots[int(input("Type Index of Root: "))]

#Step 2
def fandg():
    f = 1 - mu / (2 * rTwo ** 3) * tau ** 2
    g = tau - mu / (6 * rTwo ** 3) * tau ** 3
    return f, g
#determining r2 and r2dot
def mog1Cont():
    f1, g1 = fg2(tau1, rTwo)
    f3, g3 = fg2(tau3, rTwo)
    #Step 3 
    c1 = g3 / (f1 * g3 - g1 * f3)
    c2 = -1
    c3 = -g1 / (f1 * g3 - g1 * f3)
    #Step 4
    p1 = (c1 * D11 + c2 * D12 + c3 * D13) / (c1 * D0)
    p2 = (c1 * D21 + c2 * D22 + c3 * D23) / (c2 * D0)
    p3 = (c1 * D31 + c2 * D32 + c3 * D33) / (c3 * D0)
    #Step 5
    r1v = p1 * np.array(p[0]) - R[0]
    r2 = p2 * np.array(p[1]) - R[1]
    r3v = p3 * np.array(p[2]) - R[2]
    
    d1 = -f3 / (f1 * g3 - f3 * g1)
    d3 = f1 / (f1 * g3 - f3 * g1)
    #Step 6
    r2dot = d1 * r1v + d3 * r3v
    return r2, p1, p2, p3, r2dot
r2, p1, p2, p3, r2dot = mog1Cont()
or2dot = rho2[0]

#Mog 2: new r2 and r2dot values
def newValues(r2, r2dot, p1, p2, p3):
    nt1 = t0 - p1 / cAU
    nt2 = t1 - p2 / cAU
    nt3 = t2 - p3 / cAU
    ntau3 = k * (nt3 - nt2)
    ntau1 = k * (nt1 - nt2)
    ntau = ntau3 - ntau1

    nf1, ng1 = fg(ntau1, r2, r2dot, "funct")
    nf3, ng3 = fg(ntau3, r2, r2dot, "funct")

    nc1 = ng3 / (nf1 * ng3 - ng1 * nf3)
    nc2 = -1
    nc3 = -ng1 / (nf1 * ng3 - ng1 * nf3)

    np1 = (nc1 * D11 + nc2 * D12 + nc3 * D13) / (nc1 * D0)
    np2 = (nc1 * D21 + nc2 * D22 + nc3 * D23) / (nc2 * D0)
    np3 = (nc1 * D31 + nc2 * D32 + nc3 * D33) / (nc3 * D0)

    nr1v = np1 * np.array(p[0]) - R[0]
    nr2v = np2 * np.array(p[1]) - R[1]
    nr3v = np3 * np.array(p[2]) - R[2]

    nd1 = -nf3 / (nf1 * ng3 - nf3 * ng1)
    nd3 = nf1 / (nf1 * ng3 - nf3 * ng1)

    nr2dot = nd1 * nr1v + nd3 * nr3v
    return nr2v, nr2dot, np1, np2, np3
nr2, nr2dot, np1, np2, np3 = newValues(r2, r2dot, p1, p2, p3)

#main iteration loop
def Mog2(nr2, nr2dot, np1, np2, np3, p2):
    while abs(np2 - p2) > 1E-12:
        p2 = np2
        nr2, nr2dot, np1, np2, np3 = newValues(nr2,nr2dot, np1, np2, np3)
    return nr2, nr2dot, np1, np2, np3
nr2, nr2dot, np1, np2, np3 = Mog2(nr2, nr2dot, np1, np2, np3, p2)

#ecliptic to equatorial converison
def ectoeq():
    obliquity = radians(23.4374)
    matrix1 = np.array([[1,0,0],[0,cos(obliquity),sin(obliquity)],[0,-sin(obliquity),cos(obliquity)]])
    eqr2 = matrix1 @ nr2
    eqr2dot = matrix1 @ nr2dot
    return eqr2, eqr2dot
eqr2,eqr2dot = ectoeq()
a, n, e, hv, h, i, omega, fPlusW, f, omini, M0, MCenter= orbitalElements(eqr2, eqr2dot, t1, tepoch)

# print("ecr2", ecr2)
# print("ecr2dot", ecr2dot)
print("a = ", a)
print("e= ", e)
print("i= ", i)
print("omini= ", omini)
print("Omega= ", omega)
print("M0= ", M0)
print("MCenter= ", MCenter)
    