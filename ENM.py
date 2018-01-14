#ENM on 4AKE


import sys
import numpy as np
#import subprocess
#import u3b
import math
from numpy import linalg as LA


class PDB:
    """Class for read pdb file"""
    def __init__(self,**kwargs):
        self.name = "pdb"
        self.pdblines = []
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            print "no pdb name given, use pdb directly"
        if 'pdbfile' in kwargs:
            self.pdbfile = kwargs['pdbfile']
        else:
            print "No pdbfile given"

    def Read_CA(self):
        self.pdblines = []
        pdb_in = open(self.pdbfile,"r")
        lines = pdb_in.readlines()
        for line in lines:
            if line[0:6].strip() == "ATOM":
                if line[12:16].strip() == "CA":
                    self.pdblines.append(line)
        pdb_in.close()



class ATOM:
    """Class to define single atoms"""
    def __init__(self,line):
        self.firsthalf = line[:30]
        self.lasthalf = line[54:]
        #Record first and last half for outputing
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.resnumber = int(line[22:26])

class ENM:
    """Class for performing ENM"""
    def __init__(self,atoms,**kwargs):
        self.atoms = atoms
        self.name = "ENM"
        self.cutoff = 10.0
        self.n = len(self.atoms)
        self.Kirchhoff = np.zeros((self.n,self.n))
        self.Hessian = np.zeros((self.n * 3,self.n * 3))
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            print "no enm name given, use ENM directly"

        if 'cutoff' in kwargs:
            self.cutoff = float(kwargs['cutoff'])
        else:
            print "no cutoff distance given, use 10"

    def distance_ij(self,i,j):
        return ((self.atoms[i].x-self.atoms[j].x)**2 + (self.atoms[i].y-self.atoms[j].y)**2 + (self.atoms[i].z-self.atoms[j].z)**2)**0.5
    def GenerateKirchhoff(self):
        """Krichhoff matrix with only diagonal components"""
        self.Kirchhoff = np.zeros((self.n, self.n))
        for i in range(self.n ):
            for j in range(i,self.n):
                if i == j:
                    continue
                if self.distance_ij(i,j) < self.cutoff:
                    if abs(self.atoms[i].resnumber - self.atoms[j].resnumber) == 1:
                        self.Kirchhoff[i][j] = -10.0
                        self.Kirchhoff[j][i] = -10.0
                    else:
                        self.Kirchhoff[i][j] = -1.0
                        self.Kirchhoff[j][i] = -1.0

    def GenerateHessian(self):
        """Generate Hessian matrix based on Kirchhoff matrix"""
        #Firstly off diagonal ones
        self.Hessian = np.zeros((self.n * 3, self.n * 3))
        for i in range(self.n ):
            for j in range(i, self.n):
                if self.Kirchhoff[i][j] == 0.0:
                    continue
                else:
                    coeff = self.Kirchhoff[i][j]/self.distance_ij(i,j)**2
                    self.Hessian[i * 3][j * 3] = coeff * (self.atoms[i].x - self.atoms[j].x) ** 2
                    self.Hessian[j * 3][i * 3] = self.Hessian[i * 3][j * 3]

                    self.Hessian[i * 3 + 1][j * 3 + 1] = coeff * (self.atoms[i].y - self.atoms[j].y) ** 2
                    self.Hessian[j * 3 + 1][i * 3 + 1] = self.Hessian[i * 3 + 1][j * 3 + 1]

                    self.Hessian[i * 3 + 2][j * 3 + 2] = coeff * (self.atoms[i].z - self.atoms[j].z) ** 2
                    self.Hessian[j * 3 + 2][i * 3 + 2] = self.Hessian[i * 3 + 2][j * 3 + 2]

                    self.Hessian[i * 3][j * 3 + 1] = coeff * (self.atoms[i].x - self.atoms[j].x) * (self.atoms[i].y - self.atoms[j].y)
                    self.Hessian[j * 3 + 1][i * 3] = self.Hessian[i * 3][j * 3 + 1]
                    self.Hessian[i * 3 + 1][j * 3] = self.Hessian[i * 3][j * 3 + 1]
                    self.Hessian[j * 3][i * 3 + 1] = self.Hessian[i * 3][j * 3 + 1]

                    self.Hessian[i * 3][j * 3 + 2] = coeff * (self.atoms[i].x - self.atoms[j].x) * (self.atoms[i].z - self.atoms[j].z)
                    self.Hessian[j * 3 + 2][i * 3] = self.Hessian[i * 3][j * 3 + 2]
                    self.Hessian[i * 3 + 2][j * 3] = self.Hessian[i * 3][j * 3 + 2]
                    self.Hessian[j * 3][i * 3 + 2] = self.Hessian[i * 3][j * 3 + 2]

                    self.Hessian[i * 3 + 1][j * 3 + 2] = coeff * (self.atoms[i].y - self.atoms[j].y) * (self.atoms[i].z - self.atoms[j].z)
                    self.Hessian[j * 3 + 2][i * 3 + 1] = self.Hessian[i * 3 + 1][j * 3 + 2]
                    self.Hessian[i * 3 + 2][j * 3 + 1] = self.Hessian[i * 3 + 1][j * 3 + 2]
                    self.Hessian[j * 3 + 1][i * 3 + 2] = self.Hessian[i * 3 + 1][j * 3 + 2]
        #Diagonal ones
        for i in range(self.n):
            temp11 = 0.0
            temp12 = 0.0
            temp13 = 0.0
            temp22 = 0.0
            temp23 = 0.0
            temp33 = 0.0
            for j in range(self.n):
                temp11 -= self.Hessian[i * 3][j * 3]
                temp12 -= self.Hessian[i * 3][j * 3 + 1]
                temp13 -= self.Hessian[i * 3][j * 3 + 2]
                temp22 -= self.Hessian[i * 3 + 1][j * 3 + 1]
                temp23 -= self.Hessian[i * 3 + 1][j * 3 + 2]
                temp33 -= self.Hessian[i * 3 + 2][j * 3 + 2]
            self.Hessian[i * 3][i * 3] = temp11
            self.Hessian[i * 3][i * 3 + 1] = temp12
            self.Hessian[i * 3][i * 3 + 2] = temp13
            self.Hessian[i * 3 + 1][i * 3] = temp12
            self.Hessian[i * 3 + 1][i * 3 + 1] = temp22
            self.Hessian[i * 3 + 1][i * 3 + 2] = temp23
            self.Hessian[i * 3 + 2][i * 3] = temp13
            self.Hessian[i * 3 + 2][i * 3 + 1] = temp23
            self.Hessian[i * 3 + 2][i * 3 + 2] = temp33
    def printK(self):
        """Test for Hessian printing out, successful"""
        np.savetxt("K",self.Kirchhoff)

    def printHessian(self):
        """Test for Hessian printing out, successful"""
        np.savetxt("Hessian",self.Hessian)
#Test for Hessian matrix
ake4 = PDB(name = "4ake", pdbfile="4AKE_A.pdb")
ake4.Read_CA()
atoms4 = []
for line in ake4.pdblines:
    atoms4.append(ATOM(line))
doENM = ENM(atoms4,name="4akeENM")
doENM.GenerateKirchhoff()
doENM.GenerateHessian()

eigenvalues,eigenvectors = LA.eigh(doENM.Hessian)

#TEST
##############################
#doENM.printHessian()
#########################
#doENM.printK()
#########################
#a,b = LA.eigh(doENM.Hessian)
#print a
#print b[7]
#####TEST for first 6 eigenvalues, all 0 indicating code is fine.

#Now calculate overlap
ake1 = PDB(name = "1ake", pdbfile="1AKE_A.pdb")
ake1.Read_CA()
atoms1 = []
for line in ake1.pdblines:
    atoms1.append(ATOM(line))

if len(atoms1) != len(atoms4):
    print "atom number not the same"
    exit
diff_vector = []
#vector for difference between 1ake and 4ake, use 1ake-4ake
for i in range(len(atoms4)):
    if atoms1[i].resnumber != atoms4[i].resnumber:
        print "resnumbers do not match"
        exit
    diff_vector.append((atoms1[i].x - atoms4[i].x))
    diff_vector.append((atoms1[i].y - atoms4[i].y))
    diff_vector.append((atoms1[i].z - atoms4[i].z))
#Now do overlap
#check vector normalize with both SK and np, return with the same result
diff_vec_norm = diff_vector/np.linalg.norm(diff_vector)

plot_out = open("data.dat","w")
eigenvectors = np.transpose(eigenvectors)
for i in range(len(eigenvectors)):
    s = str(i) + "  " + str(np.dot(diff_vec_norm,eigenvectors[i])) + "\n"
    plot_out.write(s)
plot_out.close()

###########################TEST FOR FORTRAN RESULTS
