#Homework 1 for PHY 583
#Han Wen
#2/15/2016


#usage: python hw1.py template_pdb_file the_other_pdb_file
#Need to load module ncbi/blast-2.2.29



import sys
import numpy as np
import subprocess
import u3b
import math



#First part, read two pdbs and outputs the sequence
###################################################

pdb_file_1 = str(sys.argv[1])
pdb_in_1 = open(pdb_file_1, "r")

pdb_file_2 = str(sys.argv[2])
pdb_in_2 = open(pdb_file_2, "r")


seq_1 = open("seq1", "w")
seq_2 = open("seq2", "w")
print "Load ncbi/blast-2.2.29 before use"



#Array to restore the CA coordinates for both pdbs
pdb1_x=[]
pdb1_y=[]
pdb1_z=[]

pdb2_x=[]
pdb2_y=[]
pdb2_z=[]
###############


#Array for seq translation:
full_seq = ['ala','asn','asp','arg','cys','gln','glu','gly','his','ile','leu','lys','met','pro','phe','ser','thr','trp','tyr','val','ALA','ASN','ASP','ARG','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PRO','PHE','SER','THR','TRP','TYR','VAL']
abb_seq ='ANDRCQEGHILKMPFSTWYVANDRCQEGHILKMPFSTWYV'




lines = pdb_in_1.readlines()
for line in lines:
    if line[0:6].strip() == "ATOM":
	if line[12:16].strip() == "CA":
            resname = line[17:20]
	    pdb1_x.append(float(line[30:38]))
            pdb1_y.append(float(line[38:46]))
            pdb1_z.append(float(line[46:54]))
	    if resname not in full_seq:
	        if resname == "UNK":
                    print "containing UNK", line
	            continue
	        else:
	            print "Error in reading resname", line
	    #print resname
	    #print full_seq.index(resname)
	    seq_1.write(abb_seq[full_seq.index(resname)])
seq_1.close()




lines = pdb_in_2.readlines()
for line in lines:
    if line[0:6].strip() == "ATOM":
        if line[12:16].strip() == "CA":
            resname = line[17:20]
            pdb2_x.append(float(line[30:38]))
            pdb2_y.append(float(line[38:46]))
            pdb2_z.append(float(line[46:54]))
            if resname not in full_seq:
                if resname == "UNK":
                    print "containing UNK", line
                    continue
                else:
                    print "Error in reading resname", line

	    seq_2.write(abb_seq[full_seq.index(resname)])
seq_2.close()

##############################################################
#####Extract seq test successfully############################
##############################################################


#Second part, run blast
#Need to load module ncbi/blast-2.2.29
#This is different from bl2seq but the results are the same

subprocess.call("module load ncbi/blast-2.2.29 ",shell=True)

subprocess.call(" blastp -query seq1 -subject seq2 > seq_result ",shell=True)
##############################################################
#####Test successfully, can generate right blast result#######
##############################################################



#Third part, to read blast result and extract CAs coordinates for align and rmsd
#Finally, list of x,y,z for both pdb CAs for rmsd wanted

#List and sequence number for generate the xyz coordinates for align/rmsd

que_coor=[[],[],[]]
sub_coor=[[],[],[]]



que_num = 0
sub_num = 0
w=[]#weight list

w_sum = 0
seq_res = open("seq_result", "r")
line_seq = seq_res.readlines()
for line in line_seq:
    if line[0:6] == "Query ":
	#print line_seq[line_seq.index(line)]
        #print line_seq[line_seq.index(line)+1]
        #print line_seq[line_seq.index(line)+2]
	#test for selection, succeed

	#recognizing needed string
	que_line = line
	mid_line = line_seq[line_seq.index(line)+1]
	sub_line = line_seq[line_seq.index(line)+2]
	words = str.split(sub_line)
	num = int(words[-1]) - int(words[1]) + 1
	str_que = que_line[12:12+num]
	str_mid = mid_line[12:12+num]
	str_sub = sub_line[12:12+num]
	###########################
	#print str_que
        #print str_mid
        #print str_sub

	#Now append the needed coordinates data
	for i in range(num):
	    if str_mid[i] == str_que[i] and str_mid[i] == str_sub[i]:
		que_coor[0].append(pdb1_x[que_num])
                que_coor[1].append(pdb1_y[que_num])
                que_coor[2].append(pdb1_z[que_num])
                sub_coor[0].append(pdb2_x[sub_num])
                sub_coor[1].append(pdb2_y[sub_num])
                sub_coor[2].append(pdb2_z[sub_num])
		w.append(1)
	 	w_sum += 1
		que_num += 1
		sub_num += 1
	    else:
		if str_que[i] != "-":
		    que_num += 1
                if str_sub[i] != "-":
                    sub_num += 1

	#Finishing with getting coordinates list and weight list and w sum
#print que_coor
#print w
#print w_sum
#print sub_coor
#print len(que_coor)
#print len(sub_coor)
####A total of 723 atoms, test successfully##########
seq_res.close()



#############################################

#Fourth part, do the rmsd alignment calculation

mode = 1 


#print w_sum
#print len(w)
#print len(sub_coor[0]) 
#print u3b.u3b.__doc__
u3b_result=u3b.u3b(w,sub_coor,que_coor,mode,w_sum)

#print u3b_result
#print rms
#print math.sqrt(rms[0]/w_sum)

##############################################################
#####Test successfully, generate right rmsd###################
#u3b_result=(31362.167663666653, array([[-0.78473756,  0.03663972,  0.61874429],
#       [-0.38788093,  0.74960034, -0.53632799],
#       [-0.48346183, -0.66087583, -0.57402769]]), array([ 108.08821032,   15.09780369,   62.49975784]), 0)
ier = u3b_result[3]
if ier==-1:
    print "\nError in bestfit_u3b rms SUPERPOSITION IS NOT UNIQUE!\n"
elif ier==-2:
    print "\nError in bestfit_u3b: NO RESULT OBTAINED!\n"
    

rmsd_out = math.sqrt(u3b_result[0]/w_sum)
U = u3b_result[1]
T = u3b_result[2] 
#print U
#print T
#print U[0][1]
print "The rmsd is ",rmsd_out

################################################################
#Final step, do the transformation UX+T and out put the
#transformed protein structure for second pdb file
################################################################

pdb2_out = open("aligned_"+pdb_file_2, "w")



#lines = pdb_in_2.readlines()
for line in lines:
    if line[0:6].strip() == "ATOM":
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
	out_x = x*U[0][0] + y*U[0][1] + z*U[0][2] + T[0]
        out_y = x*U[1][0] + y*U[1][1] + z*U[1][2] + T[1]
        out_z = x*U[2][0] + y*U[2][1] + z*U[2][2] + T[2]
	outxyz='%8.3f%8.3f%8.3f'%(out_x,out_y,out_z)
	s = line[0:30] + outxyz + line[54:] 
	pdb2_out.write(s)
pdb2_out.close()
pdb_in_1.close()
pdb_in_2.close()
#Check the rmsd between aligned and original, 0
#Check in vmd, overlapped, indicating success

