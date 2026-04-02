#!/usr/bin/env python
import sys
##########################################
class ph4vector():
    def __init__(self, x0, y0, z0, x1, y1, z1, feature):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0

        self.x1 = x1
        self.y1 = y1
        self.z1 = z1

        if (feature == 'ring'):
            self.color = 'orange'
        elif (feature == 'Rring'):
            self.color = 'yellow'
        elif (feature == 'donor'):
            self.color = 'blue'
        elif (feature == 'acceptor'):
            self.color = 'red'
 


def main():
    mol_filename = sys.argv[1]
    out_filename = sys.argv[2]

    indices = []
    vectors = []

    csv = open(mol_filename, 'r')

    acceptor_flag  = False
    donor_flag     = False
    ringUp_flag    = False
    ringDown_flag  = False
    RringUp_flag   = False
    RringDown_flag = False

    for line in csv:
        spl = line.split()

        if (ringDown_flag):
            if (not(spl[1] == 'H2')):
                print ("error with input file")
                return
            #x1, y1, z1 = spl[2], spl[3], spl[4][:-2]
            x1, y1, z1 = spl[2], spl[3], spl[4]
            vectors.append( ph4vector(x0, y0, z0, x1, y1, z1, feature) )
            ringDown_flag = False

        if (ringUp_flag):
            if (not(spl[1] == 'H1')):
                print ("error with input file")
                return
            #x1, y1, z1 = spl[2], spl[3], spl[4][:-2]
            x1, y1, z1 = spl[2], spl[3], spl[4]
            vectors.append( ph4vector(x0, y0, z0, x1, y1, z1, feature) )
            ringUp_flag = False
            ringDown_flag = True


        if (RringDown_flag):
            if (not(spl[1] == 'H4')):
                print ("error with input file")
                return
            #x1, y1, z1 = spl[2], spl[3], spl[4][:-2]
            x1, y1, z1 = spl[2], spl[3], spl[4]
            vectors.append( ph4vector(x0, y0, z0, x1, y1, z1, feature) )
            RringDown_flag = False

        if (RringUp_flag):
            if (not(spl[1] == 'H3')):
                print ("error with input file")
                return
            #x1, y1, z1 = spl[2], spl[3], spl[4][:-2]
            x1, y1, z1 = spl[2], spl[3], spl[4]
            vectors.append( ph4vector(x0, y0, z0, x1, y1, z1, feature) )
            RringUp_flag = False
            RringDown_flag = True


        if (acceptor_flag):
            if (not(spl[1] == 'HA')):
                print ("error with input file")
                return
            #x1, y1, z1 = spl[2], spl[3], spl[4][:-2]
            x1, y1, z1 = spl[2], spl[3], spl[4]
            vectors.append( ph4vector(x0, y0, z0, x1, y1, z1, feature) )
            acceptor_flag = False

        if (donor_flag):
            if (not(spl[1] == 'N')):
                print ("error with input file")
                return
            #x1, y1, z1 = spl[2], spl[3], spl[4][:-2]
            x1, y1, z1 = spl[2], spl[3], spl[4]
            vectors.append( ph4vector(x0, y0, z0, x1, y1, z1, feature) )
            donor_flag = False

        if (len(spl) == 9):
            if (spl[1] == 'S'):
                ringUp_flag = True
                x0, y0, z0 = spl[2], spl[3], spl[4]
                feature = 'ring'

            elif (spl[1] == 'P'):
                RringUp_flag = True
                x0, y0, z0 = spl[2], spl[3], spl[4]
                feature = 'Rring'

            elif (spl[1] == 'O'):
                acceptor_flag = True               
                x0, y0, z0 = spl[2], spl[3], spl[4]
                feature = 'acceptor'

            elif (spl[1] == 'HD'):
                donor_flag = True               
                x0, y0, z0 = spl[2], spl[3], spl[4]
                feature = 'donor'

    csv.close()

##############################
# print out vectors here

    output = open(out_filename, 'w')
    for vec in vectors:
        #print ".color %s" %vec.color
        output.write(".color %s\n" %vec.color)
        #print ".arrow%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f 0.01 0.04 0.75"%(float(vec.x0),float(vec.y0),float(vec.z0),float(vec.x1),float(vec.y1),float(vec.z1))
        output.write(".arrow")
        output.write("%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f" %(float(vec.x0),float(vec.y0),float(vec.z0),float(vec.x1),float(vec.y1),float(vec.z1))) 
        output.write(" 0.01 0.04 0.75\n")
    output.close()

main()
