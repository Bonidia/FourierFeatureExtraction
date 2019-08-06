#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import numpy as np
import scipy.stats
import statistics
import argparse
from Bio import SeqIO
from scipy.fftpack import fft, ifft
import warnings
warnings.filterwarnings("ignore")


#############################################################################
#############################################################################
def header():
	dataset = open(foutput, 'a')
	dataset.write("average,median,maximum,minimum,peak,nonelevatedpeak,standarddeviation,standarddeviationpop,percentile15,percentile25,percentile50,percentile75,amplitude,variance,class")
	dataset.write("\n")
	return

        
def fileRecord():
	dataset = open(foutput, 'a')
	for metric in features:
		dataset.write("%s," % (metric))
		# dataset.write("{0:.4f},".format(metric))
	dataset.write(labelDataset)
	dataset.write("\n")
	print ("Sequence Analyzed!!")
	return


def featureExtraction():
	global features
	features = []
	average = sum(spectrum)/len(spectrum)
	features.append(average)
	###################################
	median = np.median(spectrum)
	features.append(median)
	###################################
	maximum = np.max(spectrum)
	features.append(maximum)
	###################################
	minimum = np.min(spectrum)
	features.append(minimum)
	###################################
	peak = (len(spectrum)/3)/(average)
	features.append(peak)
	###################################
	peakTwo = (len(spectrumTwo)/3)/(np.mean(spectrumTwo))
	features.append(peakTwo)
	###################################
	standarddeviation = np.std(spectrum) # standard deviation
	features.append(standarddeviation)
	###################################
	standarddeviationpop = statistics.stdev(spectrum) # population sample standard deviation 
	features.append(standarddeviationpop)
	###################################
	percentile15 = np.percentile(spectrum, 15)
	features.append(percentile15)
	###################################
	percentile25 = np.percentile(spectrum, 25)
	features.append(percentile25)
	###################################
	percentile50 = np.percentile(spectrum, 50)
	features.append(percentile50)
	###################################
	percentile75 = np.percentile(spectrum, 75)
	features.append(percentile75)
	###################################
	amplitude = maximum - minimum
	features.append(amplitude)
	###################################
	# mode = statistics.mode(spectrum)
	###################################
	variance = statistics.variance(spectrum)
	features.append(variance)
	return
    

def binaryFourier():
    header()
    global spectrum, spectrumTwo
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        spectrum = []
        spectrumTwo = []
        A = []
        C = []
        T = []
        G = []
        for nucle in seq:
            if nucle == "A":
                A.append(1)
            else:
                A.append(0)
            if nucle == "C":
                C.append(1)
            else:
                C.append(0)
            if nucle == "T":
                T.append(1)
            else:
                T.append(0)
            if nucle == "G":
                G.append(1)
            else:
                G.append(0)
        FA = fft(A)
        FC = fft(C)
        FT = fft(T)
        FG = fft(G)
        for i in range(len(seq)):
        	specTotal = (abs(FA[i])**2) + (abs(FC[i])**2) + (abs(FT[i])**2) + (abs(FG[i])**2)
        	specTwo = (abs(FA[i])) + (abs(FC[i])) + (abs(FT[i])) + (abs(FG[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        featureExtraction()
        fileRecord()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print (seq)
        # print ("\n")
        # print ("A -- %s" % (A))
        # print ("C -- %s" % (C))
        # print ("T -- %s" % (T))
        # print ("G -- %s" % (G))
        # print ("\n")
        # print ("A -- %s" % (abs(FA)))
        # print ("C -- %s" % (abs(FC)))
        # print ("T -- %s" % (abs(FT)))
        # print ("G -- %s" % (abs(FG)))
        # print ("\n")
        # print (spectrumTwo)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def zcurveFourier():
    header()
    global spectrum, spectrumTwo
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        spectrum = []
        spectrumTwo = []
        ###################################
        ###################################
        R = 0 # x[n] = (An + Gn) − (Cn + Tn) ≡ Rn − Yn
        Y = 0
        M = 0 # y[n] = (An + Cn) − (Gn + Tn) ≡ Mn − Kn
        K = 0
        W = 0 # z[n] = (An + Tn) − (Cn + Gn) ≡ Wn − Sn
        S = 0
        ###################################
        ###################################
        x = []
        y = []
        z = []
        for nucle in seq:
            if nucle == "A" or nucle == "G":
            	R += 1
            	x.append((R)-(Y))
            else:
            	Y += 1
            	x.append((R)-(Y))
            if nucle == "A" or nucle == "C":
                M += 1
                y.append((M)-(K))
            else:
                K += 1
                y.append((M)-(K))
            if nucle == "A" or nucle == "T":
                W += 1
                z.append((W)-(S))
            else:
                S += 1
                z.append((W)-(S))
        FX = fft(x)
        FY = fft(y)
        FZ = fft(z)
        for i in range(len(seq)):
        	specTotal = (abs(FX[i])**2) + (abs(FY[i])**2) + (abs(FZ[i])**2)
        	specTwo = (abs(FX[i])) + (abs(FY[i])) + (abs(FZ[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        featureExtraction()
        fileRecord()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print (seq)
        # print ("\n")
        # print ("X -- %s" % (x))
        # print ("Y -- %s" % (y))
        # print ("Z -- %s" % (z))
        # print ("\n")
        # print ("X -- %s" % (abs(FX)))
        # print ("Y -- %s" % (abs(FY)))
        # print ("Z -- %s" % (abs(FZ)))
        # print ("\n")
        # print (spectrum)
        # print ("\n")
        # print (spectrumTwo)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def integerFourier():
    header()
    global spectrum, spectrumTwo
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        spectrum = []
        spectrumTwo = []
        integer = []
        for nucle in seq:
            if nucle == "T":
            	integer.append(0)
            elif nucle == "C":
            	integer.append(1)
            elif nucle == "A":
            	integer.append(2)
            else:
            	integer.append(3)
        FI = fft(integer)
        for i in range(len(seq)):
        	specTotal = (abs(FI[i])**2)
        	specTwo = (abs(FI[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        featureExtraction()
        fileRecord()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print (seq)
        # print ("\n")
        # print ("I -- %s" % (integer))
        # print ("\n")
        # print ("I -- %s" % (abs(FI)))
        # print ("\n")
        # print (spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def realFourier():
    header()
    global spectrum, spectrumTwo
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        spectrum = []
        spectrumTwo = []
        real = []
        for nucle in seq:
            if nucle == "T":
            	real.append(1.5)
            elif nucle == "C":
            	real.append(0.5)
            elif nucle == "A":
            	real.append(-1.5)
            else:
            	real.append(-0.5)
        FR = fft(real)
        for i in range(len(seq)):
        	specTotal = (abs(FR[i])**2)
        	specTwo = (abs(FR[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        featureExtraction()
        fileRecord()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print (seq)
        # print ("\n")
        # print ("R -- %s" % (real))
        # print ("\n")
        # print ("R -- %s" % (abs(FR)))
        # print ("\n")
        # print (spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def eiipFourier():
    header()
    global spectrum, spectrumTwo
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        spectrum = []
        spectrumTwo = []
        eiip = []
        for nucle in seq:
            if nucle == "T":
            	eiip.append(0.1335)
            elif nucle == "C":
            	eiip.append(0.1340)
            elif nucle == "A":
            	eiip.append(0.1260)
            else:
            	eiip.append(0.0806)
        Feiip = fft(eiip)
        for i in range(len(seq)):
        	specTotal = (abs(Feiip[i])**2)
        	specTwo = (abs(Feiip[i]))
        	spectrum.append(specTotal)
        	spectrumTwo.append(specTwo)
        featureExtraction()
        fileRecord()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print (seq)
        # print ("\n")
        # print ("R -- %s" % (eiip))
        # print ("\n")
        # print ("R -- %s" % (abs(Feiip)))
        # print ("\n")
        # print (spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


        
#############################################################################    
if __name__ == "__main__":
	print("\n")
	print("###################################################################################")
	print("########## Feature Extraction: A Fourier and Numerical Mapping Approach ###########")
	print("##########     Arguments: python3.5 -i input -o output -l mRNA -r 2     ###########")
	print("#### Representations: 1 = Binary, 2 = Z-curve, 3 = Real, 4 = Integer, 5 = EIIP ####")
	print("##########                 Author: Robson Parmezan Bonidia              ###########")
	print("###################################################################################")
	print("\n")
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', help='Fasta format file | E.g., test.fasta')
	parser.add_argument('-o', '--output', help='Csv format file | E.g., test.csv')
	parser.add_argument('-l', '--label', help='Dataset Label | E.g., lncRNA, mRNA, sncRNA ...')
	parser.add_argument('-r', '--representation', help='1 = Binary, 2 = Z-curve, 3 = Real, 4 = Integer, 5 = EIIP')
	args = parser.parse_args()
	finput = str(args.input)
	foutput = str(args.output)
	labelDataset = str(args.label)
	representation = int(args.representation)
	if representation == 1:
		binaryFourier()
	elif representation == 2:
		zcurveFourier()
	elif representation == 3:
		realFourier()
	elif representation ==4:
		integerFourier()
	elif representation == 5:
		eiipFourier()
	else:
		print("This package does not contain this representation")
#############################################################################
