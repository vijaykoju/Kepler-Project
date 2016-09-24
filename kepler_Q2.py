#!/usr/bin/env python

##################################  kepler_04.py  #############################################
# File:              kepler_02.py (for data from Q2)
# Author:            Vijay Koju (vk8736@truman.edu)
# Created:           11 June 2011
# Last modified:     29 July 2011
#
# Please feel free to contact me if you have any questions regarding this program.
#
################################################################################################
# Required modules for running this program.
#        scipy, numpy, matplotlib
# scipy, and numpy can be downloaded from scipy.org
# matplotlib can be downloaded from matplotlib.sourceforge.net/index.html
################################################################################################
# Installation:
# In order to run this program from the command line, a symbolic link to this program has to be
# created in /usr/local/bin/.
# First save this file (kepler_02.py) at a certain place and then run the command given below at
# the command line.
#       $ sudo ln -s source_file(path) target_file(path)
# eg. In my Mac Book Pro (it should be the same for all the UNIX-based patforms)
#       $ sudo ln -s /Users/vijaykoju/kepler_02.py /usr/local/bin/kplrlca
# NOTE: the symbolic link is created with whatever name you give it while installation. I gave it
#       the name kplrlca(kepler light curve analyzer). You can give it any name you desire.
# NOTE: the kepler_02.py file that you saved in your computer must have executable access.
#       this can be done by running the given command at the command line, if it is not already
#       executable.
#       $ chmod a+x kepler_02.py
################################################################################################
# Required files for running this program.
# 1) kepIDs_with_period.txt (This file contains all the informations such as keplerID, light
#    curve type, period, fill out factor, magnitude, etc. All these informations can be found
#    at the website astro4.ast.villanova.edu/aprsa/kepler. This txt file basically contains
#    all the informations found in the given website.)
# 2) effTemp.txt (This file contains keplerID and effective tempereture.)
# 3) kepler_EBdata_Prsa_filenames.txt (This file contains the full names for all the light
#    curves with a serial number assigned to each (eg. 1  kplr001026032-2009259160929_llc.txt).
# all these files must be kept in the same folder where the light curve files are stored.
################################################################################################
#      kplr (kepler_02.py) is a program for analyzing the light curves of Eclipsing 
# Binaries from the KEPLER MISSION(KM). Light curve data (from the KM) can be downloaded
# from the official website of the KM (archive.stsci.edu/kepler/kepler_fov/search.php).
# The downloaded light curves are in FITS (Flexible Image Transport System) format, which
# need to be converted to txt files in order to use them in this program. Use unfits
# (unfits_09.py) program to convert FITS files into txt files.
#       This program (kplr) takes a light curve (txt format) with data for all the 
# available cycles as a whole, and then extract the data for each cycle and stores it
# temporarily in the memory. All the computations, such as fourier fitting or polynomial
# fitting, calculating delta-m max (measure of the O'Connell Effect), delta-m min, O'Connell
# Effect Ratio(OER), ight Curve Asymmetry(LCA), and errors, are done on cycle-wise basis.
# All these measurements are then saved in a txt file with their kepler id as their name
# (eg. 03321207.txt). For the light curve types OC(Over-Contact), SD(Semi-Detached),
# ELV(Ellipsoidal), and ?(Unknown) fourier fitting is done, whereas for the D(Detached)
# type polynomial fitting (only for the out of eclipse regions) is done. Due to this, for
# D type systems delta-m min, LCA, and OER are not computed. To plot the results from this
# program use keplerPlot.py.
#       It is intended to be run from the command line in a UNIX-like environment such as
# Linux, Mac OS X, Solaris, or FreeBSD. Run this program with the -h option or see the
# display_help() function below for further details.
#       You may use, adapt, and redistribute this software subject to the terms of the
# Creative Commons Attribution-ShareAlike 2.5 Generic License.
# --> See http://creativecommons.org/licenses/by-sa/2.5/ for details.
################################################################################################

from scipy import *
from numpy import *
from pylab import *
from scipy.optimize import leastsq
import sys, os

class Kepler():
    def __init__(self,filename):
        self.file_name = filename
        
    # Return the Kepler ID for the system from the filename.
    # @param kepler_filename
    # @return int keplerID.
    def getKeplerID(self):
        kepler_ID = self.file_name[5:13]
        return kepler_ID
    
    def getInfo(self,i):
        kepler_ID = self.file_name[5:13]
        infoFile = open('kepIDs_with_period.txt','r')
        id_type1 = kepler_ID + '.00'
        id_type2 = kepler_ID + '.01'
        id_type3 = kepler_ID + '.02'
        for line in infoFile:
            line = line.split()
            if id_type1 in line[0]:
                return eval(line[i])
            elif id_type2 in line[0]:
                return eval(line[i])
            elif id_type3 in line[0]:
                return eval(line[i])

    def getInfoS(self,i):
        kepler_ID = self.file_name[5:13]
        infoFile = open('kepIDs_with_period.txt','r')
        id_type1 = kepler_ID + '.00'
        id_type2 = kepler_ID + '.01'
        id_type3 = kepler_ID + '.02'
        for line in infoFile:
            line = line.split()
            if id_type1 in line[0]:
                return line[i]
            elif id_type2 in line[0]:
                return line[i]
            elif id_type3 in line[0]:
                return line[i]
            
    # Return the period of the system.
    # @param kepler_filename
    # @return int period
    def getPeriod(self):
        colnum = 3
        return self.getInfo(colnum)

    # Return BKJD0([B]arycentric [K]epler [J]ulian [D]ate) for the system at intial data acquisition.
    # @param kepler_filename
    # @return int BKJD0 value.
    def getBKJD0(self):
        colnum = 2
        return self.getInfo(colnum)
        

    # Return K-magnitude of the star.
    # @param kepler_filename
    # @return float Kmag.
    def getKmag(self):
        colnum = 4
        return self.getInfo(colnum)

    # Return type of the lightcurve.
    # @param kepler_filename
    # @return string type.
    def getLcType(self):
        colnum = 1
        return self.getInfoS(colnum)
            
    # Return type of the lightcurve.
    # @param kepler_filename
    # @return string type.
    def getTempRatio(self):
        colnum = 5
        return self.getInfoS(colnum)

    # Return type of the lightcurve.
    # @param kepler_filename
    # @return string type.
    def getSumOfFracRadii(self):
        colnum = 6
        return self.getInfoS(colnum)

    # Return type of the lightcurve.
    # @param kepler_filename
    # @return string type.
    def getMassRatio(self):
        colnum = 7
        return self.getInfoS(colnum)

    # Return type of the lightcurve.
    # @param kepler_filename
    # @return string type.
    def getSineEccentricity(self):
        colnum = 8
        return self.getInfoS(colnum)

    # Return type of the lightcurve.
    # @param kepler_filename
    # @return string type.
    def getCosEccentricity(self):
        colnum = 9
        return self.getInfoS(colnum)

    # Return type of the lightcurve.
    # @param kepler_filename
    # @return string type.
    def getSineOfInclination(self):
        colnum = 11
        return self.getInfoS(colnum)

    # Return fill out factor of the system.
    # @param kepler_filename
    # @return string FF.
    def getFillOutFactor(self):
        colnum = 10
        return self.getInfoS(colnum)

    def getEffectiveTemp(self):
        kepler_ID = self.file_name[5:13]
        infoFile1 = open('effTemp.txt','r')
        for line in infoFile1:
            line = line.split()
            if kepler_ID in line[0]:
                return line[1]


    ###################################################################################################
    #####################Start extracting data from the kepler data file.##############################
    ###################################################################################################

    # Return the number of periods between the BKJD0 given in the Prsa_Catalog and the first BKJD
    # entry in the keplerfile.
    # @param file_name
    # @return float p_gaps
    def period_gaps(self):
        period = self.getPeriod()
        BKJD0 = self.getBKJD0()
        kep_file = open(self.file_name, 'r')
        BKJD = []
        for line in kep_file:
            line = line.split()
            BKJD.append(eval(line[0]))
        p_gaps = (BKJD[0] - BKJD0)/period
        kep_file.close
        #print p_gaps
        return p_gaps

    # Return the first primary minimum data point from the keplerfile.
    # @param file_name
    # @return float frst_pmin, first primary min data.
    def first_pmin(self):
        period = self.getPeriod()
        BKJD0 = self.getBKJD0()
        p_gaps = self.period_gaps()
        if p_gaps == 0:
            frst_pmin = BKJD0
        elif p_gaps == 1.0:
            frst_pmin = BKJD0 + p_gaps*period
        else:
            ceiling = ceil(p_gaps)
            frst_pmin = BKJD0 + ceil(p_gaps)*period
        return frst_pmin
        
    # extract BKJD(column 1), ap_corr_flux(column 10), and ap_corr_err(column 11)
    # and convert the flux and error in column 10 and 11 into magnitude scale
    # using  magnitude = -2.5*log(flux) + C, where C is a constant
    def extract_data(self):
        avg_Kmag = self.getKmag()
        kep_file = open(self.file_name,'r')
        BKJD = []
        ap_corr_mag = [] #mag without adding C (differential mag)
        ap_corr_err = [] #mag without adding C (differential mag err)
        for line in kep_file:
            line = line.split()
            if line[9] != '-inf':
                BKJD.append(eval(line[0]))
                ap_corr_mag.append(-2.5*log10(eval(line[9]))) 
                ap_corr_err.append(2.5*log10((eval(line[9]) + eval(line[10]))/eval(line[9])))
        avg_all_mag = average(ap_corr_mag)
        C = avg_Kmag - avg_all_mag
        ap_corr_mag_withC = [mag + C  for mag in ap_corr_mag] #mag after adding C (absolute mag in K-Band)
        ap_corr_err_withC = [err  for err in ap_corr_err] #mag err after adding C (absolute mag err in K-Band)
        kep_file.close()
        print len(BKJD), len(ap_corr_mag_withC), len(ap_corr_err_withC)
        scatter(BKJD, ap_corr_mag_withC, s =3)
        plot(BKJD, ap_corr_mag_withC)
        ax = gca()
        ax.set_ylim(ax.get_ylim()[::-1])
        show()


    # Return only the first column from the kepler file
    # @param file_name
    # @return list BKJD_data only, excluding the invalid data
    def BKJD_data(self):
        kep_file = open(self.file_name,'r')
        BKJD = []
        for line in kep_file:
            line = line.split()
            if line[9] != '-inf':
                BKJD.append(eval(line[0]))
        kep_file.close()
        return BKJD

    def timecorr_data(self):
        kep_file = open(self.file_name,'r')
        timecorr = []
        for line in kep_file:
            line = line.split()
            if line[9] != '-inf':
                timecorr.append(eval(line[1]))
        kep_file.close()
        return timecorr

    
    # Return only the tenth column from the kepler file
    # @param file_name
    # @return list mag_data only, excluding the invalid data
    def mag_data(self):
        avg_Kmag = self.getKmag()
        kep_file = open(self.file_name,'r')
        ap_corr_mag = [] #mag without adding C (differential mag)
        for line in kep_file:
            line = line.split()
            if line[9] != '-inf':
                ap_corr_mag.append(-2.5*log10(eval(line[9]))) 
        avg_all_mag = average(ap_corr_mag)
        C = avg_Kmag - avg_all_mag
        ap_corr_mag_withC = [mag + C  for mag in ap_corr_mag] #mag after adding C (absolute mag in K-Band)
        kep_file.close()
        return ap_corr_mag_withC

    def flux_data(self):
        kep_file = open(self.file_name,'r')
        ap_corr_flux = [] 
        for line in kep_file:
            line = line.split()
            if line[9] != '-inf':
                ap_corr_flux.append(eval(line[9])) 
        kep_file.close()
        return ap_corr_flux

    # Return only the eleventh column from the kepler file
    # @param file_name
    # @return list magErr_data only, excluding the invalid data
    def magErr_data(self):
        avg_Kmag = self.getKmag()
        kep_file = open(self.file_name,'r')
        ap_corr_mag = [] #mag without adding C (differential mag)
        ap_corr_err = [] #mag without adding C (differential mag err)
        for line in kep_file:
            line = line.split()
            if line[9] != '-inf':
                ap_corr_mag.append(-2.5*log10(eval(line[9])))
                ap_corr_err.append(2.5*log10((eval(line[9]) + eval(line[10]))/eval(line[9])))
        avg_all_mag = average(ap_corr_mag)
        C = avg_Kmag - avg_all_mag
        ap_corr_err_withC = [err  for err in ap_corr_err] #mag err after adding C (absolute mag err in K-Band)
        kep_file.close()
        return ap_corr_err_withC

    def fluxErr_data(self):
        kep_file = open(self.file_name,'r')
        ap_corr_fluxErr = [] 
        for line in kep_file:
            line = line.split()
            if line[9] != '-inf':
                ap_corr_fluxErr.append(eval(line[10])) 
        kep_file.close()
        return ap_corr_fluxErr

    # @param self
    # @return [a float] number of data cycles available in the llc file.
    def numberOfCycles(self):
        period = self.getPeriod()
        firstBKJD = self.first_pmin()
        lastBKJD = self.BKJD_data()[-1]
        num_of_cycles = (lastBKJD - firstBKJD)/period
        return floor(num_of_cycles)

    # @param self, number of cycle in interest
    # @return [a list] [index_lower, index_upper] gives the index of the first
    # and the last data points of the cycle in interest.
    def get_index(self,cycle):
        fpmin = self.first_pmin()
        period = self.getPeriod()
        bkjd_list = self.BKJD_data()
        lower_min = fpmin + (cycle - 1)*period
        upper_min = fpmin + cycle*period
        alist = []
        for item in bkjd_list:
            if item >= lower_min:
                if item < lower_min + period/4:
                    alist.append(item - lower_min)
                else:break
        desired_bkjd_lower = alist[0] + lower_min
        index_lower = bkjd_list.index(desired_bkjd_lower)
        blist = []
        for item in bkjd_list:
            if item > upper_min - period/4:
                if item <= upper_min:
                    blist.append(upper_min - item)
                else:break
        desired_bkjd_upper = upper_min - blist[-1]
        index_upper = bkjd_list.index(desired_bkjd_upper)
        return [index_lower , index_upper + 1]

    # @param self, [an int] cycle number of interest
    # @return [a list] [BKJD] of the given cycle
    def data_perCycle_BKJD(self,cycle):
        bkjd_list = self.BKJD_data()
        BKJD_perCycle = []
        for i in range(self.get_index(cycle)[0],self.get_index(cycle)[1]):
            BKJD_perCycle.append(bkjd_list[i])
        return BKJD_perCycle

    # @param self, [an int] cycle number of interest
    # @return [a list] [mag] of the given cycle
    def data_perCycle_mag(self,cycle):
        maglist = self.mag_data()
        mag_perCycle = []
        for i in range(self.get_index(cycle)[0],self.get_index(cycle)[1]):
            mag_perCycle.append(maglist[i])
        return mag_perCycle

    def data_perCycle_flux(self,cycle):
        fluxlist = self.flux_data()
        flux_perCycle = []
        for i in range(self.get_index(cycle)[0],self.get_index(cycle)[1]):
            flux_perCycle.append(fluxlist[i])
        return flux_perCycle
        
    # @param self, [an int] cycle number of interest
    # @return [a list] [phased BKJD] of the given cycle
    def phased_BKJD(self,cycle):
        fpmin = self.first_pmin()
        period = self.getPeriod()
        mag_list = self.data_perCycle_mag(cycle)
        min_mag = max(mag_list)
        min_mag_index = mag_list.index(min_mag)
    
        BKJD_list = self.data_perCycle_BKJD(cycle)
        phase = []
        for item in BKJD_list:
            phase.append(math.fmod((item-fpmin)/period,1))
        return phase
    
    # @param self, [an int] cycle number of interest
    # @return [a list] [magErr] of the given cycle
    def data_perCycle_magErr(self,cycle):
        magErrlist = self.magErr_data()
        magErr_perCycle = []
        for i in range(self.get_index(cycle)[0],self.get_index(cycle)[1]):
            magErr_perCycle.append(magErrlist[i])
        return magErr_perCycle

    def data_perCycle_fluxErr(self,cycle):
        fluxErrlist = self.fluxErr_data()
        fluxErr_perCycle = []
        for i in range(self.get_index(cycle)[0],self.get_index(cycle)[1]):
            fluxErr_perCycle.append(fluxErrlist[i])
        return fluxErr_perCycle

    # @param self, [an int] cycle number of interest
    # @return [a plot] of the given cycle
    def plot_cycle(self,cycle):
        x = self.phased_BKJD(cycle)
        y = self.data_perCycle_flux(cycle)
        y_err = self.data_perCycle_fluxErr(cycle)
        err = .9
        errorbar(x,y,y_err)
        xlim(-0.2,1.2)
        ax = gca()
        ax.set_ylim(ax.get_ylim()[::-1])
        xlabel('Phase')
        ylabel('K-magnitude')
        title( '    KeplerID : ' + self.getKeplerID() + '  Period : ' +str(self.getPeriod())+'d'+'  Cycle : '+ str(cycle))
        show()



#######################################Fitting Functions############################################       
def fourier_func(x,p):
    return p[0] + p[1]*cos(1*2*pi*x) + p[2]*cos(2*2*pi*x) + p[3]*cos(3*2*pi*x) +\
           p[4]*cos(4*2*pi*x) + p[5]*cos(5*2*pi*x) + p[6]*cos(6*2*pi*x) + p[7]*cos(7*2*pi*x) +\
           p[8]*cos(8*2*pi*x) + p[9]*sin(1*2*pi*x) + p[10]*sin(2*2*pi*x) + p[11]*sin(3*2*pi*x)+\
           p[12]*sin(4*2*pi*x) + p[13]*sin(5*2*pi*x) + p[14]*sin(6*2*pi*x) + p[15]*sin(7*2*pi*x) +\
           p[16]*sin(8*2*pi*x)

def fourier_func1(x,p):
    return p[0] + p[1]*cos(1*2*pi*x) + p[2]*cos(2*2*pi*x) + p[3]*cos(3*2*pi*x) +\
           p[4]*cos(4*2*pi*x) + p[5]*sin(1*2*pi*x) + p[6]*sin(2*2*pi*x) + p[7]*sin(3*2*pi*x)+\
           p[8]*sin(4*2*pi*x) 

def polynomial_func(x,p):
    return p[0] + p[1]*x + p[2]*x**2 + p[3]*x**3
####################################################################################################

#Plot the light curve of the specified cycle number of the given kepler file.
#@param filenum(kepler file S.N.), cyclenum(desired cycle number)
def plotCycle(filenum,cyclenum):
    opnfile = open('kepler_EBdata_Prsa_filenames.txt','r')
    fylnum = []
    fylname = []
    for line in opnfile:
        line = line.split()
        fylnum.append(line[0])
        fylname.append(line[1])
        
    kepler_file = fylname[filenum - 1]
    
    kep = Kepler(kepler_file) #opens the kepler file of the given file number.
    num_cycle = int(kep.numberOfCycles())
    lcType = kep.getLcType()
    opnfile.close()
    if lcType == 'OC' or lcType == 'SD' or lcType == 'ELV' or lcType == '?':
        try:
            x = np.array(kep.phased_BKJD(cyclenum))
            data = np.array(kep.data_perCycle_mag(cyclenum))
            yerr = np.array(kep.data_perCycle_magErr(cyclenum))

            ###If number of data points in each cycle is less than 17(no. of coeffients in the fitting function,
            ###use a fitting function with less no. of coefficients.
            if len(x)<17:
                ########## for magnitude ###########
                def residuals(p,data,x):
                    err = data - fourier_func1(x,p)
                    return err
                

                p01 = [1.37980075e+01,6.74224644e-02,-1.41086438e-01,-1.42023067e-01, \
                       -2.56529287e-02,-1.83357748e-01,-1.73075867e-01,4.10602345e-02,\
                       1.57625025e-01]
                pbest = leastsq(residuals,p01,args=(data,x),full_output = 1)#leastsq optimization for fitting
        
                bestparams = pbest[0]  #a list of best fitting coefficients.
        
                x_new = arange(0,1,0.001)
                datafit = fourier_func1(x_new,bestparams) #fitted curve
            else:
                ##### when no. of data points in each cycle is more than 17, use a 17-term fourier fitting function.

                ########## for magnitude ###########
                def residuals(p,data,x):
                    err = data - fourier_func(x,p)
                    return err
                

                p0 = [1.37980075e+01,6.74224644e-02,-1.41086438e-01,-1.42023067e-01, \
                      -2.56529287e-02,9.47168094e-02,9.46511622e-02,-9.29974207e-03,\
                      -6.01895557e-02,-1.83357748e-01,-1.73075867e-01,4.10602345e-02,\
                      1.57625025e-01,6.22460398e-02,-5.78687258e-02,0-7.23577220e-02,-1.95020746e-02]
                
                
                pbest = leastsq(residuals,p0,args=(data,x),full_output = 1)
            
                bestparams = pbest[0]
            
                x_new = arange(0,1,0.001)
                datafit = fourier_func(x_new,bestparams)
            clf()
            errorbar(x,data,yerr,fmt='bx')
            plot(x_new,datafit,'r')
            ax = gca()
            ax.set_ylim(ax.get_ylim()[::-1])
            xlabel('Phase')
            ylabel('K-magnitude')
            title( '    S.N.: ' + str(cyclenum) + '  KeplerID : ' + str(kep.getKeplerID()) + '  Period : ' + str(kep.getPeriod()) +'d'+\
                       '  LC-Type : ' + str(kep.getLcType()) + '   EffT : ' +str(kep.getEffectiveTemp())+'K', fontsize =12)
            grid(True)
            show()
            #savefig('/Destination path/'+str(cyclenum)+'.png', format = 'png')
            
        except (IndexError):pass
    if lcType == 'D':
        try:
            ##########################polynomial fitting for phase<0.5################################
            x = kep.phased_BKJD(cyclenum)
            y = kep.data_perCycle_mag(cyclenum)
            yerr = kep.data_perCycle_magErr(cyclenum)
            x_trimmed = []
            for item in x:
                if item >= .15:
                    if item <= .4:
                        x_trimmed.append(item)
            ind_first = x.index(x_trimmed[0])
            ind_last = x.index(x_trimmed[-1])
            y_trimmed = []

            for i in range(ind_first,ind_last+1):
                y_trimmed.append(y[i])
            x_array = array(x_trimmed)
            y_array = array(y_trimmed)
            def residuals(p,y_array,x_array):
                err = y_array - polynomial_func(x_array,p)
                return err
            p0 = [.23,-.02,-.5,.8]
            pbest = leastsq(residuals,p0,args=(y_array,x_array),full_output = 1)
            bestparams = pbest[0]
            x_new = arange(.15,.4,.001)
            datafit = polynomial_func(x_new,bestparams)
            
            ##############################polynomial fitting for phase>0.5#############################
            x_up_trimmed = []
            for item in x:
                if item >= .6:
                    if item <= .85:
                        x_up_trimmed.append(item)
            ind_up_first = x.index(x_up_trimmed[0])
            ind_up_last = x.index(x_up_trimmed[-1])
            y_up_trimmed = []
            for i in range(ind_up_first,ind_up_last+1):
                y_up_trimmed.append(y[i])
            x_up_array = array(x_up_trimmed)
            y_up_array = array(y_up_trimmed)
            def residuals1(p,y_up_array,x_up_array):
                err = y_up_array - polynomial_func(x_up_array,p)
                return err
            p0 = [.23,-.02,-.5,.8]
            pbest_up = leastsq(residuals1,p0,args=(y_up_array,x_up_array),full_output = 1)
            bestparams_up = pbest_up[0]
            x_up_new = arange(.6,.85,.001)
            datafit_up = polynomial_func(x_up_new,bestparams_up)
            clf()
            errorbar(x,y,yerr,fmt='bx')
            plot(x_new,datafit,'r',x_up_new,datafit_up,'r')
            ax = gca()
            ax.set_ylim(ax.get_ylim()[::-1])
            xlabel('Phase')
            ylabel('K-magnitude(mag)')
            title( '    S.N.: ' + str(cyclenum) + '  KeplerID : ' + str(kep.getKeplerID()) + '  Period : ' + str(kep.getPeriod()) +'d'+\
                       '  LC-Type : ' + str(kep.getLcType()) + '   EffT : ' +str(kep.getEffectiveTemp())+ 'K', fontsize =12)
            grid(True)
            show()
            #savefig('/Destination path/'+str(cyclenum)+'.png', format = 'png')
            
        except IndexError, ValueError:pass

#This is the process which is carried out when the the program is called from the command line.
#@param file number associated to each kepler file, which is given listed in kepler_EBdata_Prsa_filenames.txt
def process(filenum):
    opnfile = open('kepler_EBdata_Prsa_filenames.txt','r')
    fylnum = []
    fylname = []
    for line in opnfile:
        line = line.split()
        fylnum.append(line[0])
        fylname.append(line[1])

    kepler_file = fylname[filenum - 1]
    
    kep = Kepler(kepler_file) #opens the kepler file of the given file number.
    num_cycle = int(kep.numberOfCycles())
    lcType = kep.getLcType()
    opnfile.close()
    
    #################################################################################################
    ###################################begin fourier fitting#########################################
    #################################################################################################
    
    
    if lcType == 'OC' or lcType == 'SD' or lcType == 'ELV' or lcType == '?':
        print ''
        print '#####################################################################################'
        print '######################## Kepler Light Curve Information #############################'
        print ''
        print 'KeplerID : ' + kep.getKeplerID() + '    Period(days) : ' + str(kep.getPeriod()) + '    Effective Temp(K) : ' + kep.getEffectiveTemp() + '    S.N : ' + str(filenum)
        print ''
        print 'T2/T1 : ' + kep.getTempRatio() + '    r1 + r2 : ' + kep.getSumOfFracRadii() + '    q : ' + kep.getMassRatio() + '    FF : ' + kep.getFillOutFactor()
        print ''
        print '[e sin(w), e cos(w)] : ' + '[' + kep.getSineEccentricity(), kep.getCosEccentricity() + ']' + '    sin(i) : ' + kep.getSineOfInclination()
        print ''
        print 'LC-type : ' + lcType + '    K-magnitude(mag) : ' + str(kep.getKmag()) + '    No. of cycles : ' + str(kep.numberOfCycles())
        print ''
        print '#####################################################################################'
        print ''

        destination = open(str(kep.getKeplerID())+'.txt','a')#creates a new file to write in the results.
        cycle_list = []
        deltam_list = []
        deltam_listerr = []
        for i in range(1,num_cycle +1 ):
            cycle = i
            try:
                x = np.array(kep.phased_BKJD(i))
                data = np.array(kep.data_perCycle_mag(i))
                yerr = np.array(kep.data_perCycle_magErr(i))
                fluxdata = np.array(kep.data_perCycle_flux(i))
                fluxerr = np.array(kep.data_perCycle_fluxErr(i))
                min_flux = min(fluxdata)
                max_yerr = max(yerr)
                max_fluxerr = max(fluxerr)

                ###If number of data points in each cycle is less than 17(no. of coeffients in the fitting function,
                ###use a fitting function with less no. of coefficients.
                
                if len(x)<17:
                    ########## for magnitude ###########
                    def residuals(p,data,x):
                        err = data - fourier_func1(x,p)
                        return err
                    

                    p0 = [1.37980075e+01,6.74224644e-02,-1.41086438e-01,-1.42023067e-01, \
                          -2.56529287e-02,-1.83357748e-01,-1.73075867e-01,4.10602345e-02,\
                          1.57625025e-01]
                    pbest = leastsq(residuals,p0,args=(data,x),full_output = 1)#leastsq optimization for fitting
            
                    bestparams = pbest[0]  #a list of best fitting coefficients.
            
                    x_new = arange(0,1,0.001)
                    datafit = fourier_func1(x_new,bestparams) #fitted curve

                    ########## for flux ################
                    def fresiduals(p,fluxdata,x):
                        err = fluxdata- fourier_func1(x,p)
                        return err

                    p02 = [7.60968946e-01,-1.18870563e-02,-2.08312226e-01,-6.96007158e-03, \
                           -4.04757077e-02,6.28220746e-03, 2.02293971e-02,-1.59931637e-03, \
                           4.45437189e-03]
                    pfbest = leastsq(fresiduals,p02,args=(fluxdata,x),full_output = 1)
            
                    fbestparams = pfbest[0]
            
                    fdatafit = fourier_func1(x_new,fbestparams)
                    
                else:
                    ##### when no. of data points in each cycle is more than 17, use a 17-term fourier fitting function.

                    ########## for magnitude ###########
                    def residuals(p,data,x):
                        err = data - fourier_func(x,p)
                        return err
                    

                    p0 = [1.37980075e+01,6.74224644e-02,-1.41086438e-01,-1.42023067e-01, \
                          -2.56529287e-02,9.47168094e-02,9.46511622e-02,-9.29974207e-03,\
                          -6.01895557e-02,-1.83357748e-01,-1.73075867e-01,4.10602345e-02,\
                          1.57625025e-01,6.22460398e-02,-5.78687258e-02,0-7.23577220e-02,-1.95020746e-02]
                    
                    
                    pbest = leastsq(residuals,p0,args=(data,x),full_output = 1)
                
                    bestparams = pbest[0]
                
                    x_new = arange(0,1,0.001)
                    datafit = fourier_func(x_new,bestparams)

                    ########## for flux ################
                    def fresiduals(p,fluxdata,x):
                        err = fluxdata- fourier_func1(x,p)
                        return err

                    p03 = [7.60968946e-01,  -1.18870563e-02,  -2.08312226e-01,  -6.96007158e-03, \
                           -4.04757077e-02,  -4.05274560e-03,  -1.80543213e-02,  -5.79266951e-04, \
                           -7.20603774e-03,   6.28220746e-03,   2.02293971e-02,  -1.59931637e-03, \
                           4.45437189e-03,  -1.10370228e-03,   4.79737499e-03,  1.48869172e-03, \
                           1.10595143e-03]
                    pfbest = leastsq(fresiduals,p03,args=(fluxdata,x),full_output = 1)
            
                    fbestparams = pfbest[0]
            
                    fdatafit = fourier_func1(x_new,fbestparams)
            
            
                #####################################################################################
                ##############################calculate deltammax####################################
                #####################################################################################
                firsthalf = []
                for i in range(len(datafit)/2):
                    firsthalf.append(datafit[i])
                secondhalf = []
                for i in range(len(datafit)/2,len(datafit)):
                    secondhalf.append(datafit[i])
                q1 =firsthalf[:len(firsthalf)/2]+secondhalf[len(secondhalf)/2:]
                q2 =firsthalf[len(firsthalf)/2:]+secondhalf[:len(secondhalf)/2]
                min_q1 = max(q1)
                minq1_round = '%4.1f' % min_q1
                min_q2 = max(q2)
                minq2_round = '%4.1f' % min_q2
                max_firsthalf = min(firsthalf)
                maxq1_round = '%4.1f' % max_firsthalf
                
                max_secondhalf = min(secondhalf)
                maxq2_round = '%4.1f' % max_secondhalf
                deltammax = -(max_firsthalf - max_secondhalf)*1000
                deltammax_err = sqrt(max_yerr**2 + max_yerr**2)*1000 # error in deltammax
                deltammin = (min_q1 - min_q2)*1000
                deltammin_err = sqrt(max_yerr**2 + max_yerr**2)*1000
                dmaxerr_round = '%1.1f' % deltammax_err
                dmax_round = '%4.1f' % deltammax
                dminerr_round = '%1.1f' % deltammin_err
                dmin_round =  '%4.1f' % deltammin
                cycle_list.append(cycle)
                deltam_list.append(eval(dmax_round))
                deltam_listerr.append(eval(dmaxerr_round))

                #####################################################################################
                ###### calculate O'Connell Effect Ration(OER) and Light Curve Asymmetry(lca) ########
                #####################################################################################
                bins = 36
                datanum = len(x_new)/bins
                x_avg = []
                flux_avg_minusmin = []
                flx_avg = []
                for i in range(bins):
                    x_binned = [x_new[j] for j in range(i*datanum, (i+1)*datanum)]
                    flux_binned = [ fdatafit[j]  for j in range(i*datanum,(i+1)*datanum)]
                    xavg = average(x_binned)
                    fluxavg = average(flux_binned)
                    x_avg.append(xavg)
                    flux_avg_minusmin.append(fluxavg-min_flux)
                    flx_avg.append(fluxavg)
                    
                halflength = bins/2
                fluxSum_firsthalf = math.fsum(flux_avg_minusmin[:halflength])
                fluxSum_Secondhalf = math.fsum(flux_avg_minusmin[halflength:])
                OER = fluxSum_firsthalf/fluxSum_Secondhalf
                OER_round = '%1.2f' % OER
                OERerr = OER * (18*(max_fluxerr/sqrt(27)+max_fluxerr)/fluxSum_firsthalf +\
                                18*(max_fluxerr/sqrt(27)+max_fluxerr)/fluxSum_Secondhalf)
                OERerr_round = '%1.2f' % OERerr
                numerator_lcalist = [(flx_avg[i] - flx_avg[-(i+1)])**2/(flx_avg[i])**2   for i in range(halflength)]
                numerator_lca = math.fsum(numerator_lcalist)
                lca = sqrt((1./bins)*numerator_lca)*100
                plca = [ (((4*max_fluxerr/sqrt(27))/(flx_avg[i]+flx_avg[-(i+1)]))+((2*max_fluxerr/sqrt(27))/flx_avg[i])) for i in range(halflength)]
                plca_sum = math.fsum(plca)
                lcaerr = lca*(2./bins)*plca_sum*100
                lca_round = '%4.2f' % lca
                lcaerr_round = '%1.2f' % lcaerr
                
                if len(x)>=17:
                    destination.write(str(cycle) +' ' + str(dmax_round) + ' ' + str(dmaxerr_round) + ' ' + str(dmin_round) + ' ' +str(dminerr_round) + ' ' +\
                                      str(OER_round) + ' ' +str(OERerr_round) + ' ' + str(lca_round) + ' ' +str(lcaerr_round) + ' ' + str(max_firsthalf)+" " +\
                                      str(max_secondhalf) + ' ' + str(min_q1) + ' ' + str(min_q2) + ' ' + str(bestparams[0]) + ' ' + str(bestparams[1]) + ' ' +\
                                      str(bestparams[2]) + ' ' + str(bestparams[3]) + ' ' + str(bestparams[4]) + ' ' + str(bestparams[5]) + ' ' + str(bestparams[6]) + ' ' +\
                                      str(bestparams[7]) + ' ' + str(bestparams[8]) + ' ' + str(bestparams[9]) + ' ' + str(bestparams[10]) + ' ' + str(bestparams[11]) + ' ' +\
                                      str(bestparams[12]) + ' ' + str(bestparams[13]) + ' ' + str(bestparams[14]) + ' ' + str(bestparams[15]) + ' ' + str(bestparams[16]) + '\n')
                if len(x)<17:
                    destination.write(str(cycle) +' ' + str(dmax_round) + ' ' + str(dmaxerr_round) + ' ' + str(dmin_round) + ' ' +str(dminerr_round) + ' ' +\
                                      str(OER_round) + ' ' +str(OERerr_round) + ' ' + str(lca_round) + ' ' +str(lcaerr_round) + ' ' + str(max_firsthalf) + ' ' +\
                                      str(max_secondhalf) + ' ' + str(min_q1) + ' ' + str(min_q2) + ' ' + str(bestparams[0]) + ' ' + str(bestparams[1]) + ' ' +\
                                      str(bestparams[2]) + ' ' + str(bestparams[3]) + ' ' + str(bestparams[4]) + ' ' + str(bestparams[5]) + ' ' + str(bestparams[6]) + ' ' +\
                                      str(bestparams[7]) + ' ' + str(bestparams[8]) + '\n')

                #####################################################################################
                ############################### Display on the screen ###############################
                #####################################################################################
                print 'S.N. = ' +str(cycle) +'   deltammax = ' +str(dmax_round) + ' +/- ' +str(dmaxerr_round) + '   deltammin = ' +str(dmin_round) + ' +/- ' +str(dminerr_round) +\
                      '   OER = ' +str(OER_round) + ' +/- ' +str(OERerr_round) + '   LCA = ' +str(lca_round) + ' +/- ' +str(lcaerr_round)  
                
            except (IndexError):pass

    ######################################################################################################
    ##########################################begin polynomial fit########################################
    ######################################################################################################

    if lcType == 'D':
        print ''
        print '#####################################################################################'
        print '######################## Kepler Light Curve Information #############################'
        print ''
        print 'KeplerID : ' + kep.getKeplerID() + '    Period(days) : ' + str(kep.getPeriod()) + '    Effective Temp(K) : ' + kep.getEffectiveTemp() + '    S.N : ' + str(filenum)
        print ''
        print 'T2/T1 : ' + kep.getTempRatio() + '    r1 + r2 : ' + kep.getSumOfFracRadii() + '    q : ' + kep.getMassRatio() + '    FF : ' + kep.getFillOutFactor()
        print ''
        print '[e sin(w), e cos(w)] : ' + '[' + kep.getSineEccentricity(), kep.getCosEccentricity() + ']' + '    sin(i) : ' + kep.getSineOfInclination()
        print ''
        print 'LC-type : ' + lcType + '    K-magnitude(mag) : ' + str(kep.getKmag()) + '    No. of cycles : ' + str(kep.numberOfCycles())
        print ''
        print '#####################################################################################'
        print ''
        destination = open(str(kep.getKeplerID())+'.txt','a')
        cycle_list = []
        deltam_list = []
        deltam_listerr = []
        for i in range(1,num_cycle + 1):
            cycle = i
            try:
                ##########################polynomial fitting for phase<0.5################################
                x = kep.phased_BKJD(i)
                y = kep.data_perCycle_mag(i)
                yerr = kep.data_perCycle_magErr(i)
                max_yerr = max(yerr)
                x_trimmed = []
                for item in x:
                    if item >= .15:
                        if item <= .4:
                            x_trimmed.append(item)
                ind_first = x.index(x_trimmed[0])
                ind_last = x.index(x_trimmed[-1])
                y_trimmed = []

                for i in range(ind_first,ind_last+1):
                    y_trimmed.append(y[i])
                x_array = array(x_trimmed)
                y_array = array(y_trimmed)
                def residuals(p,y_array,x_array):
                    err = y_array - polynomial_func(x_array,p)
                    return err
                p0 = [.23,-.02,-.5,.8]
                pbest = leastsq(residuals,p0,args=(y_array,x_array),full_output = 1)
                bestparams = pbest[0]
                x_new = arange(.15,.4,.001)
                datafit = polynomial_func(x_new,bestparams)
                max1 = min(datafit)
                max1_round = '%4.1f' % max1
                ##############################polynomial fitting for phase>0.5#############################
                x_up_trimmed = []
                for item in x:
                    if item >= .6:
                        if item <= .85:
                            x_up_trimmed.append(item)
                ind_up_first = x.index(x_up_trimmed[0])
                ind_up_last = x.index(x_up_trimmed[-1])
                y_up_trimmed = []
                for i in range(ind_up_first,ind_up_last+1):
                    y_up_trimmed.append(y[i])
                x_up_array = array(x_up_trimmed)
                y_up_array = array(y_up_trimmed)
                def residuals1(p,y_up_array,x_up_array):
                    err = y_up_array - polynomial_func(x_up_array,p)
                    return err
                p0 = [.23,-.02,-.5,.8]
                pbest_up = leastsq(residuals1,p0,args=(y_up_array,x_up_array),full_output = 1)
                bestparams_up = pbest_up[0]
                x_up_new = arange(.6,.85,.001)
                datafit_up = polynomial_func(x_up_new,bestparams_up)
                max2 = min(datafit_up)
                max2_round = '%4.1f' % max2
                deltammax = -(max1 - max2)*1000
                deltammax_err = sqrt(max_yerr**2 + max_yerr**2)*1000
                derr_round = '%1.1f' % deltammax_err
                d_round = '%4.1f' % deltammax
                cycle_list.append(cycle)
                deltam_list.append(eval(d_round))
                deltam_listerr.append(eval(derr_round))
                destination.write(str(cycle) +' ' + str(d_round) + ' ' + str(derr_round) + ' ' + str(max1) + ' ' + str(max2) + ' ' +\
                                  str(bestparams[0]) + ' ' + str(bestparams[1]) + ' ' + str(bestparams[2]) + ' ' + str(bestparams[3]) + ' ' +\
                                  str(bestparams_up[0]) + ' ' + str(bestparams_up[1]) + ' ' + str(bestparams_up[2]) + ' ' + str(bestparams_up[3]) + '\n')
                #####################################################################################
                ################################Display on the screen################################
                #####################################################################################
                print 'S.N. = ' + str(cycle) + '    deltammax = ' + str(d_round) + ' +/- ' + str(derr_round)

            except IndexError, ValueError:pass

        

def info_gen(filenum):
    opnfile = open('kepler_EBdata_Prsa_filenames.txt','r')
    fylnum = []
    fylname = []
    for line in opnfile:
        line = line.split()
        fylnum.append(line[0])
        fylname.append(line[1])
    kepler_file = fylname[filenum - 1]
    
    kep = Kepler(kepler_file)
    num_cycle = int(kep.numberOfCycles())
    lcType = kep.getLcType()
    opnfile.close()
    print ''
    print '#####################################################################################'
    print '######################## Kepler Light Curve Information #############################'
    print ''
    print 'KeplerID : ' + kep.getKeplerID() + '    Period(days) : ' + str(kep.getPeriod()) + '    Effective Temp(K) : ' + str(kep.getEffectiveTemp()) + '    S.N : ' + str(filenum)
    print ''
    print 'T2/T1 : ' + kep.getTempRatio() + '    r1 + r2 : ' + kep.getSumOfFracRadii() + '    q : ' + kep.getMassRatio() + '    FF : ' + kep.getFillOutFactor()
    print ''
    print '[e sin(w), e cos(w)] : ' + '[' + kep.getSineEccentricity(), kep.getCosEccentricity() + ']' + '    sin(i) : ' + kep.getSineOfInclination()
    print ''
    print 'LC-type : ' + lcType + '    K-magnitude(mag) : ' + str(kep.getKmag()) + '    No. of cycles : ' + str(kep.numberOfCycles())
    print ''
    print '#####################################################################################'
    print ''

# Display detailed help
def display_help():
  print ""
  print "   ##################################  kepler_04.py  #############################################"
  print "   # File:              kepler_04.py"
  print "   # Author:            Vijay Koju (vk8736@truman.edu)"
  print "   # Created:           11 June 2011"
  print "   # Last modified:     29 July 2011"
  print "   #"
  print "   # Please feel free to contact me if you have any question regarding this program."
  print "   #"
  print "   ################################################################################################"
  print "   # Required modules for running this program."
  print "   # scipy, numpy, matplotlib"
  print "   ################################################################################################"
  print "   # Required files for running this program."
  print "   # 1) kepIDs_with_period.txt (This file contains all the informations such as keplerID, light"
  print "   #    curve type, period, fill out factor, magnitude, etc. All these informations can be found"
  print "   #    at the website astro4.ast.villanova.edu/aprsa/kepler. This txt file basically contains"
  print "   #    all the informations found in the given website.)"
  print "   # 2) effTemp.txt (This file contains keplerID and effective tempereture.)"
  print "   # 3) kepler_EBdata_Prsa_filenames.txt (This file contains the full names for all the light"
  print "   #    curves with a serial number assigned to each (eg. 1  kplr001026032-2009259160929_llc.txt)."
  print "   # all these files must be kept in the same folder where the light curve files are stored."
  print "   ################################################################################################"
  print "   #      kplr (kepler_02.py) is a program for analyzing the light curves of Eclipsing"
  print "   # Binaries from the KEPLER MISSION(KM). Light curve data (from the KM) can be downloaded"
  print "   # from the official website of the KM (archive.stsci.edu/kepler/kepler_fov/search.php)."
  print "   # The downloaded light curves are in FITS (Flexible Image Transport System) format, which"
  print "   # need to be converted to txt files in order to use them in this program. Use unfits"
  print "   # (unfits_09.py) program to convert FITS files into txt files."
  print "   #       This program (kplr) takes a light curve (txt format) with data for all the"
  print "   # available cycles as a whole, and then extract the data for each cycle and stores it"
  print "   # temporarily in the memory. All the computations, such as fourier fitting or polynomial"
  print "   # fitting, calculating delta-m max (measure of the O'Connell Effect), delta-m min, O'Connell"
  print "   # Effect Ratio(OER), ight Curve Asymmetry(LCA), and errors, are done on cycle-wise basis."
  print "   # All these measurements are then saved in a txt file with their kepler id as their name"
  print "   # (eg. 03321207.txt). For the light curve types OC(Over-Contact), SD(Semi-Detached),"
  print "   # ELV(Ellipsoidal), and ?(Unknown) fourier fitting is done, whereas for the D(Detached)"
  print "   # type polynomial fitting (only for the out of eclipse regions) is done. Due to this, for"
  print "   # D type systems delta-m min, LCA, and OER are not computed.To plot the results from this"
  print "   # program use keplerPlot.py"
  print "   #       It is intended to be run from the command line in a UNIX-like environment such as"
  print "   # Linux, Mac OS X, Solaris, or FreeBSD."
  print "   #       You may use, adapt, and redistribute this software subject to the terms of the"
  print "   # Creative Commons Attribution-ShareAlike 2.5 Generic License."
  print "   # --> See http://creativecommons.org/licenses/by-sa/2.5/ for details."
  print "   ################################################################################################"
  print "USAGE: kplr [OPTIONS] [FILES]"
  print ""
  print "Available options:"
  print "   -h              Print this page."
  print "   -i              Print information about the kepler light curve."
  print "   -p              Display the light curve of the specified cycle for the kepler file."
  print ""
  print "   Ex: kplr -i 21"
  print "       > Print information about the kepler light curve of S.N. 21."
  print "       > Here is what is prints for S.N. 21."
  print "       #####################################################################################"
  print "       ######################## Kepler Light Curve Information #############################"
  print ""
  print "       KeplerID : 02305277    Period(days) : 0.378757    Effective Temp(K) : 5694    S.N : 21"
  print ""
  print "       T2/T1 : 0.99961    r1 + r2 : N/A    q : 1.11416    FF : 0.88203"
  print ""
  print "       [e sin(w), e cos(w)] : [N/A N/A]    sin(i) : 0.42018"
  print ""
  print "       LC-type : OC    K-magnitude(mag) : 15.101    No. of cycles : 233.0"
  print ""
  print "       #####################################################################################"
  print ""
  print "   Ex: kplr 21"
  print "       > Print information about the kepler light curve of S.N. 1."
  print "       > Display the measuments for each cycle."
  print "       > Here is what is prints for S.N. 21. "
  print "       #####################################################################################"
  print "       ######################## Kepler Light Curve Information #############################"
  print ""
  print "       KeplerID : 02305277    Period(days) : 0.378757    Effective Temp(K) : 5694    S.N : 21"
  print ""
  print "       T2/T1 : 0.99961    r1 + r2 : N/A    q : 1.11416    FF : 0.88203"
  print ""
  print "       [e sin(w), e cos(w)] : [N/A N/A]    sin(i) : 0.42018"
  print ""
  print "       LC-type : OC    K-magnitude(mag) : 15.101    No. of cycles : 233.0"
  print ""
  print "       #####################################################################################"
  print "       S.N. = 1   deltammax = -6.8 +/- 0.8   deltammin = -15.3 +/- 0.8   OER = 0.91 +/- 0.03   LCA = 0.83 +/- 0.03"
  print "       S.N. = 2   deltammax = -6.1 +/- 0.8   deltammin = -15.2 +/- 0.8   OER = 0.94 +/- 0.03   LCA = 0.82 +/- 0.03"
  print "       S.N. = 3   deltammax = -4.4 +/- 0.8   deltammin = -15.0 +/- 0.8   OER = 0.96 +/- 0.03   LCA = 0.79 +/- 0.03"
  print "       S.N. = 4   deltammax = -4.7 +/- 0.8   deltammin = -10.8 +/- 0.8   OER = 0.97 +/- 0.03   LCA = 0.78 +/- 0.03"
  print "       S.N. = 5   deltammax = -0.8 +/- 0.8   deltammin = -12.6 +/- 0.8   OER = 1.00 +/- 0.03   LCA = 0.76 +/- 0.03"
  print "       and so on..."
  print ""
  print "   Ex. kplr -p 21 3"
  print "       > Displays the light curve (cycle number 3) of the file with S.N. 21."
  print ""

# Display usage information
def display_usage():
  print ''
  print "USAGE: kplr [OPTIONS] [FILES]"
  print "Type \'kplr -h\' for detailed help."
  print ''

# Initialize flags for parsing arguments...
title_flag = True
info_flag = False
empty = True

# Parse arguments...
for argument in sys.argv:
  if title_flag:
    title_flag = False
    continue
  if argument == '-h':
    empty = False
    display_help()
    break
  if argument == '-i':
    info_flag = True
    continue
  if info_flag:
    empty = False
    info_gen(eval(argument))
    continue
  if argument == '-p':
      empty = False
      plotCycle(eval(sys.argv[2]),eval(sys.argv[3]))
      break
  empty = False

  process(eval(argument))

if empty:
  display_usage()


###################################################################################################################
###################################################################################################################
# if you need to analyze many kepler light curves and you want to automate the program, so that you don't have to
# run the program again and again after each file, then uncomment the following commands and run the program from
# the python interactive shell (not from the terminal command line), or run the python interactive shell through
# the terminal. Using the interactive shell from the terminal is a better option as it is much faster than the
# interactive shell embedded in python. Please install ipython for better performance.
# To run python interactive shell from the termianl:
# $ ipython -pylab , or $ python  if you don't have ipython installed
# >> run kepler_02.py
#
# NOTE: OClistOutput.txt file is required, which has three columns
# Column 1: Serial No. (Just to know how many files are there in total)
# Column 2: KeplerID (eight digits)
# Cloumn 3: S.N. of the kepler file of the given keplerID (look up in kepler_EBdata_Prsa_filenames.txt)
###################################################################################################################

##if __name__ == '__main__':
##    f1= open('OClistOutput.txt','r')
##    snnum = []
##    for line in f1:
##        line = line.split()
##        snnum.append(eval(line[2]))
##    for item in snnum:
##        process(item)

###################################################################################################################


###################################################################################################################
###################################################################################################################
# if you want to make plots of light curves for all the cyles of any kepler file, then uncomment the following
# commands and run the program from the python interactive shell (not from the terminal command line), or run the
# python interactive shell through the terminal. Using the interactive shell from the terminal is a better option
# as it is much faster than the interactive shell embedded in python. Please install ipython for better performance.
# To run python interactive shell from the termianl:
# $ ipython -pylab , or $ python  if you don't have ipython installed
# >> run kepler_02.py
#
# NOTE: uncomment the savefig() command in the function plotCycle(filenum,cyclenum) if you want to save the graphs, and give
#       the path to the destination where you want to save the graphs.
###################################################################################################################

##if __name__=='__main__':
##    os.mkdir('/Destination Directory path/')
##    for i in range(1,lastcyclenum):  # lastcyclenum is the total no. of cycles present in the kepler file. this can be
##                                     # running kplrlca -i S.N. in the command line.
##        plotCycle(S.N.,i)            # S.N. is the serial number given to the kepler file
##                                     # (look up in kepler_EBdata_Prsa_filenames.txt)

###################################################################################################################
