#################################  kepler_04.py  #############################################
 File:              kepler_04.py

 Author:            Vijay Koju (vjk8736@gmail.com)

 Created:           11 June 2011

 Please feel free to contact me if you have any question regarding this program.

###############################################################################################

 Required modules for running this program.

 scipy, numpy, matplotlib"

 scipy, and numpy can be downloaded from scipy.org

 matplotlib can be downloaded from matplotlib.sourceforge.net/index.html

###############################################################################################

 Installation:

 In order to run this program from the command line, a symbolic link to this program has to be
 created in /usr/local/bin/.

 First save this file (kepler_04.py) at a certain place and then run the command given below at
 the command line.

       $ sudo ln -s source_file(path) target_file(path)

 eg. In my Mac Book Pro (it should be the same for all the UNIX-based patforms)

       $ sudo ln -s /Users/vijaykoju/kepler_04.py /usr/local/bin/kplrlca

 NOTE: the symbolic link is created with whatever name you give it while installation. I gave it
       the name kplrlca(kepler light curve analyzer). You can give it any name you desire.

 NOTE: the kepler_04.py file that you saved in your computer must have executable access.
       this can be done by running the given command at the command line, if it is not already
       executable.

       $ chmod a+x kepler_04.py

###############################################################################################

 Required files for running this program."

 1) kepIDs_with_period.txt (This file contains all the informations such as keplerID, light
    curve type, period, fill out factor, magnitude, etc. All these informations can be found
    at the website astro4.ast.villanova.edu/aprsa/kepler. This txt file basically contains
    all the informations found in the given website.)

 2) effTemp.txt (This file contains keplerID and effective tempereture.)

 3) kepler_EBdata_Prsa_filenames.txt (This file contains the full names for all the light
    curves with a serial number assigned to each (eg. 1  kplr001026032-2009259160929_llc.txt).

 All these files must be kept in the same folder where the light curve files are stored.

###############################################################################################

 			kplrlca (kepler_04.py) is a program for analyzing the light curves of Eclipsing
 Binaries from the KEPLER MISSION(KM). Light curve data (from the KM) can be downloaded
 from the official website of the KM (archive.stsci.edu/kepler/kepler_fov/search.php).
 The downloaded light curves are in FITS (Flexible Image Transport System) format, which
 need to be converted to txt files in order to use them in this program. Use unfits
 (unfits_09.py) program to convert FITS files into txt files.

       This program (kplrlca) takes a light curve (txt format) with data for all the
 available cycles as a whole, and then extract the data for each cycle and stores it
 temporarily in the memory. All the computations, such as fourier fitting or polynomial
 fitting, calculating delta-m max (measure of the O'Connell Effect), delta-m min, O'Connell
 Effect Ratio(OER), ight Curve Asymmetry(LCA), and errors, are done on cycle-wise basis.
 All these measurements are then saved in a txt file with their kepler id as their name
 (eg. 03321207.txt). For the light curve types OC(Over-Contact), SD(Semi-Detached),
 ELV(Ellipsoidal), and ?(Unknown) fourier fitting is done, whereas for the D(Detached)
 type polynomial fitting (only for the out of eclipse regions) is done. Due to this, for
 D type systems delta-m min, LCA, and OER are not computed.To plot the results from this
 program use keplerPlot.py

       It is intended to be run from the command line in a UNIX-like environment such as
 Linux, Mac OS X, Solaris, or FreeBSD.

       You may use, adapt, and redistribute this software subject to the terms of the
 Creative Commons Attribution-ShareAlike 2.5 Generic License.

 --> See http://creativecommons.org/licenses/by-sa/2.5/ for details.

###############################################################################################

  USAGE: kplrlca [OPTIONS] [FILES]
  
  Available options:

     -h              Print this page.

     -i              Print information about the kepler light curve.

     -p              Display the light curve of the specified cycle for the kepler file.
  
     Ex: kplrlca -i 21
         > Print information about the kepler light curve of S.N. 21.
         > Here is what is prints for S.N. 21.
         #####################################################################################
         ######################## Kepler Light Curve Information #############################
    
         KeplerID : 02305277    Period(days) : 0.378757    Effective Temp(K) : 5694    S.N : 21
    
         T2/T1 : 0.99961    r1 + r2 : N/A    q : 1.11416    FF : 0.88203
    
         [e sin(w), e cos(w)] : [N/A N/A]    sin(i) : 0.42018
    
         LC-type : OC    K-magnitude(mag) : 15.101    No. of cycles : 233.0
    
         #####################################################################################
    
     Ex: kplrlca 21
         > Print information about the kepler light curve of S.N. 1.
         > Display the measuments for each cycle.
         > Here is what is prints for S.N. 21. 
         #####################################################################################
         ######################## Kepler Light Curve Information #############################
    
         KeplerID : 02305277    Period(days) : 0.378757    Effective Temp(K) : 5694    S.N : 21
    
         T2/T1 : 0.99961    r1 + r2 : N/A    q : 1.11416    FF : 0.88203
    
         [e sin(w), e cos(w)] : [N/A N/A]    sin(i) : 0.42018
    
         LC-type : OC    K-magnitude(mag) : 15.101    No. of cycles : 233.0
    
         #####################################################################################
         S.N. = 1   deltammax = -6.8 +/- 0.8   deltammin = -15.3 +/- 0.8   OER = 0.91 +/- 0.03   LCA = 0.83 +/- 0.03
         S.N. = 2   deltammax = -6.1 +/- 0.8   deltammin = -15.2 +/- 0.8   OER = 0.94 +/- 0.03   LCA = 0.82 +/- 0.03
         S.N. = 3   deltammax = -4.4 +/- 0.8   deltammin = -15.0 +/- 0.8   OER = 0.96 +/- 0.03   LCA = 0.79 +/- 0.03
         S.N. = 4   deltammax = -4.7 +/- 0.8   deltammin = -10.8 +/- 0.8   OER = 0.97 +/- 0.03   LCA = 0.78 +/- 0.03
         S.N. = 5   deltammax = -0.8 +/- 0.8   deltammin = -12.6 +/- 0.8   OER = 1.00 +/- 0.03   LCA = 0.76 +/- 0.03
         and so on...
    
     Ex. kplrlca -p 21 3
         > Displays the light curve (cycle number 3) of the file with S.N. 21.
    
