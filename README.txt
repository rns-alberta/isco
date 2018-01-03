ISCO

This code computes the structure of rapidly rotating neutron star and outputs the orbital frequency for a particle orbiting at the ISCO.

1. INPUT: Provide the code with an equation of state file, a range of central energy densities and a spin frequency.

2. OUTPUT: Mass, Radius, and the orbital frequencies at the co-rotating and counter-rotating ISCOs. (Generally only the co-rotating ISCO is of interest)

3. COMPILATION: At the command line type "make isco" to create the executable isco
> make isco
The warnings produced should not be problematic. It's possible that you may have to change the name of your compiler if you are not using the gnu compiler. 

4. RUN: To run, you will need to choose an equation of state and decide on a useful range of central energy densities. The code RNS is useful for this task. There is also a code called "maxmass" in this bundle that could also be useful.

EOS FILES: I keep my equation of state files in a directory called "eos" that has the relative path given below. The string "../eos/eosAPR" could be replace with just "eosAPR" if your eos files are in the same directory as the executable. 

EXECUTION: To create a series of 21 models using eos APR spinning at 300 Hz type at the command line (or run the simple bash script go.sh)

> isco -f "../eos/eosAPR" -e 8e14 -l 22e14 -n 20 -s 300
This command should take about 2-3 minutes to run. 

Sample Output:
spin=300
../eos/eosAPR,  MDIVxSDIV=151x301
e_c 	 Mass 	 Radius	 Spin 	 Freq+ 	 Freq- 
e15 	 Msun 	 km    	 Hz   	 Hz    	 Hz 
0.8 	 1.067 	 11.53 	 300.0 	 1523.4 	-1550.5
0.87 	 1.211 	 11.52 	 300.1 	 1623.6 	-1580.3
0.94 	 1.346 	 11.50 	 300.1 	 1713.3 	-1442.0
1.01 	 1.469 	 11.48 	 300.0 	 1643.0 	-1333.9
1.08 	 1.577 	 11.44 	 299.9 	 1530.1 	-1251.4
1.15 	 1.673 	 11.40 	 300.0 	 1441.6 	-1186.2
1.22 	 1.759 	 11.35 	 299.9 	 1369.7 	-1133.2
1.29 	 1.834 	 11.30 	 300.1 	 1312.3 	-1090.6
1.36 	 1.899 	 11.24 	 300.1 	 1265.9 	-1056.2
1.43 	 1.955 	 11.17 	 299.9 	 1228.1 	-1028.3
1.5 	 2.004 	 11.10 	 299.9 	 1197.4 	-1005.7
1.57 	 2.045 	 11.04 	 299.9 	 1172.0 	-987.1
1.64 	 2.080 	 10.97 	 299.9 	 1151.1 	-972.0
1.71 	 2.110 	 10.90 	 300.1 	 1134.0 	-959.6
1.78 	 2.134 	 10.82 	 300.0 	 1119.8 	-949.6
1.85 	 2.155 	 10.75 	 299.9 	 1107.9 	-941.4
1.92 	 2.173 	 10.69 	 299.9 	 1098.2 	-934.9
1.99 	 2.187 	 10.62 	 299.8 	 1090.2 	-929.6
2.06 	 2.199 	 10.55 	 300.1 	 1083.7 	-925.4
2.13 	 2.208 	 10.49 	 300.1 	 1078.4 	-922.2

NOTES:
(a) To make changes to the format of the output, edit the file isco.c
(b) This code won't compute zero spin stars. But if you ask it to compute stars spinning at 1 Hz, that will be close enough to zero. In this case the co- and counter- rotating isco should be equal (with opposite sign).



