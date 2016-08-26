Directions on how to get correlation functions and powerspectra from halo catalogues

FILES NEEDED:

1.) ./powerspectra - source code is located at /home/gstein/src/powerspectra to make the binary yourself 
2.) pk_comparison.py
3.) doall_xi_pk.sh


powerspectra is the main fortran code. It reads in a halo catalogue or field, and writes out the power spectrum and correlation function. Usage: powerspectra <mergedfile> <outfile> <Lbox in [Mpc]> <nmesh> <fmt=: 0=peaks, 1=field>

pk_comparison.py is a python wrapper that takes the full halo catalogue, cuts it to only keep the N most massive halos, and saves a copy of the cut version. It then feeds in the appropriate files to ./powerspectra

doall_xi_pk.sh is a bash script that takes input parameters, sets up the proper sub directories, and runs pk_comparison.py with these parameters. EVERYTHING YOU SHOULD NEED TO CHANGE IS IN HERE

HOW TO RUN:
Change "folder" in doall_xi_pk.sh to the location of your halo catalogues and set the rest of the parameters. Then run ./doall_xi_pk.sh from the command line
