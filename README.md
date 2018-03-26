Directions on how to get correlation functions and powerspectra from halo catalogues

# HOW TO RUN:

  1.) Make the binary - source code is located on scinet at /home/r/rbond/gstein/src/twopoint_codes  or on gitlab to make the binary yourself 

  2.) Run with:
   ./powerspectra 

     usage: powerspectra <pksfile1> <outfile> <Lbox in [Mpc]> <nmesh> <fmt=: 0=peaks, 1=field, 2=peaks x peaks, 3=field x peaks> <Mmin halo> <Nhalo cut> <pksfile2>

     <pksfile1>          file you with to take power spectrum of
     <outfile>           filename to save output
     <Lbox in [Mpc]>    
     <nmesh>             number of cells of FT grid. For fmt=field this must be the same resolution as the field you are trying to read in 
     <fmt>               0=peaks, 1=field, 2=peaks x peaks, 3=field x peaks 

     <Mmin halo>         Minimum mass of halos you want to keep. Only needed if peaks are used
     <Nhalo cut>         Number of halos to use. Assumes the halo catalogue is already ranked in mass.  Only needed if peaks are used
     <pksfile2>          if you want a cross correlation (fmt=2 or fmt=3) then give the location of the second file











#  OLD CODE, only useful for monte carlo runs 
2.) pk_comparison.py
3.) doall_xi_pk.sh


powerspectra is the main fortran code. It reads in a halo catalogue or field, and writes out the power spectrum and correlation function. Usage: powerspectra <mergedfile> <outfile> <Lbox in [Mpc]> <nmesh> <fmt=: 0=peaks, 1=field>

pk_comparison.py is a python wrapper that takes the full halo catalogue, cuts it to only keep the N most massive halos, and saves a copy of the cut version. It then feeds in the appropriate files to ./powerspectra

doall_xi_pk.sh is a bash script that takes input parameters, sets up the proper sub directories, and runs pk_comparison.py with these parameters. EVERYTHING YOU SHOULD NEED TO CHANGE IS IN HERE

HOW TO RUN:
Change "folder" in doall_xi_pk.sh to the location of your halo catalogues and set the rest of the parameters. Then run ./doall_xi_pk.sh from the command line
