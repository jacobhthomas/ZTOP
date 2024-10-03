                      ====== ======= ====== ======
                          /     |    |    | |    |
                         /      |    |    | |    |
                        /       |    |    | |==== 
                       /        |    |    | |
                      ======    |    ====== |

                   Copyright 2001, 2004, 2016 by Zack Sullivan

                                 ZTOP
    NLO differential s-channel and t-channel single-top-quark production

----------------------------------------------------------------------------

ZTOP is a FORTRAN program that calculates the fully differential spin-
averaged production of a top quark plus one jet at next-to-leading order
(i.e. through t + 2 jets).  The program is a coding of the phase space
slicing with two cutoffs calculation presented in the papers below.

Please refer to these papers:

1. B.W. Harris, E. Laenen, L. Phaf, Z. Sullivan, and S. Weinzierl,
Phys. Rev. D 66, 054024 (2002) [arXiv:hep-ph/0207055].

2. Zack Sullivan, Phys. Rev. D 70, 114012 (2004) [arXiv:hep-ph/0408049].

----------------------------------------------------------------------------

Contents:

0. License
1. Compiling
2. s-channel input file details
3. t-channel input file details
4. Various notes

----------------------------------------------------------------------------

0. License:

All code provided here is written and copyrighted by Zack Sullivan (C) 2004,
except where noted in the source codes.  These source codes may be used and
redistributed without permission of the author as long as the original code
is included in an unaltered form.  This code is intended for academic use.
While extensive testing has been done to confirm that there are no bugs
under current FORTRAN standards, the author makes no guarantee of fitness
for any use or correctness.

----------------------------------------------------------------------------

2. Compiling

What you will need/want to edit:

makefile         - choose compiler options
cteqpdffastpdf.f    - set explicit path to CTEQ PDF tables (I've included 5L, 5M1)
                   NOTE: This file is a VERY fast interface to the CTEQ PDF
                    tables, which speeds up PDF evaluation by a factor of 10.  
                    cf. Zack Sullivan, ``Fast evaluation of CTEQ parton
                         distributions in Monte Carlos,'' arXiv:hep-ph/0407033.
(NOTE: LHAPDF version link is set in makefile.)

The included makefile includes good choices of compiler flags for the g77
 and Intel compilers.

To get the executables ZTOPschan.x and ZTOPtchan.x, just type:
make

You should probably run some tests by changing to the examples/ subdirectory,
 and executing the file test.bat.  This may take a few minutes, but you should
 verify that the distributions placed in the subdirectories test/schan/ and
 test/tchan/ agree within round-off errors to the ones in schan/ and tchan/.
 These numbers were generated using g77, but the Intel compiler produces much
 faster code.  These files are produced with about 1/10 the number of points
 used for publication.  For your own evaluations, use the higher statistics -
 the programs are fast enough to handle it.

----------------------------------------------------------------------------

2. s-channel production

Executable: ZTOPschan.x
Sample input files: ls_nc.in, ls_dc.in, ls_jc.in, ls_tc.in
                     s_nc.in,  s_dc.in,  s_jc.in,  s_tc.in

WARNING: There is some error checking, but using parameters other than the
recommended ones (indicated with a *) can lead to peculiar behavior.

Structure of input: (sample input files are documented)

Line 1: Pick order of calculation
* 1 : LO
  2 : NLO correction only
* 3 : Full NLO

Line 2: Base factorization and renormalization scale.  Note there is not an
         option in this release to treat the scales differently.
  0 : PT of top quark
  1 : top quark transverse mass
  2 : PT of highest-PT jet
* 3 : top quark-b jet invariant mass (just like Drell-Yan)
* 4 : top-quark mass (not a great choice for LO)
  5 : s^hat
  6 : W boson mass

Line 3: modifies line 2.
* 0 : Do not modify the scales above
  1 : 2*(line 2)
  2 : (line 2)/2

Lines 4, 6, 8: Set the initial number of iterations for vegas for LO,
                NLO 2 body, and NLO 3 body (use something like 5, 5, 6).
               For graphing, this number is doubled.

Lines 5, 7, 9: Set the initial number of points for vegas for LO,
                NLO 2 body, and NLO 3 body.  A good ratio to use is
                1 : 1.2: 10, e.g. 100000, 120000, 1000000.  The error will
                be dominated by the number of points in line 9.
               For graphing, this number is multiplied by 10.

Line 10: top-quark mass in GeV
175d0    is typical as of June 2004

Line 11: W boson mass in GeV
80.4d0   is typical as of June 2004

Line 12: sqrt(S) in TeV
1.96d0   is the energy of the Tevatron Run II, 1.8d0 was Run I, 14d0 is LHC

Line 13: collider type
* 1,2: p-pbar for Tevatron
* 3  : p-p    for LHC

Line 14: (Original version)
          CTEQ PDF set.  Currently the main function is written to used CTEQ
          PDFs through the PFTOPDG interface.  See accompanying PDF interface
          file for all numbers.
* 53 : CTEQ 5L  (use for LO cross sections)
* 58 : CTEQ 5M1 (use for NLO cross sections)
  64 : CTEQ 6L1 (use for LO cross sections)
  61 : CTEQ 6M  (use for NLO cross sections)
   # : Other choices 41-49 are CTEQ 4, 1-9
                     51-59 are CTEQ 5, 1-9
                     61-64 are CTEQ 6, 1-4
                     6100-6140 are CTEQ 6M,  100-140
                     6200-6240 are CTEQ 61,   00- 40

         (LHAPDF version)
Use the format "pdfname set", e.g. "CT14nlo 0"

Line 15: Calculate t production?  For comparison it is more useful to check
          just t or t_bar, but not both.
* 1 : Yes
  0 : No

Line 16: Calculate t_bar production?  For comparison it is more useful to
          check just t or t_bar, but not both.
* 1 : Yes
  0 : No

Line 17: k_T cone size.  Fixed cones are roughly k_T cones/1.35, e.g. 1d0 is
          like a fixed cone of 0.74.  WARNING: If you use a cone size that is
          different and try to compare, your answer can be off by more than
          all theoretical uncertainties combined.
* 1d0 : Other reasonable choices 0.54d0, 1.35d0 (like fixed cones 0.4, 1.0)

Line 18: k_T cone r_sep.  Currently ignored in the code.  
* Set it equal to line 17 for now.

Line 19: Cuts.  Also defines the suffix for output files.
* 0 : No cuts (nc)
* 1 : Require 2 jets (not including top) passing cuts on lines 26, 27. (dc)
* 2 : Require at least 1 jet (not including top) passing first cut on lines
       26, 27, and any cuts on top quark. (jc)
* 3 : Require exactly 1 jet (not including top) passing first cut on lines
       26, 27, and any cuts on top quark.  Other jets must fail third cuts
       on lines 26, 27. (tc)

Line 20: Smearing
* 0 : none
  1 : Call simple Gaussian smearing routine on jets 80%/sqrt(E)+5%
  2 : 1 + smear top quark (bad choice)

Line 21: Vary soft and collinear cut parameters.  Used to test for numerical
          independence of cuts.
  0 : Vary delta_soft
  1 : Vary delta_collinear
* 2 : Vary delta_collinear = delta_soft/(line 24)

Lines 22-23: Loop over minimum (line 22) to maximum (line 23) exponents.
* 6 : Do not change this unless you understand the code

Line 24: Set collinear/soft factor
* 300d0 : Do not change this unless you understand the code

Line 25: Prefix string (prepended to all output files)
* Put whatever you want.

Line 26: PT_min cuts (GeV) for jet_1, top quark, jet_2.  Cut on jet_2 must be
          <= cut on jet_1 (since jets are PT_ordered).  Top quark cut is
          only applied if eta cut on line 27 is greater than 10^-5.
* Typical values might be 15d0 0d0 15d0

Line 27: Eta_max cuts for jet_1, top quark, jet_2.  Cut on jet_2 should be
          >= cut on jet_1.  Top quark cut is only applied if the second value
          is greater than 10^-5.
* Typical values might be 2.5d0 0d0 2.5d0


----------------------------------------------------------------------------

3. t-channel production

Executable: ZTOPtchan.x
Sample input files: lt_nc.in, lt_dc.in, lt_jc.in, lt_tc.in
                     t_nc.in,  t_dc.in,  t_jc.in,  t_tc.in

WARNING: There is some error checking, but using parameters other than the
recommended ones (indicated with a *) can lead to peculiar behavior.

Structure of input: (sample input files are documented)

Line 1: Pick order of calculation
* 1 : LO
  2 : NLO correction only
* 3 : Full NLO

Line 2: Base factorization scale for light-quark line
Line 4: Base factorization scale for heavy-quark line
Line 6: Base renormalization scale for light-quark line
Line 8: Base renormalization scale for heavy-quark line
  0 : PT of top quark
  1 : top quark transverse mass
  2 : PT of highest-PT jet
  3 : top quark-highest-PT jet invariant mass
* 4 : top-quark mass (not a great choice for LO)
  5 : s^hat
* 6 : Q^2      (best choice for light-quark line, massless DIS scale)
* 7 : Q^2+mt^2 (best choice for heavy-quark line, massive DIS scale)
The rest are nonsensical unless looking for t+2 jets
  8 : M_jj dijet invariant mass
  9 : Q^2+M_jj^2
 10 : MT_jj dijet transverse mass
 11 : Q^2+MT_jj^2 (fac)
 12 : scalar sum of all PT
 13 : scalar sum of all transverse masses

Lines 3, 5, 7, 9 modify lines 2, 4, 6, 8, respectively.
* 0 : Do not modify the scales above
  # : Any other number is assumed to be the fraction of the scale/100,
       e.g. 200 (50) on line 3 would multiply the scale in line 2 by 2 (0.5).

Lines 10, 12, 14: Set the initial number of iterations for vegas for LO,
                    NLO 2 body, and NLO 3 body (use something like 5, 5, 6).
                  For graphing, this number is doubled.

Lines 11, 13, 15: Set the initial number of points for vegas for LO,
                    NLO 2 body, and NLO 3 body.  A good ratio to use is
                    1 : 1: 10, e.g. 100000, 100000, 1000000.  The error will
                    be dominated by the number of points in line 15.
                  For graphing, this number is multiplied by 10.

Line 16: top-quark mass in GeV
175d0    is typical as of June 2004

Line 17: W boson mass in GeV
80.4d0   is typical as of June 2004

Line 18: sqrt(S) in TeV
1.96d0   is the energy of the Tevatron Run II, 1.8d0 was Run I, 14d0 is LHC

Line 19: collider type
* 1,2: p-pbar for Tevatron
* 3  : p-p    for LHC

Line 20: (Original version)
         CTEQ PDF set.  Currently the main function is written to used CTEQ
          PDFs through the PFTOPDG interface.  See accompanying PDF interface
          file for all numbers.
* 53 : CTEQ 5L  (use for LO cross sections)
* 58 : CTEQ 5M1 (use for NLO cross sections)
  64 : CTEQ 6L1 (use for LO cross sections)
  61 : CTEQ 6M  (use for NLO cross sections)
   # : Other choices 41-49 are CTEQ 4, 1-9
                     51-59 are CTEQ 5, 1-9
                     61-64 are CTEQ 6, 1-4
                     6100-6140 are CTEQ 6M,  100-140
                     6200-6240 are CTEQ 61,   00- 40

         (LHAPDF version)
Use the format "pdfname set", e.g. "CT14nlo 0"

Line 21: Calculate t production?  For comparison it is more useful to check
          just t or t_bar, but not both.
* 1 : Yes
  0 : No

Line 22: Calculate t_bar production?  For comparison it is more useful to
          check just t or t_bar, but not both.
* 1 : Yes
  0 : No

Line 23: k_T cone size.  Fixed cones are roughly k_T cones/1.35, e.g. 1d0 is
          like a fixed cone of 0.74.  WARNING: If you use a cone size that is
          different and try to compare, your answer can be off by more than
          all theoretical uncertainties combined.
* 1d0 : Other reasonable choices 0.54d0, 1.35d0 (like fixed cones 0.4, 1.0)

Line 24: k_T cone r_sep.  Currently ignored in the code.  
* Set it equal to line 23 for now.

Line 25: Cuts.  Also defines the suffix for output files.
* 0 : No cuts (nc)
* 1 : Require 2 jets (not including top) passing cuts on lines 32, 33. (dc)
* 2 : Require at least 1 jet (not including top) passing first cut on lines
       32, 33, and any cuts on top quark. (jc)
* 3 : Require exactly 1 jet (not including top) passing first cut on lines
       32, 33, and any cuts on top quark.  Other jets must fail third cuts
       on lines 32, 33. (tc)

Line 26: Smearing
* 0 : none
  1 : Call simple Gaussian smearing routine on jets 80%/sqrt(E)+5%
  2 : 1 + smear top quark (bad choice)

Line 27: Vary soft and collinear cut parameters.  Used to test for numerical
          independence of cuts.
  0 : Vary delta_soft
  1 : Vary delta_collinear
* 2 : Vary delta_collinear = delta_soft/(line 30)

Lines 28-29: Loop over minimum (line 28) to maximum (line 29) exponents.
* 6 : Do not change this unless you understand the code

Line 30: Set collinear/soft factor
* 300d0 : Do not change this unless you understand the code

Line 31: Prefix string (prepended to all output files)
* Put whatever you want.

Line 32: PT_min cuts (GeV) for jet_1, top quark, jet_2.  Cut on jet_2 must be
          <= cut on jet_1 (since jets are PT_ordered).  Top quark cut is
          only applied if eta cut on line 33 is greater than 10^-5.
* Typical values might be 15d0 0d0 15d0

Line 33: Eta_max cuts for jet_1, top quark, jet_2.  Cut on jet_2 should be
          >= cut on jet_1.  Top quark cut is only applied if the second value
          is greater than 10^-5.
* Typical values might be 2.5d0 0d0 2.5d0

----------------------------------------------------------------------------

4. Common structural notes (things you might want to look at):

If you wish to use alternate cut schemes, the logic for cuts is in scuts.f
 and tcuts.f.

If you wish to plot alternate quantities, e.g. delta_phi_jj, delta_eta_jj,
 etc., the histogramming is performed in shist.f and thist.f.  There are two
 things to be careful about:
 1. The graphing routine may be entered with 2 or 3 jets defined, so some
     checking must be done.
 2. s-channel ps(0,3) is a b jet.
    t-channel ps(0,3) or ps(0,5) could be a b jet.  The fractional probability
      is passed in.  See thist.f to understand how to deal with this.

The CKM matrix is defined in the smain.f and tmain.f files.  Currently
 V_td = V_ts = 0, V_tb = 1, with other values commented out.

The programs are fairly modular.  To change how PDFs are set up, edit smain.f
 or tmain.f, and use any set that provides a PFTOPDG interface.

Phase space generation is separated from the momentum-dependent part of the
 matrix element code (up to overall factors).  So even VEGAS could be swapped
 out fairly transparently.

----------------------------------------------------------------------------
Jacob's Guide

Hi learner! Here's a semi-brief introduction to ZTOP, a Fortran program for single-top quark production
By the end of this, you should understand the structure and functionality of ZTOP, making it easier for you to get started with simulations, modifications, and optimizations.



0. Overview

ZTOP calculates the fully differential spin-averaged production of a top quark plus one jet at next-to-leading order (NLO). The program uses phase space slicing with two cutoffs, as detailed in the following papers:

B.W. Harris, E. Laenen, L. Phaf, Z. Sullivan, and S. Weinzierl, Phys. Rev. D 66, 054024 (2002) [arXiv:hep-ph/0207055].
Zack Sullivan, Phys. Rev. D 70, 114012 (2004) [arXiv:hep-ph/0408049].

1. Compiling
In setting up the compilation for the ZTOP program, you will need to edit the following files:

makefile: This is where compiler options are chosen.
  Purpose: Compiles the project and links the necessary object files to create the executables ZTOPschan.x and ZTOPtchan.x.
  Key Components:
    Compiler and Flags: Specifies the Fortran compiler (gfortran or ifort) and compilation flags.
    Object Files: Lists the object files (.o) required for different executables.
    Targets: Defines the targets for the executables and the rules for compiling the source files.
cteqpdflibfast.f: Set the explicit path to CTEQ PDF tables.

While we're here, this is a good time to preview all the source files used in this program:

  vsup.f: Contains subroutines for supporting functions, including Monte Carlo integration using the VEGAS algorithm.
  dilog.f: Implements the dilogarithm function, which is used in various calculations.
  histlib.f: Provides a simple histogramming package for data analysis.
  cteqpdflibfast.f: Contains routines for accessing and interpolating Parton Distribution Functions (PDFs) from the CTEQ collaboration.
  misc.f, lum.f, smearing.f, jetfinder1.f: Additional utility functions and subroutines for various tasks in the simulation.
  smain.f, scuts.f, shist.f, sevent_lo2.f, sevent_nlo2.f, sevent_nlo3.f, ssig_lo2.f, ssig_nlo2.f, ssig_nlo3.f: Source files for simulating single-top quark production in the s-channel.
  tmain.f, tcuts.f, thist.f, tevent_lo2.f, tevent_nlo2.f, tevent_nlo3.f, tsig_lo2.f, tsig_nlo2.f, tsig_nlo3.f: Source files for simulating single-top quark production in the t-channel.

  Recently, changes have been made to use the central PDF values, resulting in the following files and their respective changes:

pftopdgn.f: Interface for the CTEQ PDF routines, forwarding calls to structm and partonx12.
partonx12structmMulti.f: Contains the partonx12 function and other routines for accessing PDF values.
tmainpdf.f: Main program for t-channel single-top-quark production, setting up the simulation and calling various subroutines.
(which are associated with the executables ZTOPtchanPDF.x and ZTOPschanPDF.x)

The included makefile has good choices of compiler flags for the g77 and Intel compilers. To get the executables ZTOPschan.x and ZTOPtchan.x, simply type: make

Running Tests
You should run some tests by changing to the examples/ subdirectory and executing the file test.bat. This will verify that the distributions placed in the subdirectories test/schan/ and test/tchan/ agree within round-off errors to the ones in schan/ and tchan/.

3. s-channel Production

Executable: ZTOPschan.x
Sample Input Files: ls_nc.in, ls_dc.in, ls_jc.in, ls_tc.in, s_nc.in, s_dc.in, s_jc.in, s_tc.in

Input File Structure:

Line 1: Pick order of calculation

  1: LO
  2: NLO correction only
  3: Full NLO

Line 2: Base factorization and renormalization scale

  0: PT of top quark
  1: Top quark transverse mass
  2: PT of highest-PT jet
  3: Top quark-b jet invariant mass
  4: Top-quark mass
  5: s^hat
  6: W boson mass

Line 3: Modifies line 2

  0: Do not modify the scales above
  1: 2*(line 2)
  2: (line 2)/2

Lines 4, 6, 8: Set the initial number of iterations for VEGAS for LO, NLO 2 body, and NLO 3 body.

Lines 5, 7, 9: Set the initial number of points for VEGAS for LO, NLO 2 body, and NLO 3 body.

Line 10: Top-quark mass in GeV (e.g., 175d0).

Line 11: W boson mass in GeV (e.g., 80.4d0).

Line 12: sqrt(S) in TeV (e.g., 1.96d0 for Tevatron Run II).

Line 13: Collider type

  1,2: p-pbar for Tevatron
  3: p-p for LHC

Line 14: CTEQ PDF set

  53: CTEQ 5L (LO)
  58: CTEQ 5M1 (NLO)
  64: CTEQ 6L1 (LO)
  61: CTEQ 6M (NLO)

Line 15: Calculate t production?

  1: Yes
  0: No

Line 16: Calculate t_bar production?

  1: Yes
  0: No

Line 17: k_T cone size (e.g., 1d0).

Line 18: k_T cone r_sep (set equal to line 17).

Line 19: Cuts

  0: No cuts (nc)
  1: Require 2 jets (dc)
  2: Require at least 1 jet (jc)
  3: Require exactly 1 jet (tc)

Line 20: Smearing

  0: None
  1: Gaussian smearing on jets
  2: Smear top quark (not recommended)

Line 21: Vary soft and collinear cut parameters

  0: Vary delta_soft
  1: Vary delta_collinear
  2: Vary delta_collinear = delta_soft/(line 24)

Lines 22-23: Loop over minimum to maximum exponents (e.g., 6).

Line 24: Set collinear/soft factor (e.g., 300d0).

Line 25: Prefix string for output files.

Line 26: PT_min cuts for jet_1, top quark, jet_2 (e.g., 15d0 0d0 15d0).

Line 27: Eta_max cuts for jet_1, top quark, jet_2 (e.g., 2.5d0 0d0 2.5d0).

4. t-channel Production

Executable: ZTOPtchan.x
Sample Input Files: lt_nc.in, lt_dc.in, lt_jc.in, lt_tc.in, t_nc.in, t_dc.in, t_jc.in, t_tc.in

Input File Structure:

Line 1: Pick order of calculation

  1: LO
  2: NLO correction only
  3: Full NLO

Lines 2, 4, 6, 8: Base factorization and renormalization scales for light-quark and heavy-quark lines.

Lines 3, 5, 7, 9: Modify lines 2, 4, 6, 8.

Lines 10, 12, 14: Set the initial number of iterations for VEGAS for LO, NLO 2 body, and NLO 3 body.

Lines 11, 13, 15: Set the initial number of points for VEGAS for LO, NLO 2 body, and NLO 3 body.

Line 16: Top-quark mass in GeV (e.g., 175d0).

Line 17: W boson mass in GeV (e.g., 80.4d0).

Line 18: sqrt(S) in TeV (e.g., 1.96d0 for Tevatron Run II).

Line 19: Collider type

  1,2: p-pbar for Tevatron
  3: p-p for LHC

Line 20: CTEQ PDF set

  53: CTEQ 5L (LO)
  58: CTEQ 5M1 (NLO)
  64: CTEQ 6L1 (LO)
  61: CTEQ 6M (NLO)

Line 21: Calculate t production?

  1: Yes
  0: No

Line 22: Calculate t_bar production?

  1: Yes
  0: No

Line 23: k_T cone size (e.g., 1d0).

Line 24: k_T cone r_sep (set equal to line 23).

Line 25: Cuts

  0: No cuts (nc)
  1: Require 2 jets (dc)
  2: Require at least 1 jet (jc)
  3: Require exactly 1 jet (tc)

Line 26: Smearing

  0: None
  1: Gaussian smearing on jets
  2: Smear top quark (not recommended)

Line 27: Vary soft and collinear cut parameters

  0: Vary delta_soft
  1: Vary delta_collinear
  2: Vary delta_collinear = delta_soft/(line 30)

Lines 28-29: Loop over minimum to maximum exponents (e.g., 6).

Line 30: Set collinear/soft factor (e.g., 300d0).

Line 31: Prefix string for output files.

Line 32: PT_min cuts for jet_1, top quark, jet_2 (e.g., 15d0 0d0 15d0).

Line 33: Eta_max cuts for jet_1, top quark, jet_2 (e.g., 2.5d0 0d0 2.5d0).

5. Various Notes
  Cut Schemes: The logic for cuts is in scuts.f and tcuts.f.
  Histogramming: Performed in shist.f and thist.f.
  CKM Matrix: Defined in smain.f and tmain.f.
  PDF Setup: Edit smain.f or tmain.f to use any set that provides a PFTOPDG interface.
  Phase Space Generation: Separated from the momentum-dependent part of the matrix element code.

Summary
  Compilation: Managed by the Makefile.
  Initialization: Main programs (smain.f, tmain.f) set up the simulation.
  PDF Access: Handled by pftopdgn.f and partonx12structmMulti.f.
  Monte Carlo Integration: Performed by vsup.f.
  Event Generation: Implemented in files like tevent_lo2.f.
  Histogramming: Managed by histlib.f.
  Output: Controlled by the main programs and input files.
This structure allows for a modular and organized approach to simulating particle interactions and analyzing the results. Happy simulating!
