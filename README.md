#CMCsubPack
CMCsubPack is collection of Fortran subroutines aimed at providing submodels for the Probability Density Function (PDF), the Conditional Velocity (CV) and the Conditional Scalar Dissipation Rate (CSDR), all of which appear unclosed in the  Conditional Moment Closure (CMC) turbulent combustion model. CMCsubPack provides CV and CSDR closures that are fully consistent with the moments of the PDF.

##PDF
Three PDFs are implented. Each PDF is accompanied by consistent CV and CSDR submodels. The considered PDFs are:

1. Beta PDF. This PDF is valid for binary (two-srtream) mixing.
2. Binary Presumed Mapping Function (PMF) PDF proposed by Mortensen and Andersson [2]. This PDF is valid for binary mixing.
3. Trinary PMF PDF proposed by Mortensen et al. [3-5]. The closure is valid for trinary (three-stream) mixing.

##CSDR
Homogeneous and inhomogeneous CSDR closures are implemented based on three PDFs:

1. Beta PDF-based closre following the model proposed by Mortensen [1]. The closure is valid for binary mixing.
2. Binary PMF PDF-based closure following the model proposed by Mortensen and Andersson [2]. The closure is valid for binary mixing.
3. Trinary PMF PDF-based closure following the model proposed by Mortensen et al. [3-5] . The closure is valid for trinary mixing.

##CV
The CV is modelled using the PDF gradient diffusion model of Pope [6] based on the three presumed PDFs:

1. Beta PDF-based cloure. The closure is valid for binary mixing.
2. Binary PMF PDF-based following the model proposed by Mortensen and Andersson [2]. The closure is valid for binary mixing.
3. Trinary PMF PDF approach following the model proposed by Mortensen [3]. The closure is valid for trinary mixing.

##Other submodels
Additionally, CMCsubPack provides the following commonly used CMC submodels:

1. The delta PDF
2. The clipped Gaussian PDF [7].
3. The linear CV model [8] based on the beta PDF.
4. The Amplitude Mapping Closure CSDR model [9] based on the beta and binary PMF PDFs.
5. Girimaji's CSDR model [10] based on the beta PDF.

###References

[1] M. Mortesen Consistent modeling of scalar mixing for presumed, multiple parameter probability density functions. Phys. Fluids, 17(1):018106, 2005.

[2] M. Mortensen and B. Andersson. Presumed mapping functions for eulerian modelling of turbulent mixing. Flow
Turb. Combust., 76(2):199-219, 2006.

[3] M. Mortensen. Mathematical modelling of turbulent reacting ﬂows. PhD thesis, Chalmers University of Technology, Goteborg, Sweden, 2005.

[4] M. Mortensen, S. M. de Bruyn Kops, and C. M. Cha. Direct numerical simulations of the double scalar mixing layer: Part II: Reactive scalars. Combust. Flame, 149(4):392-408, 2007.

[5] C. M. Cha, S. M. de Bruyn Kops, , and M. Mortensen. Direct numerical simulations of the double scalar mixing
layer. Part I: Passive scalar mixing and dissipation. Phys. Fluids, 18(6):067106, 2006.

[6] S. B. Pope. The probability approach to the modelling of turbulent reacting ﬂows. Combust. Flame, 27(3):299-312, 1976.

[7] F. C. Lockwood and A. S. Naguib. The prediction of the fluctuations in the properties of free, round-jet, turbulent, diffusion flames. Combust. Flame, 24:1090-124, 1975.

[8] V. R. Kuznetsov and V. A. Sabel'nikov. Turbulence and combustion. English Edition Editor: P. A. Libby. Hemisphere Publishing Corporation, New York, U.S.A., 1990.

[9] E. E. O'Brian and T. L. Jiang. The conditional dissipation rate of an initially binary scalar in homogeneous turbulence. Phys. Fluids A, 3(12):3121-3123, 1991.

[10] S. S. Girimaji. On the modeling of scalar diffusion in isotropic flows. Phys. Fluids A, 4(11):2529-2537, 1992.
