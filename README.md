# SMPE
Two spectral methods for solving range-independent parabolic equation model in ocean acoustics

**`SMPE_readme`, Jul. 20, 2021, Houwang Tu, NUDT**

The programs `CSMPE.m` and `CCMPE.m` compute the range-independent acoustic field in
Fig.1 using the Chebyshev-Tau spectral method and Chebyshev collocation method, respectively. 
The methods are described in the articles (H. Tu, Y. Wang, X. Ma et al., Applying the 
Chebyshev-Tau spectral method to solve the  parabolic equation model of wide-angle rational 
approximation in ocean acoustics, Journal of Theoretical and Computational Acoustics, 2021, 
https://doi.org/10.1142/S2591728521500134 and Y. Wang, H. Tu, W. Liu et al., Application of 
a Chebyshev Collocation Method to Solve a Parabolic Equation Model of Underwater Acoustic 
Propagation, Acoustics Australia, 2021, https://doi.org/10.1007/s40857-021-00218-5).
Both of the programs use the same input file "`input_SMPE.txt`", '`ReadEnvParameter`' function 
is used to read "`input_SMPE.txt`" file. User can make changes to "`input_SMPE.txt`" for the
desired calculation. 

The "`input_SMPE.txt`" file contains the parameters defining the calculation. 
See the following example:

```
example2                                 ! casename
50                                       ! N (truncation order of spectral methods)
8                                        ! np (term of pade approximation)
20                                       ! freq (frequency of source)
40                                       ! zs (depth of source)
40                                       ! zr (depth of special receiver)                                                                              
2500.0                                   ! rmax (receiver ranges(m))
5                                        ! dr (discrete step in horizontal direction)                                      
150                                      ! H (depth of ocean)
0.5                                      ! dz (discrete step in depth direction)
40                                       ! tlmin (minimum value of TL in colorbar)
70                                       ! tlmax (maximum value of TL in colorbar)
30                                       ! n (profiles' points in ocean)
0.00  1480.65   1.0   0.0                ! dep (depth m) c (sound speed m/s) rho (density gm/cm$^3$) alpha (attenuation dB/wavelength)
5.79  1480.80   1.0   0.0
12.06 1481.01   1.0   0.0
19.78 1480.94   1.0   0.0
27.01 1481.16   1.0   0.0
35.21 1481.30   1.0   0.0
37.62 1480.36   1.0   0.0
40.51 1479.06   1.0   0.0
45.82 1475.15   1.0   0.0
47.75 1473.70   1.0   0.0
49.68 1472.68   1.0   0.0
53.54 1471.88   1.0   0.0
56.91 1470.80   1.0   0.0
60.29 1469.78   1.0   0.0
65.60 1468.91   1.0   0.0
72.35 1468.91   1.0   0.0
78.62 1468.62   1.0   0.0
84.41 1467.61   1.0   0.0
89.23 1467.97   1.0   0.0
93.09 1468.84   1.0   0.0
96.95 1468.26   1.0   0.0
99.84 1467.54   1.0   0.0
104.18 1466.81  1.0   0.0
109.49 1466.23  1.0   0.0
113.83 1466.45  1.0   0.0
119.61 1466.59  1.0   0.0
125.40 1466.67  1.0   0.0
130.71 1466.74  1.0   0.0
138.42 1466.67  1.0   0.0
150.00 1466.88  1.0   0.0

```

The "`input_SMPE.txt`" file include:

*  `casename` is the name of current example,

* `N` (the number to truncated
  order of the spectral methods),  

  Generally speaking, the more complicated the shape of the sound speed profile, 
  the more `N` and is needed to accurately fit.

* `np` is the number of items used for the rational approximation; the choice will also affect
   the accuracy of the approximation. More precisely, it is determined by the user according
   to the complexity of the research environment and the characteristics of the sound source. 

* `freq` (frequency of sound source, Hz), 

* `zs` (the depth of source, m), 

* `zr` (depth of a special
  receiver, user used to specify to draw the transmission loss curve of
  arbitrary depth, m), 

* `rmax` (the maximum range of horizontal direction, m), 

* `dr` (horizontal discrete step, m),

*  `H` (depth of ocean, m),

* `dz` (discrete step size in depth direction, m),

* `tlmin`
  and `tlmax` are the minmum and maximum value transmission loss,
  respectively, which used to determine the color range of the output
  transmission loss graph, `tlmin` must less than `tlmax`.

* `n` is the amount of environmental profile data in ocean. 

  There is a table of environmental parameter: the units are depth(m), speed(m/s),
  density(gm/cm$^3$) and attenuation (dB/wavelength)
