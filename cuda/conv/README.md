# CONV

## NAME

gen, ana.py, interp.py, fig.py – grid and convergence study of cns3d and ns2d

## SYNOPSIS

**gen** *case factor*[,*factor*[,*factor*]] *dir*

**gen tbendg** *factor dir inlet* [ *streamwise-factor* ]

**gen tbendd** *factor dir inlet bend-factor*

**gen cavityb** *factor dir nodes*

**python3 ana.py** [ **-s** *suffix* ] *dir* ...

**python3 ana.py** [ **-s** *suffix* ] **-d** *ref dir* ...

**python3 interp.py** [ **-x** ] *coarse fine*

**python3 fig.py** *bends steps*

## DESCRIPTION

*Gen* writes *dir*/control and *dir*/field for one case, with the base grid refined by *factor* in each direction: one factor for all directions, or one for i, one for j and one for k (in two-dimensional cases k is not refined). A direction of *n* points becomes (*n* − 1)·*factor* + 1 points; the base points are kept and the new ones placed by monotone cubic interpolation of the base coordinates over the point index. Refined grids are nested.

Cases of the report:

    cavity  driven cavity, Re 400, lid u = 1 at j = m, uniform 31 × 31
    cavityb cavity on the report's nonorthogonal grid (Figure 5b), read
            from *nodes* (../../redraw/fig5b_nodes.csv); refined over
            the point index along each grid line
    bend    laminar 90-deg bend, square duct, Re 790, 20 × 17 × 9
    lstep   laminar step, 1:2 expansion, Re 300 in control, 45 × 45
    tstep   turbulent step, 2:3 expansion, Re 45000, 45 × 43
    tstepw  tstep with the first spacing at each wall kept
    tstepf  tstep on the grid and inlet that reproduce the report's
            Figures 13 to 15 (../../TESTS, case C), 45 × 42
    tstepfw tstepf with the first spacing at each wall kept
    tstepfp tstepf with the inlet u = 8/7 (distance to the nearer
            wall)^(1/7)

The cavities are run with flag 16 (../README.md). The bend has inner radius 1.8, outer 2.8, four stations before the bend, ten in it and six after, and a developed laminar inlet profile. The laminar step has a 30h channel, x points stretched as exp(0.5 s), uniform y and a parabolic inlet over the upper half; Re 300 on the inlet height is the report's Re 600 on twice that height. The turbulent step has 9 points over the 4h before the step, 36 over 30h after it, stretched as exp(2.5 s), uniform y, a uniform inlet, inlet k from 4 percent intensity and length scale 0.14. *Tstepf* has 3 points before the step, 42 over 15.4h after it, stretched as exp(2.698 s), y spacing reduced towards the walls and the step corner, inlet k from 2.97 percent intensity and length scale 3.741; it is run with hybrid factor 0.792 and flag 8 (no wall function at the step corner). Relaxation is 0.5, except 0.3 for the pressure of *tstep* and 0.7 for the viscosity of the turbulent cases. A grid marked w below is a *tstepw* grid.

A check case, not in the report, with no measurements to compare with:

    tbend   the bend turbulent, Re 40000, initial k 0.005, ε 0.001162
    tbendw  tbend with the first spacing at each wall kept
    tduct   straight duct, 60 long, 121 points, at 45 deg to x, the
            cross-section of tbendg, uniform inlet, k 0.005, ε 0.001162
    tbendg  tbend with a graded cross-section and the inlet profile, k
            and ε taken from the last interior section of the tduct
            result inlet on the same factor
    tbendd  tbendg with designed streamwise points

The cross-section of *tduct*, *tbendg* and *tbendd* is graded from each wall with first spacing 0.024 and ratio at most 1.2^(1/*factor*): 20, 26 and 30 cells across and 10, 13 and 15 in the half height for factors 1, 2 and 4. The streamwise points of *tbendg* are refined by *streamwise-factor*, default *factor*. *Tbendd* has 9·*bend-factor* equal intervals over the bend; in the legs, 2.74 before and 1.93 after it, the spacing starts from that of the bend and grows by 1.1 per point up to 0.3. The duct is at 45 deg because the inflow is the flux of u and the outflow that of v (../README.md).

*Ana.py* prints for each result *dir*/field*suffix* (default .o):

- step: *xr*, where u on the first row above the lower wall, plane k = 2, changes from negative to non-negative downstream of the step, interpolated linearly;
- bend: *dp*, the area-weighted mean pressure of the section 1.03 before the bend less that 1.09 after it; *db*, the same from the bend's 0-deg section to its 90-deg section; *ur* and *w*, the largest |u| and |w| on the 90-deg section (u is radial there). Stations are found by centreline arc length. For laminar bends, the RMS difference of the streamwise velocity from the digitized curves of the report's Figure 19 (present), Humphrey's measurements and Rhie's results at 0, 60 and 90 deg, z/d 0 and 0.245 (../../redraw/fig19*). For turbulent bends, the median y+ of the points next to the walls.

*Interp.py* replaces the initial field of *fine*/field, except the inlet section i = 1, by the result *coarse*/field.o interpolated linearly: over the point index when *fine* is a refinement of *coarse*, or with **-x** over the x and y coordinates of a step grid.

*Fig.py* draws, from the bend results in *bends* and the others in *steps*: fig7.svg, the cavity centreline velocity on the report's grids and on the finest grids over the report's Figure 7 and Burggraf's results; fig10.svg, xr of the laminar step against Re, hybrid and central (hybrid factor 0.1) differencing, on the 45×45 and 1409×353 grids, converged runs only, over the report's Figure 10 (../../redraw/fig10*); fig19.svg, the laminar bend profiles of the 20×17×9 and 153×129×65 grids over the report's Figure 19 and Humphrey's measurements; step.svg, xr of the laminar step against the x spacing.

With **-d** it prints the RMS difference of each of the six profiles of each *dir* from those of *ref*, on 49 points across the duct. The profiles are named by angle and plane: 60a is 60 deg, z/d 0; 60b is 60 deg, z/d 0.245.

## RESULTS

Bend results are from *cns3d*, single precision, TOL 1.E-5 on the sum of the changes (../README.md), hybrid factor 0.5, unless given. Step and cavity results are from *ns2d* with multigrid, TOL 1.E-6 on the largest residual ratio, an iteration limit of 10000, double precision, or single on grids wider than 1036 points, hybrid factor 0.5, unless given. Every *ns2d* run starts from the field *gen* writes, unless given. A run is converged when the program printed done; (*) marks a run stopped at its iteration limit, which the program prints as not converged.

### Laminar bend

    grid          iterations   dp      db      ur     w      present  humphrey
    20×17×9        3370        0.3828  0.1900  0.217  0.165  0.146    0.281
    39×33×17       1070        0.3991  0.2265  0.245  0.219  0.127    0.241
    58×49×25       1640        0.3962  0.2331  0.297  0.255  0.134    0.230
    77×65×33       2330        0.3940  0.2356  0.339  0.279  0.140    0.227
    115×97×49      4060        0.3925  0.2379  0.385  0.301  0.147    0.224
    153×129×65     5800        0.3917  0.2384  0.399  0.310  0.151    0.224

    RMS from 153×129×65   0a     0b     60a    60b    90a    90b
    77×65×33              0.003  0.003  0.055  0.062  0.051  0.063
    115×97×49             0.001  0.001  0.021  0.025  0.020  0.022

Fits X = X∞ − C h^p to the 4×, 6× and 8× grids: db p = 3.31, db∞ = 0.2387; ur p = 2.39, ur∞ = 0.413; w p = 1.56, w∞ = 0.326. The report's own Figure 19 curves differ from Humphrey's measurements by RMS 0.258.

One direction refined 4 times:

    grid          refined      iterations   dp      db      ur     w      humphrey
    77×17×9       streamwise    1110        0.3475  0.1764  0.180  0.143  0.253
    20×65×9       radial        6230        0.3977  0.2160  0.219  0.165  0.259
    20×17×33      z             1350        0.4195  0.2233  0.276  0.234  0.278

### Laminar step

    grid       relax   iterations    xr
    45×45      0.5       4700        8.127
    89×89      0.5       1884        9.074
    133×133    0.5       1944        9.692
    133×133    0.3      10000 (*)    9.691
    177×177    0.5       1664       10.001
    177×177    0.3       3812       10.000

(*) the residual ratio stays near 9.E-4. A fit to the 89, 133 and 177 grids at relaxation 0.5 gives p = 1.00, xr∞ = 10.93.

x and y refined apart, relaxation 0.3; rows are the x factor, columns the y factor:

    xr      y 1      y 2      y 4      y 8
    x 1     8.127             8.452    8.492
    x 4     9.659    9.929   10.000
    x 8    10.071            10.405
    x 16   10.270            10.596   10.635
    x 32                     10.679   10.687
    x 40                     10.686

The x 1, y 1 result is the 45×45 run at relaxation 0.5. Iterations: 4700, 4872, 5188; 3216, 3348, 3812; 2872, 3232; 2988, 3148, 3312; 3812, 2548; 2548. The x 32 and x 40 runs are single precision; their residual ratios stop near 6.E-7. The 1409×353 and 1761×177 runs start from the 1409×177 and 353×177 results interpolated (*interp.py*). The report gives 8.44 with the hybrid scheme and 10.53 with central differencing; Armaly's measurement is 11.45 (../../redraw/fig10*).

Against Re, relaxation 0.5 on 45×45 and 0.3 on 1409×353, hybrid factor 0.5 (hybrid) and 0.1 (central); the 1409×353 runs start from the x 32, y 8 result above, which is the hybrid Re 600 run:

    Re     hybrid 45×45     central 45×45     hybrid 1409×353    central 1409×353
    100    3.102    920     3.122    920      3.215   1268       3.215   1176
    150    4.123    930     4.188   1150      4.318   1360       4.319   1200
    200    5.047    980     5.157   1240      5.329   1900       5.329   1744
    300    6.688   1360     6.846   1320      7.131   2336       7.129   2100
    400    7.922   2120     8.194   1900      8.637   2576       8.637   2276
    500    8.458   3540     9.117   2580      9.817   2760       9.826   2420
    600    8.127   4700     9.743   2990     10.689    148      10.724   1000
    700    7.377   4320    10.226   3340     11.421   3308      11.484   2928
    800    6.873   3920    10.603   3500     12.077   3728      12.171   3508

On 1409×353 the two schemes differ by at most 0.009 up to Re 500 and by 0.094 at Re 800. Figure 10 of the report gives, at Re 800, about 7.7 hybrid and 12.8 central.

### Driven cavity

Re 400, relaxation 0.5. umin is the smallest u on the vertical centreline and y its height:

    grid        iterations    umin      y
    31×31           530      -0.2376   0.300
    61×61          1430      -0.2992   0.285
    121×121         620      -0.3188   0.283
    241×241         700      -0.3242   0.283
    31×31 b         500      -0.2198   0.295
    61×61 b        1380      -0.2882   0.295
    121×121 b       608      -0.3138   0.285
    241×241 b       660      -0.3223   0.285

b marks the report's nonorthogonal grid (*cavityb*). Its digitized boundary nodes lie up to 0.001 off the walls; *gen* moves them onto the walls before refining. A fit to the 61, 121 and 241 grids gives p = 1.86, umin∞ = −0.3263, and on the b grids p = 1.59, umin∞ = −0.3265. The report's Figure 7 gives about −0.25 at y 0.30 on both grids; Burggraf gives −0.26.

### Turbulent step

    grid          relax     iterations     xr
    45×43         0.5/0.3      760        4.302
    89×85         0.5/0.3     1616        5.203
    133×127       0.5/0.3     3284        5.444
    133×127       0.3/0.2     7604        5.444
    177×169       0.5/0.3     5560        5.542
    177×169       0.3/0.2    10000 (*)    5.542
    88×81 w       0.5/0.3     1488        5.178
    174×157 w     0.3/0.2    10000 (*)    5.764

(*) the residual ratio is 7.E-6 (177×169) and 3.E-6 (174×157 w). A fit to the 89, 133 and 177 grids at relaxation 0.5/0.3 gives p = 1.58, xr∞ = 5.71.

x and y refined apart, relaxation 0.3/0.2:

    xr      y 1          y 2          y 4          y 8
    x 1     4.302                     4.186 (*)
    x 4     5.267        5.507        5.542 (*)
    x 8     5.285                     5.679 (*)    5.920 (*)
    x 16    5.290                     5.734 (*)    6.013 (*)
    x 32                              5.827 (*)
    x 16 w  5.399 (*)    5.502 (*)    5.771 (*)    6.130 (*)
    x 32 w                            5.492 (*)

The x 1, y 1 result is the 45×43 run at relaxation 0.5/0.3. Iterations of the converged runs: 1264, 3756; 2004; 7884. (*) stopped at 10000 iterations with the residual ratio between 6.E-6 and 4.E-5, except x 16 w, y 8 (2.E-4) and x 32 w (5.E-4, single precision). The x 32, y 8 run was not made: it starts from the x 32, y 4 result.

*Tstepf* 45×42 at relaxation 0.862/0.556 converges in 290 iterations: xr 5.201. The report gives 5.4h. Chen and Kim (NASA CR-179204, 1987), with second-order upwind differencing on a 51×41 grid, grid independent within 1 percent, give 5.2h for the standard k–ε model on this step. The measurement is about 7h. *Tstepf* 45×42 at relaxation 0.3/0.2 and 0.1/0.05:

    hybrid factor   flags   xr              iterations
    0.5             0       4.987  4.987    2460  10000 (*)
    0.5             8       5.269  5.269    2400  10000 (*)
    0.792           0       4.943  4.943    2590  10000 (*)
    0.792           8       5.201  5.201    2540  10000 (*)

(*) the residual ratio is 1.0–1.5E-5. The relaxation does not change xr.

*Tstepf* refined, relaxation 0.3/0.2: 89×83 converges in 6904 iterations, xr 5.606; 177×165 and 353×329 stop at 10000 with residual ratios 7.E-5 and 4.E-4, xr 5.735 and 6.090; 705×657 stops at 5.E-4 with no reattachment found.

### Turbulent bend (check case)

    grid          iterations    dp      db      ur     w      y+
    20×17×9        2430         0.0827  0.0989  0.280  0.194  73.7
    39×33×17       1530         0.0204  0.0515  0.258  0.195  30.6
    58×49×25       2480         0.0053  0.0399  0.250  0.188  17.5
    77×65×33       4010         0.0028  0.0376  0.247  0.186  12.6
    39×31×16 w     20000 (*)    0.0264  0.0488  0.266  0.199  61.0
    39×31×16 w     40000 (*)    0.0264  0.0488  0.266  0.199  61.0

    RMS from 77×65×33   0a     0b     60a    60b    90a    90b
    20×17×9             0.073  0.064  0.067  0.073  0.068  0.074
    39×33×17            0.023  0.016  0.022  0.018  0.010  0.008
    58×49×25            0.011  0.007  0.010  0.008  0.006  0.003

(*) relaxation 0.3, the largest change per iteration near 1.E-3. The 58×45×23 and 77×59×30 *tbendw* grids, with relaxation 0.3, were stopped at 16420 and 8050 iterations with the largest change near 1.E-3. The wall function switches to its laminar form below y+ 11.63.

Graded cross-section, inlet from the duct:

    duct     iterations   dp/ds    mean speed
    21×11    3270         0.0146   0.989
    27×14    4330         0.0139   0.989
    31×16    5170         0.0136   0.989

dp/ds is the gradient of the area-weighted mean pressure between 40 and 59 along the duct; the Blasius friction factor gives 0.0112. The mean speed is the area-weighted streamwise speed at 50. The streamwise speed at 59.5, the section *tbendg* takes, differs from that at 49.5 by at most 0.0065 (factor 1) and 0.0035 (factor 2).

    grid          factor  streamwise  iterations   dp      db      ur     w      y+
    20×21×11      1       1           1580         0.1848  0.1515  0.180  0.075  39.7
    39×21×11      1       2           1120         0.1328  0.1037  0.167  0.085  39.6
    77×21×11      1       4           1760         0.1051  0.0785  0.163  0.090  39.4
    153×21×11     1       8           3110         0.0913  0.0661  0.162  0.092  39.4
    305×21×11     1       16          6080         0.0844  0.0599  0.163  0.093  39.4
    609×21×11     1       32         12110         0.0813  0.0571  0.163  0.094  39.4
    39×27×14      2       2           1160         0.1312  0.1040  0.166  0.081  35.0
    305×27×14     2       16          7000         0.0799  0.0564  0.161  0.088  34.8
    609×27×14     2       32         13280         0.0758  0.0526  0.161  0.088  34.8
    20×31×16      4       1           1430         0.1853  0.1557  0.178  0.070  33.0
    77×31×16      4       4           2360         0.1008  0.0762  0.160  0.083  32.8
    305×31×16     4       16          7730         0.0789  0.0557  0.159  0.086  32.8

Designed streamwise points (*tbendd*):

    grid          factor  bend-factor  iterations   dp      db      ur     w      y+
    194×21×11     1       16            4140        0.0838  0.0599  0.163  0.092  39.0
    351×21×11     1       32            7580        0.0796  0.0563  0.164  0.093  38.9
    653×21×11     1       64           13970        0.0798  0.0567  0.163  0.092  38.9
    351×27×14     2       32            8320        0.0770  0.0541  0.162  0.088  34.6
    653×27×14     2       64           15030        0.0741  0.0514  0.161  0.087  34.5
    351×31×16     4       32            8880        0.0740  0.0516  0.160  0.085  32.7
    653×31×16     4       64           15170        0.0721  0.0499  0.159  0.085  32.7

    RMS from 653×31×16   0a     0b     60a    60b    90a    90b
    351×31×16            0.000  0.000  0.001  0.001  0.001  0.001
    653×27×14            0.002  0.003  0.003  0.003  0.003  0.003

On 351×21×11 the hybrid factor 0.25 converges in 7570 iterations: dp 0.0789, db 0.0554, ur 0.163, w 0.096. Factors 0.1, 0.03, 0.01 and 0.001 diverge in 10 iterations.

### Tolerance, base grids

TOL is the sum of the changes for the bends (*cns3d*) and the largest residual ratio for the steps (*ns2d*); the steps run in single precision at TOL 1.E-5 and above, double below.

    TOL      bend dp  iterations   tbend dp  iterations   lstep xr  iterations   tstep xr  iterations
    1.E-3    0.3806      220       0.0858       200        8.718      1180         4.303       430
    1.E-4    0.3830      520       0.0828       620        8.186      2410         4.301       540
    1.E-5    0.3828     3370       0.0827      2430        8.132      3560         4.302       610
    1.E-6    0.3830  1539580       0.0829     27490        8.127      4700         4.302       760
    1.E-7    not reached by the bends: bend stopped at 1785090, tbend at 1388350; lstep 8.126 in 5870, tstep 4.302 in 870

Over the same TOL ur, w and the profile RMS of *bend* and *tbend* change by at most 0.001. The 153×129×65 bend result, restarted in double precision at TOL 1.E-7, converges in 5920 more iterations: dp 0.3918, db 0.2385, ur 0.399, w 0.310, the six profiles within 0.0005 of the TOL 1.E-5 result.

## EXAMPLES

Build, make the 2× bend and run it:

    cc -O2 -o gen gen.c -lm
    gen bend 2 b2
    echo "$PWD/b2 0.5" > list
    ITER=500000 cns3d 40 34 18 list .o
    python3 ana.py b2

The laminar step refined 16 times in x and 4 in y:

    gen lstep 16,4 s
    echo "$PWD/s 0.5 0.3 0.3" > list
    TOL=1e-6 ns2d 706 178 list .o
    python3 ana.py s

Profiles of two grids against a third:

    python3 ana.py -d b4 b1 b2

The turbulent bend, designed streamwise points, 32 intervals per base bend interval:

    gen tduct 1 d1
    echo "$PWD/d1 0.5" > list
    ITER=500000 cns3d 122 22 12 list .o
    gen tbendd 1 g d1/field.o 32
    echo "$PWD/g 0.5" > list
    ITER=500000 cns3d 352 22 12 list .o
    python3 ana.py g

## FIGURES

fig7.svg, fig10.svg, fig19.svg, step.svg (*fig.py*), on a white background.

## SOURCE

`gen.c`, `ana.py`, `interp.py`, `fig.py`

## BUGS

Under a CUDA MPS server a fatal error in one run stopped all the runs on that server.

In *ns2d* k and ε are solved on the finest grid only, so turbulent runs gain less from multigrid than laminar ones: most refined turbulent steps stop at 10000 iterations a few times above TOL. *Tstepf* refined in y, where the first node lies below y+ 5, stops far above TOL.

*Lstep* 133×133 at relaxation 0.3 stops with the residual ratio near 9.E-4, while relaxation 0.5 converges.

At TOL 1.E-5 *ns2d* stops too early on fine grids: on 1409×177 xr is 10.754 at 1.E-5 and 10.679 at 1.E-6, and umin of the 241×241 cavity −0.3201 and −0.3242.

The Humphrey data of ../../redraw/fig19* hold stray points, at 0 deg (0.43, 1.23) and at 60 deg, z/d 0 four points near U/U_ref 0; the RMS values include them.
