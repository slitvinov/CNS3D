# CNS3D(1)

## NAME

cns3d, cns2d, ns2d – many steady incompressible flows at once on a GPU

## SYNOPSIS

**cns3d** *di dj dk list suffix*

**cns2d** *di dj dk list suffix*

**ns2d** *di dj list suffix*

## DESCRIPTION

The programs solve the steady incompressible Navier–Stokes equations on body-fitted grids, laminar or with the k–ε model and wall functions. They use finite differences, the hybrid scheme and SIMPLE-C. All cases of *list* run together.

*Cns3d* is three-dimensional: symmetry at the first z plane, a wall at the last. *Cns2d* is two-dimensional: four z planes (n = 4), symmetry at both ends; it solves w and links the planes. *Ns2d* is two-dimensional with one plane of unknowns and no w: its equations are those of *cns3d* for a flow that does not depend on z. The metrics are those of the plane, and the source of the pressure correction weighs its three families of control volumes as at an interior plane. It reads the field file of *cns2d*, takes plane k = 2, and writes every plane with the same values and w = 0.

The pressure is changed by the pressure correction less its value at a reference node: in *cns3d* node (2, 2, 2); in *cns2d* and *ns2d* the node at i = l − 1 or l − 2, j = (m − 1)/2 (integer part), whichever has the parity of i + j of the node changed. The correction can hold a part that alternates from node to node and moves no velocity; the two references keep it out of the pressure.

    i = 1      inlet, fixed; the inflow is the flux of u through it
    i = l      outlet: copies its neighbour, plus a uniform correction to
               v in cns3d, u in cns2d and ns2d, that makes the flux of that
               component equal the inflow
    j = 1, m   walls
    block      solid, given in the field file

At a wall the diffusion uses the wall node's own viscosity. In turbulent flow a wall function sets the shear on the link to each wall and the production of k next to it. After its equation k is kept positive, since the wall value of ε is proportional to k to the power 1.5.

*Di*, *dj*, *dk* are the array sizes, at least one more than the largest grid. *Di* ≤ 2075, 1037 in double precision; *dk* ≤ 128 in *cns3d*, 16 in *cns2d*.

Each line of *list* is a case:

    dir hyb [ relax relaxp flags c1 ]

Each case reads *dir*/control and *dir*/field and writes *dir*/field*suffix*. The program prints, for each case,

    dir iterations n done | diverged | not converged

The error of a variable is its largest change in an iteration; that of the pressure is its largest correction. A case is done when the errors of u, v, w, k, ε and the pressure add up to 1.E-5 or less (*ns2d*: no w). It has diverged when the sum passes 1.E3 after 20 iterations, or a field is not finite. The errors are tested every ten iterations, so *n* is a multiple of ten or the limit.

The cases must share the grid size, laminar or turbulent, the iteration limit and the number of pressure sweeps.

## OPTIONS

    hyb      hybrid factor in (0, 1]; 0.5 is the hybrid scheme, a smaller
             factor keeps central differencing up to cell Peclet number 1/hyb
    relax    relaxation of u, v, w, k, ε (−1: from control)
    relaxp   relaxation of p (−1: from control)
    flags    sum of
               1  k production from shear only
               2  near-wall k production from the velocity along the wall
               4  wall distance to the wall node
               8  no wall function at the step corner
              16  closed: i = l keeps its field values, as i = 1 does, and
                  j = m is not a wall, so a velocity given there drives the
                  flow (a lid; cns2d and ns2d only)
    c1       C1 of the ε equation (1.43)

## ENVIRONMENT

    ITER   iteration limit, overriding control
    TOL    convergence limit (1.E-5)
    TIMES  print the time of each stage
    VERB   print the errors of the first case every ten iterations:
           u, v, w, p; in ns2d u, v, k, ε, p

## FILES

*Control*, three lines:

    columns 1–5, 6–10, 11–15   iteration limit; 1 laminar, 2 turbulent; pressure sweeps
    7 columns each             relaxation of u, v, w, p, k, ε, viscosity, and one unused
    columns 1–7                Reynolds number; the viscosity is its inverse

*Field*, the grid and the initial field. The first line has nine integers of five columns: *l m n*, the grid size, and *i1 i2 j1 j2 k1 k2*, the solid block (zeros for none). Then one line per node, k outer, i, then j inner:

    x y z      8 columns each
    u v w      12 columns each, then 2 blank
    p t k ε    11 columns each; further columns are ignored

Lines shorter than 100 characters are skipped. The inlet keeps its initial values.

*Field suffix*, the result: the same layout, with the viscosity and the density added; it can start the next run.

## EXAMPLES

Build the programs:

    nvcc -O3 -arch=sm_90 -Xcompiler -fopenmp -o cns3d cns3d.cu
    nvcc -O3 -arch=sm_90 -Xcompiler -fopenmp -o cns2d cns2d.cu
    nvcc -O3 -arch=sm_90 -Xcompiler -fopenmp -o ns2d ns2d.cu

The same with -DDOUBLE computes in double precision.

Run a batch of bends, 20 × 17 × 9 points, and a batch of steps, 45 × 42 × 4:

    cns3d 21 18 10 bends .g
    cns2d 46 46 8 steps .g

The same steps with one plane:

    ns2d 46 46 steps .g

Compare two step results and print their reattachment lengths; `cmp.py` is for laminar results, `cmpk.py` for turbulent ones, both for *cns2d* only:

    python3 cmpk.py a/field.g b/field.g 45 42

A laminar bend of 20 × 17 × 9 points converges in 3380 iterations: 4.5 s for one case, 20 s for 512 (GH200).

## SOURCE

`cns3d.cu`, `cns2d.cu`, `ns2d.cu`

## SEE ALSO

[conv/](conv/README.md): grid and tolerance study

## DIAGNOSTICS

    cns3d: file: reason

on standard error, exit status 1.

## BUGS

*Cns2d* and *ns2d* sweep red and black; in three dimensions that stalls near 3.E-5, so *cns3d* sweeps in order.

In *cns2d* the pressure on planes 2 and 4 can differ from that on plane 3; the difference moves no u or v but drives w, and with the y spacing small next to the z spacing w grows until the case diverges. Its links to planes 1 and 4, which hold the previous iteration, act as an added relaxation, so the z spacing changes its convergence. *Ns2d* has neither.

The temperature t is not solved for.

In *cns3d* node (2, 2, 2) must be in the fluid; inside a block the pressure drifts by the same amount at each iteration.

With flag 4 a turbulent bend does not converge: the errors of u and v oscillate near 1.E-2.
