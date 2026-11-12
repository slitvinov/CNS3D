# The report's figures redrawn

Each data figure of the report is here as data (`.csv`, and `.dat` as read by the programs) and a
Fortran program that draws it with DISLIN 11. The drawing settings (lettering sizes and widths,
legend layout, symbol radii, pen widths, arrow scales, the place of each panel) are in
`figN_style.dat`, fitted to the scan with `fitfig.py N`, which maximizes the normalized
cross-correlation of the blurred ink, coarse to fine.

The scores compare the drawing as fitted, laid on the scan at the same pixel size with no further
registration (`imatch.py` also reports scores after an affine registration of its own, which are
higher). On a dense drawing a score has a floor: Figure 5b shifted by 4, 6 and 10 px scores 56, 28
and 19 percent.

| Figure | data | from | scan ink within 2 px of the redraw / redraw ink within 2 px of the scan |
|---|---|---|---|
| 5a, 5b cavity meshes | `fig5{a,b}_nodes.csv`, 31x31 nodes | digitized, `fig5.py` | 100% / 94%, 100% / 99% |
| 6 cavity velocity vectors | `fig6{a,b}.dat` | CNS3D on the Figure 5 grids (TESTS, case A) | 96% / 89% |
| 7 cavity centreline | `fig7_*.csv` | digitized, `fig7.py` | 88% / 89% |
| 8 convergence history | `fig8_{u,v,p}.csv` | digitized, `fig8.py` | 89% / 90% |
| 10 reattachment length | `fig10_*.csv` | digitized, `fig10.py` | 83% / 85% |
| 11 streamlines, 1:2 step, Re 600 | `fig11{a,b}.dat`, stream function | CNS3D, case B on a 45x45 grid: (a) as printed, (b) wallvis, pecap=20 | 78% / 71% |
| 13 flow reversal | `fig13_*.csv` | digitized, `fig13.py` | 93% / 91% |
| 14 streamlines, 2:3 step | `fig14{u,d}.dat` | CNS3D, case C on a 45x42 grid | 82% / 76% |
| 15 k contours, 2:3 step | `fig15{u,d}.dat` | CNS3D, case C on a 45x42 grid | 85% / 70% |
| 16b bend mesh | `fig16_grid.json`, `fig16_nodes.csv` | measured, `fig16.py` | 97% / 92% |
| 17 bend primary velocity | `fig17{a,b,c}.dat` | CNS3D, case D | 90% / 87% |
| 18 bend secondary velocity | `fig18{a,b,c,d}.dat` | CNS3D, case D | 72% / 76% |
| 19 bend profiles | `fig19{a,b}_theta*_*.csv` | digitized, `fig19.py` | 85% / 85% |

In `fig5b_nodes.csv` the lines across of columns i = 25 to 30 were found again after `fig5.py`: on
columns 27 and 28 it took a mark near the top for a line, which shifted the rows below by one and gave
a zigzag. Each was searched again within 9 px of the line through its neighbours.

The scans used: `../fig/p27-1.png` (Figure 10); the others are cropped from report.pdf rendered at
300 dpi (`pdftoppm -r 300 -gray`), with a margin, turned upright where the page is sideways, and
rotated to undo the skew measured from the ticks or axes: Figure 7 by 0.52 deg, 8 by -0.405, 13 by
-0.44, 16 by -0.745, 19 by -0.31 (`figN_deskew.png`, `figN_page.png`). `fig16_bscan.png` is the
Figure 16 scan with the sketch (a) blanked, as it is not redrawn.
The scans of Figures 5a, 5b and 16(b) were used tilted until `score.py` measured each against
its redraw (a rigid fit of rotation and shift); they are now turned by -0.32, -0.24 and +0.11
deg (counter-clockwise positive, `fig5a_deskew.png`, `fig5b_deskew.png`, `fig16b_deskew.png`) and
the drawings fitted again. Figure 11 is turned by -0.43 deg (`fig11_deskew.png`), which levels its
four walls to within 0.04 deg; the -0.26 deg first measured against the redraw left them sloping by
3 to 4 px over the width, as the redraw had been fitted to the tilted scan. The rigid fit still asks
for -0.12 deg and a 1.9 px shift, from the bubbles, where the drawings differ; the walls lie on the
scan's.

The lettering is Helvetica-Bold condensed (DISLIN `CHAWTH`), the symbol of the measured points a ring
with a dot (`MARK2` in `rlib.f`), the legends laid out by `MYLEG` (`rlib.f`).

## Sketches

Figures 1 to 4, 9, 12, 16(a), A-1 and A-2 are drawings, not data. Each is in `figN.sk`: polylines
(each with its pen width), filled polygons (arrow heads, dots, hatching) and labels (text, left end,
right end, baseline, capital height; text starting with `@` is set in the Symbol font). `sketch.f`
draws any of them:

    gfortran -std=legacy -o sketchd sketch.f -L$DISLIN -ldislin
    ./sketchd fig1.sk fig1_dislin.pdf "$(cat fig1_sketch.par)"

The three numbers in `figN_sketch.par` scale the capital height, shift the baselines and scale the
pen widths; `skfit.py` fits them to the scan and then moves each label by up to 6 px to its place.

`sketch.py SCAN OUT.sk [LABELS]` makes a `.sk` file: the ink outside the labels is split into filled
parts (wider than the strokes) and strokes, the strokes are thinned to their centre lines, traced
between junctions and simplified to 0.8 px. The labels of Figure A-1 are read with Tesseract (3x,
`--psm 3`); those of the others, which OCR does not read reliably, are given by hand in
`figN_labels.json` (text, approximate centre and capital height) and placed by `sklabel.py`, which
matches each label set in Nimbus Sans Bold (Helvetica) or Standard Symbols against the scan
(`figN_labelbox.json`). Figure 16(a) is cropped from `../fig/p30-1.png` (`fig16a_scan.png`).

| Figure | scan | scan ink within 2 px of the redraw / redraw ink within 2 px of the scan |
|---|---|---|
| 1 grid node labels | `../fig/p16-1.png` | 99.1% / 99.1% |
| 2 variable locations | `../fig/p18-1.png` | 95.3% / 98.2% |
| 3 control volume | `../fig/p21-1.png` | 99.4% / 98.6% |
| 4 cavity geometry | `../fig/p22-1.png` | 97.5% / 99.0% |
| 9 1:2 step geometry | `../fig/p26-1.png` | 98.7% / 98.5% |
| 12 2:3 step geometry | `../fig/p28-1.png` | 98.7% / 99.0% |
| 16(a) bend geometry | `fig16a_scan.png` | 97.4% / 98.6% |
| A-1 program structure | `../fig/p38-1.png` | 89.1% / 93.3% |
| A-2 grid indices | `../fig/p39-1.png` | 95.5% / 96.7% |

As in the report, zero in the boundary conditions is printed as the letter o.

Figures 11, 14 and 15 are contours (`CPANEL` in `clib.f`) of the stream function over the flow rate
and of k from the runs; `figstep.py` makes the `.dat` files from the unit 7 files of the runs
(`step_hyb600.f7`, `step_pc20600.f7`, `tstep.f7`; TESTS, "Figures 11, 14 and 15"). The main stream
function levels (0.1 to 0.9 of the flow rate, near-wall levels) are those of the scan; the levels
inside the recirculation bubbles are not printed and are fitted. The runs are converged to EREXT =
1.E-5: with the program's 1.E-3 a run continued from a lower Re stops after a few hundred
iterations short of the steady state. They reattach at 8.3h and 10.3h (report about 8.6h and 11h)
and 5.3h (report 5.4h), and k reaches 0.063. Figures 14 and 15 are drawn as the report draws them,
with the channel downstream of the step at 2/3 of the scale upstream (`SXD` = 88.8 px/h, fixed by the
report's 5.4h label, not fitted: fitted, it took up our shorter bubble). Their scans are turned
clockwise by 0.31 and 0.29 deg (`fig14_deskew.png`, `fig15_deskew.png`, from the top wall and the
step face); the run's grid, hybrid factor and relaxation were searched against them (TESTS).
Figure 11 differs from the report inside the bubbles: the scan's contours reach into the corner
at the foot of the step, and in (a) its vortex centre is near 4h, ours near 6h. A grid clustered
towards the step face (432 runs, TESTS) draws contours into the corner but shortens the bubble and
matches the scan worse; the grids used fit best. Both panels of the scan have a bubble contour
along the floor and the step face, 0.05h to 0.14h from them; it is drawn at -0.003 of the flow rate
(`ZM3`, the best of -0.002 to -0.013, fitted from 72/67 to 74/68). Our innermost contour in (a),
-0.084 near 6h, has no counterpart there; the scan's innermost one lies at 3h to 5h, where our
vortex is not: drawn without it, the match is the same. The upper-wall bubble of (b) is the
scan's 1.0024 contour; with (b) run with the diffusion cap |C|/20 (TESTS) the panel matches 81/75,
the figure 78/71.
Figure 15 has a 0.005 contour along the top wall upstream that the scan does not show.
Turned back as the scan is (`compare/fig14_tilt.png`, `compare/fig15_tilt.png`), the redraws match
the scans as printed as well as they match the deskewed ones, 82/76 and 85/70 percent; drawn level
against the printed scans they match 50/44 and 71/56.

Run, with DISLIN installed (`DISLIN`, `LD_LIBRARY_PATH`):

    gfortran -std=legacy -o fig8d fig8.f rlib.f -L$DISLIN -ldislin && ./fig8d     # fig8_dislin.pdf
    gfortran -std=legacy -o fig5ad fig5a.f fig5.f -L$DISLIN -ldislin && ./fig5ad
    python3 fig8.py                  # digitize again from the scan
    python3 fitfig.py 8              # fit the settings again (needs pdftoppm, or PDFTOPPM=path)
    python3 imatch.py fig8_deskew.png fig8_dislin.png m_fig8   # registration, scores, overlay
