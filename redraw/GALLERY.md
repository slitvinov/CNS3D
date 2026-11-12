# Appendix: the report's figures redrawn

Each figure of NASA CR-178818 with the redraw laid over the scan: grey where only the scan has ink, red
where only the redraw has, black where both have. The score is the share of the scan's ink within 2 px
of the redraw / the share of the redraw's ink within 2 px of the scan, in percent (one number where only
the first was measured), for the drawing as fitted, with no further registration. How each was made is in [README.md](README.md).

## Sketches

Line art traced from the scan into polylines, filled shapes and labels (`figN.sk`, drawn by `sketch.f`).

### Figure 1. Three-dimensional grid structure and labeling around a grid node P

![Figure 1, redraw over the scan](compare/fig1.png)

Match: 99.1 / 99.1

### Figure 2. Locations where the variables are stored

![Figure 2, redraw over the scan](compare/fig2.png)

Match: 95.3 / 98.2

### Figure 3. Control volumes where the mass conservation is evaluated

![Figure 3, redraw over the scan](compare/fig3.png)

Match: 99.4 / 98.6

### Figure 4. Physical geometry and wall boundary conditions of the driven cavity

![Figure 4, redraw over the scan](compare/fig4.png)

Match: 97.5 / 99.0

### Figure 9. Geometry and boundary conditions of the laminar backward-facing step (1:2)

![Figure 9, redraw over the scan](compare/fig9.png)

Match: 98.7 / 98.5

### Figure 12. Geometry and boundary conditions of the turbulent backward-facing step (2:3)

![Figure 12, redraw over the scan](compare/fig12.png)

Match: 98.7 / 99.0

### Figure 16a. Geometry of the 90-deg bend

![Figure 16a, redraw over the scan](compare/fig16a.png)

Match: 97.4 / 98.6

### Figure A-1. Global structure of the program CNS3D

![Figure A-1, redraw over the scan](compare/figA-1.png)

Match: 89.1 / 93.3

### Figure A-2. Grid mesh structure for 3-D calculations

![Figure A-2, redraw over the scan](compare/figA-2.png)

Match: 95.5 / 96.7

## Figures digitized from the scan

Curves, points and meshes digitized from the scan into data files, drawn by a DISLIN program with settings fitted to the scan.

### Figure 5a. Mesh for the driven cavity: uniform and orthogonal

![Figure 5a, redraw over the scan](compare/fig5a.png)

Match: 100 / 94

### Figure 5b. Mesh for the driven cavity: nonuniform and nonorthogonal

![Figure 5b, redraw over the scan](compare/fig5b.png)

Match: 100 / 99

### Figure 7. Velocity profiles along the mid-section of the cavity

![Figure 7, redraw over the scan](compare/fig7.png)

Match: 88 / 89

### Figure 8. Convergence history for the driven cavity, Re = 400

![Figure 8, redraw over the scan](compare/fig8.png)

Match: 89 / 90

### Figure 10. Reattachment length versus Reynolds number, laminar step

![Figure 10, redraw over the scan](compare/fig10.png)

Match: 83 / 85

### Figure 13. Locus of flow reversal, turbulent step

![Figure 13, redraw over the scan](compare/fig13.png)

Match: 93 / 91

### Figure 16b. Mesh of the 90-deg bend

![Figure 16b, redraw over the scan](compare/fig16b.png)

Match: 97 / 92

### Figure 19. Primary velocity profiles in the bend

![Figure 19, redraw over the scan](compare/fig19.png)

Match: 85 / 85

## Figures computed again with CNS3D

Computed again with CNS3D (TESTS) and drawn by a DISLIN program with settings fitted to the scan; differences show where the runs differ from the report's.

### Figure 6. Velocity vector plots in the driven cavity

![Figure 6, redraw over the scan](compare/fig6.png)

Match: 96 / 89

### Figure 11. Streamlines of the laminar step at Re = 600

![Figure 11, redraw over the scan](compare/fig11.png)

Match: 78 / 71. Our runs of case B reattach at 8.3h and 10.3h, the report at about 8.6h and 11h; in (a) our vortex centre lies near 6h, the report's near 4.75h; the scan is turned by -0.43 deg.

### Figure 14. Streamlines of the turbulent step

![Figure 14, redraw over the scan](compare/fig14.png)

Match: 82 / 76. Our run of case C reattaches at 5.3h, the report at 5.4h; the downstream scale is fixed by the report's label; the scan is turned clockwise by 0.31 deg.

Turned as the scan is, against the scan as printed: 82 / 76 (drawn level: 50 / 44).

![Figure 14, redraw turned as the scan, over the scan as printed](compare/fig14_tilt.png)

### Figure 15. Contours of turbulent kinetic energy, turbulent step

![Figure 15, redraw over the scan](compare/fig15.png)

Match: 85 / 70. k in our run peaks at 0.063; the report's highest contour is 0.045. The 0.005 contour along the top wall upstream is not in the scan.

Turned as the scan is, against the scan as printed: 85 / 70 (drawn level: 71 / 56).

![Figure 15, redraw turned as the scan, over the scan as printed](compare/fig15_tilt.png)

### Figure 17. Primary velocity in the bend

![Figure 17, redraw over the scan](compare/fig17.png)

Match: 90 / 87

### Figure 18. Secondary velocity in the bend

![Figure 18, redraw over the scan](compare/fig18.png)

Match: 72 / 76
