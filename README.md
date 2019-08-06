# ACCEL+

A fast ring dynamics simulation code.

## Usage

X. Huang  
10/22/2017

Setup commands are case insensitive.

Each block of commands starts with the block type and ends with 'End'. Block types include:

- Lattice
- Optics
- Track
- DAMA

The first block must be Lattice. Other than that, blocks can appear in any order and for multiple times. Each block can consist of multiple commands. 

### Lattice setup block

```
Lattice
  file_name = lat_lelt.txt #,type=ERING
  reset_first: famname=RF #SEPTUM
# reset_first: type=Cavity, famname=RF, index=1
  print_elem_value: type=StrMPole, famname=SFM, PolyB[2]
# print_elem_value: type=StrMPole, famname=SFM, Index,Length,PolyB[2]
# print_elem_value: type=Quad, famname=QF, Index, Length, K
  export_element: type=Cavity, famname=RF
  modify_element: type=Cavity, famname=RF, index=1, field=Voltage, value=3.0
# modify_element: type=Cavity, famname=RF, index=1, field=Voltage, delta=-0.2	
  set_flags: wake=off, cavity=on,radiation=off,crabcav=off
  correct_tune: type=Quad, famname=QF QD [0.13, 0.22]
  correct_chrom: type=StrMPole, famname=SF SD [2, 2]
  print_lat_summary
  set_max_num_threads = 6
End
```

- `file_name = FILE[,type=TYPE]` 

  Specify lattice file name. It will load the lattice. If loading fails, it will print out an error message “Exception opening/reading file: xxxx” where xxxx is the file name. Initially radiation is turned off, RF cavities are turned off, wake elements are turned off, crab cavities are turned off.

  The only lattice type recognized is ERING (electron storage ring). For this type, the RF phase is set to match the radiation energy loss. 

- `reset_first: type=TYPE, [famname=FAMNAME][,index=INDEX]`

  `reset_first: famname=FAMNAME`
  
  Shift the lattice elements to make the ring starts with the specified element. 

- `print_elem_value: type=TYPE, [famname=FAMNAME], FIELD1, [FIELD2],...`

  Print out specified field values for selected elements. Only implemented for two types so far: `StrMPole` and `Quad`. If famname is not specified, all elements of the type will be selected. Field types can be:
  ```
  StrMPole: Index, Length, PolyB[1], PolyB[2], PolyB[3], PolyA[1], PolyA[2], PolyA[3]
  Quad: Index, Length, K
  ```  
  The output block starts with:

  `$$$Print-Element_Values***********************`

- `print_lat_summary` (use this only if the lattice is a ring)

  Print out some lattice parameters, including the number of elements, the type and famname of the first element, tunes, chromaticities, closed orbit and Courant-Snyder parameters at the ring entrance (i.e., the entrance face of the first element). 

- `export_element :  type=TYPE[, famname=FAMNAME]`

  Export selected elements in lattice print out format. If type is RING or LINE, all elements are selected. 
  Output block starts with

  `$$$Export-Elements***********************`

- `correct_tune: type=TYPE, famname=FAMNAME1 FAMNAME2 [nux, nuy] [REPEAT]`

  Correct the betatron tunes, with target tunes specified in the square bracket. Two famnames are needed.  Type is usually Quad. If REPEAT appears, it will correct a second time.

  The output block starts with:

  `$$$Correct Tunes***********************`

- `correct_chrom: type=TYPE, famname=FAMNAME1 FAMNAME2 [chromx, chromy] [REPEAT]`

  Correct the chromaticities, with target chromaticities specified in the square bracket. Two famnames are needed. Type is usually StrMPole. If REPEAT appears, a second correction will be done.

  The output block starts with:

  `$$$Correct chromaticities***********************`

- `set_max_num_threads=NUMBER`

  Set the maximum number of threads in tracking. 

- `modify_element: type=TYPE[, famname=FAMNAME][, index=INDEX], field=FIELD, value=VALUE | delta=DELTA`

  Modify the fields of lattice elements. When INDEX is 0, all elements with specified type and famname are chosen. FIELD needs to be the same fields recognized by the setfield function of the corresponding element type. 

- `set_flags: [RADIATION=ON|OFF]|RADQE|LUMP][,CAVITY=ON|OFF][,CRABCAV=ON|OFF][,WAKE=ON|OFF]`

  Speficify flags for radiation, rf cavity, crab cavity, and wake types.

### Optics print out block

- `print_optics: type=TYPE, famname = FAMNAME, index=INDEX, dpp=DPP, delta=DELTA`

  Print CS functions and phase advances at selected elements. When type is ALL, all elements are selected. “dpp” is fixed momentum error at which the optics functions are calcualated. “delta” is the small step in numerical difference calculation. By default, dpp=0, delta=1.0e-8.

  The output block starts with:

  `$$$Optics-Twiss, dp=??***************************`, where `??` is the dpp value

- `print_matrix: type=TYPE, famname = FAMNAME, index=INDEX, dpp=DPP, delta=DELTA, line_style=ONELINE|4by4|6by6, print_option`

  Print out transfer matrix. Default style is 6 by 6.

  The output block starts with:

  `$$$OPTICS-MATRIX, dpp=??"******************* ???? points`

  where `????` is the number of output points. If in ‘oneline’ style, the matrix elements will be printed out in row order, i.e., first row followed by second row, etc. 

- `print_tune: dpp=DPP, delta=DELTA`
  
  `print_chrom: dpp=DPP, delta=DELTA`

  The two commands have the same effects. Both will print out the betatron tunes and the chromaticities. Note that the radiation is turned off for all elements. 

### Tracking block

In the tracking block one can specify the initial distribution, set up the monitor options, and launch tracking with options. After tracking the particle distribution can be transported to other locations and be printed out.

```
Track
#initial distribution: load from a file, random distribution, specify initial coordinate here
#	Distribution : 
#	Distribution : Load r_init.dat, format=plain  #accel-export
  Distribution : Load test0.out, format=accel-export
#	Distribution : Zero NPar = 10
#	Distribution  : Coordinate r=(0, 0.0001, 0, 0.0001, 0, 0; 0, 0.0005, 0, 0.0005, 0, 0)
#	Distribution  : Coordinate r=(0, 0.0001, 0, 0.0001, 0, 0) Print
#	Distribution : Gaussian emit_x=10E-9, emit_y=10E-12, sigz=0.006, \
#		sigdp=0.001 NPar = 10, shift=(0,0,0.001, 0, 0, 0) Print
#	Distribution : Gaussian emit_x=10E-9, betx=10, emit_y=10E-12, bety=5, sigz=0.006, \
#		sigdp=0.001 NPar = 10, shift=(0,0,0.001, 0, 0, 0) Print
#	Distribution : Gaussian sigx=0.3e-3, sigxp=0.3e-4, sigy=0.05e-3,sigyp=0.01e-3, sigz=0.006, sigdp=0.001, NPar = 1000
#	Distribution : Linear X=0:-0.001:-0.01 shift=(0,0,0.0, 0, 0, 0) print
#	Distribution : Linear z=-0.01:0.001:0.01 shift=(0,0,0.0, 0, 0, 0) print
  
#monitor setup
  Monitor : Interval = 1, Fraction = 1.0, PrintFlag=Stat2  
#default is to print out coordinate, or use Coordinate. Use Stat1 to print average and sigma, Stat2 in addition print out sigxx, sigxxp, sigxdelta, sigxpdelta, sigyz, sigypz
  #Do_Track: NTurn = 512, Radiation=Off, RFCavity=Off
  #Do_Track: NTurn = 12, Radiation=Rad, RFCavity=On, PrintFinal=yes #Full
  Do_Track: NTurn = 20, Radiation=RadQE, RFCavity=On PrintFinal=Yes

  Transport: type=Monitor, famname=FirstMonitor  ##this export the final distribution
  #Transport: type=Aperture, famname=AP, index=1
  #Transport: type=Aperture, famname=AP, index=[2, 3]
End
```

- `Distribution`

  There are multiple ways to specify the initial particle distribution.

  1. If no option is given, it uses one particle with all coordinates at zero.
  2. Load distribution from a data file. Format can be ‘plain’ or ‘accel-export’. In the case of ‘plain’, the text data file contains N rows, each row with 6 numbers for [x, px, y, py, z, δ]. 
