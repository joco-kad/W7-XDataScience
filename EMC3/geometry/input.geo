* input file for EMC3 code (1D parallel test)
*** 1. fein grid mesh representing magnetic coordinates
1                     !zone number
101 2 4     !core    !R-,P-,T_surfaces
*** 2. non default non transparent, plate surface
*** 2.1 non default surface
* radial
0                     !number of surfaces
* poloidal
* Type of surface=1:  periodic
*                 2:  up/down symmetric
*                 3:  Mapping
2                     !number of surfaces
0 0 1                 !IP  IZONE  Type of surface
0 99 0 2             !IR1-2, IT1-2
1 0 1                !IP  IZONE  Type of surface
0 99 0 2             !IR1-2, IT1-2
* toroidal
2                     !number of surfaces
0 0 1                 !IT  IZONE  Type of surface
0 99 0 0             !IR1-2, IP1-2
3 0 1                 !IT  IZONE  Type of surface
0 99 0 0             !IR1-2, IP1-2
*** 2.2 non transparent (Boundary condition must be defined)
* radial
2                     !number of surfaces
0 0 1                 !IR  IZONE  ISIDE
0 0 0 2             !IP1-2, IT1-2
100 0 -1               !IR  IZONE  ISIDE
0 0 0 2             !IP1-2, IT1-2
* POLOIDAL
0                     !number of surfaces
* TOROIDAL
0                     !number of surfaces
*** 2.3 plate surface (Bohm Boundary condition)
* radial
0                     !number of surfaces
* POLOIDAL
0                     !number of surfaces
* TOROIDAL
1                     !number of surfaces
3 0 -1               !IT  IZONE  ISIDE
0 99 0 0             !IR1-2, IP1-2
*** 3. physical cell
*<0: user defines cell
* 0: default;
*>0: read from file
1 3
* check cell ?
T