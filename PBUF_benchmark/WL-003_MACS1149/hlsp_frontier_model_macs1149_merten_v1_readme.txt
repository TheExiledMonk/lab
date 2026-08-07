Cluster: MACS J1149
Method: SaWLens 
Author: Julian Merten

Description: The "map" file contains a single reconstruction of the cluster field and shows as a FITS image the convergence, both components of the shear, the determinant of the lens mapping Jacobian and the magnification. All FITS images contain a WCS. Several versions for different redshifts are provided and the FITS headers contain the lens redshifts under the keyword "Z_L" and the source redshift under "Z_S".  The file is a multi-extension FITS file with the following image content:

primary: convergence
ext.: shear1
ext.: shear2
ext.: jacdet
ext.: magnification

Map details:
Centre: 177.39877; 22.398532 [RA/DEC/in deg]
x-dim [pixels]: 168
y-dim [pixels]: 168
Field size x [arcsec]: 1200 
Field size y [arcsec]: 1200
pixel scale [arcsec]: 7.14

Input data WL:
-Subaru Rc-band shear catalogue, 6137 sources at an effective lensing redshift of z_s = 1.28, background selected with a colour-colour method. This catalogue cover the full 1500 arcsec field. Source: CLASH team 
-HST/ACS multi-band shear catalogue, 844 sources at an effective lensing redshift of z_s = 0.99, background selected by photo-zs. This catalogue covers the inner ~300 arcsec of the field. Source: CLASH team

Input data SL:
-The following multiple image system following the nomenclature of the Frontier Field arcs Google spreadsheet used by the map makers. Source: CLASH team, J. Merten and the ST FF map makers:

ID   RA	       DEC	z     delta_z

  1.1  177.397   22.396007    1.48	0.01    
  1.2  177.39941 22.397438
  1.3  177.40341 22.402426
  2.1  177.40243 22.389739    1.894	0.01	
  2.2  177.40607 22.392484
  2.3  177.40657 22.392881
  3.1  177.39087 22.39989     2.497	0.01	
  3.2  177.39272 22.403074
  3.3  177.40129 22.40718 
  4.1  177.39301 22.396826    2.95	0.2
  4.2  177.3944  22.400729
  4.3  177.40419 22.40612 
  5.1  177.39976 22.393062    2.7	0.25
  5.2  177.40111 22.393824
  5.3  177.40794 22.403538
  6.1  177.39972 22.392545    3.0	0.3
  6.2  177.40181 22.393858
  6.3  177.40804 22.402505
  7.1  177.39895 22.391332    2.85	0.3
  7.2  177.40339 22.394269
  7.3  177.40759 22.401243
  8.1  177.39849 22.394351    2.7	0.5	
  8.2  177.39978 22.395055
  8.3  177.40706 22.405552
  9.1  177.40515 22.426221    1.6	0.5
  9.2  177.40387 22.427217
  9.3  177.40323 22.427221
  10.1 177.40447 22.425508    1.25	0.5
  10.2 177.40362 22.425629
  10.3 177.4022  22.426611
  12.1 177.39857 22.389356    1.2	0.2
  12.2 177.40375 22.392345
  12.3 177.40822 22.398801
  13.1 177.4037  22.397787    1.25	0.2
  13.2 177.40282 22.396656
  13.3 177.40003 22.393857
  14.1 177.39166 22.403504    3.6	0.4
  14.2 177.39084 22.402624


Map production details:
The SaWLens method (Merten et al. 2009, 2011) was used to produce these maps with a three-level adaptive mesh scheme. The first level is a full-field run using all WL and SL data on relatively low resolution, which then serves as a template for the run on the med-resolution regimes focussing on the area which is covered by the VLT and HST WL catalogue. This reconstruction is then again used as a template for the final reconstruction of the cluster core on high resolution and which is mostly dominated by the CL constraints.

lowres:
-600.0 -- 600.0; -600.0 -- 600.0  [x/y/in arcsec] around centre 
42x42 [pixels] with 28.57 arcsec resolution per pixel
medres:
-57.14 -- 85.71; -57.14 -- 85.71  [x/y/in arcsec] around centre
16x16 pixels with 8.93 arcsec resolution per pixel
highres:
-42.86 -- 42.82; -42.86 -- 42.82  [x/y/in arcsec] around centre
12x12 pixel with 8.33 arcsec resolution per pixel

Error maps info:
Also attached are error maps for convergence, shear and magnification, following the same FITS extension scheme as for the actual maps. All error maps contain the a WCS and are available for the same redshifts as the actual maps.

Map details:
Centre: 177.39877; 22.398532  [RA/DEC/in deg]  
x-dim [pixels]: 180
y-dim [pixels]: 180
Field size x [arcsec]: 1200 
Field size y [arcsec]: 1200
pixel scale [arcsec]: 6.66


Error maps production:
All errors are derived from bootstraps realisations of the three different regimes highlighted above. The maps show a simple 1 SD error in each pixel, as they are derived from the bootstraps using the full sample and not only certain quartiles, as discussed within the map maker teams. The bootstraps of the different regimes are presented in a separate readme file.  
    
Bootstraps:

Description: Because of the way how SaWLens works and to keep runtime at acceptable levels, we produce bootstrap realisations for three different regimes of the cluster field separately. Each bootstrap realisation is a multi-extension FITS file with the following images

primary: convergence
1.ext.: shear1
2.ext.: shear2 
3.ext.: jacdet
4.ext.: magnification
5.ext.: field mask

The last extension shows masked pixels in the reconstructed field by a pixel value of "not 0" instead of 0 for unmasked pixels. All FITS images contain in their header the redshift of the lens under the keyword "Z_L" and the source redshift the map is scaled to under keyword "Z_S".

The bootstraps were derived by bootstrap-resampling the input WL catalogues and by randomly sampling the allowed redshift range of SL features.  

Bootstrap regimes:

lowres
Filename: */rec1_BS<N>.fits
Number of bootstraps: 250/250/250/250/250 [z=1/2/4/9/20000]
Field size: -600.0 -- 600.0; -600.0 -- 600.0  [x;y in arcsec]
Field centre: 177.39877; 22.398532 [RA/DEC/in deg]  
Map size: 30x30 [pixel]

medres
Filename: */rec2_BS<N>.fits
Number of bootstraps: 250/250/250/250/250 [z=1/2/4/9/20000]
Field size: -57.14 -- 85.71; -57.14 -- 85.71 [x;y in arcsec]
Field centre: 177.39877; 22.398532 [RA/DEC/in deg]   
Map size: 10x10 [pixel]

highres
Filename: */rec3_BS<N>.fits
Number of bootstraps: 822/733/704/691/679 [z=1/2/4/9/20000]
Field size: -57.14 -- 85.71; -57.14 -- 85.71 [x;y in arcsec]
Field centre: 177.39877; 22.398532 [RA/DEC/in deg]  
Map size: 20x20 [pixel]
    
