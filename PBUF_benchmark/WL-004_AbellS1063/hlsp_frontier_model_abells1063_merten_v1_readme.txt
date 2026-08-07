Cluster: Abell S1063
Method: SaWLens 
Author: Julian Merten

Description: The "map" file contains a single reconstruction of the cluster field and shows as a FITS image the convergence, both components of the shear, the determinant of the lens mapping Jacobian and the magnification. All FITS images contain a WCS. Several versions for different redshifts are provided and the FITS headers contain the lens redshifts under the keyword "Z_L" and the source redshift under "Z_S".  The file is a multi-extension FITS file with the following image content:

primary: convergence
ext.: shear1
ext.: shear2
ext.: jacdet
ext.: magnification

Map details:
Centre: 342.18322, -44.530908 [RA/DEC/in deg]
x-dim [pixels]: 132
y-dim [pixels]: 132
Field size x [arcsec]: 1500 
Field size y [arcsec]: 1500
pixel scale [arcsec]: 11.36

Input data WL:
-ESO/WFI shear catalogue,  4008 sources at an effective lensing redshift of z_s = 1.049, background selected with a colour-colour method. This catalogue cover the full 1500 arcsec field. Source: CLASH team 
-HST/ACS multi-band shear catalogue,  sources at an effective lensing redshift of z_s = 1.12, background selected by photo-zs. This catalogue covers the inner ~300 arcsec of the field. Source: CLASH team

Input data SL:
-The following multiple image system following the nomenclature of the Frontier Field arcs Google spreadsheet used by the map makers. Source: A. Monna, J. Merten and the ST FF map makers:

ID   RA	       DEC	z     delta_z

1.1	342.19447	-44.527053	1.28	0.1
1.2	342.19582	-44.528889
1.3	342.18653	-44.521278
2.1	342.19482	-44.5274	1.28	0.1
2.2	342.19556	-44.528408
2.3	342.18624	-44.521053
3.1	342.19256	-44.530764	1.25	0.1
3.2	342.19247	-44.530475
3.3	342.17975	-44.521565
4.1	342.18783	-44.527286	1.46	0.3
4.2	342.17922	-44.523606
4.3	342.19316	-44.536536
5.1	342.18772	-44.527572	1.46	0.3
5.2	342.17881	-44.52365
5.3	342.19285	-44.536656
6.1	342.17418	-44.528325	1.34	0.2
6.2	342.17583	-44.532561
6.3	342.18842	-44.539997
7.1	342.1801	-44.538456	1.0	0.2
7.2	342.17545	-44.536008
7.3	342.1719	-44.530267
8.1	342.18185	-44.540567	2.05	0.2
8.2	342.17425	-44.537122
8.3	342.16936	-44.527286
9.1	342.18023	-44.540808	2.7	0.2
9.2	342.17491	-44.538686
9.3	342.16774	-44.526292
10.1	342.18084	-44.540894	2.9	0.2
10.2	342.17455	-44.538336
10.3	342.16792	-44.5262
11.1	342.19085	-44.537472	6.0	0.2
11.2	342.18103	-44.534636
11.3	342.17121	-44.519847

Map production details:
The SaWLens method (Merten et al. 2009, 2011) was used to produce these maps with a three-level adaptive mesh scheme. The first level is a full-field run using all WL and SL data on relatively low resolution, which then serves as a template for the run on the med-resolution regimes focussing on the area which is covered by the VLT and HST WL catalogue. This reconstruction is then again used as a template for the final reconstruction of the cluster core on high resolution and which is mostly dominated by the CL constraints.

lowres:
-750.0 -- 750.0; -750.0 -- 750.0  [x/y/in arcsec] around centre 
44x44 [pixels] with 34.1 arcsec resolution per pixel
medres:
-68.2 -- 102.3; -68.2 -- 102.3  [x/y/in arcsec] around centre
14x14 pixels with 12.18 arcsec resolution per pixel
highres:
-56.84 -- 56.76; -56.84 -- 56.76  [x/y/in arcsec] around centre
10x10 pixel with 11.36 arcsec resolution per pixel

Error maps info:
Also attached are error maps for convergence, shear and magnification, following the same FITS extension scheme as for the actual maps. All error maps contain the a WCS and are available for the same redshifts as the actual maps.

Map details:
Centre: 342.18322, -44.530908 [RA/DEC/in deg]  
x-dim [pixels]: 270
y-dim [pixels]: 270
Field size x [arcsec]: 1500 
Field size y [arcsec]: 1500
pixel scale [arcsec]: 5.55


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
Field size: -750.0 -- 750.0; -750.0 -- 750.0  [x;y in arcsec]
Field centre: 342.18322, -44.530908 [RA/DEC/in deg]  
Map size: 30x30 [pixel]

medres
Filename: */rec2_BS<N>.fits
Number of bootstraps: 250/250/250/250/250 [z=1/2/4/9/20000]
Field size: -68.2 -- 102.3; -68.2 -- 102.3  [x;y in arcsec]
Field centre: 342.18322, -44.530908 [RA/DEC/in deg]   
Map size: 10x10 [pixel]

highres
Filename: */rec3_BS<N>.fits
Number of bootstraps: 784/802/852/784/779  [z=1/2/4/9/20000]
Field size: -68.2 -- 102.3; -68.2 -- 102.3  [x;y in arcsec]
Field centre: 342.18322, -44.530908 [RA/DEC/in deg]  
Map size: 15x15 [pixel]
    
