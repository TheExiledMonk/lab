Cluster: Abell 370
Method: SaWLens 
Author: Julian Merten

Description: The "map" file contains a single reconstruction of the cluster field and shows as a FITS image the convergence, both components of the shear, the determinant of the lens mapping Jacobian and the magnification. All FITS images contain a WCS. Several versions for different redshifts are provided and the FITS headers contain the lens redshifts under the keyword "Z_L" and the source redshift under "Z_S".  The file is a multi-extension FITS file with the following image content:

primary: convergence
ext.: shear1
ext.: shear2
ext.: jacdet
ext.: magnification

Map details:
Centre: 39.971145; -1.582251  [RA/DEC/in deg]  
x-dim [pixels]: 240
y-dim [pixels]: 240
Field size x [arcsec]: 1500 
Field size y [arcsec]: 1500
pixel scale [arcsec]: 6.25

Input data WL:
-Subaru Rc-band shear catalogue, 16394 sources at an effective lensing redshift of z_s = 1.11, background selected with a colour-colour method. This catalogue cover the full 1500 arcsec field. Source: Umtesu et al. 2011 (ApJ 729, 127) 
-HST/ACS F814W shear catalogue,  sources at an effective lensing redshift of z_s = 0.8, background selected by photo-zs. This catalogue covers the inner ~300 arcsec of the field. Source: Tim Schrabback, private comm. via Marusa Bradac 

Input data SL:
-The following multiple image system following the nomenclature of the Frontier Field arcs Google spreadsheet used by the map makers. Source: J. Merten and the ST FF map makers:

ID   RA	       DEC	z     delta_z

1.1	39.96686	-1.57691	0.806	0.01	
1.2	39.97609	-1.57603
1.3	39.96847	-1.57661
2.1	39.97364	-1.58421	0.725	0.01
2.2	39.97077	-1.58507
2.3	39.96855	-1.58450
2.4	39.96919	-1.58471
2.5	39.96942	-1.58483
3.1	39.96546	-1.56686	1.421	0.1
3.2	39.96833	-1.56582
3.3	39.97708	-1.56718
4.1	39.97944	-1.57632	1.275	0.1
4.2	39.97055	-1.57626
4.3	39.96174	-1.57793
5.1	39.97328	-1.58907	1.28	0.3
5.2	39.97095	-1.58924
5.3	39.96862	-1.58903
6.1	39.96927	-1.57719	1.063	0.01
6.2	39.96432	-1.57825
6.3	39.97944	-1.57714
7.1	39.96958	-1.58041	1.78	0.3
7.2	39.96967	-1.58074
7.3	39.96861	-1.58561
7.4	39.96137	-1.58001
7.5	39.98360	-1.57797
8.1	39.96429	-1.56979	2.41	0.3
8.2	39.96168	-1.57368
8.3	39.98391	-1.57335
9.1	39.96222	-1.57790	1.54	0.4
9.2	39.96929	-1.57626
9.3	39.98183	-1.57655
11.1	39.96363	-1.56936	5.5	0.6
11.2	39.96058	-1.57415
12.1	39.96953	-1.56667	3.1	0.4
12.2	39.95901	-1.57530
12.3	39.98392	-1.57091
    


Map production details:
The SaWLens method (Merten et al. 2009, 2011) was used to produce these maps with a three-level adaptive mesh scheme. The first level is a full-field run using all WL and SL data on relatively low resolution, which then serves as a template for the run on the med-resolution regimes focussing on the area which is covered by the VLT and HST WL catalogue. This reconstruction is then again used as a template for the final reconstruction of the cluster core on high resolution and which is mostly dominated by the CL constraints.

lowres:
-750.0 -- 750.0; -750.0 -- 750.0  [x/y/in arcsec] around centre 
60x60 [pixels] with 25 arcsec resolution per pixel
medres:
-50.0 -- 75.0; -50.0 -- 75.0  [x/y/in arcsec] around centre
13x13 pixels with 9.62 arcsec resolution per pixel
highres:
-31.3 -- 75.0; -31.3 -- 75.0  [x/y/in arcsec] around centre
17x17 pixel with 6.25 arcsec resolution per pixel

Error maps info:
Also attached are error maps for convergence, shear and magnification, following the same FITS extension scheme as for the actual maps. All error maps contain the a WCS and are available for the same redshifts as the actual maps.

Map details:
Centre: 39.971145; -1.582251  [RA/DEC/in deg]  
x-dim [pixels]: 240
y-dim [pixels]: 240
Field size x [arcsec]: 1500 
Field size y [arcsec]: 1500
pixel scale [arcsec]: 6.25


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
Field centre: 39.971145; -1.582251  [RA/DEC/in deg]  
Map size: 30x30 [pixel]

medres
Filename: */rec2_BS<N>.fits
Number of bootstraps: 250/250/250/250/250 [z=1/2/4/9/20000]
Field size: -50.0 -- 75.0; -50.0 -- 75.0  [x;y in arcsec]
Field centre: 39.971145; -1.582251  [RA/DEC/in deg]  
Map size: 12x12 [pixel]

highres
Filename: */rec3_BS<N>.fits
Number of bootstraps: 698/783/773/769/763  [z=1/2/4/9/20000]
Field size: -50.0 -- 75.0; -50.0 -- 75.0  [x;y in arcsec]
Field centre: 39.971145; -1.582251  [RA/DEC/in deg]  
Map size: 20x20 [pixel]
    
