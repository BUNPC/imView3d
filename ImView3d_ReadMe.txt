Graphing Steps:

I) Set parameters first

	1) Auto-seed 
	Set Hxy = 40, Hz = 20 (box to mask around each seed)
	I thresh = 0.05 (check first that all vessels are exposed for this 	 
 	threshold)
	Ntracks = 25

	2) Seed 1 to...
	BIDIRECTIONAL
	#Trcks = 25 (same as above)
	Hxy = 4; Hz = 4; (box in which to calculate Monte Carlo PDF)
	Function = 0.05 0.5 (Threshold and power-square root)

	3) Seed conectivity to
	Hxy = 5; Hz = 1; Nthresh = 2

II) Click Auto-seed and choose fast (large RAM) option

III) Save seeds after auto-seeding

IV) Click on Path... to directory with seyshelles scripts
	There you must have launcher.sh running on seyshelles
	log directory must be created in advance as well
	image file will be saved on first Path command
(Expanded): 
	Start terminal window (tc shell terminal)
	ssh seychelles.nmr.mgh.harvard.edu
	change directory (use cd command) to Seed directory
	check that in Seed directory you have log subdirectory
	in command window (from seychelles) run launcher.sh
	go back to imView3d program and click on Path and select Seed dir

	
V) Click on Run to start Monte Carlo seed propagations...

VI) Click on Process Seeds - this will ask you to import all the
	processed seed files (tracks). After all seed tracks are 
	loaded, save the imView3d file under a new name 
	(with name like AfterSeedProcesing.seed)

VII) Click on Seed Conectivity to fill conectivity matrix. 
	Before this, set Hxy=5, Hz=2, Nthresh=2

VIII) Check 'Save' at 'Image Tracks' and click on Image Tracks to plot all tracks. Save checkbox will allow you to save tiff file with tracks at the end

VIII) Save data as AfterSeedConnectivity.seed (or something like this); later you can open imView3d from this point and continue...

IX) Start Graph tracks
 
