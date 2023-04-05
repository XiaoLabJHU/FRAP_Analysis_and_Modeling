For RNAP promoter search and transcription kinetics in E. coli cells Bettridge and Harris et al 2023:

FRAP_batch_process_concatenate_files
	Concatenates prebleach and postbleach tifs

FRAP_ring_circleFH.m - Initial Image Analysis For FRAP
	Uses FRAP_ring_circleFH.fig to initialize GUI.
	Upload brightfield and concatenated FRAP image from ImageJ macro
	Select top and bottom of the cell (longwise) for rotation
	Select photobleached area
	Select Cell ROI
	Select Background ROI (representative region away from cell)
	Save Structure

FRAP_analysis
	Expects folders for each day with subfolders named "Cell_01" etc
	Will go through and produce normalized FRAP signal for each cell
	Next section combines cells from multiple days
	Should produce a .mat file with avg_FRAP, std_FRAP, and time
	Also a .mat file with a cfit object for exponential fitting of avg_FRAP

FRAP_modeling 
	One round of what was used for bootstrapping
	Has one step modeling and two step modeling

CreateBootStrappedFrapData
	Creates bootstrapped structures. Sampling numbers (40 or 5 for example) and indexes should be adjusted based on sample size

FRAPbootstrapping
	Two Step Model

FRAPbootstrappingOneStep
	One Step Model

Contact fharri19@jhu.edu or jhu.xiaolab@gmail.com with questions.

