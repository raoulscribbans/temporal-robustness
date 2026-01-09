# temporal-robustness
code for: Phenologically Explicit Robustness Metric Reveals Increased Vulnerabilities in Temporal Plant-Pollinator Networks


The following scripts are used to produce the results and figures shown in the main text

	extinction_curve.m is used to generate the extinction curve that is a subfigure of figure 1
	donana_analysis_figures.m generates figures 2 and 4
	uniform_model_fig_3.m generates figure 3
	supplementary_analysis.m generates all figures in the supplemenray material

Note that uniform_model_fig_3.m and supplementary_analysis.m were evaluated using parallel computing on a computer cluster 

The scripts above are dependent on (some or all of) the following files, which contain class definitions that are used throughout

	Species.m
	Pollinator.m
	Interaction.m
	TemporalNetwork.m
	plotHelper.m

Additionally the following data files are needed for supplementary_analysis.m

	Interaction_data.csv - this contains interaction data from EuPPollNet, and can be found at https://zenodo.org/records/15183272
   - Note that due to the large file size, we have only included entries with study ID "14_Dupont" here, however the code is written to run on the full data file
	greenlandData.xlsx - this contains data extracted and formatted from Olesen et. al 2008 and Pradal et. al. 2009

The following MATLAB versions and packages were used, and may need to be installed before running

MATLAB                                                Version 24.1        (R2024a)
MATLAB Compiler                                       Version 24.1        (R2024a)
MATLAB Compiler SDK                                   Version 24.1        (R2024a)
Parallel Computing Toolbox                            Version 24.1        (R2024a) (?)
Statistics and Machine Learning Toolbox               Version 24.1        (R2024a)
