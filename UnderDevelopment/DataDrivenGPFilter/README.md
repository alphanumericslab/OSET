Data-drive GP filter four pager article corresponding implementation
The Functions & Scripts used for results reported in the *A Data-Driven Gaussian Process Filter for Electrocardiogram Denoising* four page article.

The main function, i.e. the GP filter:
*GaussianProcessFilterInPhaseDiagFast.py*
-corresponds to the diagonal case, as reported in the article;
-corresponds to the Algorithm 1 from the article

The script generating the CSV files corresponding to the GP filtered signals: 
*GenerateCSVRawNoisyAndGPPhaseFilteredSignals.ipynb*
-Basically a two for loops applying the GP filter function for selected levels of noise and repetitions;
-Contains the full pipe-line, i.e. the baseline wander removal, the R-peak detection;
-It generates also the pure/noisy corresponding signals (used as inputs for when applying the wavelet function in Matlab)
-The GP filtered signals, the noise signals, the noisy signals are all saved as .csv files containing the GP filtered signals or noise signals or noisy signals for all records in the database (in this case QT database) for all leads (in this case 2 leads).

The script generating the CSV files corresponding to the Wavelet filtered (I.e. the benchmark) : 
*WaveletFilter.m*
-using the specific wavelet denoiser reported by Prof. Sameni as being the most efficient (hence the use as benchmark);
-producing a corresponding .csv file with the wavelet filtered signals;

The script generating the CSV files corresponding to the Prof. Sameni GP implementation: 
*RezaFilter.m*
- used for testing purposes, in order to compare the Python implementation;
- the implementation does not account for the presence of gramian in the model (the offset created by this is clearly visible for when the MDMA baseline wander removal is used when looking at the signals, it is also numerically visible for the case when "double LP filter" baseline wander method is used;
-producing a corresponding .csv file with the wavelet filtered signals;

The script generating the CSV files corresponding to the QT estimations:
*QTEstimationForGPFilter.m*
- implemented to execute over all levels of noise and all repetitions, for GP filtered signals, wavelet filtered signals and pure ECG signals (in order to compute the differences);
 - it is using as the core function  Prof. Li’s QT detection function: [QT_analysis_single_lead.m](QT_analysis_single_lead.m)

FIG3: the script generating the individual ECG comparisons (measurements vs. filter) for Prior, Posterior and Wavelet filtered signals
*ArticleECGFigsGeneration2*
- produces the comparison for any selected recording from the database, for any selected level of noise;

FIG4: the script generating the SNR curves comparing the GP filter and wavelet filter performances for all selected levels of noise & all repetitions
*AnalyseGPPhaseFiltAndWaveletFiltCSVFiles2.ipynb*

FIG5: the script generating the QT IQR comparing the GP filter & wavelet filter performances for all selected levels of noise & all repetitions
*QTIQRFigureGenerator.ipynb*


ADDENDUM
The script allowing to see how the GP filter behaves on one signal, i.e. visualize the filtered compared with the pure/noisy signals and the SNR improvement for any chosen signal from the database, use: 
*testGPPhaseModels4.ipynb*
-uses the main *GaussianProcessFilterInPhaseDiagFast.py* function;

The script allowing to see how how Prof. Sameni's implementation behaves depending on the baseline wander removal method:
*test_MAPFilter_QT_analysisMircea.m*

MENTIONS:
-all functions used by the Python GP filter implementation are using .py functions (i.e. the baseline wander removal, the R-peak, the transformation matrices, etc);

-there are lots of flavours for this GP filter that can be considered: filtering the phase domaine mean, using a filtered version of the transformation matrices, etc; they are also implemented and can be in the repository;
