# Noninvasive Fetal ECG Dataset
The data in this folder contains real fetal cardiac signals recorded from 8 channels from the abdomen of pregnant women. This database has been donated to [OSET](https://github.com/alphanumericslab/OSET) by Dr. A. Tokarev from the [Biomedical Signal Processing Laboratory of National Aerospace University, Kharkov, Ukraine](http://www.xai-medica.com/), and may be used under the GNU General Public License, with proper reference to the authors.

The data have variable qualities, ranging from low to high. An analog bandpass filter with a passband of 0.05-100 Hz has been applied during preprocessing. For research interests, the raw data (after the A/D) and the data after preprocessing (using "preprocessing.m" that uses basic [OSET](https://github.com/alphanumericslab/OSET) tools) are provided. The preprocessed data have the same name as the original data with a "_preprocessed" extension. In order to allow advanced processing, the preprocessing is rather minimal, consisting of baseline wander removal and out-of band noise rejection. More advanced processing tools for preprocessing, filtering, and extracting the fetal ECG from this data, can be found in the [Open-Source Electrophysiological Toolbox](https://github.com/alphanumericslab/OSET).

File index | File name | Mother age (years) | Pregnancy term (weeks) | Sampling rate (Hz) | Number of channels ADC (bits) | ADC full scale (mV) | Quantization (µV/bit)
------------ | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | -------------
1	| signal_01.mat	| 25 |	33 |	1000	| 8	| 16	| ± 8	| 0.25
2	| signal_02.mat	| 25	| 35	| 1000	| 8	| 16	| ± 8	| 0.25
3	| signal_03.mat	| 27	| 36	| 1000	| 8	| 16	| ± 8	| 0.25
4	| signal_04.mat	| 31	| 34	| 500	| 8	| 12	| ± 1	| 0.5
5	| signal_05.mat	| 30	| 25	| 500	| 8	| 12	| ± 1	| 0.5
6	| signal_06.mat	| 28	| 38	| 500	| 8	| 12	| ± 1	| 0.5
7	| signal_07.mat	| 26	| 40	| 500	| 8	| 12	| ± 1	| 0.5
8	| signal_08.mat	| 32	| 33	| 500	| 8	| 12	| ± 1	| 0.5
9	| signal_09.mat	| 28	| 37	| 500	| 8	| 12	| ± 1	| 0.5
10	| signal_10.mat	| 25	| 31	| 500	| 8	| 12	| ± 1	| 0.5
11	| signal_11.mat	| 27	| 35	| 1000	| 8	| 16	| ± 8	| 0.25
12	| signal_12.mat	| 25	| 39	| 1000	| 8	| 16	| ± 8	| 0.25
13	| signal_13.mat	| 25	| 39	| 1000	| 8	| 16	| ± 8	| 0.25
14	| signal_14.mat	| 26	| 36	| 1000	| 8	| 16	| ± 2	| 0.0625
15	| signal_15.mat	| 34	| 39	| 1000	| 8	| 16	| ± 2	| 0.0625
16	| signal_16.mat	| 34	| 39	| 1000	| 8	| 16	| ± 2	| 0.0625
17	| signal_17.mat	| 32	| 37	| 500	| 8	| 12	| ± 1	| 0.5
18	| signal_18.mat	| 31	| 38	| 500	| 8	| 16	| ± 8	| 0.25
19	| signal_19.mat	| 33	| 38	| 500	| 8	| 16	| ± 8	| 0.25
20	| signal_20.mat	| 33	| 38	| 500	| 8	| 16	| ± 8	| 0.25
