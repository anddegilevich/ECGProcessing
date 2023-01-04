# ECG Processing

## About repository

<p><b>ECG Processing</b> contains different methods for processing electrocardiosignal using MATLAB.</p>
<p>The following processing procedures are presented:</p>

 - [Digital filtering](#digital-filtering) 
 - [Spectrum analysis](#spectrum-analysis)
 - [Calculation of arterial blood oxygen saturation](#calculation-of-arterial-blood-oxygen-saturation)
 - [Correlation analysis](#correlation-analysis)

## Digital filtering
<p>The downloaded electrocardiogram signal samples (sampling rate = 250 Hz) were corrupted with noise (60 Hz). To filter the signal, a digital notch filter was developed in the <ins><b>Filter Designer</b></ins> utility, the AFC of which is shown in the figure:</p>

<p><img src="screenshots/screen1.JPG" /></p>

<p>The result of the signal filtering is presented below:</p>

<p><img src="screenshots/screen2.JPG"/> </p>

## Spectrum analysis
<p>Power spectral density and amplitude spectrum of the electrocardio signal were plotted using the <ins><b>Fast Fourier Transform</b></ins> algorithm.</p>
<p>The obtained spectra are shown in the figure below:</p>
<p><img src="screenshots/screen3.JPG" /></p>

## Calculation of arterial blood oxygen saturation
<p>To calculate blood oxygen saturation using the pulse oximetry signal, an algorithm for finding QRS complexes was developed. It is based on differentiation of ECG signal to hyperbolize signal changes and use a threshold rule with empirically selected coefficients.</p>
<p>Then we calculated pulse oximetry signal ranges of two wavelengths (600...700 nm and 810...960 nm) for obtained RR intervals. Calculating their average and using the formula the saturation value (~94%) was calculated.</p>
<p>The algorithm is illustrated by the figure below:</p>

<p><img src="screenshots/screen4.JPG" /></p>

## Correlation analysis
<p>Correlation analysis methods were used to find periodicities and correlations between signals:</p>
<p><img src="screenshots/screen5.JPG" /></p>
<p>The autocorrelation function of ECG shows that there are repeating segments in the signal after fixed intervals of time.</p>
<p>Cross-correlation between pulse oximetry signals shows that these signals have the same dynamics of change but are shifted relative to each other by a certain time interval.</p>
