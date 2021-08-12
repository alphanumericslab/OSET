# Python Implementation of The Open-Source Electrophysiological Toolbox
This sub-project is a Python implementation of the OSET package

## Tools
Function | Description
---------|------------
[AdaptiveFilter]() | Adaptive Noise Canceller/Line Enhancer
[KFNotch]() | Linear Kalman Filter and Smoother (Using Q,R,gamma)
[KFNotch2]() | Linear Kalman Filter and Smoother (Using Qbar,Wlen)
[NSCA]() | Nonstationary Component Analysis Given Time Windows
[PiCA]() | Pseudo-Periodic Component Analysis Given One or Two Peaks Inputs
[PolyFit]() | Least Squares Polynomial Fit for Signal Segment
[RLS]() | Recursive Least Squares Filter
[SynchPhaseTimes2]() | Finds Synchronous Times Instants Given R-peaks, for PiCA
[TikhonovRegularization]() | Solves Constrained Least Squares Problem per Channel

## testPrograms
Function | Description
---------|------------
[testECGAdaptiveFilter]() | Uses sample ECG data to test active noice canceling with [AdaptiveFilter]()
[TestKFNotch]() | Uses sample ECG data to test Kalman filter notch with [KFNotch]()