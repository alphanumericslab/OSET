# Python Implementation of The Open-Source Electrophysiological Toolbox
This sub-project is a Python implementation of the OSET package

## Tools
Function | Description
---------|------------
[AdaptiveFilter](Tools/AdaptiveFilter.py) | Adaptive Noise Canceller/Line Enhancer
[KFNotch](Tools/KFNotch.py) | Linear Kalman Filter and Smoother (Using Q,R,gamma)
[KFNotch2](Tools/KFNotch2.py) | Linear Kalman Filter and Smoother (Using Qbar,Wlen)
[NSCA](Tools/NSCA.py) | Nonstationary Component Analysis Given Time Windows
[PiCA](Tools/PiCA.py) | Pseudo-Periodic Component Analysis Given One or Two Peaks Inputs
[PolyFit](Tools/PolyFit.py) | Least Squares Polynomial Fit for Signal Segment
[RLS](Tools/RLS.py) | Recursive Least Squares Filter
[SynchPhaseTimes2](Tools/SynchPhaseTimes2.py) | Finds Synchronous Times Instants Given R-peaks, for PiCA
[TikhonovRegularization](Tools/TikhonovRegularization.py) | Solves Constrained Least Squares Problem per Channel

## testPrograms
Function | Description
---------|------------
[testECGAdaptiveFilter](testPrograms/testECGAdaptiveFilter.ipynb) | Uses sample ECG data to test active noice canceling with [AdaptiveFilter](Tools/AdaptiveFilter.py)
[TestKFNotch](testPrograms/testKFNotch.ipynb) | Uses sample ECG data to test Kalman filter notch with [KFNotch](Tools/KFNotch.py)