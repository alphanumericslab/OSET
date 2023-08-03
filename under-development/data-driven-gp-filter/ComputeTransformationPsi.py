# Python implementation for computing the beat-dependent transformation matrices corresponding to an ECG signal x
# Inputs: x - the ECG
#         pshb - the R-peak positions
#         TauFirseHalf - beats length in the phase domain for the interval between the first sample and R-peak
#         TauSecondHalf - beats length in the phase domain for the interval between the R-peak and the last sample
# Outputs: Psi - transformation matrices


import numpy as np
from ComputeHeartbeatPsi import ComputeHeartbeatPsi


def ComputeTransformationPsi(x, pshb, TauFirseHalf, TauSecondHalf):
    Psi = []
    for i, hb in enumerate(x):
        # print(i, len(hb[0:pshb[i]]), int(Tau/2), int(Tau/2) - len(hb[0:pshb[i]]))
        Ps1 = ComputeHeartbeatPsi(len(hb[0 : pshb[i]]), TauFirseHalf)
        # ----- ---- --- -- -
        # print(i, len(hb[pshb[i]:]), Tau-int(Tau/2), Tau - int(Tau/2) - len(hb[pshb[i]:]))
        Ps2 = ComputeHeartbeatPsi(len(hb[pshb[i] :]), TauSecondHalf)
        # ----- ---- --- -- -
        Ps = np.zeros((Ps1.shape[0] + Ps2.shape[0], Ps1.shape[1] + Ps2.shape[1]))
        Ps[0 : Ps1.shape[0], 0 : Ps1.shape[1]] = Ps1
        Ps[
            Ps1.shape[0] : Ps1.shape[0] + Ps2.shape[0],
            Ps1.shape[1] : Ps1.shape[1] + Ps2.shape[1],
        ] = Ps2
        Psi.append(Ps)
    return Psi
