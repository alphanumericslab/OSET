function phase = PhaseShifting(phasein, teta)
% PhaseShifting has been deprecated. Use phase_shifter instead.
warning('PhaseShifting has been deprecated. Use phase_shifter instead.');
phase = phase_shifter(phasein, teta);