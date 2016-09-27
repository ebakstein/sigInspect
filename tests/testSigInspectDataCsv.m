% testSigInspectDataCsv
% E. Bakstein 2015-07-07
intf=sigInspectDataCsv('artifacts/sigInspect/csvDemo/');
intf.settings.NORMALIZE_SIGNAL_PER_CHANNEL = 0;
sigInspect(intf)
