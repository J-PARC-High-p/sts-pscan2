======================================================

pscan analysis

   original by Osnan et al.(?)
   modified by K. Aoki
======================================================
pscan_samples/
	pscan files created at GSI.
pscan/
	pscan files created at KEK.

COMMENTS:
	execute() in execute.C will execute the process.
	if inputfilename.root exists, just open the root file. Do not invoke the procedure.
	filename should be without ".txt" as it will be added later in the class.


	vcnt[ch][d][ivp] : count of ADC counters.
			 ch: channel #. (0-127)
			 d: ADC count. 0-31
			 ivp: pulse height. (0-255)
		 
		 if vcnt exceeds the maximum (cut_db_pulses) vcnt is made to be cut_db_pulses.
	
	inside Analysis()
	       Soft_val(true)
	       		start with ivp=2;
			basically copy vcn into vnc_soft.
	
		d_cnt is the difference of [ivp]and[ivp-1]. vnc_soft[ch][d][ivp]-vcn_soft[ch][d][ivp-1];
		hdcnt[ch][d] is the difference (d_cnt)
		hscurve[ch][d] is the vnc_soft itself.
		
		Fit_values (d<=30)
			   Fit hdcnt
			       find maximum bin. fit range: maximum-bin +- width. width is 35.specified in execute.C
			       for d>28, fit range is (0,80) why?
			   Fit option: SWLQR
			       S: results of the fit is returned in the TFitResultPTr
			       W: ignore the bin uncertainties. skip empty bin.
			       L: log likelihood
			       Q: Quite mode
			       R: Use the range specified in the function range. (If this is true, I think R should not be used!)
			       WL: Use loglikelihood method. histogram is weighted. 
			   make hmeanf 2D histo. ch vs d.
			   make hsigef 2D histo. ch vs d.
		Fit_values_erfc (d<=30)
		Fitting_Fast (only for d = 31)
		
		
		