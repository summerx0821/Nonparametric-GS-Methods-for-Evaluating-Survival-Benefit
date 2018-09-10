# User's guide of R functions for "Nonparametric Group Sequential Methods for Evaluating Survival Benefit from Multiple Short-Term Follow-up Windows"
We provide guidelines for using the XiaMurrayTayob.gstest function that calculates critical values for group sequential monitoring of the test proposed in Nonparametric Group Sequential Methods for Evaluating Survival Benefit from Multiple Short-Term Follow-up Windows by Xia, Murray and Tayob.  A simulated data example is given so that correct use of this function can be easily verified.  As output, the function provides the calculated test statistics at each analysis time (\$standardized.test.statistics), the estimated correlation between test statistics at the different interim analyses (\$correlation.of.statistics), OF efficacy critical values (\$efficacy.critical.value) and safety critical values (\$safety.critical.value) as specified by a user supplied argument  and  indicators that the trial should be stopped for efficacy ($efficacy.stop) or safety (\$safety.stop).

### Usage:
XiaMurrayTayob.gstest(X, delta, E, group, s, Tau, alpha.efficacy=0.025, alpha.safety=0.025, safety.bound="JT", overall.alpha.safety.JT=0.2)

### Arguments:
* X: a vector of the observed follow up times at the current analysis time, s, for the individuals in the study (`$X_i(s), i=1,\dots,n(s)$` from Xia, Murray and Tayob reference) $\sum_{i=1}^n X_i$
* delta: a vector of the status indicators at the current analysis time, $s$, for the individuals in the study,  1=event, 0=otherwise ($\delta_i(s), i=1,\dots,n(s)$ from Xia, Murray and Tayob reference)
* E: a vector of individual entry times; this is on the same scale with the analysis time, s ($E_i, i=1,\dots,n(s)$ from Xia, Murray and Tayob reference)
* group: a vector of the treatment group indicators for individuals in the study, group=1 or 2 (X, delta, E and group must have the same length)
* s: a vector of interim analysis times ($s_1, s_2, \dots$ from Xia, Murray and Tayob reference)
* Tau: desired short-term follow-up window length ($\tau$ from Xia, Murray and Tayob reference)
* alpha.efficacy: false positive clinical trial rate when null hypothesis is true 
* alpha.safety: overall rate of stopping incorrectly for safety when the null hypothesis is true
* safety.bound: a choice of "OF" (O'Brien and Fleming), "Pocock" or "JT" (the modified Jennison and Turnbull bound proposed by Xia, Murray and Tayob)
* overall.alpha.safety.JT: this argument is only used when safety.bound="JT". A user-specified overall error rate for exceeding the safety boundary under the null hypothesis and incorrectly stopping the trial as a result, namely,  ($\alpha_{safety}$ described in Section 4.2 of Xia, Murray and Tayob)
