# Wang model foundation — 2026-11-12

Investigator supplied own 2016 paper and explicitly requested using that model for moving celestial patterns and replacing phenomenological turbidity depolarization with unpolarized ALAN where appropriate. Read complete main text and visually inspected formula page. Eq7 exponent verified as -T/4+.12, not (-T+.12)/4. Main text refers to supplement for spectral coefficient details; publisher access/search did not recover it. No supplementary coefficients inferred from memory.

Implemented equation-traceable kernel, explicit full-eq17 input branch, reduced eq6+eq8 monochromatic pilot, I-only ALAN addition, conditional turbidity equivalence, fixed 30deg-radius heading response, optional explicit-UTC Astropy lunar geometry. Details and run command in docs/wang2016_dynamic_spider_model.md. Neutral separation must be supplied explicitly (24deg used only in synthetic demo). No empirical timestamp presumed exact, no real image fit and no HPS/LED spectral conversion claimed.

Eight kernel checks pass; elevation-sweep pilots execute successfully both with baseline T3 and with natural amplitude=1 plus ALAN. Pilot PNG visually inspected. Astropy lunar runtime not installed/tested here; API checked against official documentation. Model source positions only above horizon. Uniform relative intensity and solar-to-lunar extension remain explicit limitations. Full spectral/empirical model is not complete until supplementary inputs and capture metadata are verified.

Outputs are separate; locked references/fits untouched. Code, tests, dependency file and guide committed to publication-figures. User can run reduced pilot now and supply supplement/original model code to complete the spectral branch.
