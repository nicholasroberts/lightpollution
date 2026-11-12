# Berry normalization audit — 12 November 2026

Source: Berry, Dennis & Lee (2004), New J Phys 6 162, equation 4.2 and paragraph immediately below it. https://michaelberryphysics.wordpress.com/wp-content/uploads/2013/07/berry373.pdf ; DOI 10.1088/1367-2630/6/1/162.

Berry eq4.2 has prefactor -4 multiplying the quartic divided by the stereographic denominator. The authors explicitly choose a unity maximum of |w|. Wang (2016) printed eq3 has no -4. Our literal transcription therefore had maximum 0.25. This is a discrepancy between printed formulas, not evidence about what code generated Wang's published figures. Original Wang simulation code has not been verified.

The active model now uses berry_normalized_pattern = -4 * source_meridian_pattern. Literal complex_pattern and source_meridian_pattern remain unchanged for printed-equation audits. The minus sign rotates raw AoP 90 degrees; the measured comparison's fitted constant basis offset compensates, so fitted AoP maps should remain the same axially. Pilot uses the same corrected source-meridian kernel. Atmospheric turbidity stays excluded; ALAN is only additive unpolarized intensity. No fitted Cottage gain is used.

Numerical checks: exact -4 relation, unity magnitude at the analytical maximum, and the small neutral-separation limit |w|=sin^2(scattering angle). That limit is polarization intensity, not normalized single-scattering Rayleigh DoLP sin^2(gamma)/(1+cos^2(gamma)). Berry explicitly distinguishes unnormalized polarization intensity from radiometric DoLP below eq2.1. Consequently, treating the unity-peak Berry pattern with Wang's horizon correction as a natural DoLP envelope remains an explicit modelling assumption. This correction fixes the ceiling; it does not establish absolute lunar radiometric validity or reproduce all Wang figures, which also used spectral and atmospheric factors excluded here.

From uploaded no-turbidity summary (2), multiplication alone predicts Cottage field mean 0.8559194388 rather than 0.2139798597, versus observed 0.3135794838; Oxford field mean 0.6612658401 rather than 0.1653164600, versus observed 0.1370901847. These are algebraic predictions, not new measured runs. Old ALAN/dilution estimates are superseded. Do not report ratios of these means as established ALAN measurements.

Validation: 22 existing Wang-related tests plus 3 new Berry tests pass. The synthetic elevation-sweep pilot completed. Measured rerun needs EcoVis TIFF-derived saved arrays. Original alignment, anchors and measurements remain unchanged.
