# Normalization audit before revising the biological sky model

## Aim and decision

The goal remains to predict photon-weighted polarization receptor contrast (PRC)
within the animal's fixed 60deg full field (30deg radius), for a horizontal body
rotating through 360deg, as lunar position/brightness changes overnight and HPS
or LED skyglow is added. The unknown species response is currently represented
by a green RH1-type alpha-band pigment with lambda_max=510nm. This is a proxy,
not a species calibration. Actual historical site radiances and empirical PRC
remain independent checks, not targets used to choose sky geometry.

**No production model, ALAN setting, atmospheric factor, PRC series, geometry or
AoP convention was changed by this audit.** An isolated ideal-limit script and
four tests were added. Earlier pilot outputs remain scenarios conditional on
an assumed natural DoLP envelope; their width mismatch is unresolved.

## Sources inspected

- Berry, Dennis & Lee (2004), *Polarization singularities in the clear sky*,
  New J Phys 6 162, DOI 10.1088/1367-2630/6/1/162.
  Author archive: https://michaelberryphysics.wordpress.com/wp-content/uploads/2013/07/berry373.pdf
- Investigator-supplied published Wang et al. (2016), J Opt 18 065601,
  DOI 10.1088/2040-8978/18/6/065601, and supplementary equations.
  Author manuscript: https://research-information.bris.ac.uk/ws/files/69774662/Revised_Manuscript_final_version.pdf

Berry equation 2.1 defines the unnormalized complex Stokes quantity Q+iU;
its magnitude is polarized intensity. The text below explicitly distinguishes
DoLP, obtained by dividing by total intensity. Equation 4.2 has maximum magnitude
one; equation 4.3 is |sin(theta_P+) sin(theta_P-)|. Unity-peak scaling is not
normalization by total intensity. Figure 7 compares polarized-intensity maps
with a horizon function fitted to one measured sky (equation 6.1), then applied
to a second. It supplies no quantitative FWHM comparison.

Wang equation 1 likewise introduces an unnormalized complex Stokes parameter,
but calls its modulus DoP. Equations 3 and 6 then use it as a phenomenological
DoP pattern. Equation 16 introduces a direction- and wavelength-dependent
inverse-spectral-power correction, involving SD, S90 and Ssun. Equation 17
combines that term with the angular horizon factor and turbidity multiplier:

DoP = |w| * [theta E(theta) + (pi/2 - theta) S(gamma,theta,lambda)] * M(T).

This is not an explicit physical Q,U solution followed by a separately supplied
I division. The paper validates the full empirical model against libRadtran;
that validation does not transfer automatically to our reduced, lunar pilot.
The supplement provides Perez/Preetham coefficients, not a separate missing
physical Stokes denominator. Previously identified coefficient/convention issues
remain separate from this normalization audit.

**Qualification to the preceding discussion:** it is premature to label this a
missing-denominator programming error. Wang intentionally reuses the geometry
as DoP. A Rayleigh denominator would change that phenomenological model; it is
not simply restoring a division that our implementation accidentally omitted.

## Equation-to-code trace

| Stage | Implementation | Meaning / audit result |
|---|---|---|
| Literal printed Wang geometry | wang2016_sky.py: complex_pattern | Archival equation form |
| Selected geometric phase/magnitude | source_meridian_pattern, berry_normalized_pattern | Explicit rotation convention and Berry -4 prefactor; magnitude interpreted as unit-peak DoLP envelope |
| Reduced natural Stokes | natural_stokes, mode eq6_horizon | p=|W| H eta, then Q+iU=I p exp(2i psi) |
| Full spectral equation audit | natural_stokes mode eq17, wang_spectral_inputs | Archived equation-17 expression, not the production nightly model |
| Natural lunar intensity | spider_spectral_lunar.py: ks91_brightness | KS91 V-band photometric brightness; normalized to a reference and transferred to fixed spectral shape |
| Spectral snapshot | plot_spider_spectral_pilot.py: sky_dop_snapshot | p I_natural times sensitivity-weighted polarized spectrum, divided by natural+ALAN photon catch |
| Field integration | model_spider_nightly_prc.py: integrated_pattern | Intensity-weighted projected Q/U, with analyzer orientation transported over field |
| Receptor contrast | response_contrast | Natural+ALAN mean response and orthogonal polarization modulation; log response ratio |

The downstream Stokes mixing is internally consistent **conditional on p being
a legitimate natural DoLP**. Natural-only snapshot DoLP does not change merely
by replacing one common I with another, because I appears in both numerator
and denominator. Adding another arbitrary I division would double-normalize
that chain. The unresolved assumption is the angular law assigned to p,
including omission of Wang's angular spectral term, rather than the final
Stokes-to-DoLP arithmetic.

The production horizon factor is Wang cos(theta)^(1/10), not Berry's fitted
polynomial in stereographic radius. Production spectra are fixed-shape lunar
and lamp templates with prescribed transfer; they do not implement equation 17.
KS91 is a brightness model, not a full polarized radiative-transfer solution.
It cannot by itself fix the missing physical relationship between the assumed
polarized amplitude and natural total intensity.

## Ideal-limit diagnostic

With coincident neutral points, Berry's geometric magnitude tends to sin²gamma.
The new script verifies this numerically using separation=.001deg. In ideal
single Rayleigh scattering, consistent relative quantities are:

I = 1+cos²gamma; |Q+iU| = sin²gamma;
p = sin²gamma/(1+cos²gamma).

Analytic FWHM values are 90deg for sin²gamma and 70.5288deg for p. These are
ideal limits without horizon attenuation; they are not predictions of the
full Wang model or fits to the three camera datasets. The observed 62–68deg
widths motivate checking this normalization, but do not prove that Rayleigh is
the appropriate real-atmosphere replacement.

For uniform absolute unpolarized addition a in these same ideal intensity
units, p_mix=sin²gamma/(1+cos²gamma+a). Its peak-normalized FWHM is
2 asin(sqrt((1+a)/(3+2a))). It **broadens**, rather than narrows, with a:
70.53deg at a=0 towards 90deg as a increases without bound. This corrects any
impression from the earlier discussion that additive ALAN necessarily narrows
a band. Direction of the effect depends on the natural intensity distribution.
A constant eta or constant ALAN/natural ratio only rescales and preserves width.

## Run the isolated check

```bash
git pull origin publication-figures
source .venv/bin/activate
python src/audit_spider_intensity_normalization.py \
  --output data/processed/spider_intensity_normalization_audit_01
```

This creates intensity_normalization_audit.png, ideal_profiles.csv and audit.json
in one folder. No measured image is read or fitted. Four tests verify the
Berry limit, Rayleigh half maximum, intensity-scale invariance and unpolarized
Stokes addition. The existing animal-model tests are also run for regression.

## Short route back to the biological aim

1. Establish one consistent natural I,Q,U prescription. Keep Wang's successful
   AoP/singularity geometry; explicitly compare its phenomenological DoLP envelope
   with a physically normalized candidate. Do not blindly divide by KS91 or
   append an arbitrary power chosen to match widths.
2. Evaluate candidate angular DoLP using the same measured sky pixels, masks,
   corrected Moon positions and binning. Check widths and amplitudes separately
   across Cottage, Skillinge and Slabhouse; reserve Oxford as a further check.
3. Only once the natural angular law is acceptable, rerun the existing spectral
   HPS/LED additions and overnight heading analysis. Add unpolarized light to I,
   retain Q/U, integrate spectral photon catches over the fixed animal field,
   and report PRC ranges across headings and time.
4. Retain the actual radiance context and explicit uncertainty in atmosphere,
   lamp mix and upward-radiance-to-ground-skyglow transfer. Do not re-fit natural
   geometry or DoLP to historical PRC to manufacture agreement.

This is a focused normalization check supporting the PRC prediction, not a new
atmospheric radiative-transfer project or a change to the empirical Europe
radiance–PRC analyses.
