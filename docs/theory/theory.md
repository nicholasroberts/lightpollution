# Theory — spectral sky polarization and receptor contrast

Technical working reference for Paper 1 and subsequent historical scenarios.
This documents the current mathematical framework and distinguishes adopted
approximations from analyses still to implement. It is not manuscript prose.

## 1. Coordinates, geometry and Stokes conventions

Represent each sky direction by a unit vector v. Source directions s_M and s_S
come from fixed lunar and solar ephemerides at the capture location/time.
Image-to-sky registration is established independently using landmarks and
observed Moon position. Camera parity, analyzer sequence sign, analyzer zero
and the publication display transform are separate decisions.

All Q and U contributions must use the same local polarization basis. With
complex linear polarization Z = Q + iU and angle psi,

\[
Z=I p\,e^{2i\psi}, \qquad p=\sqrt{Q^2+U^2}/I.
\]

A passive rotation of the reference axes through beta gives
\(Z'=Z e^{-2i\beta}\). A reflection requires the corresponding handedness
transformation; it is not merely an arbitrary constant angle adjustment.
The publication sky view is north up/east left. Changing the display alone
does not change the inferred physical source geometry.

## 2. Natural solar and lunar contributions

For incoherent sources, add spectral Stokes vectors:

\[
\mathbf S_{\rm nat}(\lambda,\mathbf v)
=\mathbf S_M(\lambda,\mathbf v)
+\mathbf S_S(\lambda,\mathbf v),\qquad \mathbf S=(I,Q,U).
\]

Do not add DoLP values or average AoP angles. Different source polarization
directions can cancel in Q/U even when each source is strongly polarized.

The present exploratory component uses \(\mu=\mathbf v\cdot\mathbf s\):

\[
I_s \propto 1+\mu^2,\qquad
|Z_s|\propto p_{\max}(1-\mu^2).
\]

The source intensity normalization in the existing fitter makes zenith I=1
for each component. Its polarization phase comes from the scattering plane
and ideal meridional transport into the camera analyzer plane. The pilot
fixes \(p_{\max}=0.8\). These ideal components are not a complete twilight
radiative-transfer model, especially with the Sun below the horizon.

The current fitted natural field is \(\mathbf S_M+r\mathbf S_S\). The source
ratio r is presently estimated with a DoLP objective and separately assessed
using spatial AoP changes. A joint AoP/DoLP objective is planned, not already
implemented. The current spectral comparison assumes a common transferred
natural spectral template for Sun and Moon; this is explicitly provisional.

## 3. Artificial-light contribution and normalization

For a lamp scenario j, define its transferred energy spectrum e_j(lambda)
over the common band B such that

\[
\int_B e_j(\lambda)\,d\lambda=1.
\]

The current unpolarized, spatially uniform skyglow approximation adds

\[
\mathbf S_j=A_j e_j(\lambda)(1,0,0).
\]

More general spatial weighting W_j(v) can replace the constant field, but is
not fitted by the present equivalence calculation. A_j is sky band radiance
relative to a specified natural reference, not electrical power or illuminance.
The pilot uses 352–700 nm and lunar zenith band radiance as its reference.

Adding unpolarized light leaves Q/U and hence AoP unchanged. It changes DoLP
and receptor contrast through I. A solar/lunar mixture can change the angular
pattern through Q/U. Keep these mechanisms distinct when discussing the band.

## 4. Spectral receptor weighting

Let R(lambda) be photon sensitivity and L_lambda an energy spectrum per nm.
Photon weighting is proportional to lambda R(lambda)/(hc). The common hc
cancels in contrast ratios; wavelength units must remain consistent.

\[
(I_R,Q_R,U_R)=\int_B {\lambda R(\lambda)\over hc}
             (I_\lambda,Q_\lambda,U_\lambda)\,d\lambda.
\]

Only after these integrations form receptor-weighted DoLP and AoP if needed.
For wavelength-dependent dichroic absorbance, compute receptor responses
inside the spectral integral, rather than applying constant absorbance
coefficients afterwards.

The investigator accepts reconstructed green-camera PRC as an approximation
to RH1 PRC. It is not a unique recovery of spectral animal Stokes from RGB.
The 510 nm alpha-band template is the current animal proxy. The measured ELF
D850 green response provides a sensitivity comparison; it excludes this
project's added polarizer transmission.

## 5. Orthogonal receptor responses and PRC

For the existing wavelength-independent absorbance approximation, let

\[
t_\parallel=10^{-a_\parallel},\quad
t_\perp=10^{-a_\perp},\quad
b=1-\frac{t_\parallel+t_\perp}{2},\quad
g=\frac{t_\perp-t_\parallel}{2}.
\]

The pilot uses \(a_\parallel=0.4\), \(a_\perp=0.04\). A receptor axis theta
has absorbed response \(bI_R+g[Q_R\cos2\theta+U_R\sin2\theta]\); the orthogonal
receptor reverses the modulation. Equal overall gain cancels in the ratio.

With a declared field support and spatial weights w(v), integrate the Stokes
components first, obtaining T_I, T_Q and T_U. Then

\[
\mathrm{PRC}(\theta)
=\ln\frac{bT_I+g(T_Q\cos2\theta+T_U\sin2\theta)}
          {bT_I-g(T_Q\cos2\theta+T_U\sin2\theta)}.
\]

This is the current common camera-plane/ideal transport proxy. It is not a
complete calibration of spider retinal analyzer transport. Equal included
pixels correspond to the adopted provisional equisolid image sampling.

For \(k=g|T_Q+iT_U|/(bT_I)\),

\[
\max_\theta|\mathrm{PRC}|=2\,\mathrm{atanh}(k),\quad
\mathrm{median}_{\theta\ {\rm uniform}}|\mathrm{PRC}|
=2\,\mathrm{atanh}(k/\sqrt2).
\]

The continuous-orientation median is not the same statistic as an arbitrary
finite set of headings. A constant analyzer-zero rotation changes phase but
not these two amplitudes. Maximum, heading median and specified compass
heading are different endpoints and must not be interchanged.

## 6. Fitting equal PRC under alternative lamp spectra

Write the integrated response mean as M_0 + A_j H_j, and the fixed natural
modulation amplitude as P_0. H_j includes the spectral response, spatial
support and baseline absorbance coefficient for unit added band energy.

\[
C_j(A_j)=2\,\mathrm{atanh}\left({P_0\over M_0+A_jH_j}\right).
\]

For a positive target C_t,

\[
A_j={P_0/\tanh(C_t/2)-M_0\over H_j}.
\]

A target above the natural maximum cannot be reached by adding unpolarized
light; report it as unreachable rather than fitting a negative amount.
Zero target is only approached asymptotically for nonzero natural modulation.

The equal-effect comparison is

\[
\rho={A_{\rm LED}\over A_{\rm HPS}},\qquad
\Delta_E=100(1-\rho).
\]

Under this fixed-natural-field/uniform-addition approximation, the common
numerator cancels, so \(\rho=H_{\rm HPS}/H_{\rm LED}\). Therefore the ratio is
determined by spectral overlaps, while the site target and geometry determine
absolute required amounts. Identical ratios across sites are not independent
empirical replication. At matched PRC the present HPS/LED RH1 maps coincide.

The current pilot ratio is 0.609388 (39.061% less LED band energy), conditional
on the selected HPS/mid-LED spectra, transferred spectral shapes and RH1 proxy.
The D850-green sensitivity calculation gives 10.459% at the same numerical
target; neither number is asserted to apply universally to LEDs or sites.

## 7. Equal human luminance — planned calculation

For spectral sky radiance L_lambda in W m^-2 sr^-1 nm^-1, photopic luminance is

\[
L_v=683\int_B V(\lambda)L_\lambda\,d\lambda
\]

in cd m^-2 when the integral and wavelength convention are consistent.
Use a documented V(lambda) dataset; no dataset has yet been selected in this
workflow setup. Record any restriction to a finite wavelength band.

For unit-energy lamp spectrum e_j, define
\(K_j=683\int_B V(\lambda)e_j(\lambda)d\lambda\). Equal added sky luminance
L_t requires \(A_j=L_t/K_j\) in physical radiance units. If A_j is dimensionless
relative to natural sky radiance, carry that reference conversion explicitly.
Then calculate the animal receptor Stokes/responses and PRC; do not equate
equal luminance with equal biological response.

Equal sky luminance is not equal road illuminance, lumens, fixture power or
unchanged upward flux. The intervention case requires separate evidence.

## 8. Historical interpretation

VIIRS-band radiance and animal receptor response use different spectral
weightings. Without source-spectrum information a measured VIIRS trajectory
does not uniquely specify a biological PRC trajectory. HPS-only and LED-only
histories are conditional scenarios, not proven extrema over all spectra and
atmospheres. Stable-spectrum dimming can still support useful biological
comparisons, with starting brightness affecting the nonlinear PRC change.

## 9. Scope and source pointers

This setup does not establish a universal twilight elevation threshold,
behavioral navigability threshold, absolute analyzer zero or a calibrated
site-specific HPS/LED mixture. Code and numerical assumptions take precedence
over an attractive figure layout; disagreements must be documented.

- [Paper 1 decisions](../paper1/decisions.md)
- [Existing lamp-equivalence methods and calibration attribution](../prc_lamp_equivalence.md)
- [Existing PRC implementation](../../src/compare_spider_mixed_prc.py)
- [Natural component implementation](../../src/fit_spider_twilight_moon.py)
- [Spectral/lunar helper](../../src/spider_spectral_lunar.py)
- [Two-paper research plan](../lab_notebook/2026-11-23-two-paper-plan.md)
