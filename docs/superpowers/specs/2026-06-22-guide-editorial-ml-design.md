# Irrigant Editorial Guides and Machine Learning Page

## Goal

Make the guides feel like a natural extension of the Irrigant landing page, then replace the potato page with a technical Machine Learning page for investors and crop scientists.

## Guide layout

- Remove the boxed-card treatment from every guide.
- Use an editorial flow: eyebrow, large heading, lead, full-width content sections, thin dividers, and deliberate vertical spacing.
- Preserve the existing Geist typography, off-white canvas, hairline rules, black pill CTA, and responsive navigation.
- Keep the grower guides plain-language and practical. They remain focused on irrigation scheduling, soil-moisture forecasting, and how Helios helps.

## Machine Learning page

- Replace `potato-irrigation-guide.html` with `machine-learning.html`.
- Update sitemap, homepage guide links, guide navigation, titles, descriptions, canonical URLs, and Open Graph metadata.
- Target investors and crop scientists. Use technical language where it explains the system accurately.
- Describe the current Helios stack: soil-moisture, crop, weather, ET, and operating-constraint inputs; XGBoost multi-output forecasts at 24, 48, and 72 hours; OpenET monthly ET enrichment; NOAA weather backfill; cache and fallback behavior; rule-based irrigation timing/depth recommendations; FastAPI, static frontend, API endpoints, and nearby-feedback loop.
- Do not name unsupported metrics, field-validation results, or commercial outcomes.
- Do not foreground prototype limitations on this page. The copy must stay accurate without positioning the page as a limitations report.

## Verification

- Extend the static SEO test so the new Machine Learning URL is indexed and the removed potato URL is absent from sitemap and navigation.
- Verify every guide has no boxed section class, retains its canonical and social metadata, and links to the waitlist.
- Run local desktop/mobile browser QA, then live production verification after deployment.
