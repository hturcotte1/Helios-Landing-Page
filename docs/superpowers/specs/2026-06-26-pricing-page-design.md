# Pricing Page — Cost-Arbitrage Reframe

**Date:** 2026-06-26
**Status:** Approved design, ready for implementation plan
**Scope:** New standalone `pricing.html`, then rework the existing `#pricing` section on `index.html`. No other homepage copy changes. Guides untouched.

## Goal

Sell Helios on resource arbitrage instead of leading with the price. Anchor first on what a grower already spends pumping water (and the share of it wasted by over-application), then make Helios at ~$10/acre read as a rounding error against that. Reinforce with speed: Helios forecasts 24-72 hours ahead and answers instantly, versus waiting on a consultant.

Two framings run through the page and the homepage pricing section:
1. **Cost arbitrage** — anchor on the old/ongoing cost (wasted water + pump energy), then reveal the price as trivial by comparison.
2. **Speed** — lead time (forecast ahead) and instant answer (no waiting on a human).

**Hard constraint:** never use em dashes anywhere in copy.

## Verified anchor numbers (cited)

All figures are conservative and sourced. The page states the basis of any dollar figure so nothing is asserted bare.

- Pumping cost: ~$48/acre (small farms) to ~$68/acre (large). Source: USDA ERS, Western Irrigated Agriculture.
- Purchased water: $26-$71/acre depending on farm size; plus ~$20/acre maintenance. Same source.
- Total water-related spend lands ~$50-90/acre/year.
- Average Idaho farm ~505 acres; 3.17M irrigated acres statewide. Source: USDA NASS Idaho Overview.
- Waste is real: a 4-state study found ~25% of pumping energy wasted to poor efficiency; a Kansas study found systems burn ~40% more fuel than properly tuned. Over-application of ~20-30% is the standard extension figure. Source: ATTRA, Energy Saving Tips for Irrigators.

**Headline number:** "$20,000+ spent on water each year" is conservative for a representative ~300-acre irrigated operation (300 × ~$65/acre ≈ $19,500); a 505-acre Idaho farm runs well past $30k. The page labels this basis next to the number.

Sources:
- https://www.ers.usda.gov/amber-waves/2013/september/western-irrigated-agriculture-production-value-water-use-costs-and-technology-vary-by-farm-size
- https://www.nass.usda.gov/Quick_Stats/Ag_Overview/stateOverview.php?state=IDAHO
- https://attra.ncat.org/publication/energy-saving-tips-for-irrigators/

## Approach

Approach A: static anchor narrative, with a small set of scroll-triggered animated widgets in the anchor section. No calculator, no user input. Fully static single file (inline CSS + JS), no build step, matching `index.html`.

## Design system (match the existing site exactly)

The new page is built from the same tokens and components as `index.html` so it is visually indistinguishable from the rest of the site.

- Fonts: `--font` = Geist; `--mono` = Geist Mono (already loaded via Google Fonts). All numeric/data values render in `--mono`.
- Color tokens (reuse `:root` from `index.html`): `--ink #171717`, `--body #4d4d4d`, `--mute #888`, `--hairline #ebebeb`, `--canvas-soft #fafafa` (page bg), `--link`/`--water #0070f3`, `--wait #f5a623`. Dark band uses `--app-bg #0a0b0d` / `--app-text #ededee`, same surface as the homepage product mockup and `.dark-panel`.
- Layout: `.wrap` (max 1100px), section rhythm via `--pad`, `.eyebrow` mono uppercase labels, `.h-hero`/`.h-lg`/`.h-md`/`.lead` type scale, `.btn`/`.btn-primary`/`.btn-secondary`, `.reveal` scroll-in animation, fixed `nav` with `.scrolled` state and `#progress` bar.
- Reused markup, copied verbatim from `index.html`: the `<nav>` block and the footer block (so brand logo, links, and behavior stay identical).
- Motion honors `prefers-reduced-motion`: animated counters and bars render their final values statically when reduced motion is requested (same guard the hero canvas already uses).

## Page structure: `pricing.html`

1. **Nav** — reused from `index.html`, plus one added "Pricing" link (navigation only).
2. **Anchor hero (pain first).** Eyebrow "Pricing". Headline that names the status quo cost, not the price (e.g. "You are already paying for water the crop never used."). Subhead states the cited waste reality. The price is deliberately absent here.
3. **Animated arbitrage block.**
   - Spin-up anchor counter: on scroll into view, a Geist Mono counter animates 0 -> "$20,000+", with a small caption stating the basis (~300 irrigated acres at USDA per-acre water costs, cited).
   - Waste reveal: a second smaller counter animates to the wasted slice (~$4,000-6,000, the 20-30% over-application), labeled "watering ground the crop never needed".
   - Arbitrage bar: an animated bar comparison, a tall water-spend bar versus a barely-visible Helios sliver (~$10/acre), driving the "rounding error" read. A one-line caption makes the point explicit.
4. **Speed band (dark panel, `--app-bg`).** Two beats: forecast 24-72 hours ahead (act before the field forces it) and an instant Water/Wait call the moment you open it (no waiting on a consultant's visit). Framed as: cheaper and faster.
5. **What's included.** The existing `.price-feats` list (24-72 hr forecasts, confidence-scored Water/Wait calls, dashboard on any device, email and chat support), ~$10/acre, per growing season, no contracts.
6. **CTA.** Join the waitlist, free access at launch, pricing locked in for good. Links to the homepage `#waitlist` form (single source of truth for the form).
7. **Footer** — reused from `index.html`.

### Animated widgets (technical)

- Vanilla JS, no dependencies. `IntersectionObserver` triggers each widget when it scrolls into view (same pattern as `.reveal`).
- Count-up: `requestAnimationFrame` easing from 0 to target over ~1.2s, formatting with thousands separators and the `$`/`+` affixes; value sits in a `--mono` element.
- Arbitrage bar: CSS width/height transition kicked by an `.in` class toggled by the observer.
- `prefers-reduced-motion: reduce` short-circuits all of the above to final values.
- All animation code is inline in `pricing.html`, consistent with the single-file pattern.

## Homepage change: `index.html` `#pricing` section

Rework the existing pricing section to a condensed version of the same framing:
- Lead with the anchor line (the status quo cost) instead of "Simple, honest pricing."
- A compact version of the arbitrage comparison (no full animated suite required; may reuse the bar or a static comparison line).
- Keep ~$10/acre and the feature list.
- Add a "See full pricing ->" link to `pricing.html`.
- Keep the existing waitlist CTA.

The hero and every other section of `index.html` stay exactly as they are.

## Plumbing

- Add `pricing.html` to `sitemap.xml`.
- Add the "Pricing" nav link on `index.html` and `pricing.html`.
- Page `<title>`, meta description, canonical, and OG/Twitter tags follow the same format as `index.html`, pointed at `https://irrigant.xyz/pricing.html`.

## Out of scope

- No calculator or user input.
- No changes to guide pages, the hero, or other homepage copy.
- No new dependencies or build step.

## Verification

- Serve locally; screenshot `pricing.html` and the homepage `#pricing` section at desktop and mobile widths.
- Confirm: no em dashes anywhere; all numbers render in Geist Mono; counters animate on scroll and settle on correct values; reduced-motion shows static finals; nav/footer match the rest of the site; the "Pricing" link and "See full pricing" link resolve; cited claims match the sources above.
- Confirm visual parity with `DESIGN.md` and `index.html` (color, type scale, spacing, motion).
