# Irrigant SEO foundation and grower content

## Goal

Make irrigant.xyz discoverable for U.S. growers looking for irrigation scheduling and soil-moisture guidance, while positioning Helios as a confident, practical way to make better water decisions.

## Positioning

- Audience: U.S. growers and farm operators evaluating irrigation decision support.
- Primary intents: irrigation scheduling, soil-moisture forecasting, and crop-water decisions.
- Tone: direct, practical, and grower-facing. Lead with the outcome, then explain enough to earn trust.
- Conversion: waitlist sign-up.

## Claim boundaries

- Lead with the useful promise: Helios gives growers a clearer read on when to water, how much to apply, and what is driving the call.
- Support that promise with the current Helios facts: it considers field conditions, crop context, near-term weather, and crop-water demand to make a recommendation.
- Do not claim measured water savings, calibrated confidence, universal live coverage, or field-validated accuracy that the business has not earned yet.
- Do not name algorithms, feature engineering, data schemas, model versions, or private data on public pages.
- State that Helios supports grower judgment and does not autonomously control irrigation. Put this in the explanatory content, not in the hero message.

## Sales posture

- The hero and page titles answer the buyer's question first: what decision gets easier, and why does it matter now.
- Use decisive language where it is true: clearer calls, earlier warning, and better visibility into the drivers behind an irrigation decision.
- Keep caveats specific and local. A visitor should not have to read a disclaimer before understanding the value.
- Frame pilot access as an advantage: early access to a new decision-support tool, not an apology for an unfinished product.
- Avoid weak filler such as "may help," "intended to," and "potentially," unless it is necessary to avoid a claim the evidence cannot support.

## Public pages

1. Homepage: retain the visual system and waitlist flow. Update the title, main heading, and supporting copy so the page clearly says what Irrigant helps with in plain language.
2. Irrigation scheduling guide: explain the decisions a grower needs to make, common signals to consider, and where Helios is intended to help.
3. Soil-moisture guide: explain why soil moisture matters and how forecasts can support a decision without promising a result.
4. How Helios helps: a short, non-technical explanation of the decision-support workflow, with its pilot status and non-autonomous role stated clearly but without undercutting the sales message.
5. Potato irrigation guide: publish only source-backed, general guidance. It must not claim field-specific prescription capability or universal application rates.

Each guide links to the waitlist and to the other relevant guides. No state or city landing pages, no keyword-stuffed pages, and no programmatic content expansion.

## Technical SEO

- Add a canonical URL to every indexable page.
- Add a root `robots.txt` and `sitemap.xml` listing each indexable public page.
- Add concise, page-specific titles and descriptions.
- Add complete Open Graph and Twitter metadata, including a social image and canonical URL.
- Add only `Organization` and `WebSite` JSON-LD. Do not add review, rating, FAQ, LocalBusiness, or software-performance markup.
- Preserve the existing HTTPS and host redirect behavior.
- Optimize the current logo asset without changing the visual design.

## Measurement and indexing

- Create or connect Google Search Console for `irrigant.xyz`, submit the sitemap, and request indexing for the homepage and guides.
- Create or connect GA4 and record a successful waitlist form submission as a conversion.
- Review search queries, indexed pages, and conversion activity monthly. Expand content only from observed demand or defensible original material.

## Verification

- Add a test-first static check for canonical URLs, metadata, JSON-LD, `robots.txt`, sitemap coverage, and internal links.
- Verify the rendered site at desktop and mobile widths, including the waitlist flow.
- Verify the deployed URLs return `200`, the canonical host resolves correctly, and the new crawl assets are live.
- Verify Search Console sitemap acceptance and the GA4 waitlist conversion after account access is available.

## Out of scope

- Product dashboard changes.
- New model features or revised model claims.
- Local SEO or a Google Business Profile.
- Backlink campaigns, paid search, and fabricated case studies.
