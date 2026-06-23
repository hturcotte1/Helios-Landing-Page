# Irrigant SEO Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Make irrigant.xyz crawlable, persuasive for U.S. growers, and measurable without unsupported product claims.

**Architecture:** Keep the static site. Add a Bun static-content test, root crawl assets, one shared stylesheet for four focused guides, and stronger homepage metadata and links. The homepage remains the conversion page. Guides answer grower questions and point into the waitlist.

**Tech Stack:** Static HTML, CSS, vanilla JavaScript, Bun test runner, Vercel, Google Search Console, Google Analytics 4.

## Global Constraints

- Use plain-language, grower-facing copy. Lead with the outcome, then explain enough to earn trust.
- Lead with clearer irrigation decisions, earlier warning, and visibility into the call drivers.
- Do not claim measured water savings, calibrated confidence, universal live coverage, or field-validated accuracy.
- Do not name algorithms, feature engineering, private data, or model versions on public pages.
- State that Helios supports grower judgment and does not autonomously control irrigation only in explanatory content.
- Preserve the existing visual system and Formspree waitlist endpoint.
- Do not add city/state pages, review/rating/FAQ/LocalBusiness markup, or performance markup.

---

## File Structure

- Create: tests/seo-static.test.mjs, static SEO regression tests.
- Create: robots.txt and sitemap.xml, crawler discovery assets.
- Create: guides.css, shared responsive guide styling.
- Create: irrigation-scheduling.html, soil-moisture-forecasting.html, how-helios-helps.html, and potato-irrigation-guide.html.
- Create: assets/irrigant-social-card.png, 1200 × 630 social-preview image.
- Modify: index.html, homepage metadata, structured data, sales copy, and guide links.

## Task 1: Define the static SEO contract

**Files:**
- Create: tests/seo-static.test.mjs

**Interfaces:**
- Consumes: root-level static HTML files, robots.txt, and sitemap.xml.
- Produces: bun test tests/seo-static.test.mjs.

- [ ] **Step 1: Write the failing test**

~~~js
import { describe, expect, test } from "bun:test";
import { existsSync, readFileSync } from "node:fs";

const pages = [
  "index.html",
  "irrigation-scheduling.html",
  "soil-moisture-forecasting.html",
  "how-helios-helps.html",
  "potato-irrigation-guide.html",
];

const read = (path) => existsSync(path) ? readFileSync(path, "utf8") : "";

describe("public SEO assets", () => {
  test("publishes crawler discovery for every public page", () => {
    expect(read("robots.txt")).toContain("Sitemap: https://irrigant.xyz/sitemap.xml");
    const sitemap = read("sitemap.xml");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/irrigation-scheduling.html</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/soil-moisture-forecasting.html</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/how-helios-helps.html</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/potato-irrigation-guide.html</loc>");
  });

  test("gives every public page canonical and social metadata", () => {
    for (const page of pages) {
      const html = read(page);
      expect(html).toMatch(/<link rel="canonical" href="https:\/\/irrigant\.xyz\/.+">/);
      expect(html).toMatch(/<meta name="description" content="[^"]+">/);
      expect(html).toContain("https://irrigant.xyz/assets/irrigant-social-card.png");
      expect(html).toContain('<meta name="twitter:card" content="summary_large_image">');
    }
  });

  test("keeps unsupported technical claims out of public copy", () => {
    const copy = pages.map(read).join("\n");
    expect(copy).toContain("Join the waitlist");
    expect(copy).not.toMatch(/XGBoost|LightGBM|calibrated confidence|guaranteed water savings/i);
  });
});
~~~

- [ ] **Step 2: Run the test to verify it fails**

Run: bun test tests/seo-static.test.mjs

Expected: FAIL because the crawl assets and guide pages do not exist.

- [ ] **Step 3: Commit the failing test**

~~~bash
git add tests/seo-static.test.mjs
git commit -m "test: define Irrigant SEO crawl checks"
~~~

## Task 2: Build the crawl foundation and homepage positioning

**Files:**
- Create: robots.txt
- Create: sitemap.xml
- Create: assets/irrigant-social-card.png
- Modify: index.html lines 1-20, hero at lines 297-306, and nav at lines 283-295

**Interfaces:**
- Consumes: Task 1 test.
- Produces: metadata, crawl assets, Organization/WebSite JSON-LD, and a guide entry point.

- [ ] **Step 1: Extend the failing test**

Add this assertion:

~~~js
test("positions the homepage around irrigation decisions", () => {
  const html = read("index.html");
  expect(html).toContain("Irrigation Scheduling & Soil-Moisture Forecasting | Irrigant");
  expect(html).toContain("Know when to water. And how much.");
  expect(html).toContain('"@type":"Organization"');
  expect(html).toContain('"@type":"WebSite"');
  expect(html).toContain('href="irrigation-scheduling.html"');
});
~~~

- [ ] **Step 2: Run the test to verify it fails**

Run: bun test tests/seo-static.test.mjs

Expected: FAIL because the homepage still has its old title, copy, JSON-LD, and navigation.

- [ ] **Step 3: Implement the minimal crawl and homepage changes**

Create robots.txt:

~~~text
User-agent: *
Allow: /

Sitemap: https://irrigant.xyz/sitemap.xml
~~~

Create sitemap.xml with exactly these URLs:

~~~xml
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
  <url><loc>https://irrigant.xyz/</loc></url>
  <url><loc>https://irrigant.xyz/irrigation-scheduling.html</loc></url>
  <url><loc>https://irrigant.xyz/soil-moisture-forecasting.html</loc></url>
  <url><loc>https://irrigant.xyz/how-helios-helps.html</loc></url>
  <url><loc>https://irrigant.xyz/potato-irrigation-guide.html</loc></url>
</urlset>
~~~

Generate assets/irrigant-social-card.png at 1200 × 630 using the current off-white, black, blue, cyan, and green visual system. Include the Irrigant wordmark and only: “Know when to water. And how much.”

Replace homepage metadata with:

~~~html
<title>Irrigation Scheduling & Soil-Moisture Forecasting | Irrigant</title>
<meta name="description" content="Know when to water and how much to apply. Irrigant gives growers a clearer read on soil moisture, weather, and crop water needs.">
<link rel="canonical" href="https://irrigant.xyz/">
<meta property="og:title" content="Know when to water. And how much. | Irrigant">
<meta property="og:description" content="Clearer irrigation decisions, built for growers.">
<meta property="og:type" content="website">
<meta property="og:url" content="https://irrigant.xyz/">
<meta property="og:image" content="https://irrigant.xyz/assets/irrigant-social-card.png">
<meta name="twitter:card" content="summary_large_image">
<meta name="twitter:title" content="Know when to water. And how much. | Irrigant">
<meta name="twitter:description" content="Clearer irrigation decisions, built for growers.">
<meta name="twitter:image" content="https://irrigant.xyz/assets/irrigant-social-card.png">
~~~

Add this JSON-LD before closing head:

~~~html
<script type="application/ld+json">{"@context":"https://schema.org","@graph":[{"@type":"Organization","name":"Irrigant","url":"https://irrigant.xyz/","logo":"https://irrigant.xyz/assets/logos/IrrigantLIGHT.png"},{"@type":"WebSite","name":"Irrigant","url":"https://irrigant.xyz/"}]}</script>
~~~

Use this hero copy:

~~~html
<h1 class="h-hero">Know when to water.<br><span class="em">And how much.</span></h1>
<p class="lead">Irrigant gives growers a clearer read on soil moisture, weather, and crop water needs, so every irrigation decision starts with better information.</p>
~~~

Add a Guides nav link to irrigation-scheduling.html and a four-link guide row after the existing feature section.

- [ ] **Step 4: Run the test to verify it passes as far as possible**

Run: bun test tests/seo-static.test.mjs

Expected: homepage test passes; the full suite still fails until Task 3 creates all guide files.

- [ ] **Step 5: Commit**

~~~bash
git add index.html robots.txt sitemap.xml assets/irrigant-social-card.png tests/seo-static.test.mjs
git commit -m "feat: add Irrigant SEO foundation"
~~~

## Task 3: Publish the four grower guides

**Files:**
- Create: guides.css
- Create: irrigation-scheduling.html
- Create: soil-moisture-forecasting.html
- Create: how-helios-helps.html
- Create: potato-irrigation-guide.html

**Interfaces:**
- Consumes: the social image, the homepage waitlist URL /#waitlist, and shared guides.css.
- Produces: four indexable pages with unique metadata, cross-links, and conversion CTA.

- [ ] **Step 1: Write the failing guide test**

Append:

~~~js
test("gives each guide a unique grower question and a waitlist path", () => {
  const expectations = {
    "irrigation-scheduling.html": "How to make a better irrigation call",
    "soil-moisture-forecasting.html": "See the trend before it becomes stress",
    "how-helios-helps.html": "Better information before the water decision",
    "potato-irrigation-guide.html": "Potato irrigation starts with the field, not the calendar",
  };

  for (const [page, heading] of Object.entries(expectations)) {
    const html = read(page);
    expect(html).toContain("<h1>" + heading + "</h1>");
    expect(html).toContain('href="/#waitlist"');
    expect(html).toContain('href="guides.css"');
  }
});
~~~

- [ ] **Step 2: Run the test to verify it fails**

Run: bun test tests/seo-static.test.mjs

Expected: FAIL because the guide files and shared stylesheet do not exist.

- [ ] **Step 3: Implement the shared shell and guide content**

Create guides.css using the existing Geist type scale, off-white canvas, 1100 px content width, black pill CTA, blue text accent, and one-column mobile layout below 720 px. Do not copy the homepage product-widget styles.

Every guide must include:
- Page-specific title, description, canonical, Open Graph, and Twitter metadata, using the Task 2 social image.
- Wordmark, a Back to Irrigant link, four-guide navigation, a Join the waitlist CTA to /#waitlist, and the privacy footer.
- These titles, H1s, and section headings:

| File | Title | H1 | Sections |
| --- | --- | --- | --- |
| irrigation-scheduling.html | Irrigation Scheduling for Growers | Irrigant | How to make a better irrigation call | The decision is bigger than a date on the calendar; What to watch before you irrigate; Make the next call with more confidence |
| soil-moisture-forecasting.html | Soil-Moisture Forecasting for Growers | Irrigant | See the trend before it becomes stress | Soil moisture tells you what the field needs now; A forecast gives you time to act; Know more before the next pass |
| how-helios-helps.html | How Helios Helps Growers Make Irrigation Decisions | Irrigant | Better information before the water decision | Start with the conditions in front of you; See what is driving the call; Your judgment stays in control |
| potato-irrigation-guide.html | Potato Irrigation Guide for Growers | Irrigant | Potato irrigation starts with the field, not the calendar | Potatoes leave little room for guesswork; Read the signals before you apply water; Build the next decision from the field |

Guide body copy may describe only soil moisture, crop context, near-term weather, and crop-water demand. The potato guide must not give numerical application rates, crop-stage prescriptions, disease claims, or yield claims.

- [ ] **Step 4: Run the test to verify it passes**

Run: bun test tests/seo-static.test.mjs

Expected: PASS with four tests.

- [ ] **Step 5: Commit**

~~~bash
git add guides.css irrigation-scheduling.html soil-moisture-forecasting.html how-helios-helps.html potato-irrigation-guide.html tests/seo-static.test.mjs
git commit -m "feat: add Irrigant grower guides"
~~~

## Task 4: Verify and deploy the static site

**Files:**
- Verify: all files from Tasks 1-3.

**Interfaces:**
- Consumes: passing static tests and the connected Vercel project.
- Produces: live crawl assets and public guide URLs.

- [ ] **Step 1: Run the full static test suite**

Run: bun test tests/seo-static.test.mjs

Expected: PASS with no failures.

- [ ] **Step 2: Inspect locally**

Run: python3 -m http.server 4173.

Use a browser to check the homepage and all four guide URLs at desktop and 375 px. Confirm navigation, waitlist CTA, canonical tags, social-image URL, and no console errors.

- [ ] **Step 3: Commit browser fixes only if needed**

~~~bash
git add index.html guides.css irrigation-scheduling.html soil-moisture-forecasting.html how-helios-helps.html potato-irrigation-guide.html robots.txt sitemap.xml tests/seo-static.test.mjs
git commit -m "fix: polish Irrigant SEO pages"
~~~

- [ ] **Step 4: Deploy through Vercel**

Use the existing connected Vercel project. Request confirmation immediately before the final production deployment action because it changes the public website.

- [ ] **Step 5: Verify live**

Check the homepage, robots.txt, sitemap.xml, and all four guide URLs for a 200 response. Re-run metadata, canonical, JSON-LD, mobile, console-error, and network checks.

## Task 5: Connect search and conversion measurement

**Files:**
- Modify: index.html, only if the selected GA4 setup requires a site tag.

**Interfaces:**
- Consumes: the deployed sitemap plus owner access to Google Search Console and Analytics.
- Produces: sitemap submission, indexing requests, and a verified waitlist conversion.

- [ ] **Step 1: Inspect existing Google setup without changing state**

Open Search Console and Analytics in separate browser tabs. Check for an existing irrigant.xyz Domain property and a GA4 web data stream.

- [ ] **Step 2: Hand off only when necessary**

Ask for login, MFA, CAPTCHA completion, or DNS verification only when the screen requires it. Never enter credentials, solve a CAPTCHA, or change DNS records.

- [ ] **Step 3: Create or connect the Google properties**

At the final creation step, request confirmation before creating a Search Console Domain property or GA4 property. Prefer Domain-property DNS verification.

- [ ] **Step 4: Add GA4 measurement only after a Measurement ID exists**

Write a failing static test requiring the approved Measurement ID and this event inside the existing successful Formspree-submit branch:

~~~js
window.gtag?.("event", "generate_lead", {
  form_id: "waitlist",
  method: "formspree",
});
~~~

Add the approved gtag snippet and event, run the test to green, and deploy through Task 4.

- [ ] **Step 5: Submit and verify**

With immediate confirmation before each external change, submit the sitemap, request indexing for all five indexable pages, and verify the generate_lead event in GA4 Realtime or DebugView. Use a safe test identity only if the user supplies one.

- [ ] **Step 6: Start the monthly review loop**

Record indexed-page count, search queries, impressions, clicks, waitlist conversions, and the next content decision. Do not expand content without demand evidence or source-backed material.
