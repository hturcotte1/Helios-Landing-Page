# Editorial Guides and Machine Learning Page Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Turn the guides into a continuous Irrigant editorial experience and replace Potato Irrigation with a technical Machine Learning page.

**Architecture:** Keep the shared static guide shell. Replace card rules with a vertical article layout and thin dividers. Swap the potato file, URL, navigation, and sitemap entry for machine-learning.html.

**Tech Stack:** Static HTML, CSS, Bun tests.

## Global Constraints

- Grower guides remain plain-language.
- The Machine Learning page describes the documented Helios stack accurately without unsupported metrics or outcomes.
- Remove section-card styling. Preserve Geist, off-white canvas, hairline rules, black CTA, and responsive navigation.

---

### Task 1: Define the replacement and layout contract

**Files:**
- Modify: tests/seo-static.test.mjs

- [ ] **Step 1: Write the failing test**

~~~js
test("indexes the Machine Learning page and no longer publishes the potato guide", () => {
  const sitemap = read("sitemap.xml");
  const machineLearning = read("machine-learning.html");
  const guideStyles = read("guides.css");

  expect(sitemap).toContain("<loc>https://irrigant.xyz/machine-learning.html</loc>");
  expect(sitemap).not.toContain("potato-irrigation-guide.html");
  expect(machineLearning).toContain("<h1>Machine learning for better irrigation decisions</h1>");
  expect(machineLearning).toMatch(/XGBoost multi-output regression/);
  expect(machineLearning).toMatch(/OpenET monthly ET enrichment/);
  expect(guideStyles).not.toContain(".section{");
});
~~~

- [ ] **Step 2: Run the test to verify it fails**

Run: bun test tests/seo-static.test.mjs --test-name-pattern 'Machine Learning page'

Expected: FAIL because the potato URL remains and machine-learning.html does not exist.

- [ ] **Step 3: Commit the failing test**

~~~bash
git add tests/seo-static.test.mjs
git commit -m "test: define editorial guide contract"
~~~

### Task 2: Replace cards with editorial flow

**Files:**
- Modify: guides.css
- Modify: irrigation-scheduling.html
- Modify: soil-moisture-forecasting.html
- Modify: how-helios-helps.html

- [ ] **Step 1: Replace the section CSS**

~~~css
.sections{margin-top:80px}.section{padding:42px 0;border-top:1px solid var(--line)}.section:last-child{border-bottom:1px solid var(--line)}.section p{max-width:720px;margin:0;color:var(--body);font-size:18px}.section h2{max-width:760px;margin:0 0 14px}
~~~

- [ ] **Step 2: Keep existing guide content semantic**

Retain each existing H1, three H2 sections, waitlist CTA, metadata, guide navigation, and privacy footer. Do not add background fills, rounded corners, or borders around individual sections.

- [ ] **Step 3: Run the static suite**

Run: bun test tests/seo-static.test.mjs

Expected: existing metadata and waitlist tests remain green; the Machine Learning test remains red until Task 3.

### Task 3: Add the Machine Learning page

**Files:**
- Create: machine-learning.html
- Delete: potato-irrigation-guide.html
- Modify: index.html
- Modify: sitemap.xml
- Modify: tests/seo-static.test.mjs

- [ ] **Step 1: Create the page**

Use title Machine Learning for Irrigation Decisions | Irrigant, canonical https://irrigant.xyz/machine-learning.html, and the existing social image. Use H1 Machine learning for better irrigation decisions.

The sections are: System inputs, Forecasting layer, Runtime enrichment, Decision layer, and Product architecture.

Describe: XGBoost multi-output regression at 24, 48, and 72 hours; soil-moisture, crop, weather, ET, and operating-constraint inputs; OpenET monthly ET enrichment; NOAA weather backfill; cache and fallback behavior; pump, budget, infiltration, and water-window constraints; FastAPI, static frontend, API endpoints, and nearby feedback.

- [ ] **Step 2: Replace Potato Irrigation references**

Replace potato-irrigation-guide.html with machine-learning.html in index.html, all guide navigation blocks, and sitemap.xml. Delete the potato guide.

- [ ] **Step 3: Run tests to green**

Run: bun test tests/seo-static.test.mjs

Expected: PASS with the updated five-page indexable inventory.

- [ ] **Step 4: Commit**

~~~bash
git add guides.css index.html irrigation-scheduling.html soil-moisture-forecasting.html how-helios-helps.html machine-learning.html sitemap.xml tests/seo-static.test.mjs
git rm potato-irrigation-guide.html
git commit -m "feat: add Irrigant Machine Learning page"
~~~

### Task 4: Verify locally and publish

- [ ] **Step 1: Run all tests**

Run: bun test tests/seo-static.test.mjs && node --test tests/privacy-policy.test.mjs

Expected: both suites pass.

- [ ] **Step 2: Browser QA**

Serve the repository with python3 -m http.server 4173. Check every guide at desktop and 375 px. Verify readable article flow, no boxed sections, working nav/CTA, and no console errors.

- [ ] **Step 3: Deploy only with confirmation**

Merge and push the verified branch to main, then confirm the live Machine Learning page, sitemap, and homepage link return 200.

