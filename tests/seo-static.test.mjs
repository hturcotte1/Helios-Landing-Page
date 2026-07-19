import { describe, expect, test } from "bun:test";
import { existsSync, readFileSync } from "node:fs";

const pages = [
  "index.html",
  "irrigation-scheduling.html",
  "soil-moisture-forecasting.html",
  "how-helios-helps.html",
  "machine-learning.html",
];

const read = (path) => existsSync(path) ? readFileSync(path, "utf8") : "";

describe("public SEO assets", () => {
  test("positions the homepage around irrigation decisions", () => {
    const html = read("index.html");

    expect(html).toContain("Irrigation Scheduling &amp; Soil-Moisture Forecasting | Irrigant");
    expect(html).toContain("Know when to water. And how much.");
    expect(html).toContain('"@type":"Organization"');
    expect(html).toContain('"@type":"WebSite"');
    expect(html).toContain('href="irrigation-scheduling.html"');
  });

  test("publishes crawler discovery for every public page", () => {
    expect(read("robots.txt")).toContain("Sitemap: https://irrigant.xyz/sitemap.xml");
    const sitemap = read("sitemap.xml");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/irrigation-scheduling.html</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/soil-moisture-forecasting.html</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/how-helios-helps.html</loc>");
    expect(sitemap).toContain("<loc>https://irrigant.xyz/machine-learning.html</loc>");
  });

  test("gives every public page canonical and social metadata", () => {
    expect(existsSync("assets/irrigant-social-card.png")).toBe(true);

    for (const page of pages) {
      const html = read(page);
      expect(html).toMatch(/<link rel="canonical" href="https:\/\/irrigant\.xyz\/.*">/);
      expect(html).toMatch(/<meta name="description" content="[^"]+">/);
      expect(html).toContain("https://irrigant.xyz/assets/irrigant-social-card.png");
      expect(html).toContain('<meta name="twitter:card" content="summary_large_image">');
    }
  });

  test("keeps unsupported technical claims out of public copy", () => {
    const copy = pages.map(read).join("\n");
    expect(copy).toContain("Join the waitlist");
    expect(copy).not.toMatch(/LightGBM|calibrated confidence|guaranteed water savings/i);
  });

  test("gives each guide a unique grower question and a waitlist path", () => {
    const expectations = {
      "irrigation-scheduling.html": "How to make a better irrigation call",
      "soil-moisture-forecasting.html": "See the trend before it becomes stress",
      "how-helios-helps.html": "Better information before the water decision",
      "machine-learning.html": "Machine learning for better irrigation decisions",
    };

    for (const [page, heading] of Object.entries(expectations)) {
      const html = read(page);
      expect(html).toContain("<h1>" + heading + "</h1>");
      expect(html).toContain('href="/#waitlist"');
      expect(html).toContain('href="guides.css"');
    }
  });

  test("indexes the Machine Learning page and no longer publishes the potato guide", () => {
    const sitemap = read("sitemap.xml");
    const machineLearning = read("machine-learning.html");
    const guideStyles = read("guides.css");

    expect(sitemap).toContain("<loc>https://irrigant.xyz/machine-learning.html</loc>");
    expect(sitemap).not.toContain("potato-irrigation-guide.html");
    expect(machineLearning).toContain("<h1>Machine learning for better irrigation decisions</h1>");
    expect(machineLearning).toMatch(/XGBoost multi-output regression/);
    expect(machineLearning).toMatch(/OpenET monthly ET enrichment/);
    expect(guideStyles).not.toContain("border-radius:16px");
  });
});
