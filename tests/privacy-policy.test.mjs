import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

const readSiteFile = (name) => readFileSync(fileURLToPath(new URL(`../${name}`, import.meta.url)), 'utf8');

test('links the privacy policy from the collection notice and footer', () => {
  const landing = readSiteFile('index.html');
  const privacyLinks = landing.match(/href="privacy\.html"/g) ?? [];

  assert.equal(privacyLinks.length, 2);
  assert.match(landing, /By joining the waitlist, you agree to our <a href="privacy\.html">Privacy Policy<\/a>\./);
  assert.match(landing, /<footer>[\s\S]*?<a\b[^>]*href="privacy\.html"[^>]*>Privacy<\/a>/);
});

test('states the policy scope and privacy request route', () => {
  const policy = readSiteFile('privacy.html');

  assert.match(policy, /<title>Privacy Policy \| Irrigant<\/title>/);
  assert.match(policy, /Effective date: June 21, 2026/);
  assert.match(policy, /first name, last name, email address, farm-size range, and location/i);
  assert.match(policy, /Formspree/i);
  assert.match(policy, /Google Fonts/i);
  assert.match(policy, /do not sell your personal information or share it for cross-context behavioral advertising/i);
  assert.match(policy, /marcotrotta909@gmail\.com/);
  assert.doesNotMatch(policy, /<script\b/i);
});
