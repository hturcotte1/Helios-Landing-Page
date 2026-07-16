# Privacy Policy Design

## Goal

Add a minimal, accurate U.S. privacy policy for Irrigant's Helios waitlist landing page. The policy must describe the information the page actually collects and give visitors a clear way to make privacy requests.

## Scope

- Add a standalone `privacy.html` page styled with the existing Geist-based light design system.
- Add a `Privacy` link in the landing-page footer.
- Add one short notice beneath the waitlist submit button that links to the policy.
- Do not add a cookie banner, consent manager, analytics, tracking, or new data collection.

## Policy Content

The policy will identify Irrigant as the operator, use an effective date of June 21, 2026, and state that the site collects first name, last name, email address, farm-size range, and location when a visitor joins the waitlist.

It will explain that Irrigant uses this information to manage the waitlist and communicate about Helios; Formspree processes form submissions; the site loads Google Fonts; and Irrigant does not sell personal information or share it for cross-context behavioral advertising. It will not claim CCPA applicability. Instead, it will provide a single privacy-request email route for applicable access, correction, deletion, and opt-out requests.

The page will include concise sections covering data use, sharing, retention, security, children, policy changes, and contact details. It will not make guarantees that cannot be supported by the current site implementation.

## Architecture and Data Flow

`index.html` remains the landing page. Its footer and waitlist note link to the new sibling page, `privacy.html`. The policy page is informational only and introduces no JavaScript or new network requests beyond the fonts already used by the landing page.

## Error Handling and Accessibility

Both links use ordinary anchors so the policy remains available without JavaScript. The policy page uses semantic headings, a main landmark, readable line length, and a visible route back to the landing page.

## Verification

1. Automated content assertions confirm both `index.html` links target `privacy.html` and that the policy contains the required disclosures, rights contact email, and effective date.
2. Browser checks confirm the policy opens from the footer and the waitlist notice, contains no console errors, and remains readable at mobile and desktop widths.

## Assumptions

- "At the bottom" means the primary permanent link belongs in the landing-page footer.
- The one-line notice below the waitlist button remains necessary because the page collects personal information at that point.
- The existing contact email, `marcotrotta909@gmail.com`, handles privacy requests until Irrigant provides a dedicated privacy address.
