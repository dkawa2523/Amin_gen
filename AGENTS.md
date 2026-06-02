# AGENTS.md

## Project goal

This repository implements a Semiconductor Gas Candidate Screening MVP.
The system generates and screens candidate molecules for semiconductor process gas and precursor selection.

## Core principles

1. Local-first
   Always use local libraries and local CSV tables before remote APIs.

2. API conservation
   Do not call external APIs unless:
   - the candidate passed prefilter,
   - the candidate has sufficient API priority score,
   - the value is missing locally,
   - cache and negative cache do not already answer the query.

3. Safe defaults
   Remote APIs must be disabled by default.
   Exploration mode must have candidate and API limits.

4. Evidence separation
   Store every raw candidate value in Evidence.
   Show only selected values and engineering decisions in Summary.

5. Semiconductor engineering focus
   Optimize for:
   - phase at 25C / 1 atm,
   - vapor pressure at 25/40/60C,
   - supply class,
   - PFAS screening,
   - persistence screening,
   - reactivity flags,
   - kinetics coverage.
   Do not expand scope to full EHS approval, SDS parsing, or process emission modeling.

6. Chemistry caution
   CAS is not the primary key.
   Use InChIKey when available.
   Treat salts, mixtures, radicals, formula-only molecules, and unsupported inorganic species as review-required unless resolved.

7. Testing
   Unit tests must never require live external API calls.
   Mock remote APIs or use cache fixtures.
   Always run pytest before final response.

8. Remote providers
   PubChem PUG-REST is for identity and structure.
   PubChem PUG-View is for shortlist GHS/annotation enrichment only.
   NIST, Chemeo, CompTox, CAS, LXCat, and QDB integrations must be config-gated.

9. Output
   Keep Summary, Evidence, Coverage, Review Required, Rejected, and Run Stats sheets stable.
