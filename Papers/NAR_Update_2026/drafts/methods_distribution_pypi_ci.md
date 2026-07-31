# Methods — Distribution + Continuous Integration (draft)

**Target sections:** Materials and Methods, "Distribution: PyPi + API" and the community-contributions / CI paragraph (see `PAPER_2026_SKELETON.md`).
**Guide references:** `PAPER_2026_GUIDE.md` §12, §13.
**Status:** first pass. Both are short infrastructure paragraphs, not novel contributions.

---

## Distribution: PyPi + API

The 2020 release published the database as a Git repository and a Solr REST endpoint; both remain available. The 2026 release adds an installable PyPi package (`modelseed-biochem`, [TBD confirm exact package name]) that ships the database contents plus a documented Python API for programmatic access. Typical use cases include fetching a compound record by ID, walking a reaction's reagents, querying the thermodynamics ledger, and retrieving atom mappings. Installation and usage examples are documented in the repository README and mirror the shape of the existing Solr access patterns, so existing users can migrate incrementally.

Version pinning: the PyPi package version tracks the ModelSEED release tag (i.e. installing `modelseed-biochem==2.0.0` pulls the database contents as of the `v2.0.0` release tag). This lets downstream tools depend on a specific database state rather than the moving `master` head.

## Continuous integration: GitHub Actions

The 2020 release described a Travis CI setup that validated pull requests to `dev` at submission time. Travis is no longer used; the same intent — validating any newly curated entries at PR time — is now implemented via a GitHub Actions workflow. The workflow runs on every PR to `dev` and checks: (i) that the compound and reaction JSON shards parse and match their expected schema; (ii) that newly added or modified reactions remain mass- and charge-balanced against the (possibly modified) compound set; (iii) that newly added override or exclusion file rows conform to the appropriate schema; (iv) that every curated pick carries the required provenance fields. Failures block the PR from being merged until the contributor addresses them.

The workflow definition is under `.github/workflows/` in the repository; contributors can run the same checks locally before opening a PR by invoking the validation script directly (`Scripts/Validation/`, [TBD confirm script path]).

---

## Open loose ends flagged during drafting

- The PyPi package doesn't yet exist. Guide §12 tracks the packaging + publication work; Methods paragraph needs the confirmed package name once published.
- The GitHub Actions workflow file(s) may not yet exist under `.github/workflows/`. Guide §13 tracks that work; Methods paragraph will link to the workflow file once landed.
- Validation-script path (`Scripts/Validation/`?) needs confirming — the 2020-era Travis validation logic needs finding, refactoring where necessary, and pointing at explicitly.
