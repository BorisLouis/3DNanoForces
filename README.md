# Repository release templates

These templates are for software repositories. Copy `LICENSE.template` to
`LICENSE` and `CITATION.cff.template` to `CITATION.cff` in each repository's root,
and replace the bracketed placeholders. The `.template` suffix prevents these
drafts from being mistaken for the current workspace's license or citation.

## License choices

The license is a custom draft, not a standard approved license or legal advice.
Have your institution's legal or technology-transfer team review it before use,
including who owns the copyright and can authorize commercial licensing.

It allows noncommercial use, modification, and redistribution, requires credit
and software citation, and reserves commercial use for separate written
permission. Its definition treats paid consulting, business operations, and
commercial product development as commercial use, even when undertaken at a
university. Check that this matches your intended boundary. Citation of an
associated paper is optional unless you explicitly add that requirement.

This is source-available software with a noncommercial restriction, rather than
open-source software under the [Open Source Definition](https://opensource.org/osd).
Creative Commons [recommends against CC licenses for software](https://creativecommons.org/faq/#can-i-apply-a-creative-commons-license-to-software),
which is why this draft is tailored to software. Separately identify any datasets,
figures, or other materials requiring different terms.

Check existing licenses and contributor rights first. A new license cannot
withdraw permissions already granted for earlier copies, or override third-party
licenses. Do not replace an existing license unless you have the necessary rights.

## Release sequence

1. Fill in the copyright holder(s), years, project, repository, and licensing
   contact in `LICENSE`. Commit it before creating the release tag.
2. Preferably fill in and commit `CITATION.cff` at the same time, leaving the DOI
   fields commented out. Set the version, actual release date, repository URL,
   and license URL. Optional fields must contain real values or be removed.
3. Enable the repository's Zenodo integration, then create the GitHub release
   from the commit containing the intended license and citation metadata.
4. After archiving succeeds, copy that release's version DOI into the top-level
   `doi` field and commit the updated `CITATION.cff`. Optionally add the concept
   DOI under `identifiers` for all versions. You may instead first add the
   citation file at this stage, as originally planned.
5. An edit after release does not change the already archived snapshot. Keep
   its version/date/DOI together; for a later release, update the version/date
   and remove the old version DOI before archiving, then add the new DOI.
6. Check the published Zenodo record's authors, citation, and license. Do not
   select MIT, Apache, or a CC license as a substitute for this custom license;
   use Zenodo's custom-license mechanism and provide the actual terms.

Zenodo parses supported CFF metadata on a best-effort basis. If `.zenodo.json`
also exists, its metadata takes precedence. See [Zenodo's CFF guidance](https://help.zenodo.org/docs/github/describe-software/citation-file/)
and [metadata precedence](https://support.zenodo.org/help/en-gb/24-github-integration/96-how-does-a-citation-cff-file-affect-metadata-of-my-github-release).

Validate the completed file against [CFF 1.2.0](https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md),
for example with `cffconvert --validate` from the repository root after installing
`cffconvert`. This template intentionally contains unfinished metadata and is
not a publication-ready citation until filled in.
