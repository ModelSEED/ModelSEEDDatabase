# Applied manifests

Apply_Manifest.py archives every successfully-applied manifest here
with a `YYYY-MM-DD-` prefix on the original filename. The archive is
the audit trail for curated changes — `git log` of this directory
gives a date-ordered narrative of every non-trivial database edit.

Don't edit files in this directory. To "redo" an applied change,
copy the manifest back up to `Biochemistry/Updates/<filename>.yaml`,
edit, and re-apply.
