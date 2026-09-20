# Project working agreements

- Commit each significant change after completing the relevant verification. Use a clear commit message and include only files belonging to that change.
- Push once per working day, at the end of that day's work, rather than after every commit. Check the conversation and available Git evidence to avoid duplicate daily pushes. Explicit user requests to push take precedence.
- This push cadence is a working convention, not an automatic schedule. Do not claim that unattended daily pushes are configured.
- Keep Python dependencies isolated in this project's `.venv`; use `.venv/bin/python` and `.venv/bin/python -m pip`. Do not modify the user's conda environment.
- Keep virtual environments, caches, credentials, and the local course reference PDF out of Git, following `.gitignore`.
