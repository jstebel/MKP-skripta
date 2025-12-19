# Repository Guidelines

## Development Focus (marimo)

- Primary development happens in `marimo/`. Treat `prednaska/`, `cviceni/`, and `web/` as reference/context unless a change is explicitly requested.
- Keep changes scoped: prefer updating notebooks and environment files over touching LaTeX/PDF/web assets.

## Project Structure

- `prednaska/`: LaTeX sources for lecture notes (main: `prednaska/mkp_prednaska.tex`) plus images in `prednaska/img/`.
- `cviceni/`: exercise materials (LaTeX: `cviceni/mkp_cvika.ltx`, MATLAB scripts in `cviceni/matlab/` and `cviceni/matlab-Pavel/`, assignments/tests in `cviceni/testy/`).
- `marimo/`: Python “notebook-as-code” demos and an optional Conda environment (`marimo/conda-requirements.yml`).
- `web/`: simple static site (`web/index.html`, `web/style.css`) and published PDFs in `web/pdf/`.

## Build, Test, and Development Commands

- Build PDFs and publish to `web/pdf/`:
  - `./web_update.sh` (runs `pdflatex` twice for both `cviceni` and `prednaska`, then copies PDFs).
- Note: generated PDFs in `prednaska/` and `cviceni/` are gitignored; the published copies live in `web/pdf/`.
- Build individually (useful for quicker iteration):
  - `cd prednaska && pdflatex mkp_prednaska.tex && pdflatex mkp_prednaska.tex`
  - `cd cviceni && pdflatex mkp_cvika.ltx && pdflatex mkp_cvika.ltx`
- (Optional) Work with marimo notebooks:
  - `bash marimo/mamba_env.sh edit marimo/demo.ntb.py`
  - `bash marimo/mamba_env.sh run marimo/demo.ntb.py`
  - Key notebook (recreate env + run): `bash marimo/mamba_env.sh --force run marimo/amr_fenicsx_adaptivity.py`
  - Headless (don’t auto-open a browser): `bash marimo/mamba_env.sh --force run --headless marimo/amr_fenicsx_adaptivity.py`

## Coding Style & Naming Conventions

- marimo: keep notebooks as `.ntb.py` and preserve the header metadata block (dependencies/versions) when present.
- marimo: functions must not use local variables; return global variables with a single `return` statement on the last line.
- LaTeX in marimo/Markdown: display math blocks `\\[ ... \\]` must be preceded by a blank line, otherwise they may not render.
- Python: use 4-space indentation, explicit imports, and descriptive names; keep cells small and focused.
- Other folders: avoid drive-by reformatting; if changes are required, keep diffs minimal.

## Testing Guidelines

- No automated test suite is configured. Use “build as test”:
  - Ensure `pdflatex` completes without errors and produces expected PDFs.
  - If you change LaTeX sources, regenerate and verify the published artifacts in `web/pdf/`.
  - For `marimo/`, run the notebook (`marimo run …`) and sanity-check outputs.

## Commit & Pull Request Guidelines

- Commits in history are short, descriptive, and often scoped to content (e.g., “Update…”, “Add…”, “Fix…”). Follow the same pattern.
- Suggested convention: prefix with area when helpful (`prednaska: …`, `cviceni: …`, `web: …`, `marimo: …`).
- PRs should include: a brief summary, what was rebuilt (e.g., “ran `./web_update.sh`”), and any render-impacting changes (screenshots or notes for major formatting changes).

## Security & Configuration Notes

- `marimo/mamba_env.sh` may install Miniconda and uses `sudo` + network downloads. Review before running, and avoid invoking it in restricted environments/CI without explicit approval.
