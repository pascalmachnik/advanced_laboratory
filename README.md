# Advanced Laboratory

This repository contains LaTeX templates and source files for reports written
for the Master's laboratory course in physics at TU Dortmund University
(Masterpraktikum / Advanced Laboratory). Reports can be prepared in English or
German.

## Repository structure

```text
advanced_laboratory/
├── english/
│   ├── Makefile
│   ├── header.tex
│   ├── header-matplotlib.tex
│   ├── lit.bib
│   ├── programme.bib
│   └── vXXX/
│       ├── Makefile
│       ├── vXXX.tex
│       ├── plot.py
│       └── content/
└── german/
    └── ... (same layout as english/)
```

- `english/` and `german/` hold language-specific report templates and shared
  formatting settings.
- Each `vXXX/` directory is a report template. Its main LaTeX document includes
  the individual section files from `content/`.
- `header.tex` defines the report-wide LaTeX setup and loads the bibliography
  files. `header-matplotlib.tex` and `matplotlibrc` provide matching settings
  for plots made with Matplotlib.
- `lit.bib` and `programme.bib` are shared bibliography databases for reports
  in that language directory.
- Each report's `Makefile` builds its PDF into `build/`. The `build/` directory
  contains generated files and can be recreated by building the report.

## Producing a report

1. Choose `english/` or `german/` and copy its `vXXX/` directory to a
   directory named for the experiment, for example `v123/`.
2. Rename `vXXX.tex` to match the new directory (for example `v123.tex`) and
   update the corresponding names and output paths in that report's
   `Makefile`. The template currently uses `vXXX` in its build targets, so
   these names need to be updated together.
3. Edit the main `.tex` file: set the experiment number, title, and dates, and
   include the report sections you need. Replace the template text in the
   files under `content/`. The English template has objective, theory,
   procedure, analysis, and discussion sections; the German template has the
   equivalent sections and an optional `messunsicherheiten.tex` section.
4. If needed, update the author details in the language directory's
   `header.tex`. Add references to `lit.bib` or `programme.bib`, or add another
   bibliography resource there if the report needs its own database.
5. Put analysis or plotting code in the report directory. The supplied
   `plot.py` is a starting point for Matplotlib plots; add generated plots and
   other inputs as dependencies in the report's `Makefile` where indicated.
6. Build and preview the PDF:

   ```sh
   make -C english/v123
   make -C english/v123 preview
   ```

   Substitute `german` or your report directory as appropriate. The PDF is
   written to `build/v123.pdf`. `preview` watches for changes and rebuilds;
   stop it with Ctrl+C. To remove generated build files, run
   `make -C english/v123 clean`.

The report Makefile uses `latexmk` with LuaLaTeX and expects a LaTeX
installation with the packages used by the headers, including Biber for
bibliographies. The supplied Makefile also runs `plot.py` as part of the
build, so Python with Matplotlib and NumPy is required unless you remove or
adjust the plotting target for a report that does not use it.

The language-level Makefiles currently target the `vXXX` template directory.
If you want `make -C english` or `make -C german` to build a newly named report,
update that language-level Makefile to point to the new directory.
