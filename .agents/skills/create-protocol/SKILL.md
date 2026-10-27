---
name: create-protocol
description: Create a new advanced laboratory protocol from the vXXX template. Use when the user asks to create a protocol and provides or wants to choose an experiment number and language.
---

# Create a new protocol

Create a protocol by copying the language-specific `vXXX` template and replacing every `vXXX` occurrence with the requested experiment identifier (for example, `v61`).

## Required details

A new protocol requires an experiment number and a language (`english` or `german`). If either is missing, ask the user before creating files. If both are provided, proceed. Add the `v` prefix to the number: experiment `61` becomes `v61`.

## Select the template

- English: copy `english/vXXX/` to `english/v<NUMBER>/`.
- German: copy `german/vXXX/` to `german/v<NUMBER>/`.

Create the protocol in only the requested language unless the user explicitly asks for both. Preserve the complete template directory, including its subdirectories and files.

## Create the protocol

1. Copy the selected `vXXX` directory to `v<NUMBER>` in the same language folder.
2. Rename `vXXX.tex` in the new directory to `v<NUMBER>.tex`.
3. Replace every literal `vXXX` occurrence in every file in the new directory with `v<NUMBER>`. This includes the main TeX document and the protocol `Makefile`. Search the entire copied directory to catch all occurrences.
4. Preserve other template content. Do not invent experiment-specific titles, authors, dates, or scientific content unless supplied by the user or requested.

## Verify and report

Before reporting completion, confirm that:

- The new directory is in the requested language folder and named `v<NUMBER>`.
- The main TeX file is named `v<NUMBER>.tex`.
- No `vXXX` occurrences remain anywhere in the new directory.
- The Makefile uses the new identifier wherever the template used `vXXX`.

Report the created path and language, and mention any files that could not be copied or updated.
