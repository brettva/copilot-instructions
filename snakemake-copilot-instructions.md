# 📘 Snakemake Workflow Copilot Instructions

> **Copilot MUST read and adhere to this file at the beginning of every prompt.**  
> This file defines required behavior for Copilot's operation across all sessions.  
> If these instructions conflict with model defaults, these take precedence.

---

## 1. 📚 Use the Most Current and Authoritative Resources

- For any question involving a software or tool **library**, **tool**, **API**, or **framework**:
  1. Always **Execute `resolve-library-id`** to identify the specific entity being referenced .
  2. Always **Execute `get-library-docs`** to retrieve up-to-date documentation from **Context7**.
  3. If Context7 documentation is **unavailable or incomplete**, fall back to internal model knowledge or other trusted sources **only if necessary**, and clearly state that this fallback occurred.

- When encountering ambiguity or inconsistency, **prioritize official sources** and documentation.

---

## 2. 🛠️ Write Clean, Modular Snakemake Rules

- Snakemake rules must:
  - Follow the **single responsibility principle** — each rule should do **one thing well**.
  - Avoid embedding multiple steps, conditional logic, or long shell pipelines in a single rule.
  - Be split into multiple rules if complexity increases, improving clarity and maintainability.

- Rules should be readable and follow best practices, including:
  - Clear `input`, `output`, and `params` sections.
  - Descriptive naming.
  - Avoiding unnecessary shell logic inline.
  - Maintain linear rule order when possible: when the output of one rule is used as the input for another, place those rules next to each other in the Snakefile in the same order the data are processed (producer immediately before consumer). This sequencing improves organization and readability, makes data flow easier to follow, and reduces maintenance overhead. If a rule has multiple downstream consumers, place its primary or most closely related consumers immediately after it and use clear section comments to separate stages.

---

## 3. 📝 Snakemake Script Interface Syntax (R & Python)

When generating Snakemake rules that use the `script:` directive to call an external R or Python script, always use the correct Snakemake-provided variable interface inside the script:

- **For R scripts:**
  - Use `snakemake@input`, `snakemake@output`, `snakemake@params`, etc.
  - Example: `snakemake@input[[1]]`, `snakemake@output[[1]]`

- **For Python scripts:**
  - Use `snakemake.input`, `snakemake.output`, `snakemake.params`, etc.
  - Example: `snakemake.input[0]`, `snakemake.output[0]`

Never use hardcoded filenames or incorrect variable names for Snakemake script interfaces. Always follow the above conventions for maximum compatibility and reproducibility.

---

## 4. 🧬 Snakemake R and Python Script Placement

- **R code in Snakemake rules:**
  - Always place R code in a separate script file.
  - Use the `script:` directive in the rule to call the R script.
  - The R script must use the correct Snakemake variable interface (e.g., `snakemake@input`, `snakemake@output`).

- **Python code in Snakemake rules:**
  - If the logic is short/simple, it may be placed directly in the rule using the `run:` directive.
  - If the logic is lengthy or complex, place it in a dedicated Python script and call it with `script:`.
  - The Python script must use the correct Snakemake variable interface (e.g., `snakemake.input`, `snakemake.output`).

This ensures modular, maintainable, and reproducible workflows. Follow this approach for all Snakemake-related code unless otherwise specified.

---

## 5. 🗂️ Project-based Output Organization and Naming

- **Always use a `project` key in the Snakemake config file.**
- All output files and directories must be organized within a folder called results under a project-specific subfolder, and filenames should include the project name at the very beginning of the name for clarity and reproducibility. 
- **When constructing output paths in Python f-strings:**
  - The `project` key should be referenced directly as `config['project']` (not as a wildcard).
  - Any Snakemake wildcards should be escaped with curly braces `{}` inside the f-string, e.g., `"{wildcard}"`.
- Example pattern:
  - `f"results/{config['project']}/plots/{config['project']}.{{wildcard}}.my_plot.png"`
  - `f"results/{config['project']}/data/{config['project']}.{{wildcard}}.my_data.txt"`

- Example rule:

rule example_rule:
    input:
        "data/{sample}.input.txt"
    output:
        f"results/{config['project']}/data/{config['project']}.{{sample}}.output.txt"
    shell:
        "process_data {input} {output}"

- This ensures outputs from different runs/configurations are kept separate and are easily identifiable, even if files are moved outside their original folders.

- The exception to this rule is that files in the resources dir should not include the project name in their filenames, as they are shared resources across multiple projects. For example, `resources/bravo.{chrom}.norm.bcf`.

---

### 5.1 🗃️ Output Directory Structure: Data, Results, and Plots

- **Project outputs must be organized by type within the project-specific results subfolder:**
  - **Intermediate data files** (e.g., temporary or working files used as inputs for downstream steps) must be placed in a `data` subdirectory:  
    `results/{project}/data/`
  - **Final results** (e.g., summary tables, files intended for review or reporting) must be placed in a `results` subdirectory:  
    `results/{project}/results/`
  - **Plots and figures** must be placed in a `plots` subdirectory:  
    `results/{project}/plots/`
- **Filenames must still begin with the project name** as described above.
- **Examples:**
  - Intermediate data: `results/{config['project']}/data/{config['project']}.{wildcard}.intermediate.txt`
  - Final result: `results/{config['project']}/results/{config['project']}.{wildcard}.summary.tsv`
  - Plot: `results/{config['project']}/plots/{config['project']}.{wildcard}.qc.png`
- This ensures clear separation of intermediate files, final results, and visualizations, improving workflow clarity and reproducibility.

---

## 6. 🧪 Always Keep `rule all` Up-to-Date

- **When creating or editing Snakemake rules:**
  - Always update the `rule all` outputs to include the most downstream output(s) of any new or modified rules, so that the workflow is immediately testable after changes.
  - If new rules make some previous `rule all` outputs obsolete (i.e., they are no longer produced by any rule), remove those outputs from `rule all` to prevent bloat and confusion.
  - This ensures that `rule all` always reflects the current expected outputs of the workflow and supports efficient, up-to-date testing.

---

## 🧩 Notes

- Additional domain-specific or tool-specific rules may be added in future revisions of this document.
- These guidelines are designed to maximize maintainability, reproducibility, and compatibility with evolving infrastructure.

