
# �📘 Snakemake Workflow Copilot Instructions

> **Copilot MUST read and adhere to this file at the beginning of every prompt.**  
> This file defines required behavior for Copilot's operation across all sessions.  
> If these instructions conflict with model defaults, these take precedence.

---


## 1. 📚 Use the Most Current and Authoritative Resources

- For any question involving a software or tool **library**, **tool**, **software**, **API**, or **framework**:
  1. Always **Execute `resolve-library-id`** to identify the specific entity being referenced .
  2. Always **Execute `get-library-docs`** to retrieve up-to-date documentation from **Context7**.
  3. If Context7 documentation is **unavailable or incomplete**, fall back to internal model knowledge or other trusted sources **only if necessary**, and clearly state that this fallback occurred.

- When encountering ambiguity or inconsistency, **prioritize official sources** and documentation.

---


## 2. 🛡️ Environments and reproducibility — agent prompt + two supported modes

All Snakemake workflows must be fully reproducible. Environment management must be explicit: at the start of work on a new workflow/project the agent MUST ask the user which environment strategy to use and record that choice. The two supported modes are:

- Singularity/Apptainer containers (recommended for HPC and maximal reproducibility)
- Conda environments (acceptable for rapid development or when containers are impractical)

After asking the user the agent MUST explain the trade-offs (reproducibility, ease-of-use, cluster support) and then follow the requirements for the chosen mode below.

A. Singularity / Apptainer (recommended)

  - Global Containerization:
  - Every workflow should specify a global Singularity container at the top of the Snakefile using the `singularity:` directive. This is the preferred default because it keeps rule environments consistent and reproducible.

 - Per-rule Containerization:
  - If a rule cannot run inside the global container (for example conflicting native libraries), use the `singularity:` directive on that rule to point to a per-rule container image.

 - Container Build & Distribution Requirements (remote-first policy):
  - The agent MUST attempt a remote/user-mode build first. After loading the system Singularity/Apptainer module (see HPC Module Handling below), the agent shall run the remote build command (for example `singularity build --remote` or `apptainer build --remote`) as the first attempt.
  - If the remote build command fails because the installed Singularity/Apptainer flavor does not support `--remote`, the agent MUST detect that specific failure and report it. The agent may then offer a documented fallback (for example requesting an older module that supports `--remote`, using a remote build service, or asking the user whether a local `--fakeroot` build is acceptable). The agent must NOT automatically attempt a `--fakeroot` or local root build without explicit user approval.
  - When building locally is explicitly approved by the user, document the exact build commands and any external requirements.

- HPC Module Handling and Fallbacks:
  - On HPC systems that provide environment modules, the agent MUST first attempt to load the system-provided Singularity/Apptainer module (for example `ml load singularity` or `module load apptainer`). After loading the module the agent MUST check whether the binary supports remote builds (for example by running `singularity --version` and testing `singularity build --help` or invoking `singularity build --remote` as a probe). If the loaded module supports remote builds, use it to perform the remote build.
  - If a suitable module isn't available or the loaded binary does not support remote builds, the agent MUST report this and present alternatives to the user (request loading an older module that supports remote builds, use a remote build service, or ask for permission to perform a local `--fakeroot` build). The agent may only install a user-local Apptainer/Singularity via Conda after the user approves that fallback.
  - If installing via Conda, create a named environment (for example `singularity-build`) and document the exact conda commands used so others can reproduce or remove the installation.

By following these Singularity requirements workflows will be portable and robust across shared HPC environments.

B. Conda environments (when chosen)

- Agent prompt & repo requirements:
  - If the user chooses Conda environments, the project MUST include and maintain a Conda environment YAML file in the repository (for example `envs/{config['project']}-py310.yml` or `environment-{config['project']}.yml` rather than a generic `environment.yml`).
  - The agent MUST ensure the YAML exists and is up-to-date. If missing, the agent should create a minimal, documented YAML and add it to the repo, explaining the package choices.

- YAML contents & documentation:
  - The YAML must list channels and packages, and include any pip dependencies under a `pip:` section when needed. Include a top-line comment or README entry that states the intended environment name and exact recreate command.
  - Recommended commands to document in project docs (using the project-prefixed filename above):

    - `conda env create -f envs/{config['project']}-py310.yml` (or `mamba env create -f envs/{config['project']}-py310.yml` when mamba is available)
    - `conda activate {config['project']}-py310`

  - For reproducibility pin major/minor versions where appropriate (for example `python=3.10`) and include channels, for example:
    ```yaml
    channels:
      - conda-forge
      - defaults
    ```
  - For projects requiring stronger reproducibility the agent should recommend or optionally generate a lock file (for example using `conda-lock`) and document how to use it.

- Naming guidance:
  - The Conda environment name and YAML filename must be informative but succinct. Avoid overly generic names such as `env`, `environment`, or `env.yml` that don't convey the project or purpose. Prefer including the project key and, optionally, the Python version or purpose. Examples:

    - YAML: `envs/{config['project']}-py310.yml` or `environment-{config['project']}.yml`
    - environment name: `{config['project']}-py310` or `{config['project']}-analysis`

- Agent responsibilities when using Conda:
  - Ensure the YAML is present and documented in the repo.
  - Create a minimal YAML if requested or if missing, and explain package choices.
  - Keep `rule all` and other workflow outputs consistent with the chosen environment strategy.

These two subsections make the supported choices explicit and separate their requirements so the agent and collaborators can confidently follow the chosen approach.

---

## 3. 🛠️ Write Clean, Modular Snakemake Rules

- Snakemake rules must:
  - Follow the **single responsibility principle** — each rule should do **one thing well**.
  - Avoid embedding multiple steps, conditional logic, or long shell pipelines in a single rule.
  - Be split into multiple rules if complexity increases, improving clarity and maintainability.

- Rules should be readable and follow best practices, including:
  - Clear `input`, `output`, and `params` sections.
  - Descriptive naming.
  - Avoiding unnecessary shell logic inline.

---

## 4. 📝 Snakemake Script Interface Syntax (R & Python)

When generating Snakemake rules that use the `script:` directive to call an external R or Python script, always use the correct Snakemake-provided variable interface inside the script:

- **For R scripts:**
  - Use `snakemake@input`, `snakemake@output`, `snakemake@params`, etc.
  - Example: `snakemake@input[[1]]`, `snakemake@output[[1]]`

- **For Python scripts:**
  - Use `snakemake.input`, `snakemake.output`, `snakemake.params`, etc.
  - Example: `snakemake.input[0]`, `snakemake.output[0]`

Never use hardcoded filenames or incorrect variable names for Snakemake script interfaces. Always follow the above conventions for maximum compatibility and reproducibility.

---

## 5. 🧬 Snakemake R and Python Script Placement

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


## 6. 🗂️ Project-based Output Organization and Naming

**Unless otherwise specified, always write Snakemake pipelines with the expectation that a YAML config file will be passed at runtime (e.g., using `snakemake --configfile config.yaml`).** This approach maintains a clear separation of concerns and maximizes flexibility of inputs. All rules, scripts, and output paths should reference configuration values from the config file, not hardcoded values or wildcards, unless explicitly instructed otherwise.


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

### 6.1 🗃️ Output Directory Structure: Data, Results, and Plots

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

## 7. 🧪 Always Keep `rule all` Up-to-Date

- **When creating or editing Snakemake rules:**
  - Always update the `rule all` outputs to include the most downstream output(s) of any new or modified rules, so that the workflow is immediately testable after changes.
  - If new rules make some previous `rule all` outputs obsolete (i.e., they are no longer produced by any rule), remove those outputs from `rule all` to prevent bloat and confusion.
  - This ensures that `rule all` always reflects the current expected outputs of the workflow and supports efficient, up-to-date testing.

---


## 8. 📝 README: required and recommended contents

- Every project MUST include a top-level README developed alongside the code and workflow, and it MUST be updated when major implementations or revisions are made. Treat the README as the single source of truth for how to use, reproduce, and contribute to the project.

- The README must include at a minimum the following sections:
  - Project Title & Tagline
  - Name + one-line what it does
  - Description / Overview
    - A short paragraph describing the problem it solves and why it exists.
  - Installation / Setup
    - Clear steps to get it running (dependencies, commands). If the project uses Conda or Singularity, include exact recreate/build commands and point to environment files or Singularity definitions.
  - Usage / Examples
    - Concrete examples showing how to run the pipeline or key scripts (commands, example config, expected inputs/outputs, and short snippets). Prefer copy-pasteable commands and one small end-to-end example.

- Additional recommended sections (add as relevant):
  - Quick Start (a minimal, one-command example to run the most common use case)
  - Configuration (explain `config.yaml` keys used by the workflow; reference `project` key expectations)
  - Inputs & Outputs (describe expected input files and layout, and where final/intermediate outputs live; reference the `results/{project}/data|results|plots/` structure from these instructions)
  - Reproducibility / Environment (explain whether Conda or Singularity is used, include environment YAMLs or Singularity build instructions, and document how to rebuild environments)
  - Examples / Screenshots (command examples, example configs, and small screenshots or sample output files for clarity)
  - Development & Testing (how to run unit tests or smoke runs, linting, and how to add new rules)
  - CHANGELOG or Release Notes (link to `CHANGELOG.md` or include a brief history of major revisions; keep the README updated on major changes and link to more detailed changelog when appropriate)
  - Contributing (how to submit changes, code style, PR process, and who to contact)
  - License and Maintainers / Contact
  - Status badges (CI, license, coverage) and short usage examples for common CI checks if present

- README maintenance rules:
  - Update the README for any major feature, interface, or workflow change. Each major update should include a one-line summary and a date either directly in the README or in a `CHANGELOG.md` linked from the README.
  - Keep examples and commands up-to-date with the canonical environment/setup instructions (if environment files change, update both the env file and README simultaneously).

- Minimal README template (suggested order):
  1. Project Title & Tagline
  2. One-line name/description
  3. Table of Contents (for longer READMEs)
  4. Description / Overview
  5. Quick Start
  6. Installation / Setup
  7. Configuration
  8. Usage / Examples
  9. Inputs & Outputs
 10. Development & Testing
 11. CHANGELOG (or link)
 12. Contributing
 13. License & Maintainers / Contact

This README guidance complements the project reproducibility and environment rules in section 2; follow both when documenting environment and run instructions.

## 9. 🧩 Notes

- Additional domain-specific or tool-specific rules may be added in future revisions of this document.
- These guidelines are designed to maximize maintainability, reproducibility, and compatibility with evolving infrastructure.
