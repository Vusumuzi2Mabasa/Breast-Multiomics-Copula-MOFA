"""Generate the research README from audited numerical result tables."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
URL='https://github.com/Vusumuzi2Mabasa/Breast-Multiomics-Copula-MOFA'
SCENARIOS=['regular','skewed','technical','missing','heavy_tail','dependence_only']
LABELS=['Regular','Skewed','Technical variation','Missing proteins','Heavy tails','Dependence change only']


def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(str(x) for x in r)+' |' for r in rows)


def main():
    s=pd.read_csv(ROOT/'results/summary.csv');c=pd.read_csv(ROOT/'results/paired_contrasts.csv')
    d=pd.read_csv(ROOT/'results/dependence_summary.csv');diag=json.loads((ROOT/'results/diagnostics.json').read_text())
    audit=json.loads((ROOT/'results/audit.json').read_text())
    def val(scenario,n,method,metric='auc_mean'):
        return float(s[(s.scenario==scenario)&(s.n==n)&(s.method==method)].iloc[0][metric])
    performance=[];paired=[]
    for scenario,label in zip(SCENARIOS,LABELS):
        performance.append([label]+[f'{val(scenario,360,m):.3f}' for m in ['MOFA+','MOFA+ Gaussian','Copula','Rank PCA','Oracle factors']])
        g=c[(c.scenario==scenario)&(c.n==360)&(c.metric=='auc')&(c.comparator=='MOFA+')].iloc[0]
        paired.append([label,f'{g.difference:+.3f}',f'{g.lower:+.3f} to {g.upper:+.3f}',f'{val(scenario,360,"MOFA+","brier_mean"):.3f}',f'{val(scenario,360,"Copula","brier_mean"):.3f}'])
    diagrows=[[a['method'],a['successful'],a['flagged'],f"{a['median_seconds']:.2f}"] for a in diag]
    deprows=[]
    for _,r in d.iterrows():
        deprows.append([int(r.n),'Null' if r.scenario=='null' else 'Dependence change',r.method,f'{r.rate:.3f}',f'{r.lower:.3f} to {r.upper:.3f}'])
    copdiag=next(a for a in diag if a['method']=='Copula')
    avg=c[(c.metric=='auc')&(c.comparator=='MOFA+')&(c.scenario!='dependence_only')].difference.mean()
    mean_recovery=s[(s.method=='MOFA+')&(s.scenario!='dependence_only')].recovery_mean.mean()
    cop_recovery=s[(s.method=='Copula')&(s.scenario!='dependence_only')].recovery_mean.mean()
    main=rf'''# Copula and MOFA Integration of Simulated Blood Multiomics

**Vusumuzi Mabasa**

**Simulation study.** All molecular measurements, participants, countries and subtype labels are synthetic. No TACA, Clin-MAP, patient or controlled-access data were used. This repository evaluates methods for a proposed blood-based breast cancer multi-omics study; it does not report biomarkers or results from that study.

[Word report](report/Copula_MOFA_Simulation_Report.docx) · [Quarto manuscript](report/manuscript.qmd) · [Numerical results](results/summary.csv) · [Workflow](workflow.txt)

<!-- manuscript -->
## Abstract

Integration of blood-derived genomic, transcriptomic and proteomic measurements requires methods that accommodate different measurement scales without mistaking shared technical variation for biology. We compared an extended-rank-likelihood Gaussian copula factor model with the official MOFA+ engine, a Gaussian MOFA+ sensitivity analysis and rank PCA. The benchmark comprised 360 independently generated training cohorts, each paired with 200 independent test participants, across six mechanisms and two training sizes. Each panel contained 48 features. A separate experiment used 1,000 calibration and 2,000 evaluation datasets to assess a prespecified RNA–protein dependence test. Across the five scenarios containing subtype mean signals, the unweighted mean paired ROC-AUC difference between copula integration and mixed-likelihood MOFA+ was {avg:+.3f}; this aggregate is descriptive and does not establish superiority. The analysis separates factor recovery, subtype prediction and detection of altered dependence. The small panels and synthetic mechanisms support methodological assessment rather than clinical inference.

## Introduction

Patient-matched genomic, RNA and protein data can reveal molecular programmes that are difficult to identify within one assay. In blood-based breast cancer studies, these programmes may reflect immune composition, host response, technical processing or disease-associated biology. Integrating measurements from the same participant does not establish tumour origin or a causal sequence from genomic variation to protein abundance.

MOFA+ represents shared and view-specific variation through latent factors and supports Gaussian, Bernoulli and Poisson observation models [1]. Copula factor models provide an alternative by separating marginal distributions from dependence [2,3]. This separation may be useful for skewed abundance measurements and mixed discrete/continuous data. It does not guarantee better recovery or prediction, and a Gaussian copula does not represent arbitrary nonlinear or asymmetric tail dependence.

We asked whether a manageable copula integration pipeline improves recovery of known molecular programmes and held-out subtype prediction relative to established factor integration. We also examined a distinct hypothesis: whether a prespecified RNA–protein relationship changes between simulated subtypes even when their univariate distributions remain unchanged. These questions correspond to different estimands and are evaluated separately.

## Methods

### Design and measurement structure

The numerical design and source hashes were frozen before evaluation in [the protocol](protocol/design.json) and [source lock](protocol/source_lock.json). Each of six settings was evaluated at training sizes 180 and 360 with 30 independent replicates per cell. Every replicate had a separate test cohort of 200 participants. Thus 360 training/test pairs contributed 2,160 method records; methods within a pair used identical participants and outcomes. Pilot and unit-test seeds were separate from evaluation seeds.

Four matched views contained six germline dosages, six somatic indicators, 24 RNA measurements and 12 proteins. The two genomic views were retained separately. Germline dosages used binomial quantile thresholds with minor-allele frequencies 0.12–0.40, and somatic prevalences ranged from 0.04 to 0.20. These are illustrative marginal frequencies; the panels do not reproduce empirical African linkage disequilibrium, tumour evolution or a clinical sequencing assay.

Three latent programmes generated the dependence backbone, comprising two stronger molecular programmes and one weaker shared programme:

$$
H_i\sim N_3(0,I_3),\qquad
Z_{{ij}}=\frac{{\lambda_j^T H_i+\epsilon_{{ij}}}}{{\sqrt{{1+\|\lambda_j\|^2}}}},\qquad
\epsilon_{{ij}}\sim N(0,1).
$$

The first two loadings alternated across features, with selected opposing signs. The weaker third loading was 0.25 throughout. Under the regular mechanism, RNA and protein were affine Gaussian measurements. Under the skewed mechanism, RNA counts followed negative-binomial margins with dispersion parameter 8 and means spanning approximately 20–400 before library scaling. Library sizes were lognormal around 20 million. RNA was represented by log2(CPM+1), and lognormal protein abundance by its log2 value. This reduced panel is not a whole-transcriptome count matrix; its known simulated library size is used rather than estimating a library size by summing the retained genes.

The technical mechanism added blood-cell composition, country and batch shifts to RNA and protein. Five unnamed synthetic countries had proportions 0.35, 0.25, 0.20, 0.12 and 0.08. They represent imbalanced sampling strata, not named populations or ancestry categories. The missing-protein mechanism additionally removed values with an abundance-dependent probability:

$$
P(M_{{ij}}=1\mid Y_{{ij}})=\operatorname{{expit}}\{{-2-0.8(Y_{{ij}}-\bar Y_j)\}}.
$$

Here $M=1$ indicates a missing measurement. This is non-ignorable missingness; neither fitted integration model explicitly models its detection process. The heavy-tail mechanism used a common chi-square scale mixture with four degrees of freedom, creating a multivariate t dependence structure that misspecifies both fitted Gaussian dependence models.

For five mechanisms, the simulated subtype followed

$$
\operatorname{{logit}}P(T_i=1)=
-1+0.9H_{{i1}}+0.5H_{{i2}}+0.15A_i+0.12(C_i-1.4).
$$

In the dependence-only setting, $T_i\sim\mathrm{{Bernoulli}}(0.30)$ independently of the backbone factors. One prespecified RNA/protein latent pair instead had correlation 0 in non-TNBC and 0.65 in TNBC. Its marginal standard-normal distributions were unchanged. The resulting labels are synthetic subtype proxies, not a mechanistic model of breast cancer.

### Gaussian copula factor integration

For continuous margins, a copula represents a joint distribution as

$$
F(y_1,\ldots,y_p)=C_R\{{F_1(y_1),\ldots,F_p(y_p)\}},\qquad
C_R(u)=\Phi_p\{{\Phi^{{-1}}(u_1),\ldots,\Phi^{{-1}}(u_p);R\}}.
$$

Discrete margins do not identify a unique unrestricted copula. We imposed a latent Gaussian factor model and used the extended rank likelihood, preserving order constraints rather than treating tied genotypes or binary alterations as continuous ranks during fitting [2,3]:

$$
Z_i=\Lambda H_i+e_i,\qquad e_i\sim N_p(0,\Psi),\qquad
\mathcal D(Y)=\{{Z:Y_{{ij}}<Y_{{i'j}}\Rightarrow Z_{{ij}}<Z_{{i'j}}\}}.
$$

Loadings had independent $N(0,0.5)$ priors, residual variances had inverse-gamma priors with shape 2 and scale 0.5, and factors had standard-normal priors. Proper priors regularise latent scales; reported dependence is the scale-free correlation derived from $\Lambda\Lambda^T+\Psi$. Gibbs updates sampled loadings, factors, residual variances and truncated latent observations. Alternating nonadjacent value groups were updated together to preserve strict order constraints, while observations tied at the same value shared bounds without sharing latent values. Missing observations had unconstrained latent updates.

Two independent chains each used 600 burn-in iterations and 900 subsequent iterations, retaining every third draw. A maximum split Rhat of at least 1.10 across 16 fixed, rotation-invariant correlation entries triggered a prespecified rerun with three times these lengths. Flags remaining after extension were retained and reported. This limited diagnostic does not establish convergence of every posterior quantity. The implementation is project-specific Python code, not the published bfa software; its validation includes recovery of a known bivariate correlation and explicit checks of rank constraints.

Prediction used a computationally inexpensive plug-in projection from the posterior mean correlation. Its top three eigenvectors defined a low-rank representation; residual variances were floored at 0.1. Empirical marginal maps were learned on training participants only. Tied categories used midpoint empirical CDF scores at projection, and missing coordinates were omitted. This projection is an approximation and is not exact posterior prediction under the rank likelihood.

### MOFA and comparison methods

MOFA+ was fitted with the official Python engine, mofapy2 0.7.5 [1,4]. The primary configuration used Gaussian germline dosages, Bernoulli somatic indicators and Gaussian transformed RNA/protein measurements. A sensitivity configuration used Gaussian likelihoods for all standardised views. Both used three factors, automatic relevance determination and spike-and-slab weights, with up to 500 variational iterations and the engine's medium convergence criterion. Official convergence messages were recorded. Default engine factor selection was constrained to the same initial factor budget as the copula model; the true budget was known by design and is an optimistic simplification.

Continuous feature means and standard deviations were estimated from training data only. Gaussian-view missing values remained missing. New-participant MOFA scores were conditional MAP estimates using fitted loadings, Gaussian noise precisions and the Bernoulli likelihood, with a standard-normal factor prior. The same MAP convention was used for training scores. This extension is implemented in the benchmark rather than attributed to an official MOFA prediction API.

Rank PCA used training empirical normal scores, with zero imputation on the normal-score scale, followed by three principal components. A covariate-only reference included age, country, blood-cell composition and batch; it is labelled `Clinical only` in raw records for compatibility. An oracle reference used the three true latent factors plus the same covariates. The oracle is a reference for accessible signal and is not a clinical method or guaranteed finite-sample upper bound.

All representations entered an otherwise identical ridge logistic classifier, with prespecified scikit-learn $C=1$ and training-only standardisation. No classifier hyperparameters were selected from test outcomes. The classifier has additive effects and no interaction terms; it is therefore not designed to exploit a pure covariance change. This makes the separate dependence test necessary.

### Evaluation and uncertainty

Backbone factor recovery was the mean squared cosine of the principal angles between centred estimated scores and the three generating backbone factors:

$$
S=\frac1{{3}}\sum_{{k=1}}^3\cos^2(\theta_k).
$$

This score is invariant to factor sign, permutation and rotation. In the dependence-only setting it describes recovery of the original backbone, not recovery of the additional subtype-specific pair relationship. Independent test outcomes were evaluated using ROC AUC, average precision and Brier score,

$$
\mathrm{{Brier}}=\frac1{{n_{{test}}}}\sum_i(\widehat p_i-T_i)^2.
$$

Paired differences were calculated within each replicate. Percentile intervals used 2,000 paired bootstrap resamples of the 30 replicate differences. Plot error bars for mean performance use 1.96 times the Monte Carlo standard error across independent replicates. These are descriptive simulation intervals, not multiplicity-adjusted clinical evidence. Failures and diagnostic flags are reported separately.

### Prespecified dependence test

A companion experiment generated one RNA/protein pair with negative-binomial RNA counts and log-protein measurements. Subtype affected the pair's correlation only. Pearson correlation on transformed abundance and correlation after empirical normal scoring were compared using the absolute Fisher-z difference:

$$
Q=\frac{{|\operatorname{{atanh}}(r_1)-\operatorname{{atanh}}(r_0)|}}
{{\sqrt{{1/(n_1-3)+1/(n_0-3)}}}}.
$$

Rather than relying on a normal approximation, each method and sample size used 500 independent null calibration datasets. The strict rejection threshold was the 476th ascending calibration statistic. Under exchangeability with a new null dataset, the resulting averaged rejection bound is $25/501\approx0.0499$. Fresh evaluation used 500 null and 500 alternative datasets per sample size. Wilson intervals describe evaluation-binomial uncertainty conditional on the realised calibration bank; they do not include all calibration-bank uncertainty. This is a targeted pair test without feature search, not a test of causal regulation or an output of the factor model.

## Results

### Measurement structure and recovery of molecular programmes

The generated panels combine discrete genomic observations with continuous, transformed abundance data. Missing proteins depend on abundance and therefore challenge an ignorable-missingness analysis. The illustration uses a fixed synthetic dataset independent of the evaluation bank.

![Distributions of example synthetic measurements](figures/measurements.png)

Figure 1. One feature from each view in an independent synthetic missing-protein dataset. Protein histograms include observed values only. The regular setting uses continuous transformed-scale measurements, while the remaining settings include simulated RNA counts.

Across the five mean-signal mechanisms, the unweighted mean backbone recovery was {mean_recovery:.3f} for mixed-likelihood MOFA+ and {cop_recovery:.3f} for the copula pipeline. Recovery assesses the complete integration-and-projection pipelines, so differences cannot be attributed solely to the choice of likelihood. Technical shifts and abundance-dependent missingness test whether the estimated factors remain aligned with the generating biological backbone.

![Recovery of the generating factor subspace](figures/factor_recovery.png)

Figure 2. Mean squared principal-angle recovery, with Monte Carlo error bars, from 30 replicates per cell. Higher values indicate closer agreement with the three generating backbone factors. The subtype-specific pair perturbation is outside this backbone in the dependence-only setting.

### Independent subtype prediction

Table 1 reports performance at the larger training size. All methods used the same independent test participants within a replicate. Oracle values quantify performance with access to the true factors, while rank PCA tests whether a simpler transformation-based representation is sufficient.

Table 1. Mean independent-test ROC AUC at training N = 360 across 30 replicates.

'''+table(['Mechanism','MOFA+','MOFA+ Gaussian','Copula','Rank PCA','Oracle factors'],performance)+rf'''

![Independent test prediction across all settings](figures/prediction.png)

Figure 3. Mean ROC AUC and Monte Carlo error bars. The covariate-only reference includes clinical and technical covariates. Chance-level discrimination is 0.5. An additive classifier is deliberately not expected to identify the pure dependence-change mechanism.

Paired contrasts isolate differences between methods from variation in the generated cohorts. A positive ROC-AUC difference favours the copula pipeline. The intervals are based on only 30 independently generated pairs per setting and should be interpreted alongside their magnitude, recovery scores and computation rather than as a ranking of universal superiority.

Table 2. Copula minus mixed-likelihood MOFA+ ROC-AUC difference, paired bootstrap interval and mean Brier scores at training N = 360. Lower Brier scores are preferable.

'''+table(['Mechanism','AUC difference','95% paired interval','MOFA+ Brier','Copula Brier'],paired)+rf'''

![Paired performance differences](figures/paired_auc.png)

Figure 4. Copula minus comparator ROC AUC with 2,000-resample paired percentile intervals. No adjustment for multiple contrasts was applied. The comparison is intended to reveal where differences arise, not to select favourable settings after inspection.

### Dependence estimation and computational diagnostics

The fixed example shows the generating and estimated latent correlation patterns. Some genomic correlations are weakly informed because alteration indicators are sparse. A visually plausible correlation map alone does not establish factor recovery, posterior convergence or biological validity.

![Fixed example correlations and posterior trace](figures/dependence_example.png)

Figure 5. First regular-setting replicate at training N = 180, chosen by design order. Black lines separate germline, somatic, RNA and protein views. The trace displays one of 16 prespecified monitored correlation entries and is not a diagnostic of every model parameter.

The benchmark produced {audit['method_records']:,} method records with {audit['failed_fits']} failed fits. Of {copdiag['successful']} successful copula fits, {copdiag['extended']} required the prespecified sampling extension and {copdiag['flagged']} retained a monitored Rhat flag. Diagnostic flags are retained in the primary summaries, preventing convergence-based selection from concealing difficult settings. Successful computation is not synonymous with a validated model.

Table 3. Successful fits, unresolved diagnostic flags and median elapsed time. Timing includes fitting and projection and was measured during concurrent execution on the local workstation.

'''+table(['Method','Successful fits','Diagnostic flags','Median seconds'],diagrows)+rf'''

![Chain diagnostics and computational cost](figures/diagnostics.png)

Figure 6. Copula maximum split Rhat across the fixed monitored entries and elapsed-time distributions. The dotted diagnostic threshold is 1.10. Runtime boxplots hide outliers for readability; all runtimes remain in the released tables.

### Detecting subtype-specific coordination

The companion experiment tests a dependence change without an abundance shift. Its calibration is separate from factor fitting and subtype classification. A successful dependence test can therefore coexist with chance-level performance of the additive classifier in the dependence-only benchmark.

Table 4. Fresh evaluation of the prespecified pair test after independent calibration. Each row contains 500 evaluation datasets; intervals are conditional Wilson intervals.

'''+table(['N','Mechanism','Statistic','Rejection rate','95% interval'],deprows)+r'''

![Calibrated dependence testing](figures/dependence_test.png)

Figure 7. Null rejection and alternative detection probabilities with Wilson intervals. The alternative changes latent RNA/protein correlation from 0 to 0.65. This relatively strong effect is a feasibility benchmark and does not establish power for smaller biological differences.

## Discussion

This study evaluates whether flexible dependence modelling adds value to a compact multi-omics workflow under specified, reproducible mechanisms. The appropriate interpretation is conditional: a model can recover dependence adequately without improving phenotype prediction, and a representation that predicts well need not recover the complete generating factor space. Differences between the copula and MOFA pipelines also include priors, sparsity, inference and score projection, preventing attribution of every difference to the copula alone.

The Gaussian MOFA sensitivity analysis and rank PCA are important controls. When they perform similarly to a copula, normalisation and low-rank structure may explain much of the available gain. Conversely, a dependence test addresses information that an additive classifier may miss. It should be evaluated as a separate, prespecified scientific analysis rather than used to rescue a disappointing classification comparison.

For the proposed breast cancer Study 3, MOFA+ remains an established integration baseline. A copula analysis could complement it by examining selected cross-layer relationships, especially when continuous distributions remain non-Gaussian or dependence itself is the target. Any extension must account for actual specimen fractions, sequencing design, cell composition, treatment timing, assay batches and country imbalance. Molecular co-variation cannot establish tumour origin, causal mediation, shared causal variants or therapeutic utility.

Several limitations constrain extrapolation. The feature panel is deliberately small, genomic dependence is synthetic, and the factor count is known. The model does not reproduce real African LD, full tumour mutation processes or a whole-transcriptome assay. The missingness mechanism is not explicitly fitted. Rank likelihood reduces sensitivity to monotone marginal transformations but does not solve confounding or non-ignorable missingness. A Gaussian copula is misspecified for the heavy-tail dependence mechanism. The plug-in projection may discard posterior information and may disadvantage the copula relative to exact predictive inference.

Thirty integration replicates per cell provide a manageable benchmark, not precise evidence of small power differences. The companion test has more replicates but examines only one prespecified pair and a strong alternative. It does not address genome-wide multiplicity. No leave-one-country-out validation, empirical reference-data validation, full posterior calibration, asymmetric copula comparison or extensive factor-number sensitivity analysis was performed. These are prospective extensions rather than completed analyses.

## Conclusion

Copula modelling can be tested as a complementary integration strategy for mixed blood multi-omics measurements. Its value must be judged separately for latent programme recovery, prediction and subtype-specific coordination. This benchmark provides executable evidence for those distinctions while retaining MOFA+ and simpler baselines. The simulated results do not justify replacing an established integration method solely because the alternative is more complex.

## Reproducibility

The project uses Python 3.12, the official mofapy2 engine, a documented project-specific copula sampler and independent R checks of selected metrics. All evaluation records, compact synthetic fixtures, figures, calibration thresholds and validation summaries are released. Work directories, installed dependencies and document-rendering intermediates are excluded.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m unittest discover -s tests -v
python src/run_all.py --workers 2
python src/analyse.py
python src/validate.py
Rscript --vanilla src/check_metrics.R .
python src/manuscript.py
python src/build_word.py
```

On Windows, activate `.venv\Scripts\Activate.ps1` and set the three thread environment variables using PowerShell. The Toshiba run used a task-local dependency directory rather than changing global Python packages. See [workflow.txt](workflow.txt) for the complete Windows and Linux procedures. The committed source lock must be retained; `src/freeze.py` is only for a genuinely new, prospective design and refuses to overwrite an existing lock. Existing replicate JSON files are resumed rather than recomputed. Use a separate checkout/output folder for a fresh full rerun.

Word generation requires Quarto; independent metric verification requires R and jsonlite. The manuscript is generated from the same result tables and prose as this README, with editable equations. No GitHub access is required to reproduce the numerical analyses. GitHub publishing instructions are included in the workflow file.

## References

1. Argelaguet R, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biology. 2020;21:111. [Methods paper](https://doi.org/10.1186/s13059-020-02015-1).
2. Murray JS, Dunson DB, Carin L, Lucas JE. Bayesian Gaussian copula factor models for mixed data. Journal of the American Statistical Association. 2013;108:656–665. [Methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3753118/).
3. Hoff PD. Extending the rank likelihood for semiparametric copula estimation. Annals of Applied Statistics. 2007;1:265–283. [Author manuscript](https://arxiv.org/abs/math/0610413).
4. MOFA developers. Official implementation and installation documentation. [MOFA documentation](https://biofam.github.io/MOFA2/installation.html).

## Data and code availability

All released observations are synthetic. The repository contains no proposal document, personal contact information, clinical records, raw patient sequencing or collaborator data. The source is released under the MIT licence. Cite the repository and the original MOFA and copula-method publications when reusing the benchmark. Implementation of an existing method and its application to this simulation are not claims of methodological novelty.
'''
    (ROOT/'README.md').write_text(main,encoding='utf-8')
    print('Research README generated from audited tables.')


if __name__=='__main__':main()
