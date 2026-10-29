"""Research manuscript generated from complete audited Study 3 results."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent


def table(headers,rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in rows)


def main():
    read=lambda n:pd.read_csv(ROOT/f'results/{n}.csv',keep_default_na=False,na_values=[''])
    d=read('dependence_summary');i=read('integration_summary');raw=read('integration_evaluation');di=read('diablo_evaluation');contr=read('dependence_contrasts')
    a=json.loads((ROOT/'results/dependence_audit.json').read_text());ia=json.loads((ROOT/'results/integration_audit.json').read_text())
    rr=json.loads((ROOT/'results/validation/R_conditional.json').read_text())
    methods=['Conditional copula','Product regression','Normal-score difference','Partial Spearman','Oracle margins']
    nulls=['null','abundance_only','scale_only','technical_null','missing_null','heavy_tail_null'];nl=['Baseline conditional null','Abundance only','Scale only','Technical variation','Missing proteins','Heavy-tail stress']
    alts=['coordination_moderate','coordination_strong','both','opposing_coordination','missing_coordination'];al=['Moderate coordination','Strong coordination','Abundance and coordination','Opposing changes','Missing and coordination']
    ims=['abundance_only','coordination_moderate','both','missing_coordination'];il=['Abundance only','Coordination only','Both','Missing and both']
    def dep(n,s,m,metric='module_detection'):
        return float(d[(d.n==n)&(d.scenario==s)&(d.method==m)].iloc[0][metric])
    def integ(n,s,m,metric='auc'):
        return float(i[(i.n==n)&(i.scenario==s)&(i.method==m)].iloc[0][metric])
    nr=[[label]+[f'{dep(480,s,m,"false_module"):.3f}' for m in methods] for s,label in zip(nulls,nl)]
    ar=[[label]+[f'{dep(480,s,m):.3f}' for m in methods] for s,label in zip(alts,al)]
    ir=[[label]+[f'{integ(480,s,m):.3f}' for m in ['MOFA+ tuned','MOFA+ rank tuned','Gaussian Gibbs matched','Copula Gibbs matched','Rank PCA']] for s,label in zip(ims,il)]
    ba=[]
    for s,label in zip(ims,il):
        dg=di[(di.n==480)&(di.scenario==s)&(di.status=='ok')]
        ba.append([label,f'{dg.balanced_accuracy.mean():.3f}']+[f'{integ(480,s,m,"balanced_accuracy"):.3f}' for m in ['MOFA+ tuned','Copula Gibbs matched','Covariates only']])
    diagnostics=[]
    for m in ['MOFA+ tuned','MOFA+ rank tuned','Gaussian Gibbs matched','Copula Gibbs matched']:
        g=raw[(raw.method==m)&(raw.status=='ok')]
        diagnostics.append([m,len(g),int((g.diagnostic_pass==False).sum()),f'{g.seconds.median():.1f}',int((g.k==10).sum())])
    control=d[(d.method=='Conditional copula')&d.scenario.isin(nulls[:-1])]
    maxerr=max(v['max_absolute_error'] for v in rr['checks'])
    provisional=len(list((ROOT/'provisional_pca/integration').glob('*.json')))
    s=rf'''# Cross Omics Coordination and Integration in Simulated Breast Cancer

**Vusumuzi Mabasa**

**Simulation study.** All participants, genomic features, gene modules, measurements, countries and TNBC labels are synthetic. No patient, TACA, Clin-MAP or controlled-access data were used. The work supports methodological planning for Study 3 of a proposed blood-based breast cancer PhD; it does not identify clinical biomarkers or tumour mechanisms.

[Word report](v2/report/Study3_Coordination_Report_v2.docx) · [Quarto manuscript](v2/report/manuscript.qmd) · [Small synthetic example](v2/data/example) · [Study 3 protocol](v2/protocol/design.json) · [Complete workflow](workflow.txt) · [Preserved original study](baseline/README_v1.md)

<!-- manuscript -->
## Abstract

Matched blood-derived genomic, transcriptomic and proteomic profiles offer opportunities to investigate molecular coordination as well as abundance. We evaluated a covariate-adjusted Gaussian copula score test for subtype-specific RNA–protein coordination using biologically structured synthetic gene modules. The dependence experiment comprised {a['calibration_datasets']:,} calibration and {a['evaluation_datasets']:,} fresh evaluation cohorts, with training-independent family thresholds, simple adjusted-correlation benchmarks and latent-margin oracles. A distinct integration experiment used 96 independent replicates to compare tuned MOFA+, rank-transformed MOFA+, matched-prior Gaussian and copula factor samplers, rank PCA and an official supervised DIABLO analysis. Independent test participants included a synthetic country excluded from training. At N = 240, moderate coordination detection was {dep(240,'coordination_moderate','Conditional copula'):.1%} for the conditional copula and {dep(240,'coordination_moderate','Product regression'):.1%} for adjusted product regression; corresponding N = 480 values were {dep(480,'coordination_moderate','Conditional copula'):.1%} and {dep(480,'coordination_moderate','Product regression'):.1%}. Realised copula null module rejection ranged from {control.false_module.min():.1%} to {control.false_module.max():.1%} across the included null settings and sizes. Simpler adjusted correlation tests exceeded the copula in some settings. The experiment distinguishes detection of altered coordination from subtype prediction and latent programme recovery. It supports further evaluation of targeted dependence analysis alongside established integration, without claiming universal superiority, causal regulation or clinical validation.

## Introduction

Study 3 of the proposed PhD concerns matched genomic, RNA and protein measurements in blood from breast cancer cases, comparing TNBC with non-TNBC disease. A useful integration framework should identify reproducible molecular programmes while distinguishing disease-associated coordination from technical variation and differences in individual abundance. Blood-derived signals can reflect immune composition, host response, treatment or assay processing; matching participants does not establish tumour origin or a causal molecular chain.

The original repository benchmark found a small average prediction gain for its copula pipeline, approximately 0.003 ROC AUC. Those results remain preserved. They motivate a distinct experiment with a different primary hypothesis: whether RNA–protein dependence changes by subtype after adjustment for measured genomic and clinical features. Improved classification is not required for such a relationship to be scientifically interesting.

MOFA+ is an established unsupervised approach to shared and modality-specific molecular variation [1]. DIABLO is a supervised multiblock signature method [2]. Neither is automatically a hypothesis test of conditional dependence. Copula models separate marginal distributions from dependence [3,4], providing a framework for that additional question. The present work therefore keeps programme discovery, supervised prediction and differential coordination as related but separate tasks.

## Methods

### Experimental structure and prospective controls

The original study was preserved at commit `c56089d384bc3b447c2d6f4601c8ccdf5cd1e28c`. Original scientific sources, results, figures and the 13-page Word report remain in place. Archived copies of the original README, workflow and manifest are in [baseline](baseline). New sources and outputs are confined to [v2](v2), with versioned source locks.

The primary dependence experiment evaluated six statistics at N = 240 and 480. Five complete-null mechanisms each supplied 500 independent calibration cohorts per N. Separate evaluation used 300 cohorts per null mechanism, 200 per alternative and 200 per out-of-bank stress mechanism. Calibration thresholds were frozen before evaluation. The integration experiment comprised four mechanisms, two training sizes and 12 independent replicates per cell. Each replicate had a 240-participant independent test set from represented countries and a separate 120-participant test set from a synthetic country excluded from training.

Table 1. Independent simulation components. Different methods within a replicate share participants and outcomes.

'''+table(['Component','Design','Independent replicates'],[
['Dependence calibration','5 null mechanisms × 2 sizes × 500','5,000 cohorts'],
['Dependence evaluation','5 null × 300; 5 alternatives × 200; 2 stress × 200; each at 2 sizes','5,800 cohorts'],
['Integration','4 mechanisms × 2 sizes × 12','96 training and two-test-set replicates'],
['Supervised DIABLO','Same 96 integration replicates; internal threefold validation','96 fitted workflows']])+rf'''

### Biologically structured measurements

Each panel contained eight germline dosages, eight somatic indicators, 16 RNA measurements and 16 matched protein measurements. Four prespecified synthetic modules contained four gene–protein pairs each. Germline dosages were generated in four two-variant latent Gaussian blocks, with country-dependent allele frequencies. This supplies illustrative LD and frequency variation, not an empirical African reference panel. Somatic indicators were generated separately and could be associated with subtype in abundance-change settings.

For gene j, the RNA and protein location functions contained its mapped germline dosage, a binary somatic alteration and prespecified clinical terms. Most genetic effects had concordant directions, with selected protein effects reversed to represent discordance. These equations define simulation truth and are not causal conclusions from observed data. Genomic classification, read processing and somatic calling were not performed.

Four latent module variables generated correlated RNA residual programmes:

$$
H_i\sim N_4(0,I),\qquad
Z^R_{{ij}}=0.45H_{{i,m(j)}}+\sqrt{{1-0.45^2}}\,e_{{ij}},
$$

where m(j) denotes a prespecified module. Protein residuals were generated through

$$
Z^P_{{ij}}=\rho_{{ij}}Z^R_{{ij}}+\sqrt{{1-\rho_{{ij}}^2}}\,\xi_{{ij}},
\qquad
\rho_{{ij}}=\tanh\{{\eta_0(X_i)+\beta_jT_i\}}.
$$

The baseline predictor was $\eta_0(X_i)=\operatorname{{atanh}}(0.2)+0.25\,\mathrm{{cells}}_i+0.15\,\mathrm{{age}}_i+0.08(C_i-1.4)$. Covariates included age, cell composition, treatment, country and batch. The genotype- and covariate-dependent subtype sampling model produced approximately one-third TNBC cases. All participants were simulated breast cancer cases; the comparison does not estimate population breast cancer susceptibility.

RNA used negative-binomial quantile measurements with a known simulated whole-assay library size followed by log2(CPM+1). Optional latent-scale perturbations altered its conditional distribution. Protein was represented on a log-abundance scale. Technical settings added country and batch shifts. Missing-protein settings imposed abundance-dependent detection, creating non-ignorable missingness. The exact generator, assay library sizes and feature map are available in the source and small example data.

![Structure of the synthetic gene modules](v2/figures/biological_design.png)

Figure 1. Simulation relationships among distinct germline, somatic, RNA and protein views. Arrows describe the generator, not a causal interpretation of patient data. Module 1 has altered conditional RNA–protein coordination in alternatives; modules 2–4 remain unchanged.

### Null mechanisms and alternatives

The calibration bank included no direct subtype effect, abundance-only shifts, scale-only shifts, technical confounding and abundance shifts with missing proteins. All have generating $\beta_j=0$. Unadjusted group means can nevertheless differ because subtype sampling depends on genotype and clinical composition. This prevents a favourable result from depending solely on a clean, balanced null. Alternatives modified the first four gene–protein pairs: moderate Fisher-scale shifts $\beta=0.35$, strong shifts $\beta=0.65$, abundance plus moderate coordination, alternating positive/negative shifts, and moderate coordination with missing proteins. A Fisher-scale shift is not a constant difference in correlation: its effect depends on the baseline covariates.

Heavy-tailed null and alternative mechanisms used a shared t-scale mixture with five degrees of freedom. They were excluded from calibration and retained as out-of-bank stress tests. No threshold was retuned in response to their results.

### Covariate-adjusted dependence analysis

Each pair's marginal adjustment included subtype, clinical and technical terms, its mapped germline dosage and somatic indicator, and prespecified subtype interactions. Twofold cross-fitting estimated location and log-residual-scale models without fitting the same observations used for the residual prediction. Within-subtype empirical normal scores supplied a semiparametric approximation to the remaining margins. Ties received midranks. Complete observed pairs were used for testing; the detection process itself was not fitted. Consequently, calibration under a specified missingness mechanism is not a general solution to non-ignorable missingness.

The working conditional dependence hypothesis was

$$
\operatorname{{atanh}}\{{\rho_i\}}=D_i^T\gamma+\beta T_i,
\qquad H_0:\beta=0.
$$

Here D includes an intercept and the clinical/technical dependence covariates. This is a working observed-pair parameterisation. Under informative protein selection or misspecified margins, its coefficient need not equal the complete-data generating coefficient. Reported detection/error classifications use the generating truth; missing-protein calibration therefore addresses a specified latent-null selection mechanism rather than proving identifiability of latent dependence in real MNAR data. With normal-score residuals u and v, the Gaussian log likelihood, up to marginal constants, is

$$
\ell_i=-\tfrac12\log(1-\rho_i^2)
-\frac{{u_i^2+v_i^2-2\rho_iu_iv_i}}{{2(1-\rho_i^2)}}.
$$

The nuisance parameters were fitted under the null by bounded numerical optimisation with an analytic gradient. The Fisher-scale score contribution was

$$
s_i=\rho_i+
\frac{{(1+\rho_i^2)u_iv_i-\rho_i(u_i^2+v_i^2)}}{{1-\rho_i^2}}.
$$

The subtype indicator was residualised against D using weights $w_i=1+\rho_i^2$. The studentised efficient-score statistic was

$$
Z_\beta=\frac{{\sum_i\widetilde T_i s_i}}
{{\sqrt{{\sum_i(\widetilde T_i s_i)^2}}}},\qquad Q=Z_\beta^2.
$$

Its reference distribution was established empirically, rather than assumed to be an exact finite-sample chi-square distribution. Estimated margins, working dependence assumptions and complete-pair selection remain important qualifications.

Five comparators were evaluated on the same datasets: adjusted Pearson correlation differences; adjusted Spearman differences; normal-score correlation differences using the same location/scale preparation; HC3-adjusted product regression; and an oracle score test given the unobserved generating RNA/protein residuals. Product regression included clinical terms and quadratic age/cell terms. It is a serious conditional-moment benchmark, but its functional form differs from the copula model. The oracle is unavailable in real data and is not a guaranteed upper bound under selection or misspecification.

### Multiplicity and calibration

The primary module statistic was the sum of Q over a module's four pairs. Calibration used the maximum statistic across the four modules. A secondary gene family used the maximum Q across 16 pairs. For each method, N and null mechanism, the strict cutoff was the 476th ascending value among 500 family maxima; the released cutoff was the largest across the five null banks. Under exchangeability with an included complete-null bank, the averaged rejection bound is $25/501\approx0.0499$.

This bound averages over calibration-bank randomness. It does not force every realised bank to have exactly 5% conditional error, establish strong family-wise control for arbitrary mixed alternatives or guarantee calibration in new clinical data. False discoveries among unchanged modules and genes were evaluated explicitly in the alternatives. The gene and module families are separate analyses; there is no joint 5% guarantee across both families and no general FDR claim.

Wilson intervals describe evaluation uncertainty conditional on the realised bank. Selected method contrasts additionally used 250 paired bootstrap resamples of both calibration cohorts and evaluation cohorts, recomputing the family cutoff each time. Bootstrap intervals can be degenerate at complete detection or no detection; they are descriptive rather than evidence of equivalence.

### Fairer integration and supervised benchmarks

Integration used all four molecular views. Technical adjustment excluded subtype labels, preserving the separation between unsupervised discovery and supervised inference. Transformation parameters were learned from training participants only. MOFA+ used its official Python engine [1,5], with mixed Gaussian/Bernoulli likelihoods for standard inputs and Gaussian likelihoods for rank-transformed inputs.

For a Gaussian MOFA view, $Y^{{(m)}}=ZW^{{(m)T}}+E^{{(m)}}$; the Bernoulli view uses $P(Y_{{ij}}^{{(m)}}=1)=\operatorname{{expit}}(Z_iW_j^{{(m)T}})$. MOFA's relevance and sparsity priors are distinct from the matched-sampler priors below, so its comparison remains a workflow comparison.

Factor budgets 2, 4, 6 and 10 were evaluated by training-only masked RNA/protein reconstruction. Fifteen percent of observed RNA/protein training entries were masked; marginal preprocessing for this stage excluded those entries. Starting states were selected by training ELBO. The selected budget was refitted on the full training data. The initial pilot selected the boundary of the smaller grid, prompting the uniform prospective [factor-grid amendment](v2/protocol/integration_amendment.json).

Source inspection also showed that seed-only reruns can share MOFA's deterministic PCA start. The final comparison therefore used an official PCA start and a distinct PCA start perturbed with fixed-scale Gaussian noise in the variational factor means, with expectations refreshed before training. Initialization hashes document this distinction. The official engine performed subsequent inference. The {provisional} completed provisional seeded-PCA fits were retained in [the provisional archive](v2/provisional_pca), excluded from final evaluation, and replaced uniformly according to [the initialization amendment](v2/protocol/initialization_amendment.json).

Two additional factor samplers isolated an important source of unfairness: they used identical Gaussian loading priors, inverse-gamma residual-variance priors, factor budgets, sampling lengths and diagnostic rules, differing in Gaussian observations versus a Gaussian-copula extended rank likelihood [3,4]. The factor budget matched the raw MOFA selection. Each used two chains, 600 burn-in and 900 subsequent iterations, thinning by three. A monitored split Rhat of at least 1.10 triggered one threefold-length extension. Flags were retained, not used to select favourable datasets.

Their common factor model was $L_i=\Lambda H_i+\epsilon_i$, with $H_i\sim N_K(0,I)$, $\epsilon_i\sim N(0,\operatorname{{diag}}\psi)$, $\Lambda_{{jk}}\sim N(0,0.5)$ and $\psi_j\sim\operatorname{{InvGamma}}(2,0.5)$ in the shape/scale parameterisation. In the Gaussian arm, L is the standardised observed matrix. In the copula arm, L is latent and constrained by the observed feature orders, including tied discrete genomic observations. Both update the same factor-model blocks, with the latter additionally updating the rank-constrained latent observations.

Rank PCA and a covariate-only classifier provided simpler references. The true generating module factors were an oracle for programme capture, not an oracle classifier: subtype also changes measurements through mechanisms outside those factors. All practical representations entered the same prespecified ridge logistic classifier. Scores were projected for new participants using fixed fitted quantities. Copula projection used empirical marginal maps and a plug-in latent correlation, which remains an approximation.

Full module capture was measured by

$$
S=\frac1{{4}}\|Q_{{\widehat H}}^TQ_H\|_F^2,
$$

with absent estimated dimensions contributing zero. Cross-omics prediction withheld the complete protein view from score inference for both training and test participants, then predicted observed standardised protein values using a common ridge procedure. A covariate-only protein reference used the same train-fitted adjustment and scaling. Missing test proteins were not treated as known targets; this evaluation does not establish accuracy for unobserved MNAR values.

DIABLO was fitted through official R mixOmics [2], with the four molecular blocks and an additional covariate block. Threefold internal validation selected one or two components and four or eight retained features per block, subject to availability. Imputation and technical adjustment were repeated inside the training folds. Official weighted-centroid voting supplied class predictions. Balanced accuracy was reported; hard class labels were not presented as probabilistic ROC scores. Both supervised and unsupervised methods used the same independent seen- and unseen-country outcomes.

Integration uncertainty used 2,000 paired bootstrap resamples across the 12 independently generated replicates per cell. These are small simulation comparisons, not multiplicity-adjusted clinical claims. The held-out fifth country is synthetic distribution-shift assessment, not external validation in an African population.

## Results

### Null controls distinguish coordination from abundance

Table 2 reports realised module-family rejection at the larger sample size. The abundance-only and scale-only controls directly test whether changes in individual distributions are mistaken for altered coordination. All methods received independently calibrated thresholds; they were not compared using unmatched nominal cutoffs.

Table 2. Probability of any false module at N = 480. Each included null has 300 fresh evaluation cohorts; the heavy-tail stress setting has 200. No evaluation outcomes were used to change thresholds.

'''+table(['Mechanism','Copula','Product','Normal scores','Spearman','Oracle'],nr)+rf'''

Across included null cells and both sizes, the conditional copula's point estimates ranged from {control.false_module.min():.1%} to {control.false_module.max():.1%}. In particular, the missing-protein N = 240 cell had {dep(240,'missing_null','Conditional copula','false_module'):.1%} rejection. That result remains visible and qualifies any claim of calibration. The finite-bank guarantee is averaged rather than conditional, and the evaluation intervals are important when interpreting individual cells. Neither missingness nor model misspecification should be treated as solved by these simulations.

![Null family behaviour across all methods](v2/figures/null_control.png)

Figure 2. Any false module or false gene, with conditional 95% Wilson intervals. The dotted line is 5%. The heavy-tail setting is outside the calibration bank. The two families are reported separately.

### Altered molecular coordination can be detected without classification gains

At N = 240, moderate module detection was {dep(240,'coordination_moderate','Conditional copula'):.1%} with the conditional copula and {dep(240,'coordination_moderate','Product regression'):.1%} with product regression. At N = 480, detection increased to {dep(480,'coordination_moderate','Conditional copula'):.1%} and {dep(480,'coordination_moderate','Product regression'):.1%}. These differences concern an explicitly defined conditional-coordination target, not the power of MOFA+ or DIABLO to test that hypothesis. Strong alternatives approach a power ceiling, making the moderate and missing-data cases more informative.

The advantage over product regression did not extend to all simpler tests. For moderate coordination at N = 240, adjusted Pearson detected the altered module in {dep(240,'coordination_moderate','Partial Pearson'):.1%} of cohorts, Spearman in {dep(240,'coordination_moderate','Partial Spearman'):.1%} and the normal-score difference in {dep(240,'coordination_moderate','Normal-score difference'):.1%}, compared with {dep(240,'coordination_moderate','Conditional copula'):.1%} for the copula. All methods are shown in Figure 3 and the released tables. These results support retaining the simple adjusted-correlation benchmarks rather than selecting the copula by complexity.

Table 3. Detection probability for the altered module at N = 480, with 200 independent evaluation cohorts per mechanism.

'''+table(['Mechanism','Copula','Product','Normal scores','Spearman','Oracle'],ar)+r'''

![Detection of the prespecified altered module](v2/figures/module_detection.png)

Figure 3. Module detection with Wilson intervals. Moderate and strong effects are Fisher-scale shifts of 0.35 and 0.65; they are not fixed changes in Pearson correlation. Opposing changes alter signs across genes within the same module.

The paired contrasts resample both calibration and evaluation banks. A positive difference favours the copula score. Interpretation remains conditional on the chosen mechanisms: the copula's correlation link is correctly specified in the primary generator, whereas the product-regression benchmark uses a different conditional-mean form. The results do not establish superiority for arbitrary dependence, arbitrary marginal estimation or data-derived pathways.

![Calibration-aware paired method differences](v2/figures/paired_detection.png)

Figure 4. Copula minus comparator module-detection probability with 250-resample paired bootstrap intervals. Family cutoffs are recomputed within each draw. Boundary intervals may be degenerate, and contrasts are not multiplicity-adjusted evidence of general superiority.

### Model selection and convergence remain part of the comparison

MOFA factor budgets were selected without test outcomes or subtype-informed feature selection. Figure 5 shows the selected budgets, monitored sampler behaviour and total workflow cost. Runtime includes candidate selection and refitting, so it describes the implemented workflows rather than equal-cost microbenchmarks. Gaussian and copula Gibbs arms share priors and factor budgets; neither is automatically equivalent to MOFA's sparsity assumptions or variational approximation.

Table 4. Integration diagnostics across 96 replicates per method. Flags are unresolved convergence diagnostics after the prespecified procedure. Boundary selections at K = 10 are reported rather than followed by additional outcome-driven tuning.

'''+table(['Method','Successful','Flags','Median seconds','K = 10'],diagnostics)+rf'''

![Integration tuning and computational diagnostics](v2/figures/integration_diagnostics.png)

Figure 5. Factor-budget choices for the two MOFA configurations, monitored split Rhat and workflow runtimes. Runtime boxplots suppress outliers visually; all timings remain in the numerical tables. The monitored entries do not establish convergence of every posterior quantity.

### Integration and prediction answer complementary questions

The integration experiment produced {ia['method_records']} method records, with {ia['fit_failures']} failed fits and {ia['diagnostic_flags']} unresolved selected-fit diagnostic flags. Table 5 compares independent-test discrimination at N = 480. Programme capture, protein-view prediction and supervised balanced accuracy provide additional perspectives: a method should not be ranked from AUC alone, nor from recovery of latent factors unrelated to subtype.

All successful fits remain in the primary results. A separate [diagnostic sensitivity table](v2/results/diagnostic_sensitivity.csv) restricts each paired comparison to replicates where both selected-fit diagnostics passed and reports the retained sample counts. That subset can itself be selected and is not an unbiased replacement for the full bank. Unresolved sampler or variational diagnostics therefore continue to qualify the comparison.

Table 5. Mean independent-test ROC AUC at N = 480 across 12 replicates per cell.

'''+table(['Mechanism','MOFA raw','MOFA ranks','Gaussian Gibbs','Copula Gibbs','Rank PCA'],ir)+r'''

DIABLO is a supervised signature benchmark, with internal tuning and the same independent outcomes. Its class predictions permit balanced-accuracy comparisons but do not supply the same probability scale as the ridge classifiers. Apparent differences therefore concern complete predictive workflows, not an isolated copula-versus-supervision effect.

Table 6. Mean independent-test balanced accuracy at N = 480. The covariate reference includes the same clinical and technical inputs available to the other classifiers.

'''+table(['Mechanism','DIABLO','MOFA raw','Copula Gibbs','Covariates'],ba)+r'''

![Complementary integration endpoints](v2/figures/integration_results.png)

Figure 6. Independent-test discrimination, generating-module capture, complete protein-view prediction and balanced accuracy at N = 480. Error bars are Monte Carlo intervals across 12 replicates. Protein error is evaluated on a common observed standardised scale; lower is preferable. The purely generating-factor reference is not a classifier upper bound.

### Synthetic country transportability

The fifth country was excluded from all integration training and model selection. Figure 7 compares independent-test AUC in that stratum with AUC in represented countries. The comparison includes changes in allele frequencies, covariate composition and country effects under the same generator. It is a useful stress test of portability, but it does not reproduce the breadth of real African population, clinical or platform differences.

![Independent synthetic-country transportability](v2/figures/country_transport.png)

Figure 7. Mean unseen-minus-seen-country ROC AUC with paired bootstrap intervals for each practical integration workflow. Negative values indicate lower discrimination in the held-out stratum. The plot is descriptive, not evidence of external population validation.

## Discussion

This extension makes the connection to Study 3 more specific. MOFA+ provides an established route to discovering shared molecular variation, DIABLO evaluates supervised cross-omics signatures, and the conditional copula addresses a separate question about subtype-dependent coordination. A dependence signal can therefore add scientific information without substantially increasing additive subtype-classification AUC.

The detection gains over product regression in the specified moderate-coordination settings support further investigation of that targeted question. Higher detection by simpler adjusted correlation tests in some cells prevents a general efficiency claim for the copula. Their interpretation has limits. The dependence model has the correct link form in the main generator; alternative methods use different transformations, moments and efficiency assumptions. Matched-prior factor samplers help separate marginal/dependence modelling from prior and factor-budget differences, but the main MOFA comparisons remain comparisons of complete pipelines. Neither a favourable cell nor a strong-signal ceiling establishes a universally better method.

The negative controls are central. Abundance, heteroscedasticity, technical variation and non-ignorable protein detection can all complicate the interpretation of a coordination analysis. Independent family calibration addresses those specified simulation mechanisms; it does not certify validity for unknown missingness or incorrect covariate models in a real cohort. Realised false-positive rates, uncertainty and out-of-bank stress results remain part of the evidence, including cells exceeding the nominal point estimate.

For the PhD, a defensible analysis would prespecify a limited set of gene–protein pairs or pathways independently of the confirmation data. Modules discovered using MOFA+ or DIABLO should be evaluated in held-out participants or with a procedure that accounts for their selection. The current prespecified-module benchmark does not validate testing thousands of data-derived pairs on the same participants used to select them. Cross-country effect heterogeneity remains principally a Study 2 question; the present transport experiment supplies a complementary sensitivity analysis.

Further limitations are deliberate. Only 48 molecular features were simulated. Germline LD is synthetic, somatic alterations are simplified binary inputs, and no real tumour fractions, copy-number processes or cell-deconvolution errors are represented. The observed cell-composition covariate is idealised. Genomic view availability and biological interpretation must follow the actual assay manifest, including defensible somatic classification. Blood RNA and proteins may come from different compartments; cross-layer association does not establish tumour origin or regulatory direction.

Integration uses plug-in score projection, a limited factor grid and only 12 replicates per cell. The null score tests use estimated normal-score margins rather than an exact arbitrary conditional-CDF model. Sampler diagnostics monitor a restricted set of invariant quantities, not full posterior calibration. No real-data biomarker replication, mechanistic validation or clinical utility assessment was performed. These constraints define what the benchmark supports and what remains for future study.

## Conclusion

The more relevant contribution for Study 3 is a calibrated, explicitly conditional analysis of cross-omics coordination alongside established programme discovery and supervised signature analysis. The simulation provides evidence about this defined target under structured and challenging mechanisms, while preserving simpler comparators and the original small-gain study. It supports a complementary role for copulas rather than their automatic replacement of MOFA+.

## Reproducibility

The original study remains available through the archived README and commit. New data generation, calibration, evaluation, integration selection, DIABLO, validation, figures and manuscripts are implemented under `v2/src`. Scientific model sources were frozen before their evaluation phases. The factor-grid and initialization amendments are explicit, and provisional fits are not pooled with final results. Failed fits and diagnostic flags are retained in the released records.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m unittest discover -s tests -v
python -m unittest discover -s v2/tests -v
python v2/src/run_dependence.py --workers 2
python v2/src/run_integration.py --workers 2
Rscript --vanilla v2/src/diablo.R .
python v2/src/validate.py
Rscript --vanilla v2/src/check_conditional.R .
python v2/src/validate_integration.py
Rscript --vanilla v2/src/check_integration.R .
python v2/src/analysis.py
python v2/src/figures.py
python v2/src/manuscript.py
python v2/src/build_word.py
python v2/src/audit.py --fresh
```

DIABLO requires R and mixOmics; independent checks additionally use jsonlite. The tested software versions are recorded in the released metadata. Quarto renders the Word report from the shared manuscript, with editable equations. Existing simulation chunks are resumed. For a fresh run, use the supplied reproduction wrapper with a new output directory so original results are not overwritten. Windows/Toshiba commands and GitHub publication details are in [workflow.txt](workflow.txt).

## References

1. Argelaguet R, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biology. 2020;21:111. [Methods paper](https://doi.org/10.1186/s13059-020-02015-1).
2. Singh A, et al. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. Bioinformatics. 2019;35:3055–3062. [Methods paper](https://doi.org/10.1093/bioinformatics/bty1054).
3. Murray JS, Dunson DB, Carin L, Lucas JE. Bayesian Gaussian copula factor models for mixed data. Journal of the American Statistical Association. 2013;108:656–665. [Methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3753118/).
4. Hoff PD. Extending the rank likelihood for semiparametric copula estimation. Annals of Applied Statistics. 2007;1:265–283. [Author manuscript](https://arxiv.org/abs/math/0610413).
5. MOFA developers. [Official preprocessing and interpretation guidance](https://biofam.github.io/MOFA2/faq.html).

## Data and code availability

All released observations are synthetic. Sources, results and reports are available in [the public repository](https://github.com/Vusumuzi2Mabasa/Breast-Multiomics-Copula-MOFA). No proposal documents, clinical records, private genomic data or collaborator datasets are included. Sources are released under the MIT licence. Cite this repository and the original methodological publications when reusing the benchmark. Existing-method implementation and simulation evaluation do not by themselves establish methodological novelty.
'''
    (BASE/'README.md').write_text(s,encoding='utf-8')
    print('Study 3 research manuscript generated.')


if __name__=='__main__':main()
