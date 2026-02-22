# Limma Differential Expression

**ID**: `limma`  
**Version**: 1.0.0  
**Category**: analysis  
**Author**: CauldronGO Team

## Description

Linear models for differential expression analysis using limma

## Runtime

- **Environments**: `r`

- **Script**: `limma.R`

## Inputs

| Name | Label | Type | Required | Default | Visibility |
|------|-------|------|----------|---------|------------|
| `input_file` | Input Data File | file | Yes | - | Always visible |
| `annotation_file` | Annotation File | file | Yes | - | Always visible |
| `index_col` | Index Column | text | No | - | Always visible |
| `log2` | Apply Log2 Transformation | boolean | No | false | Always visible |
| `comparisons` | Comparisons | file | No | - | Always visible |
| `impute_order` | Imputation Order | select (before, after) | No | before | Always visible |
| `impute` | Imputation Method | select (none, knn, MinDet, MinProb, min, zero, mixed, nbavg, with, QRILC, MLE, bpca) | No | none | Always visible |
| `normalize` | Normalization Method | select (none, quantiles, quantiles.robust, vsn, center.median, center.mean) | No | none | Always visible |
| `aggregate_column` | Aggregate By Column | text | No | - | Always visible |
| `aggregate_method` | Aggregation Method | text | No | MsCoreUtils::robustSummary | Always visible |

### Input Details

#### Input Data File (`input_file`)

Proteomics or expression data file


#### Annotation File (`annotation_file`)

Sample annotation file with conditions


#### Index Column (`index_col`)

Column name to use as feature identifier


#### Apply Log2 Transformation (`log2`)

Apply log2 transformation before analysis


#### Comparisons (`comparisons`)

Comparison groups for differential analysis


#### Imputation Order (`impute_order`)

When to perform imputation relative to normalization (before or after normalization)

- **Options**: `before`, `after`

#### Imputation Method (`impute`)

Imputation method for missing values (select 'none' to skip imputation)

- **Options**: `none`, `knn`, `MinDet`, `MinProb`, `min`, `zero`, `mixed`, `nbavg`, `with`, `QRILC`, `MLE`, `bpca`

#### Normalization Method (`normalize`)

Normalization method to apply (select 'none' to skip normalization)

- **Options**: `none`, `quantiles`, `quantiles.robust`, `vsn`, `center.median`, `center.mean`

#### Aggregate By Column (`aggregate_column`)

Column name to aggregate features (e.g., aggregate peptides to proteins). Leave empty to skip aggregation.


#### Aggregation Method (`aggregate_method`)

Method to use for feature aggregation (e.g., MsCoreUtils::robustSummary, median, mean)


## Outputs

| Name | File | Type | Format | Description |
|------|------|------|--------|-------------|
| `differential_results` | `differential_analysis.txt` | data | tsv | Differential expression analysis results |
| `contrast_matrix_info` | `contrast_matrix_info.txt` | data | tsv | Contrast matrix information showing comparison directions |

## Requirements

- **R Version**: >=4.0

### R Dependencies (External File)

Dependencies are defined in: `r-packages.txt`

- `QFeatures`
- `limma`
- `MsCoreUtils`

> **Note**: When you create a custom environment for this plugin, these dependencies will be automatically installed.

## Example Data

This plugin includes example data for testing:

```yaml
  input_file: diann/imputed.data.txt
  annotation_file: differential_analysis/annotation.txt
  comparisons: differential_analysis/comparison.bca.txt
  index_col: Protein.Ids
  log2: true
```

Load example data by clicking the **Load Example** button in the UI.

## Usage

### Via UI

1. Navigate to **analysis** → **Limma Differential Expression**
2. Fill in the required inputs
3. Click **Run Analysis**

### Via Plugin System

```typescript
const jobId = await pluginService.executePlugin('limma', {
  // Add parameters here
});
```
