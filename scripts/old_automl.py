#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AutoGluon pipeline to predict binary tumor stage from merged_data.csv.

Target mapping:
  - TUMOR_STAGE in {0,1,2} -> 0 (early)
  - TUMOR_STAGE in {3,4}   -> 1 (late)

Outputs are written under: results/autogluon_stage/
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from autogluon.tabular import TabularPredictor

from sklearn.metrics import (
    accuracy_score,
    auc,
    classification_report,
    confusion_matrix,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
    f1_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.calibration import calibration_curve

import matplotlib.pyplot as plt

try:
    import shap  # type: ignore
    SHAP_AVAILABLE = True
except Exception:
    SHAP_AVAILABLE = False

try:
    from lime.lime_tabular import LimeTabularExplainer  # type: ignore
    LIME_AVAILABLE = True
except Exception:
    LIME_AVAILABLE = False


LOGGER = logging.getLogger("autogluon_stage")


DEFAULT_OUTPUT_DIR = Path("results/autogluon_stage")
DEFAULT_DATA_PATH = Path("data/merged_data.csv")
DEFAULT_SAMPLE_PATH = Path("data/data_clinical_sample.txt")

ID_COLS = ["PATIENT_ID", "SAMPLE_ID"]
TARGET_COL = "mapped_stage"
TARGET_MASK_COL = "mapped_stage_mask"
RAW_STAGE_COL = "TUMOR_STAGE"

SAMPLE_CLIN_COLS = [
    "PATIENT_ID",
    "SAMPLE_ID",
    "CANCER_TYPE",
    "CANCER_TYPE_DETAILED",
    "ER_STATUS",
    "HER2_STATUS",
    "GRADE",
    "ONCOTREE_CODE",
    "PR_STATUS",
    "SAMPLE_TYPE",
    "TUMOR_SIZE",
    "TUMOR_STAGE",
    "TMB_NONSYNONYMOUS",
]

CLINICAL_MASK_CANDIDATES = [
    "ER_STATUS",
    "HER2_STATUS",
    "PR_STATUS",
    "GRADE",
    "TUMOR_SIZE",
    "TMB_NONSYNONYMOUS",
    "TUMOR_STAGE",
]


def setup_logging(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="AutoGluon tumor stage classifier")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA_PATH)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE_PATH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--num-folds", type=int, default=10)
    # parser.add_argument("--presets", type=str, default="medium_quality")
    parser.add_argument("--presets", type=str, default="good_quality")
    parser.add_argument("--time-limit", type=int, default=0)
    parser.add_argument(
        "--secondary-stratify",
        type=str,
        nargs="*",
        default=["CANCER_TYPE_DETAILED", "ER_STATUS"],
    )
    parser.add_argument(
        "--threshold-metric",
        type=str,
        default="f1",
        choices=["f1", "roc", "pr", "accuracy"],
        help="Metric to maximize when choosing decision threshold on validation set.",
    )
    parser.add_argument("--max-shap-samples", type=int, default=200)
    parser.add_argument("--max-lime-samples", type=int, default=5)
    parser.add_argument("--permutation-subsample", type=int, default=2000)
    parser.add_argument(
        "--max-memory-ratio",
        type=float,
        default=1.0,
        help="AutoGluon ag.max_memory_usage_ratio; >1.0 may risk OOM.",
    )
    parser.add_argument(
        "--disable-ray",
        action="store_true",
        help="Disable Ray to avoid dashboard errors on Windows.",
    )
    parser.add_argument(
        "--num-bag-folds",
        type=int,
        default=0,
        help="Bagging folds inside AutoGluon (set 0 to disable).",
    )
    parser.add_argument(
        "--num-stack-levels",
        type=int,
        default=0,
        help="Stacking levels inside AutoGluon (set 0 to disable).",
    )
    parser.add_argument(
        "--stratify-min-count",
        type=int,
        default=5,
        help="Minimum samples per stratum for secondary stratification.",
    )
    return parser.parse_args()


def read_sample_clinical(sample_path: Path) -> pd.DataFrame:
    """Read clinical sample data from tab-separated file, skipping header rows."""
    if not sample_path.exists():
        raise FileNotFoundError(f"Missing sample file: {sample_path}")
    df = pd.read_csv(sample_path, sep="\t", skiprows=[0, 1, 2, 3])
    return df


def merge_sample_clinical(merged_df: pd.DataFrame, sample_df: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [c for c in SAMPLE_CLIN_COLS if c in sample_df.columns]
    sample_df = sample_df.loc[:, sample_cols].copy()

    merged = merged_df.merge(sample_df, on="PATIENT_ID", how="left", suffixes=("", "_sample"))
    # Coalesce overlapping columns
    for col in sample_cols:
        sample_col = f"{col}_sample"
        if sample_col in merged.columns:
            merged[col] = merged[col].combine_first(merged[sample_col])
            merged.drop(columns=[sample_col], inplace=True)
    return merged


def map_stage(series: pd.Series) -> Tuple[pd.Series, pd.Series]:
    """Map tumor stage to binary (early=0, late=1) and return validity mask.
    
    Early: stages 0, 1, 2 -> 0
    Late: stages 3, 4 -> 1
    """
    stage_num = pd.to_numeric(series, errors="coerce")
    mapped = stage_num.map({0: 0, 1: 0, 2: 0, 3: 1, 4: 1}).astype("float")
    mask = mapped.notna().astype(int)
    return mapped.astype("Int64"), mask


def build_missingness_masks(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[f"{col}_mask"] = df[col].isna().astype(int)
    return df


def data_profile(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "column": df.columns,
            "dtype": [str(t) for t in df.dtypes],
            "missing": df.isna().sum().values,
            "missing_pct": (df.isna().mean().values * 100).round(3),
            "unique": [df[c].nunique(dropna=True) for c in df.columns],
        }
    ).sort_values("missing_pct", ascending=False)


def get_feature_columns(df: pd.DataFrame) -> List[str]:
    drop_cols = set(ID_COLS + [RAW_STAGE_COL, TARGET_COL, TARGET_MASK_COL])
    feature_cols = [c for c in df.columns if c not in drop_cols]
    # Drop all-missing columns
    feature_cols = [c for c in feature_cols if not df[c].isna().all()]
    return feature_cols


def build_stratify_labels(
    df: pd.DataFrame,
    y: pd.Series,
    secondary_cols: Optional[List[str]],
    min_count: int = 5,
) -> pd.Series:
    if not secondary_cols:
        return y.astype(str)

    cols = [c for c in secondary_cols if c in df.columns]
    if not cols:
        return y.astype(str)

    labels = y.astype(str).copy()
    for col in cols:
        labels = labels + "|" + df[col].astype(str).fillna("NA")

    counts = labels.value_counts(dropna=False)
    if (counts < min_count).any():
        LOGGER.warning("Secondary stratification too sparse; fallback to target-only stratification.")
        return y.astype(str)
    return labels


def save_predictions(
    df: pd.DataFrame,
    y_true: pd.Series,
    y_pred: np.ndarray,
    y_proba: np.ndarray,
    out_path: Path,
) -> None:
    out_df = df.loc[:, [c for c in ID_COLS if c in df.columns]].copy()
    out_df["true_label"] = y_true.values
    out_df["predicted_label"] = y_pred
    out_df["predicted_probability"] = y_proba
    out_df.to_csv(out_path, index=False)


def safe_auc(y_true: np.ndarray, y_proba: np.ndarray) -> float:
    try:
        return float(roc_auc_score(y_true, y_proba))
    except Exception:
        return float("nan")


def plot_confusion(y_true: np.ndarray, y_pred: np.ndarray, out_path: Path) -> None:
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(6, 5))
    plt.imshow(cm, cmap="Blues")
    plt.title("Confusion Matrix")
    plt.xlabel("Predicted")
    plt.ylabel("True")
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, cm[i, j], ha="center", va="center", color="black")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_roc(y_true: np.ndarray, y_proba: np.ndarray, out_path: Path) -> None:
    if len(np.unique(y_true)) < 2:
        return
    fpr, tpr, _ = roc_curve(y_true, y_proba)
    rocA = auc(fpr, tpr)
    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, label=f"AUC={rocA:.3f}")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.title("ROC Curve")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_pr(y_true: np.ndarray, y_proba: np.ndarray, out_path: Path) -> None:
    if len(np.unique(y_true)) < 2:
        return
    prec, rec, _ = precision_recall_curve(y_true, y_proba)
    prA = auc(rec, prec) if len(rec) > 1 else float("nan")
    plt.figure(figsize=(6, 6))
    plt.plot(rec, prec, label=f"AUC={prA:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision-Recall Curve")
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_calibration(y_true: np.ndarray, y_proba: np.ndarray, out_path: Path) -> None:
    """Plot calibration curve comparing predicted probabilities to actual outcomes."""
    if len(np.unique(y_true)) < 2:
        return
    try:
        frac_pos, mean_pred = calibration_curve(y_true, y_proba, n_bins=10)
    except ValueError as e:
        LOGGER.warning("Calibration curve failed: %s", e)
        return
    plt.figure(figsize=(6, 6))
    plt.plot([0, 1], [0, 1], "k--", label="Perfect")
    plt.plot(mean_pred, frac_pos, "s-", label="Model")
    plt.xlabel("Mean Predicted Probability")
    plt.ylabel("Fraction of Positives")
    plt.title("Calibration Curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_proba: np.ndarray) -> dict:
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "roc_auc": safe_auc(y_true, y_proba),
    }


def find_best_threshold(
    y_true: np.ndarray,
    y_proba: np.ndarray,
    metric: str = "f1",
) -> float:
    if len(np.unique(y_true)) < 2:
        return 0.5

    thresholds = np.linspace(0.0, 1.0, 101)
    best_t = 0.5
    best_score = -np.inf

    if metric == "roc":
        fpr, tpr, thr = roc_curve(y_true, y_proba)
        if len(thr) == 0:
            return 0.5
        scores = tpr - fpr
        return float(thr[np.argmax(scores)])

    if metric == "pr":
        prec, rec, thr = precision_recall_curve(y_true, y_proba)
        if len(thr) == 0:
            return 0.5
        f1 = 2 * (prec[:-1] * rec[:-1]) / (prec[:-1] + rec[:-1] + 1e-12)
        return float(thr[np.argmax(f1)])

    for t in thresholds:
        y_pred = (y_proba >= t).astype(int)
        if metric == "accuracy":
            score = accuracy_score(y_true, y_pred)
        else:
            score = f1_score(y_true, y_pred, zero_division=0)
        if score > best_score:
            best_score = score
            best_t = float(t)

    return best_t


def run_shap(
    predictor: TabularPredictor,
    X: pd.DataFrame,
    out_dir: Path,
    max_samples: int,
) -> Optional[pd.DataFrame]:
    if not SHAP_AVAILABLE:
        LOGGER.warning("SHAP not available; skipping.")
        return None

    sample = X.sample(min(max_samples, len(X)), random_state=42) if len(X) else X

    def predict_proba(data: np.ndarray) -> np.ndarray:
        df = pd.DataFrame(data, columns=X.columns)
        proba = predictor.predict_proba(df)
        if isinstance(proba, pd.DataFrame):
            return proba.iloc[:, 1].to_numpy()
        return proba

    explainer = shap.KernelExplainer(predict_proba, sample)
    shap_values = explainer.shap_values(sample, nsamples=min(200, len(sample)))
    shap_values = np.array(shap_values)

    # Handle 3D array from binary classification (class x samples x features)
    if shap_values.ndim == 3:
        # Use positive class (index 1) SHAP values
        shap_values = shap_values[1] if shap_values.shape[0] == 2 else shap_values[0]

    mean_abs = np.abs(shap_values).mean(axis=0)
    shap_df = pd.DataFrame({"feature": X.columns, "mean_abs_shap": mean_abs})
    shap_df.sort_values("mean_abs_shap", ascending=False).to_csv(out_dir / "shap_importance.csv", index=False)

    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, sample, feature_names=X.columns, show=False)
    plt.tight_layout()
    plt.savefig(out_dir / "shap_summary.png", dpi=300, bbox_inches="tight")
    plt.close()

    return shap_df


def run_lime(
    predictor: TabularPredictor,
    X_train: pd.DataFrame,
    X_val: pd.DataFrame,
    out_dir: Path,
    max_samples: int,
) -> Optional[pd.DataFrame]:
    if not LIME_AVAILABLE:
        LOGGER.warning("LIME not available; skipping.")
        return None

    explainer = LimeTabularExplainer(
        training_data=X_train.values,
        feature_names=X_train.columns.tolist(),
        class_names=["early", "late"],
        mode="classification",
    )

    def predict_fn(data: np.ndarray) -> np.ndarray:
        df = pd.DataFrame(data, columns=X_train.columns)
        proba = predictor.predict_proba(df)
        if isinstance(proba, pd.DataFrame):
            return proba.values
        return np.column_stack([1 - proba, proba])

    n_explain = min(max_samples, len(X_val))
    rows = []
    lime_dir = out_dir / "lime_explanations"
    lime_dir.mkdir(parents=True, exist_ok=True)

    for i in range(n_explain):
        exp = explainer.explain_instance(X_val.values[i], predict_fn, num_features=20)
        fig = exp.as_pyplot_figure()
        fig.savefig(lime_dir / f"lime_explanation_{i}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        for feature, weight in exp.as_list():
            rows.append({"instance": i, "feature": feature, "weight": weight})

    if rows:
        lime_df = pd.DataFrame(rows)
        lime_df.to_csv(out_dir / "lime_feature_weights.csv", index=False)
        return lime_df
    return None


def main() -> None:
    args = parse_args()
    setup_logging(args.output)
    LOGGER.info("Starting AutoGluon tumor stage pipeline")

    if not args.data.exists():
        raise FileNotFoundError(f"Missing merged data: {args.data}")

    df = pd.read_csv(args.data)
    LOGGER.info("Loaded merged data: %s rows, %s columns", df.shape[0], df.shape[1])

    sample_df = read_sample_clinical(args.sample)
    df = merge_sample_clinical(df, sample_df)

    mapped, mask = map_stage(df.get(RAW_STAGE_COL, pd.Series(index=df.index, dtype=float)))
    df[TARGET_COL] = mapped
    df[TARGET_MASK_COL] = mask

    df = build_missingness_masks(df, CLINICAL_MASK_CANDIDATES)

    profile_df = data_profile(df)
    profile_df.to_csv(args.output / "data_profile.csv", index=False)

    label_map = df[[RAW_STAGE_COL, TARGET_COL]].copy()
    label_map.to_csv(args.output / "label_mapping.csv", index=False)

    run_cfg = {
        "data": str(args.data),
        "sample": str(args.sample),
        "output": str(args.output),
        "seed": args.seed,
        "num_folds": args.num_folds,
        "presets": args.presets,
        "time_limit": args.time_limit,
        "secondary_stratify": args.secondary_stratify,
        "threshold_metric": args.threshold_metric,
        "max_memory_ratio": args.max_memory_ratio,
        "disable_ray": args.disable_ray,
        "num_bag_folds": args.num_bag_folds,
        "num_stack_levels": args.num_stack_levels,
        "stratify_min_count": args.stratify_min_count,
    }
    (args.output / "run_config.json").write_text(json.dumps(run_cfg, indent=2))

    train_df = df[df[TARGET_MASK_COL] == 1].copy()
    if train_df.empty:
        raise ValueError("No rows with observed tumor stage for training.")

    feature_cols = get_feature_columns(train_df)
    LOGGER.info("Using %s feature columns", len(feature_cols))

    y = train_df[TARGET_COL].astype(int)
    strat_labels = build_stratify_labels(
        train_df, y, args.secondary_stratify, min_count=args.stratify_min_count
    )

    skf = StratifiedKFold(n_splits=args.num_folds, shuffle=True, random_state=args.seed)

    cv_rows = []
    perm_importance_rows = []
    shap_rows = []
    lime_rows = []

    for fold, (train_idx, val_idx) in enumerate(skf.split(train_df, strat_labels), start=1):
        fold_dir = args.output / f"fold_{fold}"
        fold_dir.mkdir(parents=True, exist_ok=True)

        train_split = train_df.iloc[train_idx].copy()
        val_split = train_df.iloc[val_idx].copy()

        train_data = train_split[feature_cols + [TARGET_COL]]
        val_data = val_split[feature_cols + [TARGET_COL]]

        predictor_path = fold_dir / "autogluon_models"
        predictor = TabularPredictor(
            label=TARGET_COL,
            path=str(predictor_path),
            problem_type="binary",
            eval_metric="roc_auc",
        )

        fit_kwargs = {
            "train_data": train_data,
            "presets": args.presets,
            "time_limit": args.time_limit if args.time_limit > 0 else None,
            "num_bag_folds": args.num_bag_folds,
            "num_stack_levels": args.num_stack_levels,
        }
        fit_kwargs = {k: v for k, v in fit_kwargs.items() if v is not None}
        ag_args_fit = {}
        if args.max_memory_ratio != 1.0:
            ag_args_fit["ag.max_memory_usage_ratio"] = args.max_memory_ratio
        if args.disable_ray:
            ag_args_fit["enable_ray"] = False
        if ag_args_fit:
            fit_kwargs["ag_args_fit"] = ag_args_fit
        predictor.fit(**fit_kwargs)

        leaderboard = predictor.leaderboard(val_data, silent=True)
        leaderboard.to_csv(fold_dir / "leaderboard.csv", index=False)

        pred_val_proba = predictor.predict_proba(val_split[feature_cols])
        if isinstance(pred_val_proba, pd.DataFrame):
            pred_val_proba = pred_val_proba.iloc[:, 1].values

        pred_train_proba = predictor.predict_proba(train_split[feature_cols])
        if isinstance(pred_train_proba, pd.DataFrame):
            pred_train_proba = pred_train_proba.iloc[:, 1].values

        best_threshold = find_best_threshold(
            val_split[TARGET_COL].values,
            pred_val_proba,
            metric=args.threshold_metric,
        )

        pred_val_labels = (pred_val_proba >= best_threshold).astype(int)
        pred_train_labels = (pred_train_proba >= best_threshold).astype(int)

        save_predictions(
            train_split,
            train_split[TARGET_COL],
            pred_train_labels,
            pred_train_proba,
            fold_dir / "train_predictions.csv",
        )
        save_predictions(
            val_split,
            val_split[TARGET_COL],
            pred_val_labels,
            pred_val_proba,
            fold_dir / "val_predictions.csv",
        )

        metrics = compute_metrics(val_split[TARGET_COL].values, pred_val_labels, pred_val_proba)
        cv_rows.append({"fold": fold, "best_threshold": best_threshold, **metrics})

        # Reports
        with open(fold_dir / "classification_report.txt", "w", encoding="utf-8") as f:
            f.write(classification_report(val_split[TARGET_COL], pred_val_labels, zero_division=0))

        # Plots
        plot_confusion(val_split[TARGET_COL].values, pred_val_labels, fold_dir / "confusion_matrix.png")
        plot_roc(val_split[TARGET_COL].values, pred_val_proba, fold_dir / "roc_curve.png")
        plot_pr(val_split[TARGET_COL].values, pred_val_proba, fold_dir / "pr_curve.png")
        plot_calibration(val_split[TARGET_COL].values, pred_val_proba, fold_dir / "calibration_curve.png")

        # Permutation importance
        try:
            fi = predictor.feature_importance(
                data=val_data,
                num_shuffle_sets=5,
                subsample_size=min(args.permutation_subsample, len(val_data)),
                silent=True,
            )
            fi = fi.reset_index().rename(columns={"index": "feature"})
            fi.to_csv(fold_dir / "permutation_importance.csv", index=False)
            fi["fold"] = fold
            perm_importance_rows.append(fi)
        except Exception as exc:
            LOGGER.warning("Permutation importance failed on fold %s: %s", fold, exc)

        # SHAP + LIME
        shap_df = run_shap(predictor, val_split[feature_cols], fold_dir, args.max_shap_samples)
        if shap_df is not None:
            shap_df["fold"] = fold
            shap_rows.append(shap_df)

        lime_df = run_lime(predictor, train_split[feature_cols], val_split[feature_cols], fold_dir, args.max_lime_samples)
        if lime_df is not None:
            lime_df["fold"] = fold
            lime_rows.append(lime_df)

    cv_summary = pd.DataFrame(cv_rows)
    cv_summary.to_csv(args.output / "cv_summary.csv", index=False)

    overall = {
        "mean": cv_summary.mean(numeric_only=True).to_dict(),
        "std": cv_summary.std(numeric_only=True).to_dict(),
        "n_folds": args.num_folds,
    }
    (args.output / "overall_metrics.json").write_text(json.dumps(overall, indent=2))

    feat_dir = args.output / "feature_importance"
    feat_dir.mkdir(parents=True, exist_ok=True)

    if perm_importance_rows:
        pd.concat(perm_importance_rows).to_csv(feat_dir / "permutation_importance_all_folds.csv", index=False)

    if shap_rows:
        pd.concat(shap_rows).to_csv(feat_dir / "shap_importance_all_folds.csv", index=False)

    if lime_rows:
        pd.concat(lime_rows).to_csv(feat_dir / "lime_top_features_all_folds.csv", index=False)

    LOGGER.info("Done.")


if __name__ == "__main__":
    main()
