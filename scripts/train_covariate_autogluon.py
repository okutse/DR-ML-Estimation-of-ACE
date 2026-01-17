#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoGluon Covariate Model for Tumor Stage Prediction

This script:
1. Loads and merges clinical_patient and clinical_sample data
2. Maps TUMOR_STAGE to binary (0=early [0,1,2], 1=late [3,4])
3. Performs stratified K-fold CV on TUMOR_STAGE, CANCER_TYPE_DETAILED, ER_STATUS
4. Trains AutoGluon TabularPredictor on each fold
5. Computes comprehensive metrics: AUC, ROC, sensitivity, specificity, F1, etc.
6. Extracts feature importance, SHAP values, and LIME explanations
7. Saves models and results to results/ subdirectory

Usage:
    # Quick test with 3 folds
    python scripts/train_covariate_autogluon.py --n-folds 3 --time-limit 60
    
    # Full run with 10 folds
    python scripts/train_covariate_autogluon.py --n-folds 10 --presets best_quality
"""

import argparse
import json
import logging
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import (
    accuracy_score,
    auc,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold

try:
    from autogluon.tabular import TabularPredictor
    AUTOGLUON_AVAILABLE = True
except ImportError:
    AUTOGLUON_AVAILABLE = False
    warnings.warn("AutoGluon not installed. Install with: pip install autogluon")

try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False
    warnings.warn("SHAP not installed. Install with: pip install shap")

try:
    from lime.lime_tabular import LimeTabularExplainer
    LIME_AVAILABLE = True
except ImportError:
    LIME_AVAILABLE = False
    warnings.warn("LIME not installed. Install with: pip install lime")

warnings.filterwarnings('ignore', category=FutureWarning)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s'
)
LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="AutoGluon covariate model for tumor stage prediction"
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Directory containing data files"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/covariate_autogluon"),
        help="Output directory for models and results"
    )
    parser.add_argument(
        "--n-folds",
        type=int,
        default=10,
        help="Number of folds for cross-validation (use 3 for quick test)"
    )
    parser.add_argument(
        "--presets",
        type=str,
        default="medium_quality",
        choices=["medium_quality", "good_quality", "best_quality", "optimize_for_deployment"],
        help="AutoGluon presets"
    )
    parser.add_argument(
        "--time-limit",
        type=int,
        default=3600,
        help="Time limit per fold in seconds (0=no limit)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.2,
        help="Proportion of data for final test set (default: 0.2 = 20%)"
    )
    parser.add_argument(
        "--use-nested-cv",
        action="store_true",
        help="Use nested CV: split train/test first, then K-fold on train"
    )
    parser.add_argument(
        "--max-shap-samples",
        type=int,
        default=200,
        help="Max samples for SHAP analysis"
    )
    parser.add_argument(
        "--max-lime-samples",
        type=int,
        default=5,
        help="Max samples for LIME explanations"
    )
    parser.add_argument(
        "--eval-metric",
        type=str,
        default="roc_auc",
        help="AutoGluon evaluation metric"
    )
    return parser.parse_args()


def load_and_merge_clinical_data(data_dir: Path) -> pd.DataFrame:
    """
    Load clinical_patient and clinical_sample data, then merge.
    
    Replicates the merge step from the notebook.
    """
    LOGGER.info("Loading clinical data...")
    
    patient_file = data_dir / "data_clinical_patient.txt"
    sample_file = data_dir / "data_clinical_sample.txt"
    
    if not patient_file.exists():
        raise FileNotFoundError(f"Missing file: {patient_file}")
    if not sample_file.exists():
        raise FileNotFoundError(f"Missing file: {sample_file}")
    
    # Read files (skip first 4 rows as header)
    clinical_patient = pd.read_csv(patient_file, sep='\t', header=4)
    clinical_sample = pd.read_csv(sample_file, sep='\t', header=4)
    
    LOGGER.info(f"Patient data: {clinical_patient.shape}")
    LOGGER.info(f"Sample data: {clinical_sample.shape}")
    
    # Merge on PATIENT_ID
    merged = pd.merge(
        clinical_patient,
        clinical_sample,
        on='PATIENT_ID',
        suffixes=('_patient', '_sample'),
        how='inner'
    )
    
    LOGGER.info(f"Merged data: {merged.shape}")
    return merged


def map_tumor_stage(df: pd.DataFrame) -> pd.DataFrame:
    """
    Map TUMOR_STAGE to binary: 0 (early: 0,1,2) and 1 (late: 3,4).
    """
    LOGGER.info("Mapping TUMOR_STAGE to binary...")
    
    df = df.copy()
    df['TUMOR_STAGE'] = df['TUMOR_STAGE'].replace({
        0: 0,
        1: 0,
        2: 0,
        3: 1,
        4: 1,
    })
    
    stage_counts = df['TUMOR_STAGE'].value_counts()
    LOGGER.info(f"TUMOR_STAGE distribution:\n{stage_counts}")
    
    return df


def build_stratified_labels(
    df: pd.DataFrame,
    stratify_cols: List[str],
    min_count: int = 5
) -> pd.Series:
    """
    Build combined stratification labels from multiple columns.
    
    Handles missing values and checks for minimum stratum size.
    """
    LOGGER.info(f"Building stratified labels from: {stratify_cols}")
    
    # Combine columns with | separator, filling NA with 'MISSING'
    labels = df[stratify_cols[0]].astype(str).fillna('MISSING')
    for col in stratify_cols[1:]:
        if col in df.columns:
            labels = labels + '|' + df[col].astype(str).fillna('MISSING')
    
    # Check stratum sizes
    counts = labels.value_counts()
    small_strata = counts[counts < min_count]
    
    if len(small_strata) > 0:
        LOGGER.warning(
            f"Found {len(small_strata)} strata with < {min_count} samples. "
            f"Consider using only primary stratification (TUMOR_STAGE)."
        )
        # Fallback to just target variable
        return df[stratify_cols[0]].astype(str).fillna('MISSING')
    
    return labels


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_proba: np.ndarray) -> Dict:
    """Compute comprehensive classification metrics."""
    metrics = {
        'accuracy': float(accuracy_score(y_true, y_pred)),
        'precision': float(precision_score(y_true, y_pred, zero_division=0)),
        'recall': float(recall_score(y_true, y_pred, zero_division=0)),
        'sensitivity': float(recall_score(y_true, y_pred, zero_division=0)),  # Same as recall
        'f1': float(f1_score(y_true, y_pred, zero_division=0)),
    }
    
    # Specificity (True Negative Rate)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    metrics['specificity'] = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
    
    # AUC
    try:
        metrics['roc_auc'] = float(roc_auc_score(y_true, y_proba))
    except:
        metrics['roc_auc'] = None
    
    # PR AUC
    try:
        precision, recall, _ = precision_recall_curve(y_true, y_proba)
        metrics['pr_auc'] = float(auc(recall, precision))
    except:
        metrics['pr_auc'] = None
    
    return metrics


def plot_roc_curve(y_true: np.ndarray, y_proba: np.ndarray, save_path: Path):
    """Plot and save ROC curve."""
    fpr, tpr, _ = roc_curve(y_true, y_proba)
    roc_auc = auc(fpr, tpr)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    LOGGER.info(f"Saved ROC curve to {save_path}")


def plot_confusion_matrix(y_true: np.ndarray, y_pred: np.ndarray, save_path: Path):
    """Plot and save confusion matrix."""
    cm = confusion_matrix(y_true, y_pred)
    
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    LOGGER.info(f"Saved confusion matrix to {save_path}")


def extract_feature_importance(predictor: TabularPredictor, save_path: Path, data: pd.DataFrame = None) -> pd.DataFrame:
    """Extract and save feature importance."""
    try:
        importance = predictor.feature_importance(data=data, subsample_size=1000)
        importance_df = pd.DataFrame({
            'feature': importance.index,
            'importance': importance.values
        }).sort_values('importance', ascending=False)
        
        importance_df.to_csv(save_path, index=False)
        LOGGER.info(f"Saved feature importance to {save_path}")
        
        # Plot top 20 features
        plt.figure(figsize=(10, 8))
        top_n = min(20, len(importance_df))
        plt.barh(range(top_n), importance_df['importance'].head(top_n))
        plt.yticks(range(top_n), importance_df['feature'].head(top_n))
        plt.xlabel('Importance')
        plt.title(f'Top {top_n} Feature Importance')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(save_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        return importance_df
    except Exception as e:
        LOGGER.warning(f"Could not extract feature importance: {e}")
        return pd.DataFrame()


def run_shap_analysis(
    predictor: TabularPredictor,
    X: pd.DataFrame,
    output_dir: Path,
    max_samples: int = 200
):
    """Run SHAP analysis and save results."""
    if not SHAP_AVAILABLE:
        LOGGER.warning("SHAP not available, skipping")
        return
    
    LOGGER.info("Running SHAP analysis...")
    
    # Sample data for explanation
    sample = X.sample(min(max_samples, len(X)), random_state=42) if len(X) > max_samples else X
    
    # Use shap.sample() to create a smaller background dataset (speeds up KernelExplainer)
    background_size = min(50, len(sample))  # Use 50 samples for background
    background = shap.sample(sample, background_size, random_state=42)
    
    LOGGER.info(f"Using {len(background)} background samples and {len(sample)} samples for explanation")
    
    try:
        def predict_fn(data):
            df = pd.DataFrame(data, columns=X.columns)
            return predictor.predict_proba(df, as_multiclass=False).values
        
        explainer = shap.KernelExplainer(predict_fn, background)
        shap_values = explainer.shap_values(sample, nsamples=min(100, len(sample)))
        
        # Handle 3D output (binary classification)
        if isinstance(shap_values, list):
            shap_values = shap_values[1]  # Positive class
        elif shap_values.ndim == 3:
            shap_values = shap_values[:, :, 1]
        
        # Save summary
        mean_abs_shap = np.abs(shap_values).mean(axis=0)
        shap_df = pd.DataFrame({
            'feature': X.columns,
            'mean_abs_shap': mean_abs_shap
        }).sort_values('mean_abs_shap', ascending=False)
        
        shap_df.to_csv(output_dir / 'shap_importance.csv', index=False)
        LOGGER.info(f"Saved SHAP importance to {output_dir / 'shap_importance.csv'}")
        
        # Plot summary
        plt.figure(figsize=(10, 8))
        shap.summary_plot(shap_values, sample, feature_names=X.columns, show=False)
        plt.tight_layout()
        plt.savefig(output_dir / 'shap_summary.png', dpi=300, bbox_inches='tight')
        plt.close()
        LOGGER.info(f"Saved SHAP summary plot")
        
    except Exception as e:
        LOGGER.warning(f"SHAP analysis failed: {e}")


def run_lime_analysis(
    predictor: TabularPredictor,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    output_dir: Path,
    max_samples: int = 5
):
    """Run LIME analysis and save results."""
    if not LIME_AVAILABLE:
        LOGGER.warning("LIME not available, skipping")
        return
    
    LOGGER.info("Running LIME analysis...")
    
    try:
        explainer = LimeTabularExplainer(
            X_train.values,
            feature_names=X_train.columns.tolist(),
            class_names=['Early', 'Late'],
            mode='classification'
        )
        
        def predict_fn(data):
            df = pd.DataFrame(data, columns=X_train.columns)
            proba = predictor.predict_proba(df, as_multiclass=False)
            # Return as 2D array [prob_class_0, prob_class_1]
            return np.column_stack([1 - proba.values, proba.values])
        
        lime_results = []
        n_explain = min(max_samples, len(X_test))
        
        for i in range(n_explain):
            exp = explainer.explain_instance(
                X_test.iloc[i].values,
                predict_fn,
                num_features=10
            )
            
            # Save HTML
            exp.save_to_file(output_dir / f'lime_explanation_{i}.html')
            
            # Extract feature weights
            for feat, weight in exp.as_list():
                lime_results.append({
                    'sample_idx': i,
                    'feature': feat,
                    'weight': weight
                })
        
        if lime_results:
            lime_df = pd.DataFrame(lime_results)
            lime_df.to_csv(output_dir / 'lime_results.csv', index=False)
            LOGGER.info(f"Saved LIME results to {output_dir / 'lime_results.csv'}")
        
    except Exception as e:
        LOGGER.warning(f"LIME analysis failed: {e}")


def run_nested_cv_experiment(
    df: pd.DataFrame,
    args: argparse.Namespace
) -> Dict:
    """
    Run nested CV: 
    1. Split data into train/test (80/20)
    2. K-fold CV on train set
    3. Final evaluation on hold-out test set
    """
    from sklearn.model_selection import train_test_split
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    LOGGER.info("="*80)
    LOGGER.info("NESTED CV: Train/Test Split + K-fold CV")
    LOGGER.info("="*80)
    
    # Save config
    config = {
        'n_folds': args.n_folds,
        'test_size': args.test_size,
        'presets': args.presets,
        'time_limit': args.time_limit,
        'seed': args.seed,
        'nested_cv': True,
    }
    with open(output_dir / 'config.json', 'w') as f:
        json.dump({k: str(v) for k, v in config.items()}, f, indent=2)
    
    # Prepare data
    target_col = 'TUMOR_STAGE'
    id_cols = ['PATIENT_ID', 'SAMPLE_ID']
    df_clean = df[df[target_col].notna()].copy()
    
    drop_cols = id_cols + [target_col]
    feature_cols = [c for c in df_clean.columns if c not in drop_cols]
    
    X = df_clean[feature_cols].copy()
    y = df_clean[target_col].astype(int)
    
    # Build stratification labels
    stratify_cols = [target_col, 'CANCER_TYPE_DETAILED', 'ER_STATUS']
    stratify_labels = build_stratified_labels(df_clean, stratify_cols)
    
    # 1. Split into train and final test set
    LOGGER.info(f"\n{'='*60}")
    LOGGER.info(f"Step 1: Splitting data into train ({1-args.test_size:.0%}) and test ({args.test_size:.0%})")
    LOGGER.info(f"{'='*60}")
    
    X_train_full, X_test_final, y_train_full, y_test_final, strat_train, strat_test = train_test_split(
        X, y, stratify_labels, 
        test_size=args.test_size, 
        random_state=args.seed, 
        stratify=stratify_labels
    )
    
    LOGGER.info(f"Training set: {len(X_train_full)} samples")
    LOGGER.info(f"Final test set: {len(X_test_final)} samples")
    LOGGER.info(f"Training target distribution:\n{y_train_full.value_counts()}")
    LOGGER.info(f"Test target distribution:\n{y_test_final.value_counts()}")
    
    # Save test set for later
    test_set = X_test_final.copy()
    test_set[target_col] = y_test_final
    test_set.to_csv(output_dir / 'final_test_set.csv', index=False)
    
    # 2. Train model with AutoGluon's internal validation, then refit on full training data
    LOGGER.info(f"\n{'='*60}")
    LOGGER.info(f"Step 2: Training with AutoGluon (internal validation + refit_full)")
    LOGGER.info(f"{'='*60}")
    LOGGER.info("Strategy for covariate selection:")
    LOGGER.info("  1. AutoGluon trains models with internal validation (~80/20 split)")
    LOGGER.info("  2. Finds best model ensemble and hyperparameters")
    LOGGER.info("  3. refit_full() retrains best models on 100% of training data")
    LOGGER.info("  4. Feature importance extracted from the refitted model")
    LOGGER.info("  5. Final evaluation on hold-out test set")
    
    train_data_full = X_train_full.copy()
    train_data_full[target_col] = y_train_full
    
    model_dir = output_dir / 'model'
    model_dir.mkdir(exist_ok=True)
    
    # Initial training with internal validation
    predictor = TabularPredictor(
        label=target_col,
        eval_metric=args.eval_metric,
        path=str(model_dir / 'ag_models'),
        verbosity=2
    )
    
    LOGGER.info("Training initial models with internal validation...")
    predictor.fit(
        train_data=train_data_full,
        presets=args.presets,
        time_limit=args.time_limit if args.time_limit > 0 else None,
    )
    
    # Get leaderboard before refit
    LOGGER.info(f"\n{'='*60}")
    LOGGER.info("Step 3: Model Performance (with internal validation holdout)")
    LOGGER.info(f"{'='*60}")
    
    leaderboard_before = predictor.leaderboard(silent=True)
    leaderboard_before.to_csv(output_dir / 'leaderboard_before_refit.csv', index=False)
    
    LOGGER.info(f"\n📊 Top Models (trained on ~80% with validation on ~20%):")
    LOGGER.info("="*80)
    for idx, row in leaderboard_before.head(10).iterrows():
        LOGGER.info(f"  {row['model']:30s} | score_val: {row['score_val']:.4f} | score_test: {row.get('score_test', 'N/A')}")
    LOGGER.info("="*80)
    
    best_model_before = predictor.model_best
    LOGGER.info(f"\nBest model (before refit): {best_model_before}")
    
    # Refit on full training data (no holdout)
    LOGGER.info(f"\n{'='*60}")
    LOGGER.info(f"Step 4: Refitting best models on FULL training data (refit_full)")
    LOGGER.info(f"{'='*60}")
    LOGGER.info("Using refit_full() to retrain best models on 100% of training data...")
    LOGGER.info("This ensures:")
    LOGGER.info("  - Maximum use of training data")
    LOGGER.info("  - Consistent feature importance for covariate selection")
    LOGGER.info("  - Single source of truth for feature rankings")
    
    # refit_full() trains models on full dataset
    refit_models = predictor.refit_full()
    
    LOGGER.info(f"\n✅ Successfully refitted {len(refit_models)} models")
    LOGGER.info(f"Refitted models: {refit_models}")
    
    # Get best model after refit
    best_model_after = predictor.model_best
    LOGGER.info(f"Best model (after refit): {best_model_after}")
    
    # Extract feature importance from refitted model
    LOGGER.info(f"\n{'='*60}")
    LOGGER.info("Step 5: Feature Importance (from refitted model on full training data)")
    LOGGER.info(f"{'='*60}")
    
    importance_df = extract_feature_importance(predictor, model_dir / 'feature_importance.csv', data=train_data_full)
    
    # SHAP analysis on refitted model
    if SHAP_AVAILABLE:
        shap_dir = model_dir / 'shap'
        shap_dir.mkdir(exist_ok=True)
        LOGGER.info(f"\n{'='*60}")
        LOGGER.info("Step 6: SHAP Analysis (on training data)")
        LOGGER.info(f"{'='*60}")
        run_shap_analysis(predictor, X_train_full, shap_dir, args.max_shap_samples)
        
        # Feature selection based on SHAP values
        shap_file = shap_dir / 'shap_importance.csv'
        if shap_file.exists():
            shap_importance = pd.read_csv(shap_file)
            
            LOGGER.info("\n" + "="*60)
            LOGGER.info("Feature Selection Summary")
            LOGGER.info("="*60)
            
            # Method 1: Top N features
            top_10_features = shap_importance.head(10)['feature'].tolist()
            top_20_features = shap_importance.head(20)['feature'].tolist()
            
            # Method 2: Threshold-based (SHAP > 0.01)
            threshold = 0.01
            threshold_features = shap_importance[shap_importance['mean_abs_shap'] > threshold]['feature'].tolist()
            
            # Method 3: Cumulative importance (80%)
            shap_importance['cumulative_importance'] = (
                shap_importance['mean_abs_shap'] / shap_importance['mean_abs_shap'].sum()
            ).cumsum()
            features_80pct = shap_importance[shap_importance['cumulative_importance'] <= 0.8]['feature'].tolist()
            
            # Save feature selection results
            feature_selection = {
                'method': 'single_model_refit_full',
                'top_10_features': top_10_features,
                'top_20_features': top_20_features,
                'threshold_features': {
                    'threshold': threshold,
                    'features': threshold_features,
                    'n_features': len(threshold_features)
                },
                'cumulative_80pct_features': {
                    'features': features_80pct,
                    'n_features': len(features_80pct)
                }
            }
            
            with open(output_dir / 'selected_features.json', 'w') as f:
                json.dump(feature_selection, f, indent=2)
            
            LOGGER.info(f"\n📋 Top 10 Most Important Features (by SHAP):")
            for i, feat in enumerate(top_10_features, 1):
                shap_val = shap_importance[shap_importance['feature'] == feat]['mean_abs_shap'].values[0]
                LOGGER.info(f"  {i}. {feat}: {shap_val:.4f}")
            
            LOGGER.info(f"\n📊 Feature Selection Methods:")
            LOGGER.info(f"  - Top 10: {len(top_10_features)} features")
            LOGGER.info(f"  - Top 20: {len(top_20_features)} features")
            LOGGER.info(f"  - Threshold (>{threshold}): {len(threshold_features)} features")
            LOGGER.info(f"  - Cumulative 80%: {len(features_80pct)} features")
            LOGGER.info(f"\n✅ Feature selection results saved to {output_dir / 'selected_features.json'}")
    
    # LIME analysis (optional)
    if LIME_AVAILABLE:
        lime_dir = model_dir / 'lime'
        lime_dir.mkdir(exist_ok=True)
        LOGGER.info(f"\n{'='*60}")
        LOGGER.info("Step 7: LIME Analysis (sample explanations)")
        LOGGER.info(f"{'='*60}")
        run_lime_analysis(predictor, X_train_full, X_train_full, lime_dir, args.max_lime_samples)
    
    # Final evaluation on hold-out test set
    LOGGER.info(f"\n{'='*60}")
    LOGGER.info("Step 8: Final Evaluation on Hold-out Test Set")
    LOGGER.info(f"{'='*60}")
    
    test_data_final = X_test_final.copy()
    test_data_final[target_col] = y_test_final
    
    y_pred_final = predictor.predict(test_data_final.drop(columns=[target_col]))
    y_proba_final = predictor.predict_proba(test_data_final.drop(columns=[target_col]), as_multiclass=False)
    
    final_metrics = compute_metrics(y_test_final.values, y_pred_final.values, y_proba_final.values)
    
    LOGGER.info("\n🎯 FINAL TEST SET PERFORMANCE:")
    LOGGER.info("="*80)
    for k, v in final_metrics.items():
        if v is not None:
            LOGGER.info(f"  {k}: {v:.4f}")
    LOGGER.info("="*80)
    
    # Save final results
    final_results = {
        'internal_cv_metrics': {
            'best_model': best_model_after,
            'score_val': float(leaderboard_before.iloc[0]['score_val'])
        },
        'test_metrics': final_metrics,
        'n_train': len(X_train_full),
        'n_test': len(X_test_final),
        'n_features': X.shape[1]
    }
    
    with open(output_dir / 'final_results.json', 'w') as f:
        json.dump(final_results, f, indent=2)
    
    # Final test predictions
    final_pred_df = pd.DataFrame({
        'true': y_test_final.values,
        'predicted': y_pred_final.values,
        'probability': y_proba_final.values
    })
    final_pred_df.to_csv(output_dir / 'final_test_predictions.csv', index=False)
    
    # Final plots
    plot_roc_curve(y_test_final.values, y_proba_final.values, model_dir / 'test_roc_curve.png')
    plot_confusion_matrix(y_test_final.values, y_pred_final.values, model_dir / 'test_confusion_matrix.png')
    
    # Final leaderboard on test set
    final_leaderboard = predictor.leaderboard(test_data_final, silent=True)
    final_leaderboard.to_csv(model_dir / 'test_leaderboard.csv', index=False)
    
    LOGGER.info(f"\nAll results saved to: {output_dir}")
    LOGGER.info(f"\n✅ COMPLETE: Model trained with refit_full(), features selected, test performance evaluated")
    
    return final_results


def run_cv_experiment(
    args: argparse.Namespace
) -> Dict:
    """Run K-fold cross-validation experiment."""
    
    if not AUTOGLUON_AVAILABLE:
        raise ImportError("AutoGluon is required. Install with: pip install autogluon")
    
    # Setup output directory
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save experiment config
    config = vars(args)
    with open(output_dir / 'config.json', 'w') as f:
        json.dump({k: str(v) for k, v in config.items()}, f, indent=2)
    
    # Prepare data
    target_col = 'TUMOR_STAGE'
    id_cols = ['PATIENT_ID', 'SAMPLE_ID']
    
    # Drop rows with missing target
    df_clean = df[df[target_col].notna()].copy()
    LOGGER.info(f"Samples with valid target: {len(df_clean)}")
    
    # Prepare features and target
    drop_cols = id_cols + [target_col]
    feature_cols = [c for c in df_clean.columns if c not in drop_cols]
    
    X = df_clean[feature_cols].copy()
    y = df_clean[target_col].astype(int)
    
    LOGGER.info(f"Features: {X.shape[1]}, Samples: {len(X)}")
    LOGGER.info(f"Target distribution:\n{y.value_counts()}")
    
    # Build stratification labels
    stratify_cols = [target_col, 'CANCER_TYPE_DETAILED', 'ER_STATUS']
    stratify_labels = build_stratified_labels(df_clean, stratify_cols)
    
    # Setup K-fold CV
    cv = StratifiedKFold(n_splits=args.n_folds, shuffle=True, random_state=args.seed)
    
    # Store results
    all_results = {
        'fold_metrics': [],
        'predictions': [],
        'feature_importance': [],
        'leaderboards': [],  # 儲存每個 fold 的 leaderboard
    }
    
    # Run CV
    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, stratify_labels)):
        LOGGER.info(f"\n{'='*60}")
        LOGGER.info(f"Fold {fold_idx + 1}/{args.n_folds}")
        LOGGER.info(f"{'='*60}")
        
        # Split data
        X_train, X_test = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_train, y_test = y.iloc[train_idx].copy(), y.iloc[test_idx].copy()
        
        LOGGER.info(f"Train: {len(X_train)}, Test: {len(X_test)}")
        LOGGER.info(f"Train target dist: {y_train.value_counts().to_dict()}")
        LOGGER.info(f"Test target dist: {y_test.value_counts().to_dict()}")
        
        # Prepare AutoGluon data
        train_data = X_train.copy()
        train_data[target_col] = y_train
        
        test_data = X_test.copy()
        test_data[target_col] = y_test
        
        # Setup fold directory
        fold_dir = output_dir / f'fold_{fold_idx}'
        fold_dir.mkdir(exist_ok=True)
        
        # Train AutoGluon
        LOGGER.info("Training AutoGluon...")
        predictor = TabularPredictor(
            label=target_col,
            eval_metric=args.eval_metric,
            path=str(fold_dir / 'ag_models'),
            verbosity=2
        )
        
        predictor.fit(
            train_data=train_data,
            presets=args.presets,
            time_limit=args.time_limit if args.time_limit > 0 else None,
        )
        
        # Get leaderboard (每個模型的表現)
        LOGGER.info("Getting model leaderboard...")
        leaderboard = predictor.leaderboard(test_data, silent=True)
        leaderboard.to_csv(fold_dir / 'model_leaderboard.csv', index=False)
        LOGGER.info(f"Saved leaderboard with {len(leaderboard)} models to {fold_dir / 'model_leaderboard.csv'}")
        
        # Store leaderboard for aggregation
        leaderboard['fold'] = fold_idx
        all_results['leaderboards'].append(leaderboard)
        
        # Predict
        LOGGER.info("Making predictions...")
        y_pred = predictor.predict(test_data.drop(columns=[target_col]))
        y_proba = predictor.predict_proba(test_data.drop(columns=[target_col]), as_multiclass=False)
        
        # Compute metrics
        metrics = compute_metrics(y_test.values, y_pred.values, y_proba.values)
        metrics['fold'] = fold_idx
        all_results['fold_metrics'].append(metrics)
        
        LOGGER.info(f"Fold {fold_idx + 1} Metrics:")
        for k, v in metrics.items():
            if k != 'fold' and v is not None:
                LOGGER.info(f"  {k}: {v:.4f}")
        
        # Log best model info
        best_model = leaderboard.iloc[0]['model']
        best_score = leaderboard.iloc[0]['score_test']
        LOGGER.info(f"Best model: {best_model} (score: {best_score:.4f})")
        
        # Save predictions
        pred_df = pd.DataFrame({
            'fold': fold_idx,
            'true': y_test.values,
            'predicted': y_pred.values,
            'probability': y_proba.values
        })
        all_results['predictions'].append(pred_df)
        
        # Plot ROC and confusion matrix
        plot_roc_curve(y_test.values, y_proba.values, fold_dir / 'roc_curve.png')
        plot_confusion_matrix(y_test.values, y_pred.values, fold_dir / 'confusion_matrix.png')
        
        # Feature importance (AutoGluon needs full DataFrame with target)
        importance_df = extract_feature_importance(predictor, fold_dir / 'feature_importance.csv', data=test_data)
        if not importance_df.empty:
            importance_df['fold'] = fold_idx
            all_results['feature_importance'].append(importance_df)
        
        # SHAP analysis
        if SHAP_AVAILABLE:
            shap_dir = fold_dir / 'shap'
            shap_dir.mkdir(exist_ok=True)
            LOGGER.info(f"Running SHAP analysis for fold {fold_idx}...")
            run_shap_analysis(predictor, X_test, shap_dir, args.max_shap_samples)
        
        # LIME analysis
        if LIME_AVAILABLE:
            lime_dir = fold_dir / 'lime'
            lime_dir.mkdir(exist_ok=True)
            LOGGER.info(f"Running LIME analysis for fold {fold_idx}...")
            run_lime_analysis(predictor, X_train, X_test, lime_dir, args.max_lime_samples)
    
    # Aggregate results
    LOGGER.info(f"\n{'='*60}")
    LOGGER.info("Aggregating results across folds...")
    LOGGER.info(f"{'='*60}")
    
    # Save all predictions
    all_preds = pd.concat(all_results['predictions'], ignore_index=True)
    all_preds.to_csv(output_dir / 'all_predictions.csv', index=False)
    
    # Save fold metrics
    metrics_df = pd.DataFrame(all_results['fold_metrics'])
    metrics_df.to_csv(output_dir / 'fold_metrics.csv', index=False)
    
    # Compute average metrics
    avg_metrics = metrics_df.drop(columns=['fold']).mean().to_dict()
    std_metrics = metrics_df.drop(columns=['fold']).std().to_dict()
    
    LOGGER.info("\nAverage Metrics Across Folds:")
    for k in avg_metrics:
        if avg_metrics[k] is not None:
            LOGGER.info(f"  {k}: {avg_metrics[k]:.4f} ± {std_metrics[k]:.4f}")
    
    # Save summary
    summary = {
        'avg_metrics': avg_metrics,
        'std_metrics': std_metrics,
        'n_folds': args.n_folds,
        'n_samples': len(df_clean),
        'n_features': X.shape[1]
    }
    
    with open(output_dir / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Aggregate feature importance
    if all_results['feature_importance']:
        LOGGER.info("\nAggregating feature importance across folds...")
        all_importance = pd.concat(all_results['feature_importance'], ignore_index=True)
        avg_importance_df = all_importance.groupby('feature')['importance'].agg(['mean', 'std']).reset_index()
        avg_importance_df.columns = ['feature', 'importance_mean', 'importance_std']
        avg_importance_df = avg_importance_df.sort_values('importance_mean', ascending=False)
        avg_importance_df.to_csv(output_dir / 'avg_feature_importance.csv', index=False)
        
        # Plot average importance
        plt.figure(figsize=(10, 8))
        top_n = min(20, len(avg_importance_df))
        plt.barh(range(top_n), avg_importance_df['importance_mean'].head(top_n).values)
        plt.yticks(range(top_n), avg_importance_df['feature'].head(top_n).values)
        plt.xlabel('Average Importance')
        plt.title(f'Top {top_n} Features (Average Across Folds)')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(output_dir / 'avg_feature_importance.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Aggregate SHAP importance across folds
    shap_files = list(output_dir.glob('fold_*/shap/shap_importance.csv'))
    if shap_files:
        LOGGER.info("\nAggregating SHAP importance across folds...")
        shap_dfs = [pd.read_csv(f) for f in shap_files]
        all_shap = pd.concat(shap_dfs, ignore_index=True)
        avg_shap_importance = all_shap.groupby('feature')['mean_abs_shap'].agg(['mean', 'std']).reset_index()
        avg_shap_importance.columns = ['feature', 'shap_mean', 'shap_std']
        avg_shap_importance = avg_shap_importance.sort_values('shap_mean', ascending=False)
        avg_shap_importance.to_csv(output_dir / 'avg_shap_importance.csv', index=False)
        LOGGER.info(f"Saved aggregated SHAP importance to {output_dir / 'avg_shap_importance.csv'}")
        
        # Feature selection based on SHAP values
        LOGGER.info("\n" + "="*60)
        LOGGER.info("Feature Selection Summary")
        LOGGER.info("="*60)
        
        # Method 1: Top N features
        top_10_features = avg_shap_importance.head(10)['feature'].tolist()
        top_20_features = avg_shap_importance.head(20)['feature'].tolist()
        
        # Method 2: Threshold-based (SHAP > 0.01)
        threshold = 0.01
        threshold_features = avg_shap_importance[avg_shap_importance['shap_mean'] > threshold]['feature'].tolist()
        
        # Method 3: Cumulative importance (80%)
        avg_shap_importance['cumulative_importance'] = (
            avg_shap_importance['shap_mean'] / avg_shap_importance['shap_mean'].sum()
        ).cumsum()
        features_80pct = avg_shap_importance[avg_shap_importance['cumulative_importance'] <= 0.8]['feature'].tolist()
        
        # Save feature selection results
        feature_selection = {
            'top_10_features': top_10_features,
            'top_20_features': top_20_features,
            'threshold_features': {
                'threshold': threshold,
                'features': threshold_features,
                'n_features': len(threshold_features)
            },
            'cumulative_80pct_features': {
                'features': features_80pct,
                'n_features': len(features_80pct)
            }
        }
        
        with open(output_dir / 'selected_features.json', 'w') as f:
            json.dump(feature_selection, f, indent=2)
        
        LOGGER.info(f"\nTop 10 Most Important Features (by SHAP):")
        for i, feat in enumerate(top_10_features, 1):
            shap_val = avg_shap_importance[avg_shap_importance['feature'] == feat]['shap_mean'].values[0]
            LOGGER.info(f"  {i}. {feat}: {shap_val:.4f}")
        
        LOGGER.info(f"\nFeature Selection Methods:")
        LOGGER.info(f"  - Top 10: {len(top_10_features)} features")
        LOGGER.info(f"  - Top 20: {len(top_20_features)} features")
        LOGGER.info(f"  - Threshold (>{threshold}): {len(threshold_features)} features")
        LOGGER.info(f"  - Cumulative 80%: {len(features_80pct)} features")
        LOGGER.info(f"\nFeature selection results saved to {output_dir / 'selected_features.json'}")
    
    # Aggregate model leaderboards
    if all_results['leaderboards']:
        LOGGER.info("\nAggregating model leaderboards across folds...")
        all_leaderboards = pd.concat(all_results['leaderboards'], ignore_index=True)
        
        # Average scores by model
        avg_leaderboard = all_leaderboards.groupby('model').agg({
            'score_test': ['mean', 'std'],
            'score_val': ['mean', 'std'],
            'pred_time_test': 'mean',
            'fit_time': 'mean',
        }).reset_index()
        
        # Flatten column names
        avg_leaderboard.columns = ['model', 'score_test_mean', 'score_test_std', 
                                     'score_val_mean', 'score_val_std',
                                     'pred_time_test', 'fit_time']
        
        # Sort by test score
        avg_leaderboard = avg_leaderboard.sort_values('score_test_mean', ascending=False)
        avg_leaderboard.to_csv(output_dir / 'avg_model_leaderboard.csv', index=False)
        
        LOGGER.info(f"\n📊 Average Model Leaderboard (across {args.n_folds} folds):")
        LOGGER.info("="*80)
        for idx, row in avg_leaderboard.head(10).iterrows():
            LOGGER.info(f"  {row['model']:30s} | Test: {row['score_test_mean']:.4f} ± {row['score_test_std']:.4f} | "
                       f"Val: {row['score_val_mean']:.4f} ± {row['score_val_std']:.4f}")
        LOGGER.info("="*80)
    
    LOGGER.info(f"\nAll results saved to: {output_dir}")
    
    return summary


def main():
    """Main execution function."""
    args = parse_args()
    
    LOGGER.info("Starting AutoGluon Covariate Model Training")
    LOGGER.info(f"Configuration: {vars(args)}")
    
    # Load and merge data
    df = load_and_merge_clinical_data(args.data_dir)
    
    # Map tumor stage to binary
    df = map_tumor_stage(df)
    
    # Choose CV strategy
    if args.use_nested_cv:
        LOGGER.info("\n🔄 Using refit_full() Strategy: Train/Test Split + AutoGluon Internal Validation")
        summary = run_nested_cv_experiment(df, args)
        LOGGER.info(f"\n✅ Internal Validation ROC AUC: {summary['internal_cv_metrics']['score_val']:.4f}")
        LOGGER.info(f"🎯 Final Test ROC AUC: {summary['test_metrics'].get('roc_auc', 'N/A'):.4f}")
    else:
        LOGGER.info("\n🔄 Using standard K-fold CV (no hold-out test set)")
        summary = run_cv_experiment(df, args)
        LOGGER.info(f"\n✅ Average ROC AUC: {summary['avg_metrics'].get('roc_auc', 'N/A'):.4f}")
    
    LOGGER.info("\n" + "="*60)
    LOGGER.info("Experiment Complete!")
    LOGGER.info("="*60)
    LOGGER.info(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
