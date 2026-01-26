# Test script to verify results output and visualization code
# This script checks if all necessary objects exist before running visualization

cat("=== Testing Results Output Setup ===\n\n")

# Check required objects
required_objects <- c(
  "dml_split_tbl",
  "dml_full_tbl", 
  "clinical_sensitivity_tbl",
  "selected_pathway_id",
  "selected_pathway_label",
  "selected_pathway",
  "test_df",
  "train_df",
  "gsva_rank_tbl"
)

cat("Checking required objects:\n")
for (obj in required_objects) {
  exists_flag <- exists(obj)
  status <- ifelse(exists_flag, "✓", "✗")
  cat(sprintf("  %s %s\n", status, obj))
  
  if (!exists_flag) {
    cat(sprintf("    WARNING: %s not found. Make sure to run the full analysis first.\n", obj))
  }
}

cat("\n")

# Check data structure
if (exists("dml_split_tbl")) {
  cat("dml_split_tbl structure:\n")
  print(str(dml_split_tbl))
  cat("\nPreview:\n")
  print(head(dml_split_tbl))
}

if (exists("test_df")) {
  cat("\n\ntest_df columns:\n")
  print(colnames(test_df))
  cat("\nRequired columns present:\n")
  cat("  OS_MONTHS:", "OS_MONTHS" %in% colnames(test_df), "\n")
  cat("  OS_STATUS_BINARY:", "OS_STATUS_BINARY" %in% colnames(test_df), "\n")
  cat("  pathway_score:", "pathway_score" %in% colnames(test_df), "\n")
}

cat("\n=== Test Complete ===\n")
cat("If all objects are present (✓), you can run the visualization code.\n")
cat("If any objects are missing (✗), run the full analysis pipeline first.\n")
