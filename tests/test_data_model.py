import unittest

import pandas as pd

from src.data_model import (
    build_gene_results,
    filter_lineage_hits,
    lineage_label,
    perturbation_profile,
    timepoint_label,
)


class DataModelTests(unittest.TestCase):
    def setUp(self):
        self.lineage = pd.DataFrame(
            {
                "gene": ["TF.1", "TF.1", "TF2"],
                "perturbation": ["TF1_P1", "TF1_P1", "TF2_P1"],
                "lineage": ["Amnion", "Neural_Ectoderm", "Epiblast"],
                "observed_glass_delta": [0.8, -1.2, 0.2],
                "q_value": [0.01, 0.03, 0.001],
                "n_cells": [100, 100, 80],
            }
        )

    def test_literal_search_and_effect_filters(self):
        result = filter_lineage_hits(
            self.lineage,
            search_term="TF.1",
            min_abs_delta=1.0,
            significant_only=True,
        )
        self.assertEqual(result["lineage"].tolist(), ["Neural_Ectoderm"])

    def test_gene_results_are_one_row_per_gene(self):
        knockdown = pd.DataFrame(
            {"target_gene": ["TF.1"], "knockdown_pct": [72.5]}
        )
        result = build_gene_results(self.lineage, knockdown)
        tf1 = result[result["Gene"] == "TF.1"].iloc[0]
        self.assertEqual(tf1["Strongest lineage"], "Neural ectoderm")
        self.assertEqual(tf1["Effect"], "Decreased")
        self.assertAlmostEqual(tf1["Max |Δ|"], 1.2)
        self.assertAlmostEqual(tf1["Knockdown"], 72.5)

    def test_profile_is_deterministic_and_evidence_based(self):
        result = perturbation_profile(self.lineage, "TF1_P1")
        self.assertEqual(result["n_significant_lineages"], 2)
        self.assertEqual(result["strongest_lineage"], "Neural ectoderm")
        self.assertEqual(result["significant_effects"][0]["lineage"], "Neural ectoderm")

    def test_display_labels(self):
        self.assertEqual(lineage_label("Trophoblast_Like"), "Trophoblast-like")
        self.assertEqual(timepoint_label("iPSC"), "Baseline")
        self.assertEqual(timepoint_label(4), "Day 4")


if __name__ == "__main__":
    unittest.main()
