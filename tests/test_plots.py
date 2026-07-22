import unittest

import numpy as np
import pandas as pd

from src.plots import (
    plot_compositional_bars,
    plot_deg_volcano,
    plot_lineage_effect_bars,
    plot_probability_shift_radar,
    plot_timecourse_expression,
)


class PlotTests(unittest.TestCase):
    def test_timecourse_uses_baseline_label_and_real_scale(self):
        data = pd.DataFrame(
            {
                "gene": ["TF1", "TF1", "TF1"],
                "day": ["2", "iPSC", "0"],
                "mean_expression": [1.2, 0.4, 0.8],
                "pct_detected": [50.0, 20.0, 35.0],
                "expression_scale": ["log1p_cp10k"] * 3,
            }
        )
        figure = plot_timecourse_expression(data, "TF1")
        self.assertEqual(list(figure.data[0].x), ["Baseline", "Day 0", "Day 2"])
        self.assertEqual(figure.layout.yaxis.title.text, "Mean log1p(CP10K)")

    def test_lineage_bars_distinguish_significance(self):
        data = pd.DataFrame(
            {
                "perturbation": ["TF1_P1", "TF1_P1"],
                "lineage": ["Amnion", "Epiblast"],
                "observed_glass_delta": [1.0, -0.5],
                "q_value": [0.01, 0.2],
                "n_cells": [100, 100],
            }
        )
        figure = plot_lineage_effect_bars(data, "TF1_P1")
        self.assertEqual(sorted(figure.data[0].marker.opacity), [0.38, 1.0])

    def test_current_compositional_schema_renders(self):
        data = pd.DataFrame(
            {
                "perturbation": ["TF1_P1"],
                "driving_lineage": ["Amnion"],
                "prob_effect_size": [0.3],
                "q_value": [0.01],
                "n_cells": [100],
            }
        )
        figure = plot_compositional_bars(data, "TF1")
        self.assertEqual(len(figure.data), 1)
        self.assertEqual(list(figure.data[0].y), [0.3])

    def test_probability_radar_uses_paper_log2fc_metric(self):
        data = pd.DataFrame(
            {
                "perturbation": ["TF1_P1", "TF1_P1"],
                "lineage": ["Amnion", "Epiblast"],
                "mean_prob_pert": [0.2, 0.3],
                "mean_prob_ntc": [0.1, 0.4],
                "log2fc": [1.0, -0.4],
                "ci_low": [0.5, -0.8],
                "ci_high": [1.5, 0.1],
                "q_value_boot": [0.01, 0.2],
                "n_cells": [100, 100],
            }
        )
        figure = plot_probability_shift_radar(data, "TF1_P1")
        effect_trace = figure.data[-1]
        # The minimum manuscript-style range is -1.5..+1.5. Plotly receives
        # shifted coordinates, preserving signed effects on their lineage axes.
        self.assertEqual(list(effect_trace.r), [2.5, 1.1, 2.5])
        self.assertEqual(
            [row[0] for row in effect_trace.customdata],
            [1.0, -0.4, 1.0],
        )
        self.assertEqual(
            [row[7] for row in effect_trace.customdata],
            ["Amnion", "Epiblast", "Amnion"],
        )
        self.assertEqual(list(effect_trace.theta), [0.0, 180.0, 0.0])
        self.assertTrue(
            all(
                isinstance(theta, (int, float, np.integer, np.floating))
                for trace in figure.data
                for theta in trace.theta
            )
        )
        self.assertEqual(
            list(figure.layout.polar.angularaxis.ticktext),
            ["Amnion", "Epiblast"],
        )
        self.assertEqual(figure.layout.polar.radialaxis.title.text, "NTC-centered log2FC")
        self.assertNotIn("%", figure.layout.polar.radialaxis.title.text)

    def test_deg_volcano_uses_sensible_hover_precision(self):
        data = pd.DataFrame(
            {
                "gene": ["GENE1", "GENE2"],
                "baseMean": [123.4567, 2.3456],
                "log2FoldChange": [1.23456, -0.98765],
                "padj": [1.2345e-7, 0.123456],
            }
        )
        figure = plot_deg_volcano(data, "TF1_P1")
        hover_templates = "\n".join(trace.hovertemplate for trace in figure.data)
        self.assertIn("%{x:.2f}", hover_templates)
        self.assertIn("%{y:.2f}", hover_templates)
        self.assertIn(":.1f", hover_templates)
        self.assertIn(":.3g", hover_templates)
        self.assertEqual(figure.layout.xaxis.tickformat, ".2f")
        self.assertEqual(figure.layout.yaxis.tickformat, ".1f")


if __name__ == "__main__":
    unittest.main()
