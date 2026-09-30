#!/usr/bin/env python3
"""Minimum clean-clone integrity checks for the frozen Macaque E module."""

from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

MODULE = Path(__file__).resolve().parents[1]
PARAMS = json.loads((MODULE / "config/analysis_parameters.json").read_text(encoding="utf-8"))
FEATURES = PARAMS["features"]


class ArchiveIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = pd.read_csv(MODULE / "data/02_frozen_input/final19_raw_complete_n390.csv")
        cls.z = pd.read_csv(
            MODULE / "data/03_transformed/Ca_Pu_NAC_final19_transformed_Zscore_matrix_n390.csv"
        )
        cls.pca = pd.read_csv(MODULE / "data/03_transformed/Ca_Pu_NAC_final19_PCA_scores_n390.csv")
        cls.cells = pd.read_csv(
            MODULE / "data/04_frozen_classification/frozen_cell_assignments_and_tsne_n390.csv"
        )

    def test_source_manifest(self):
        manifest = json.loads(
            (MODULE / "data/00_source_manifest/SOURCE_ARCHIVE_MANIFEST.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["allen_processed_release"]["version"], "20260228")
        for key in ("local_source_archive_not_redistributed", "expression_matrix_not_redistributed"):
            entry = manifest[key]
            self.assertGreater(entry["bytes"], 0)
            self.assertRegex(entry["sha256"], r"^[A-F0-9]{64}$")

    def test_cohort_counts(self):
        qc = pd.read_csv(MODULE / "data/01_cohort_QC/dSTR_plus_vSTR_cell_missingness_QC.csv")
        stage1 = qc.loc[
            qc.cell_qc_status.eq("Keep")
            & qc.Lib_region_of_interest_label.str.upper().isin(["CA", "PU", "NAC"])
        ]
        self.assertEqual(len(stage1), PARAMS["stage1_cells"])
        self.assertEqual(len(self.raw), PARAMS["complete_case_cells"])
        self.assertEqual(self.raw.donor_label.nunique(), PARAMS["complete_case_donors"])

    def test_features_order_ids_and_missingness(self):
        self.assertEqual(len(FEATURES), 19)
        self.assertEqual(self.raw.columns[-19:].tolist(), FEATURES)
        self.assertEqual(self.z.columns[1:].tolist(), FEATURES)
        self.assertFalse(self.raw[FEATURES].isna().any().any())
        self.assertFalse(self.z[FEATURES].isna().any().any())
        expected = self.raw.cell_label.astype(str).tolist()
        self.assertEqual(self.z.cell_label.astype(str).tolist(), expected)
        self.assertEqual(self.pca.cell_label.astype(str).tolist(), expected)
        self.assertEqual(self.cells.cell_label.astype(str).tolist(), expected)

    def test_pca_and_frozen_classes(self):
        variance = pd.read_csv(MODULE / "data/03_transformed/Ca_Pu_NAC_final19_PCA_variance.csv")
        self.assertAlmostEqual(100 * variance.explained_variance_ratio.iloc[:3].sum(), 54.6211, places=4)
        consensus = self.cells.loc[self.cells.Consensus].copy()
        consensus["E_class"] = consensus.HC_class.map(PARAMS["class_label_map"])
        self.assertEqual(len(consensus), PARAMS["consensus_cells"])
        donor_map = self.raw.set_index("cell_label").donor_label
        self.assertEqual(consensus.cell_label.map(donor_map).nunique(), PARAMS["consensus_donors"])
        self.assertEqual(consensus.E_class.value_counts().sort_index().to_dict(), PARAMS["class_counts"])
        self.assertEqual(consensus.T_class.value_counts().to_dict(), PARAMS["T_counts_consensus"])

    def test_consensus_statistics_and_rrr_cells(self):
        consensus = self.cells.loc[self.cells.Consensus].copy()
        consensus["E_class"] = consensus.HC_class.map(PARAMS["class_label_map"])
        stats_labels = pd.read_csv(MODULE / "results/01_feature_statistics/04_frozen_cell_assignments.csv")
        self.assertSetEqual(set(stats_labels.MSN_unique_ID.astype(str)), set(consensus.cell_label.astype(str)))
        stable = consensus.loc[consensus.T_class.isin(["D1", "D2"])]
        rrr = pd.read_csv(MODULE / "results/03_T_E_RRR/RRR_cell_scores_n346.csv")
        self.assertEqual(len(rrr), PARAMS["RRR"]["cells"])
        self.assertSetEqual(set(rrr.cell_label.astype(str)), set(stable.cell_label.astype(str)))
        c_to_e = PARAMS["class_label_map"]
        self.assertTrue(
            rrr.set_index("cell_label").E_class.map(c_to_e).sort_index().equals(
                stable.set_index("cell_label").E_class.sort_index()
            )
        )

    def test_ml_parameters_and_saved_results(self):
        ml = json.loads((MODULE / "results/02_ML500/core/00_parameters.json").read_text(encoding="utf-8"))
        self.assertEqual(ml["grouping_unit"], "donor_label")
        self.assertEqual(ml["repeats"], 500)
        self.assertFalse(ml["primary_transform"]["yeo_johnson"])
        repeats = pd.read_csv(MODULE / "results/02_ML500/core/08_cross_algorithm_E4_donor_resampling_500.csv")
        self.assertFalse(repeats.empty)
        self.assertFalse(repeats.isna().all(axis=1).any())
        self.assertEqual(PARAMS["features"], FEATURES)  # guards against cross-module predictors

    def test_key_tables_have_no_unexpected_missing_values(self):
        for relative in [
            "results/01_feature_statistics/01_descriptive_statistics_by_E_class.csv",
            "results/01_feature_statistics/02_omnibus_KruskalWallis_BH19.csv",
            "results/02_ML500/core/05_supervised_summary_95CI.csv",
            "results/03_T_E_RRR/RRR_cell_scores_n346.csv",
        ]:
            frame = pd.read_csv(MODULE / relative)
            self.assertFalse(frame.empty, relative)
            self.assertFalse(frame.isna().all(axis=1).any(), relative)


if __name__ == "__main__":
    unittest.main(verbosity=2)
