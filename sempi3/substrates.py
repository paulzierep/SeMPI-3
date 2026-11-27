"""Predicting the substrate of an adenylation domain from its sequence.

antiSMASH leaves 1,191 of the A domains in MIBiG without a call, and its
calls are not always right. Both cost SeMPI-3 whole residues of the backbone.

The model here is deliberately simple: character k-mers of the domain
sequence, TF-IDF weighted, and a nearest-neighbour vote under cosine
similarity. Substrate specificity is carried by a handful of binding-pocket
residues, so sequence similarity transfers it well, and a nearest-neighbour
model makes that assumption explicit instead of hiding it in a fitted model.

**Evaluation must group by cluster.** Homologous clusters share near-identical
domains, so a random split would train and test on the same sequence and
report an accuracy that cannot be reproduced on a new genome.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
)

import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.neighbors import KNeighborsClassifier

from . import monomers

#: k-mer length over the amino acid sequence.
KMER = 3

#: Neighbours consulted for a prediction.
N_NEIGHBOURS = 1


def normalise(name: Optional[str]) -> Optional[str]:
    """Canonical substrate name, so MIBiG and antiSMASH can be compared.

    ``Leu``, ``leucine`` and ``D-leucine`` all become ``leucine``; the
    stereochemistry is applied separately by the assembly engine.
    """
    if not name:
        return None
    key = name.strip().lower()
    if key in ("x", "pk", "nrp", "?"):
        return None
    key = monomers.ALIASES.get(key, key)
    for prefix in ("d-", "l-", "(2s)-", "(2r)-", "r-", "s-"):
        if key.startswith(prefix):
            stripped = key[len(prefix) :]
            key = monomers.ALIASES.get(stripped, stripped)
            break
    return key


@dataclass
class EvaluationResult:
    """Accuracy of one predictor over the held-out folds."""

    name: str
    n: int
    correct: int
    abstained: int

    @property
    def accuracy(self) -> float:
        """Share of all domains predicted correctly; abstaining counts as wrong."""
        return self.correct / self.n if self.n else 0.0

    @property
    def accuracy_when_called(self) -> float:
        """Accuracy over the domains where the predictor committed to a call."""
        called = self.n - self.abstained
        return self.correct / called if called else 0.0

    def to_dict(self) -> Dict[str, object]:
        return {
            "name": self.name,
            "n": self.n,
            "correct": self.correct,
            "abstained": self.abstained,
            "accuracy": round(self.accuracy, 4),
            "accuracy_when_called": round(self.accuracy_when_called, 4),
        }


class SubstratePredictor:
    """Nearest-neighbour substrate predictor over A-domain sequences."""

    def __init__(self, kmer: int = KMER, n_neighbours: int = N_NEIGHBOURS) -> None:
        self.vectoriser = TfidfVectorizer(
            analyzer="char",
            ngram_range=(kmer, kmer),
            lowercase=False,
            min_df=2,
        )
        self.model = KNeighborsClassifier(n_neighbors=n_neighbours, metric="cosine", weights="distance")
        self.fitted = False

    def fit(self, sequences: Sequence[str], labels: Sequence[str]) -> "SubstratePredictor":
        matrix = self.vectoriser.fit_transform(sequences)
        self.model.fit(matrix, list(labels))
        self.fitted = True
        return self

    def predict(self, sequences: Sequence[str]) -> List[str]:
        if not self.fitted:
            raise RuntimeError("the predictor must be fitted first")
        return list(self.model.predict(self.vectoriser.transform(sequences)))

    def predict_with_confidence(self, sequences: Sequence[str]) -> List[Tuple[str, float]]:
        """Predictions paired with the similarity to the nearest training domain."""
        if not self.fitted:
            raise RuntimeError("the predictor must be fitted first")
        matrix = self.vectoriser.transform(sequences)
        distances, _indices = self.model.kneighbors(matrix)
        predictions = self.model.predict(matrix)
        return [(label, float(1.0 - dist.min())) for label, dist in zip(predictions, distances)]


def cross_validate(
    records: List[dict],
    folds: int = 5,
    kmer: int = KMER,
    n_neighbours: int = N_NEIGHBOURS,
    seed: int = 0,
) -> Dict[str, EvaluationResult]:
    """Compare the model with antiSMASH over cluster-grouped folds.

    Clusters, not domains, are split, so no cluster contributes to both the
    training and the test side of a fold.
    """
    usable = [
        r
        for r in records
        if r.get("sequence") and normalise(r.get("substrate")) and len(r["sequence"]) > 50
    ]
    labels = np.array([normalise(r["substrate"]) for r in usable])
    sequences = np.array([r["sequence"] for r in usable])
    groups = np.array([r["bgc_id"] for r in usable])

    unique_groups = np.array(sorted(set(groups)))
    rng = np.random.default_rng(seed)
    rng.shuffle(unique_groups)
    assignment = {group: index % folds for index, group in enumerate(unique_groups)}
    fold_of = np.array([assignment[g] for g in groups])

    model_correct = 0
    antismash_correct = 0
    antismash_abstained = 0

    for fold in range(folds):
        test = fold_of == fold
        train = ~test
        if test.sum() == 0 or train.sum() == 0:
            continue
        # A label unseen in training can never be predicted; it still counts
        # against the model, which is the honest accounting.
        predictor = SubstratePredictor(kmer=kmer, n_neighbours=n_neighbours)
        predictor.fit(sequences[train], labels[train])
        predictions = predictor.predict(sequences[test])
        model_correct += int(sum(p == t for p, t in zip(predictions, labels[test])))

    for record, label in zip(usable, labels):
        call = normalise(record.get("antismash_call"))
        if call is None:
            antismash_abstained += 1
        elif call == label:
            antismash_correct += 1

    total = len(usable)
    return {
        "antismash": EvaluationResult("antiSMASH", total, antismash_correct, antismash_abstained),
        "sempi3": EvaluationResult("SeMPI-3 nearest neighbour", total, model_correct, 0),
    }
