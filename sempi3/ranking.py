"""Choosing between assembled structures by what the cluster predicts.

SeMPI-3 can build several structures for one cluster: the chain released four
ways, bare or decorated. A rule used to pick one. Here the cluster is asked
what its product should look like — how heavy, how many nitrogens and rings,
whether it carries a sugar or a macrocycle — and the candidate that fits that
description best is kept.

The chemistry still comes from the assembly line; the models only rank what
the chemistry produced, so every answer is a molecule that could have been
made.

Models are fitted from a table shipped with the package
(``data/cluster_properties.tsv``), not from pickles. Fitting takes a few
seconds and is cached for the life of the process, and the table is rebuilt
by ``sempi3 train`` whenever the underlying data changes. That keeps the
package free of model files pinned to one scikit-learn version.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
)

import numpy as np

from . import structure_features as sf

TRAINING_TABLE = Path(__file__).parent / "data" / "cluster_properties.tsv"

#: Properties used for ranking, with the weight each carries.
#:
#: The weights follow how well each is predicted: size and nitrogen content
#: are predicted well (R² above 0.7) and count for more than ring counts,
#: which are noisier.
WEIGHTS: Dict[str, float] = {
    "n_heavy_atoms": 1.0,
    "n_nitrogen": 1.0,
    "n_amide": 1.0,
    "n_oxygen": 0.7,
    "n_rings": 0.5,
    "has_macrocycle": 0.5,
    "has_sugar": 0.4,
    "has_aromatic_ring": 0.4,
}


@dataclass
class TrainingSet:
    """The table the property models are fitted from."""

    feature_names: List[str]
    features: np.ndarray
    targets: Dict[str, np.ndarray]

    @property
    def rows(self) -> int:
        return len(self.features)


def write_training_table(
    feature_names: Sequence[str],
    features: np.ndarray,
    targets: Dict[str, Sequence[float]],
    path: Optional[Path] = None,
) -> Path:
    """Write the table ``sempi3 train`` produces."""
    path = Path(path) if path else TRAINING_TABLE
    path.parent.mkdir(parents=True, exist_ok=True)
    target_names = sorted(targets)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(list(feature_names) + [f"target:{name}" for name in target_names])
        for index in range(len(features)):
            row = [f"{value:g}" for value in features[index]]
            row += [f"{targets[name][index]:g}" for name in target_names]
            writer.writerow(row)
    return path


@lru_cache(maxsize=1)
def load_training_table(path: Optional[str] = None) -> Optional[TrainingSet]:
    """Read the packaged training table, or ``None`` when it is absent."""
    table = Path(path) if path else TRAINING_TABLE
    if not table.exists():
        return None
    with open(table) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return None
        feature_names = [name for name in header if not name.startswith("target:")]
        target_names = [name[len("target:") :] for name in header if name.startswith("target:")]
        rows = [row for row in reader if len(row) == len(header)]
    if not rows:
        return None

    width = len(feature_names)
    features = np.array([[float(value) for value in row[:width]] for row in rows])
    targets = {
        name: np.array([float(row[width + index]) for row in rows])
        for index, name in enumerate(target_names)
    }
    return TrainingSet(feature_names=feature_names, features=features, targets=targets)


@lru_cache(maxsize=1)
def _models() -> Optional[Tuple[List[str], Dict[str, object]]]:
    """Fit one model per property, once per process."""
    data = load_training_table()
    if data is None or data.rows < 50:
        return None
    from sklearn.ensemble import (
        RandomForestClassifier,
        RandomForestRegressor,
    )

    fitted: Dict[str, object] = {}
    for name in WEIGHTS:
        target = data.targets.get(name)
        if target is None:
            continue
        if name in sf.BINARY_TARGETS:
            model = RandomForestClassifier(
                n_estimators=200, min_samples_leaf=2, class_weight="balanced", random_state=0, n_jobs=-1
            )
        else:
            model = RandomForestRegressor(
                n_estimators=200, min_samples_leaf=2, random_state=0, n_jobs=-1
            )
        model.fit(data.features, target)
        fitted[name] = model
    return data.feature_names, fitted


def available() -> bool:
    """True when ranking can be used."""
    try:
        return _models() is not None
    except Exception:  # noqa: BLE001 - a missing scikit-learn must not break prediction
        return False


def predict_properties(cluster_features) -> Dict[str, float]:
    """What the cluster says its product should look like."""
    models = _models()
    if models is None:
        return {}
    names, fitted = models
    vector = cluster_features.vector(names).reshape(1, -1)
    return {name: float(model.predict(vector)[0]) for name, model in fitted.items()}


def distance(candidate_smiles: str, wanted: Dict[str, float]) -> float:
    """How far a candidate sits from the predicted description. Lower is better."""
    profile = sf.describe(candidate_smiles)
    if profile is None:
        return float("inf")
    total = 0.0
    for name, weight in WEIGHTS.items():
        if name not in wanted:
            continue
        predicted = wanted[name]
        actual = profile.get(name)
        if name in sf.BINARY_TARGETS:
            total += weight * abs(predicted - actual)
        else:
            scale = max(abs(predicted), 1.0)
            total += weight * min(abs(predicted - actual) / scale, 2.0)
    return total


def choose(candidates: Dict[str, str], cluster_features) -> Optional[Tuple[str, str, float]]:
    """Pick the candidate that best fits the predicted description.

    Returns ``(label, smiles, distance)``, or ``None`` when ranking is
    unavailable, in which case the caller keeps its own choice.
    """
    if not candidates or not available():
        return None
    wanted = predict_properties(cluster_features)
    if not wanted:
        return None
    best: Optional[Tuple[str, str, float]] = None
    for label, smiles in candidates.items():
        if not smiles:
            continue
        score = distance(smiles, wanted)
        if best is None or score < best[2]:
            best = (label, smiles, score)
    return best
