"""Chromosome-arm reference coordinates.

The arm boundaries are read from the institutional reference assets so that
foundation-model arm assignment is IDENTICAL to the locked model's. Using a
different centromere definition would make every arm-relative position
incomparable between the two systems, which would silently invalidate the
baseline comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import pandas as pd

from .schema import normalise_chromosome


@dataclass(frozen=True)
class Arm:
    chromosome: str      # '1'
    arm: str             # 'p' | 'q'
    start: int
    end: int
    length: int

    @property
    def key(self) -> str:
        return f"{self.chromosome}{self.arm}"

    def contains(self, position: int) -> bool:
        return self.start <= position <= self.end

    def relative_position(self, position: int) -> float:
        """Arm-relative coordinate, 0.0 at the telomere and 1.0 at the centromere.

        This orientation matches the locked model's `mean_rel_tel_to_cent`
        convention, so p-arm positions run centromere-ward with increasing
        value and q-arm positions are reversed.
        """
        if self.length <= 0:
            raise ValueError(f"arm {self.key} has non-positive length")
        clamped = min(max(position, self.start), self.end)
        if self.arm == "p":
            # p-arm: telomere is at `start`, centromere at `end`.
            return (clamped - self.start) / self.length
        # q-arm: centromere is at `start`, telomere at `end`.
        return (self.end - clamped) / self.length


class ArmReference:
    """Lookup table mapping genomic coordinates to chromosome arms."""

    def __init__(self, arms: list[Arm], genome_build: str):
        self.genome_build = genome_build
        self._arms = arms
        self._by_chrom: dict[str, list[Arm]] = {}
        for a in arms:
            self._by_chrom.setdefault(a.chromosome, []).append(a)
        self._by_key = {a.key: a for a in arms}

    @classmethod
    def from_csv(cls, path: str | Path, genome_build: str) -> "ArmReference":
        df = pd.read_csv(path)
        expected = {"chromosome", "arm", "arm_length", "start", "end"}
        missing = expected - set(df.columns)
        if missing:
            raise ValueError(f"{path}: arm reference missing columns {sorted(missing)}")
        arms = [
            Arm(
                chromosome=normalise_chromosome(r.chromosome),
                arm=str(r.arm).strip(),
                start=int(r.start),
                end=int(r.end),
                length=int(r.arm_length),
            )
            for r in df.itertuples()
        ]
        return cls(arms, genome_build)

    @property
    def arm_keys(self) -> list[str]:
        return list(self._by_key)

    def get(self, arm_key: str) -> Arm:
        return self._by_key[arm_key]

    def assign(self, chromosome: str, position: int) -> Arm | None:
        """Return the arm containing `position`, or None if it falls in no arm.

        Returning None rather than raising is deliberate: acrocentric p-arms and
        centromere-spanning coordinates legitimately fall outside the reference,
        and those cases are counted in the QC report instead of crashing a run.
        """
        for a in self._by_chrom.get(chromosome, []):
            if a.contains(position):
                return a
        return None

    def assign_frame(self, chrom: pd.Series, pos: pd.Series) -> pd.Series:
        """Vectorised arm assignment returning arm keys (e.g. '1p'), NA if none."""
        out = pd.Series(pd.NA, index=chrom.index, dtype="string")
        for a in self._arms:
            hit = (chrom == a.chromosome) & (pos >= a.start) & (pos <= a.end)
            out[hit] = a.key
        return out


@lru_cache(maxsize=4)
def load_arm_reference(path: str, genome_build: str) -> ArmReference:
    return ArmReference.from_csv(path, genome_build)


def load_frozen_hotspot_windows(path: str | Path) -> pd.DataFrame:
    """Load the locked model's frozen arm-relative hotspot windows.

    These are POSITIVE CONTROLS only. They must never be hard-coded into the
    encoder or used as a training label (see prompt section 10). They exist here
    so that a rediscovered window can be scored against them.
    """
    df = pd.read_csv(path)
    required = {"arm_key", "hotspot_start_norm", "hotspot_end_norm"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: hotspot windows missing {sorted(missing)}")
    return df
