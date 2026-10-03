
from __future__ import annotations

from adit.errors import AditError
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from adit.config import Config, Profile
from adit.spec import CalculationSpec
from adit.validate_types import ValidationError

Resources = Any


@dataclass
class ReadmeNotes:

    program: str
    files: list[str] = field(default_factory=list)
    prepare: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)


class GenerationError(AditError):
    pass


class InputGenerator(ABC):
    code: str
    supports_analysis: bool = True
    cli_only: bool = False
    uses_kpoints: bool = True

    @abstractmethod
    def resolve(self, spec: CalculationSpec, cfg: Config) -> Resources:
        pass

    @abstractmethod
    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        """Code-specific mechanical checks. Missing resources are returned as errors, not raised."""

    @abstractmethod
    def generate(self, spec: CalculationSpec, res: Resources) -> dict[str, str]:
        pass

    @abstractmethod
    def files_to_copy(self, spec: CalculationSpec, res: Resources) -> dict[str, Path]:
        pass

    @abstractmethod
    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        pass

    @abstractmethod
    def readme_notes(self, spec: CalculationSpec, res: Resources, copies: dict[str, Path]) -> ReadmeNotes:
        pass

    writes_velocities: bool = False

    def version_probe(self, spec: CalculationSpec) -> tuple[str, str] | None:
        return None

    def parameter_sources(self, spec: CalculationSpec, res: Resources) -> dict[str, Path]:
        return {}


def command_values(spec: CalculationSpec, *, binary: str = "") -> dict[str, object]:
    # Values for the {...} placeholders of a run command; ntasks is the MPI process count summed over all nodes.
    r = spec.runtime
    return {"mpiprocs": r.mpiprocs, "omp_threads": r.omp_threads, "ntasks": r.nodes * r.mpiprocs, "binary": binary}


GENERATORS: dict[str, InputGenerator] = {}


def register(gen: InputGenerator) -> InputGenerator:
    GENERATORS[gen.code] = gen
    return gen
