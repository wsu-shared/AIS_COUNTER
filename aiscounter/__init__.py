"""aiscounter -- AIS detection and length measurement.

A faithful Python port of ``original/ais_auto.m``, with a browser reviewer for add / delete /
join / splice and PNG + XLSX reporting. Measuring is click-driven like the original unless
``AnalysisConfig.auto_detect`` is on, in which case every AIS in the image is found first.

The numerics are validated against MATLAB R2024b; see ``tests/test_against_matlab.py``.
"""

from .config import AnalysisConfig
from .imaging import AISImage, find_images, load_image
from .measure import AISMeasurement
from .pipeline import AISRecord, AnalysisResult, analyse_image

__version__ = "1.0.0"

__all__ = [
    "AnalysisConfig",
    "AISImage",
    "AISMeasurement",
    "AISRecord",
    "AnalysisResult",
    "analyse_image",
    "find_images",
    "load_image",
]
