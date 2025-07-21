"""Pytest configuration and fixtures for kampo tests."""
import pytest
from pathlib import Path
import tempfile
import shutil
from typing import Generator

# Test data directory
TEST_DATA_DIR = Path(__file__).parent / "test_data"


@pytest.fixture
def test_data_dir() -> Path:
    """Return the path to the test data directory."""
    return TEST_DATA_DIR


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test outputs."""
    temp_path = Path(tempfile.mkdtemp())
    yield temp_path
    # Cleanup
    shutil.rmtree(temp_path)


@pytest.fixture
def sample_pdb_content() -> str:
    """Return sample PDB content for testing."""
    return """HEADER    SAMPLE PROTEIN
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.221   2.370   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       1.970  -0.790  -1.210  1.00 20.00           C
TER
HETATM    6  C1  LIG B   1       5.000   5.000   5.000  1.00 30.00           C
HETATM    7  C2  LIG B   1       6.000   5.000   5.000  1.00 30.00           C
HETATM    8  O1  LIG B   1       7.000   5.000   5.000  1.00 30.00           O
END
"""


@pytest.fixture
def sample_sdf_content() -> str:
    """Return sample SDF content for testing."""
    return """
  Sample Molecule

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""


@pytest.fixture
def sample_pdbqt_content() -> str:
    """Return sample PDBQT content for testing."""
    return """REMARK  Sample PDBQT file
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00    -0.350 N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00     0.100 C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00     0.200 C
ATOM      4  O   ALA A   1       1.221   2.370   0.000  1.00 20.00    -0.400 O
ATOM      5  CB  ALA A   1       1.970  -0.790  -1.210  1.00 20.00     0.000 C
TER
TORSDOF 0
"""


@pytest.fixture
def sample_multiframe_pdb_content() -> str:
    """Return sample multi-frame PDB content for testing."""
    return """MODEL     1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL     2
ATOM      1  N   ALA A   1       0.100   0.100   0.100  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.558   0.100   0.100  1.00 20.00           C
ENDMDL
END
"""