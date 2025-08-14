"""Tests for pdb_reader functions."""
import pytest
from pathlib import Path
from unittest.mock import patch, Mock, mock_open
from shishin.readers.pdb_reader import read_multiframe_pdb


class TestReadMultiframePdb:
    """Test read_multiframe_pdb function."""
    
    def test_read_single_frame(self, sample_multiframe_pdb_content, temp_dir):
        """Test reading a PDB file with a single frame."""
        # Create test file
        test_file = temp_dir / "single_frame.pdb"
        test_file.write_text("""MODEL     1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ENDMDL
END
""")
        
        with patch('shishin.readers.pdb_reader.Chem.MolFromPDBBlock') as mock_mol_from_pdb:
            mock_mol = Mock()
            mock_mol_from_pdb.return_value = mock_mol
            
            # Execute
            result = read_multiframe_pdb(test_file)
            
            # Verify
            assert len(result) == 1
            assert result[0] == mock_mol
            mock_mol_from_pdb.assert_called_once()
    
    def test_read_multiple_frames(self, sample_multiframe_pdb_content, temp_dir):
        """Test reading a PDB file with multiple frames."""
        # Create test file
        test_file = temp_dir / "multi_frame.pdb"
        test_file.write_text(sample_multiframe_pdb_content)
        
        with patch('shishin.readers.pdb_reader.Chem.MolFromPDBBlock') as mock_mol_from_pdb:
            mock_mol1 = Mock()
            mock_mol2 = Mock()
            mock_mol_from_pdb.side_effect = [mock_mol1, mock_mol2]
            
            # Execute
            result = read_multiframe_pdb(test_file)
            
            # Verify
            assert len(result) == 2
            assert result[0] == mock_mol1
            assert result[1] == mock_mol2
            assert mock_mol_from_pdb.call_count == 2
    
    def test_read_empty_file(self, temp_dir):
        """Test reading an empty PDB file."""
        # Create empty test file
        test_file = temp_dir / "empty.pdb"
        test_file.write_text("")
        
        with patch('shishin.readers.pdb_reader.Chem.MolFromPDBBlock') as mock_mol_from_pdb:
            # Execute
            result = read_multiframe_pdb(test_file)
            
            # Verify
            assert len(result) == 0
            mock_mol_from_pdb.assert_not_called()
    
    def test_string_path_support(self, temp_dir):
        """Test that string paths are supported."""
        # Create test file
        test_file = temp_dir / "test.pdb"
        test_file.write_text("""MODEL     1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ENDMDL
""")
        
        with patch('shishin.readers.pdb_reader.Chem.MolFromPDBBlock') as mock_mol_from_pdb:
            mock_mol = Mock()
            mock_mol_from_pdb.return_value = mock_mol
            
            # Execute with string path
            result = read_multiframe_pdb(str(test_file))
            
            # Verify
            assert len(result) == 1
            assert result[0] == mock_mol
    
    @patch('shishin.readers.pdb_reader.tqdm')
    def test_progress_bar_usage(self, mock_tqdm, temp_dir):
        """Test that tqdm progress bar is used."""
        # Create test file
        test_file = temp_dir / "test.pdb"
        test_file.write_text("""MODEL     1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ENDMDL
""")
        
        # Setup mock
        mock_tqdm.return_value = test_file.read_text().splitlines(keepends=True)
        
        with patch('shishin.readers.pdb_reader.Chem.MolFromPDBBlock'):
            # Execute
            read_multiframe_pdb(test_file)
            
            # Verify tqdm was called
            mock_tqdm.assert_called_once()
    
    
    def test_frame_content_parsing(self, temp_dir):
        """Test correct parsing of frame content."""
        # Create test file with specific content
        test_file = temp_dir / "test.pdb"
        test_file.write_text("""MODEL     1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL     2
ATOM      1  N   GLY A   2       2.000   2.000   2.000  1.00 30.00           N
ENDMDL
""")
        
        captured_frames = []
        
        def capture_frame(pdb_block):
            captured_frames.append(pdb_block)
            return Mock()
        
        with patch('kampo.readers.pdb_reader.Chem.MolFromPDBBlock', side_effect=capture_frame):
            # Execute
            read_multiframe_pdb(test_file)
            
            # Verify correct frame content
            assert len(captured_frames) == 2
            assert "ALA" in captured_frames[0]
            assert "GLY" in captured_frames[1]
            assert "MODEL" not in captured_frames[0]  # MODEL line should be excluded
            assert "ENDMDL" not in captured_frames[0]  # ENDMDL line should be excluded