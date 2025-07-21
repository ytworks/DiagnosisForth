"""Tests for obabel converter functions."""
import pytest
from pathlib import Path
from unittest.mock import patch, call
from kampo.converters.obabel import (
    convert_pdbqt_to_pdb,
    convert_pdbqt_to_sdf,
)


class TestConvertPdbqtToPdb:
    """Test convert_pdbqt_to_pdb function."""
    
    @patch('subprocess.run')
    def test_convert_with_string_paths(self, mock_run):
        """Test conversion with string paths."""
        input_path = "input.pdbqt"
        output_path = "output.pdb"
        
        convert_pdbqt_to_pdb(input_path, output_path)
        
        expected_command = f"obabel -ipdbqt {input_path} -opdb -O{output_path}"
        mock_run.assert_called_once_with(expected_command, shell=True)
    
    @patch('subprocess.run')
    def test_convert_with_path_objects(self, mock_run):
        """Test conversion with Path objects."""
        input_path = Path("input.pdbqt")
        output_path = Path("output.pdb")
        
        convert_pdbqt_to_pdb(input_path, output_path)
        
        expected_command = f"obabel -ipdbqt {input_path} -opdb -O{output_path}"
        mock_run.assert_called_once_with(expected_command, shell=True)
    


class TestConvertPdbqtToSdf:
    """Test convert_pdbqt_to_sdf function."""
    
    @patch('subprocess.run')
    def test_convert_with_string_paths(self, mock_run):
        """Test conversion with string paths."""
        input_path = "input.pdbqt"
        output_path = "output.sdf"
        
        convert_pdbqt_to_sdf(input_path, output_path)
        
        expected_command = f"obabel -ipdbqt {input_path} -osdf -O{output_path}"
        mock_run.assert_called_once_with(expected_command, shell=True)
    
    @patch('subprocess.run')
    def test_convert_with_path_objects(self, mock_run):
        """Test conversion with Path objects."""
        input_path = Path("input.pdbqt")
        output_path = Path("output.sdf")
        
        convert_pdbqt_to_sdf(input_path, output_path)
        
        expected_command = f"obabel -ipdbqt {input_path} -osdf -O{output_path}"
        mock_run.assert_called_once_with(expected_command, shell=True)
    
