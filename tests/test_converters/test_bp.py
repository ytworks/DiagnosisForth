"""Tests for bp converter functions."""
import pytest
from pathlib import Path
from unittest.mock import patch, Mock, MagicMock, mock_open
import sys
from kampo.converters.bp import create_protein_ligand_complex


class TestCreateProteinLigandComplex:
    """Test create_protein_ligand_complex function."""
    
    @patch('kampo.converters.bp.Chem.SDMolSupplier')
    @patch('kampo.converters.bp.PDB.PDBParser')
    @patch('kampo.converters.bp.Chem.MolToPDBBlock')
    @patch('builtins.open', new_callable=mock_open)
    @patch('kampo.converters.bp.PDB.PDBIO')
    def test_successful_complex_creation(
        self, mock_pdbio_class, mock_file, mock_mol_to_pdb, 
        mock_parser_class, mock_supplier_class
    ):
        """Test successful protein-ligand complex creation."""
        # Setup mocks
        mock_ligand = Mock()
        mock_supplier = MagicMock()
        mock_supplier.__getitem__.return_value = mock_ligand
        mock_supplier_class.return_value = mock_supplier
        
        mock_parser = Mock()
        mock_structure = Mock()
        mock_parser.get_structure.return_value = mock_structure
        mock_parser_class.return_value = mock_parser
        
        mock_mol_to_pdb.return_value = "LIGAND PDB BLOCK"
        
        mock_pdbio = Mock()
        mock_pdbio_class.return_value = mock_pdbio
        
        # Execute
        ligand_path = "ligand.sdf"
        protein_path = "protein.pdb"
        output_path = "complex.pdb"
        
        create_protein_ligand_complex(ligand_path, protein_path, output_path)
        
        # Verify
        mock_supplier_class.assert_called_once_with(str(ligand_path), removeHs=False)
        mock_parser_class.assert_called_once_with(QUIET=True)
        mock_parser.get_structure.assert_called_once_with('protein', str(protein_path))
        mock_mol_to_pdb.assert_called_once_with(mock_ligand)
        
        mock_file.assert_called_once_with(output_path, 'w')
        mock_pdbio.set_structure.assert_called_once_with(mock_structure)
        mock_pdbio.save.assert_called_once()
        
        # Check file writing
        handle = mock_file()
        handle.write.assert_any_call('TER\n')
        handle.write.assert_any_call("LIGAND PDB BLOCK")
    
    @patch('kampo.converters.bp.Chem.SDMolSupplier')
    @patch('builtins.print')
    def test_failed_ligand_reading(self, mock_print, mock_supplier_class):
        """Test handling of failed ligand reading."""
        # Setup mock to return None (failed reading)
        mock_supplier = MagicMock()
        mock_supplier.__getitem__.return_value = None
        mock_supplier_class.return_value = mock_supplier
        
        # Execute and expect sys.exit
        with pytest.raises(SystemExit) as exc_info:
            create_protein_ligand_complex("ligand.sdf", "protein.pdb", "output.pdb")
        
        assert exc_info.value.code == 1
        mock_print.assert_called_once_with("Error: Couldn't read ligand.")
    
    
    def test_path_object_support(self):
        """Test that Path objects are properly converted to strings."""
        with patch('kampo.converters.bp.Chem.SDMolSupplier') as mock_supplier_class:
            mock_supplier = MagicMock()
            mock_supplier.__getitem__.return_value = None  # Will cause early exit
            mock_supplier_class.return_value = mock_supplier
            
            ligand_path = Path("ligand.sdf")
            protein_path = Path("protein.pdb")
            output_path = Path("output.pdb")
            
            with pytest.raises(SystemExit):
                create_protein_ligand_complex(ligand_path, protein_path, output_path)
            
            # Verify Path was converted to string
            mock_supplier_class.assert_called_once_with("ligand.sdf", removeHs=False)