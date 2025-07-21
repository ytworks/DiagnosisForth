"""Type definitions for kampo library."""
from typing import TypedDict, List, Tuple, Literal, Dict, Any, Optional

# Pharmacophore feature types
PharmacophoreFeature = Literal[
    "Donor", "Acceptor", "Hydrophobe", 
    "PosIonizable", "NegIonizable", "Aromatic"
]

# Interaction types supported by PLIP
InteractionType = Literal[
    "hydrophobic", "hydrogen_bond", "water_bridge",
    "salt_bridge", "pi_stacking", "pi_cation", "halogen"
]


class InteractionPoint(TypedDict):
    """Represents a point in a molecular interaction."""
    residue: str
    chain: str
    position: int
    atom: str
    coordinates: Tuple[float, float, float]


class Interaction(TypedDict):
    """Represents a molecular interaction between protein and ligand."""
    type: InteractionType
    protein: InteractionPoint
    ligand: InteractionPoint
    distance: float
    angle: Optional[float]


# Type alias for interaction data structure
InteractionData = Dict[str, List[Dict[str, Any]]]

# Type alias for pharmacophore features
PharmacophoreFeatures = Dict[PharmacophoreFeature, List[Tuple[float, float, float]]]