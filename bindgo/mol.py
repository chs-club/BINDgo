from typing import Literal

import selfies
from rdkit import Chem
from rdkit.Chem import AllChem


def clean_mol_str(
    mol_str: str,
    source_type: Literal["smiles"] = "smiles",
    dest_type: Literal["smiles", "selfies"] = "smiles",
    remove: list[str] = ["[Dy]"],
    remove_hydrogens: bool = True,
    keep_max_frag: bool = True,
) -> str:
    """
    Clean and standardize a given molecular string.

    Args:
        mol_str (str): The input molecular string representing the molecular structure.
        source_type (Literal["smiles"]): The type of the input molecular string. Default is "smiles".
        dest_type (Literal["smiles", "selfies"]): The desired output type of the molecular string. Default is "smiles".
        remove (list[str]): A list of strings to remove from the molecular string. Default is ["[Dy]"].
        remove_hydrogens (bool): Whether to remove explicit hydrogen atoms from the molecule. Default is True.
        keep_max_frag (bool): Whether to retain only the largest fragment of the molecule. Default is True.

    Returns:
        str: A cleaned and standardized molecular string in the specified output format (SMILES or SELFIES).

    Raises:
        ValueError: If the input molecular string is invalid and cannot be converted to an RDKit molecule object.

    Notes:
        - This function uses the RDKit library to handle molecule manipulation and standardization.
        - The function can remove specific substrings (e.g., "[Dy]") from the input molecular string.
        - Explicit hydrogen atoms can be removed from the molecule to simplify the structure.
        - Only the largest fragment of the molecule is retained, useful for removing salts or other small fragments.
        - 2D coordinates are computed to ensure the molecule is standardized before converting it back to the desired format.

    Example:
        ```python
        >>> clean_mol_str('CCO[Dy]CC')
        'CCOCC'

        >>> clean_mol_str('C(C)(C)C(=O)O.[Na]')
        'CC(C)C(=O)O'
        ```
    """
    for rm_str in remove:
        mol_str = mol_str.replace(rm_str, "")

    # If no post-processing requiring RDKit is requested, return the string.
    if not remove_hydrogens:
        return mol_str

    if source_type.strip().lower() == "smiles":
        mol = Chem.MolFromSmiles(mol_str)
    if mol is None:
        raise ValueError("Invalid mol string")

    if remove_hydrogens:
        mol = Chem.RemoveHs(mol)

    if keep_max_frag:
        fragments = Chem.GetMolFrags(mol, asMols=True)

        # Keep the largest fragment
        mol = max(fragments, default=mol, key=lambda m: m.GetNumAtoms())

        # Standardize the molecule
        AllChem.Compute2DCoords(mol)  # Compute 2D coordinates

    # Convert the molecule back to a canonical SMILES string
    if dest_type == "smiles":
        mol_str = Chem.MolToSmiles(mol, canonical=True)
        if dest_type == "selfies":
            mol_str = selfies.encoder(mol_str)
    return mol_str
