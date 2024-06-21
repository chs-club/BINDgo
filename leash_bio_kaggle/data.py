import pyarrow as pa
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    import torch
    from torch.utils.data import Dataset

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


class MolDataset(Dataset):

    protein_seq = {
        "sEH": "TLRAAVFDLDGVLALPAVFGVLGRTEEALALPRGLLNDAFQKGGPEGATTRLMKGEITLSQWIPLMEENCRKCSETAKVCLPKNFSIKEIFDKAISARKINRPMLQAALMLRKKGFTTAILTNTWLDDRAERDGLAQLMCELKMHFDFLIESCQVGMVKPEPQIYKFLLDTLKASPSEVVFLDDIGANLKPARDLGMVTILVQDTDTALKELEKVTGIQLLNTPAPLPTSCNPSDMSHGYVTVKPRVRLHFVELGSGPAVCLCHGFPESWYSWRYQIPALAQAGYRVLAMDMKGYGESSAPPEIEEYCMEVLCKEMVTFLDKLGLSQAVFIGHDWGGMLVWYMALFYPERVRAVASLNTPFIPANPNMSPLESIKANPVFDYQLYFQEPGVAEAELEQNLSRTFKSLFRASDESVLSMHKVCEAGGLFVNSPEEPSLSRMVTEEEIQFYVQQFKKSGFRGPLNWYRNMERNWKWACKSLGRKILIPALMVTAEKDFVLVPQMSQHMEDWIPHLKRGHIEDCGHWTQMDKPTEVNQILIKWLDSDARNPPVVSKM",
        "BRD4": "NPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEETEIMIVQAKGRGRGRKETGTAKPGVSTVPNTTQASTPPQTQTPQPNPPPVQATPHPFPAVTPDLIVQTPVMTVVPPQPLQTPPPVPPQPQPPPAPAPQPVQSHPPIIAATPQPVKTKKGVKRKADTTTPTTIDPIHEPPSLPPEPKTTKLGQRRESSRPVKPPKKDVPDSQQHPAPEKSSKVSEQLKCCSGILKEMFAKKHAAYAWPFYKPVDVEALGLHDYCDIIKHPMDMSTIKSKLEAREYRDAQEFGADVRLMFSNCYKYNPPDHEVVAMARKLQDVFEMRFAKMPDE",
        "HSA": "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL",
    }

    def __init__(self, ds: pa.dataset.Dataset, *args, **kwargs):
        self.ds = ds

    def get_protein_seq(self, key: str) -> str:
        return self.protein_seq[key]

    def get_sample(self, idx, kind):
        if kind not in ("bind", "no-bind"):
            raise ValueError("kind must be `bind` or `no-bind`")

        record = None
        if kind == "bind":
            record = self.scanner_bind.take([idx]).to_pydict()
        else:
            record = self.scanner_no_bind.take([idx]).to_pydict()

        smiles = record["molecule_smiles"][0]
        protein_seq = record["protein_name"][0]
        label = record["binds"][0]

        smiles = self.clean_smiles(smiles)

        amino_acids = self.get_protein_seq(protein_seq)
        return smiles, amino_acids, label

    def __getitem__(self, bind_idx: int) -> (str, str, int):
        bind_data = self.get_sample(bind_idx, kind="bind")
        # TODO: Fix random choice to select one with the same type of protein.
        no_bind_idx = random.choice(self.indices_no_bind)
        no_bind_data = self.get_sample(no_bind_idx, kind="no-bind")
        return bind_data, no_bind_data

    def __len__(self):
        return len(self.indices_bind)

    def featurize_ligand(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add hydrogen atoms
        AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates

        # Node features
        node_features = []
        for atom in mol.GetAtoms():
            features = [
                atom.GetAtomicNum(),
                atom.GetDegree(),
                atom.GetFormalCharge(),
                atom.GetHybridization(),
                atom.GetTotalNumHs(),
                int(atom.GetIsAromatic()),
                int(atom.IsInRing()),
                int(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED),
            ]
            node_features.append(features)

        # Edge features and edge index
        edge_features = []
        edge_index = []
        for bond in mol.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            features = [
                bond.GetBondTypeAsDouble(),
                int(bond.IsInRing()),
                bond.GetBondLength(),
                # You can add bond angle here if needed
            ]
            edge_features.append(features)
            edge_features.append(features)  # Add twice for undirected graph
            edge_index.append([start, end])
            edge_index.append([end, start])  # Add reverse edge for undirected graph

        return (
            torch.tensor(node_features, dtype=torch.float),
            torch.tensor(edge_index, dtype=torch.long).t().contiguous(),
            torch.tensor(edge_features, dtype=torch.float),
        )
