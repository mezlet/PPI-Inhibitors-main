# Complete PPI Inhibitors Pipeline: End-to-End Analysis
## Comprehensive Data Loading and Preprocessing Guide

---

# TABLE OF CONTENTS
1. Dataset File Loading & Parsing
2. Three-Stream Feature Extraction
3. Preprocessing Steps
4. Model Architecture
5. Data Flow & Pipeline Integration
6. Key Technical Details

---

# 1. DATASET FILE LOADING & PARSING

## 1.1 Source File Format

**File Name**: `WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt`

**Location**: Downloaded from GitHub during setup:
```
https://github.com/adibayaseen/PPI-Inhibitors/raw/2d6bd03422602ec19147870c487e64018b52660f/Data/WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt
```

**File Format Description**:
- Each line represents one training example
- Tab-separated values with format: `TestPosComp ComplexName Ligand_SMILES Label`
- **Example**: `2XA0_A_2_B 2XA0 CC(C)Cc1ccc(cc1)C(C)C(=O)O 1.0`

## 1.2 Data Loading Process (Cell 23)

```python
print("Loading training data...")

with open(githubpath + 'Data/WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt') as f:
    D = f.readlines()

Labels = []
Ligandnames = []
Complexs = []
TestPoscomplexes = []

for d in tqdm(D, desc="Loading examples"):
    if len(d.split()) == 4:
        # Standard case: 4 tokens
        TestPoscomp, Complexname, Ligandname, label = d.split()
    else:
        # Handle SMILES with spaces
        parts = d.split()
        TestPoscomp = parts[0]
        Complexname = parts[1]
        Ligandname = ' '.join(parts[2:-1])  # Join all middle tokens
        label = parts[-1]
    
    TestPoscomplexes.append(TestPoscomp)
    Ligandnames.append(Ligandname)
    Complexs.append(Complexname)
    Labels.append(float(label))

# Create unified dictionary
Allexamples = dict(zip(zip(TestPoscomplexes, zip(Complexs, Ligandnames)), Labels))

# Dataset statistics
print(f"Positive examples: {sum(1 for v in Labels if v == 1.0)}")
print(f"Negative examples: {sum(1 for v in Labels if v == -1.0 or v == 0.0)}")
```

## 1.3 Data Structure Created

```python
# Dictionary format:
# Key: (TestPosComp, (ComplexName, Ligand_SMILES))
# Value: Label (1.0 for active, -1.0/0.0 for inactive)

Allexamples = {
    ('2XA0_A_2_B', ('2XA0', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O')): 1.0,
    ('2XA0_A_2_B', ('2XA0', 'different_smiles_here')): -1.0,
    ...
}
```

## 1.4 Dataset Statistics

- **Format**: Binary classification (1.0 = inhibitor/binder, -1.0/0.0 = non-inhibitor)
- **Variable SMILES Length**: Can contain spaces (handled by joining middle tokens)
- **Multiple Examples per Complex**: Same complex paired with different ligands
- **Pre-filtered**: Invalid SMILES are filtered during feature extraction

---

# 2. THREE-STREAM FEATURE EXTRACTION

The pipeline extracts three independent feature streams that are later concatenated:

```
Input Data
    ↓
    ├─→ Stream 1: COMPOUND FEATURES (Ligand fingerprints)
    ├─→ Stream 2: PROTEIN FEATURES (GNN from structure)
    └─→ Stream 3: INTERFACE FEATURES (Inter-chain interaction patterns)
    ↓
Concatenate & Feed to MLP
```

## 2.1 Stream 1: COMPOUND FEATURES

### Function: `getFP(s, r=3, nBits=2048)`

**Purpose**: Generate chemical fingerprints from SMILES strings

**Code**:
```python
def getFP(s, r=3, nBits=2048):
    """
    Generate Morgan fingerprint for a compound
    
    Args:
        s: SMILES string
        r: Radius for Morgan fingerprint (default=3)
        nBits: Number of bits in fingerprint (default=2048)
    
    Returns:
        np.ndarray: 2048-dimensional binary vector (int8)
        None: If SMILES is invalid
    """
    compound = Chem.MolFromSmiles(s.strip())
    if compound is not None:
        fp = AllChem.GetMorganFingerprintAsBitVect(compound, r, nBits=nBits)
        m = np.zeros((0,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, m)
        return m
    return None  # Invalid SMILES filtered automatically
```

**Feature Details**:
- **Type**: Morgan Circular Fingerprint (Identifier-based)
- **Radius**: 3 (describes atoms up to 3 bonds away)
- **Dimensionality**: 2048 bits
- **Data Type**: Binary (int8: 0 or 1)
- **Preprocessing**: Stripped whitespace, validity checked

**Processing Pipeline for Compounds**:
```python
# 1. Extract fingerprints for all training compounds
Ctr = []
Ctrname = []
for compound_name in compound_names:
    fp = getFP(smiles_dict[compound_name])
    if fp is not None:  # Filter invalid SMILES
        Ctr.append(fp)
        Ctrname.append(compound_name)

# 2. Standardize features
Cscaler = StandardScaler()
Ctr = Cscaler.fit_transform(Ctr)

# 3. Convert to PyTorch tensors
Ctrdict = dict(zip(Ctrname, torch.FloatTensor(Ctr).cuda()))
```

## 2.2 Stream 2: PROTEIN FEATURES (GNN-based)

### Function: `processProtein(UniqueProtein, PdBloc)`

**Purpose**: Convert PDB protein structures to graph neural network-ready representations

**Input**: 
- PDB file paths for protein complexes
- Approximately 290 unique PDB files in training set

**Extraction Process**:

```python
def processProtein(UniqueProtein, PdBloc):
    """
    Process PDB files into GNN-ready graph representations
    
    Args:
        UniqueProtein: List of protein PDB filenames (without .pdb)
        PdBloc: Directory path containing PDB files
    
    Returns:
        PData_dict: Dictionary mapping protein name -> GNN input data
    """
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")
    PData_dict = {}
    
    for i in range(len(UniqueProtein)):
        UniqueProtein[i] = UniqueProtein[i].split('.pdb')[0]
        P1 = PdBloc + UniqueProtein[i] + '.pdb'
        
        # 1. Parse PDB structure
        parser = PDBParser()
        with warnings.catch_warnings(record=True):
            structure = parser.get_structure("", P1)
        
        # 2. Extract structural features
        one_hot_atom = atom1(structure)      # Atom type encoding
        one_hot_res = res1(structure)        # Residue type encoding
        neigh_same_res, neigh_diff_res = neigh1(structure)  # Connectivity
        
        # 3. Convert to PyTorch tensors
        one_hot_atom = torch.tensor(one_hot_atom, dtype=torch.float32).to(device)
        one_hot_res = torch.tensor(one_hot_res, dtype=torch.float32).to(device)
        neigh_same_res = torch.tensor(neigh_same_res).to(device).long()
        neigh_diff_res = torch.tensor(neigh_diff_res).to(device).long()
        
        # 4. Store as GNN input
        GNNData = [one_hot_atom, one_hot_res, neigh_same_res, neigh_diff_res]
        PData_dict[UniqueProtein[i]] = GNNData
    
    return PData_dict
```

**Feature Components**:

| Component | Description | Dimension |
|---|---|---|
| `one_hot_atom` | One-hot encoding of atom types (C, N, O, S, etc.) | Variable |
| `one_hot_res` | One-hot encoding of residue types (20 amino acids) | Variable |
| `neigh_same_res` | Graph edges: atoms within same residue | Sparse indices |
| `neigh_diff_res` | Graph edges: atoms across residues | Sparse indices |

**GNN Processing**:
```python
GNN_model = GNN().cuda()

# Forward pass during training:
# Input: [one_hot_atom, one_hot_res, neigh_same_res, neigh_diff_res]
# Output: 512-dimensional GNN features
gnn_features = GNN_model(protein_data[complex_id])  # Shape: [1, 512]
```

## 2.3 Stream 3: INTERFACE FEATURES

### Function: `chainLabel(Cname_T, xl_T, Cname, xl)`

**Purpose**: Extract inter-chain interaction patterns from protein complex structures

**Components of Interface Feature Extraction**:

### 2.3.1 Distance Matrix Calculation

```python
def getDist(C0, C1, thr=8.0):
    """
    Calculate distances between atoms of two protein chains
    
    Args:
        C0: List of atom coordinates for chain 1
        C1: List of atom coordinates for chain 2
        thr: Distance threshold in Angstroms (default: 8.0)
    
    Returns:
        N0, N1: Lists of (residue_i, residue_j, distance) tuples
    """
    N0 = []
    N1 = []
    for i in range(len(C0)):
        for j in range(len(C1)):
            # Calculate minimum distance between residue atoms
            d = spatial.distance.cdist(C0[i], C1[j]).min()
            
            # Record if below threshold
            if (d < thr):
                N0.append((i, j, d))      # Forward direction
                N1.append((j, i, d))      # Backward direction
    
    return (N0, N1)
```

**Key Parameter**: 
- **Distance Threshold**: 8.0 Angstroms (defines "interfacing" residues)

### 2.3.2 Amino Acid Pair Feature Generation

```python
def generate_pair_features(dist_info, xl, xr):
    """
    Count amino acid pair interactions in the interface
    
    Args:
        dist_info: List of (res_i, res_j, distance) tuples
        xl, xr: Residue objects from two chains
    
    Returns:
        prot_dic: Dictionary of AA pair counts
    """
    prot_dic = make_dic()  # Initialize amino acid pair dictionary
    
    for rec in dist_info:
        try:
            # Get amino acid one-letter code
            l_letter = three_to_one(xl[rec[0]].get_resname())
            r_letter = three_to_one(xr[rec[1]].get_resname())
            
            # Count pair occurrences (handle symmetric pairs)
            if (l_letter, r_letter) in prot_dic.keys():
                prot_dic[(l_letter, r_letter)] += 1
            elif (r_letter, l_letter) in prot_dic.keys():
                prot_dic[(r_letter, l_letter)] += 1
        except:
            # Non-standard amino acids
            prot_dic[('_','_')] += 1
    
    return prot_dic
```

**Amino Acid Dictionary Initialization**:

```python
def make_dic():
    """Create dictionary for all 20 amino acid pairs + non-standard"""
    prot_dic = {}
    letters = IUPACData.protein_letters  # ['A', 'C', 'D', ..., 'Y']
    
    # Create all unique pairs (including symmetric)
    for i in range(len(letters)):
        for j in range(i, len(letters)):
            prot_dic[(letters[i], letters[j])] = 0.0
    
    # Add entry for non-standard amino acids
    prot_dic[('_','_')] = 0.0
    
    return prot_dic  # Total: 210 + 1 = 211 features
```

**Amino Acid Pairs**:
- **Standard amino acids**: 20 (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
- **Unique pairs**: C(20,2) + 20 = 210 pairs
- **Non-standard**: 1 entry for anomalies
- **Total interface dimensions**: 211

### 2.3.3 Feature Extraction & Normalization

```python
def extract_feats(dic):
    """Extract features in consistent order"""
    feats = []
    # Load pre-defined key ordering from file
    key_list = np.load('/content/PPI-Inhibitors/Features/prote_letter_pair_keys.npy')
    
    for key in key_list:
        # Decode keys (saved as bytes in numpy)
        feats.append(dic[(key[0].decode('utf-8'), key[1].decode('utf-8'))])
    
    return feats  # Ordered list of 211 features
```

**Why Ordering Matters**: 
- Same order across all examples ensures consistent feature interpretation
- Pre-computed key ordering prevents non-deterministic dictionary iteration

**L2 Normalization**:
```python
from sklearn.preprocessing import normalize

InterfaceF = np.array(Interface)
InterfaceF = normalize(np.atleast_2d(InterfaceF), 
                       norm='l2',        # L2 norm
                       copy=True, 
                       axis=1,           # Normalize per row
                       return_norm=False)
```

**Effect**: Normalizes interaction counts to unit norm, emphasizing pattern rather than absolute count

### 2.3.4 Complete Chain Labeling Process

```python
def chainLabel(Cname_T, xl_T, Cname, xl):
    """
    Complete interface feature extraction
    
    Args:
        Cname_T: Target chain name (e.g., 'A')
        xl_T: Target chain residues
        Cname: Other chain name (e.g., 'B')
        xl: Other chain residues
    """
    # Step 1: Get coordinates
    tc = getCoords(xl_T)  # List of atom coord arrays
    nc = getCoords(xl)
    
    # Step 2: Calculate distances (8.0 Å threshold)
    D = getDist(tc, nc, thr=8.0)
    
    # Step 3: Generate and extract pair features
    feats = extract_feats(generate_pair_features(D, xl_T, xl))
    
    return feats  # Normalized 211-dimensional interface features
```

### 2.3.5 Pre-extracted Interface Features

**Important Note**: The notebook uses pre-computed interface features loaded from pickle files rather than computing them on-the-fly:

```python
print("Loading interface features...")

# Load pre-extracted features
Ubench5_dict = pickle.load(open(githubpath + 'Features/NewUbench5InterfaceandSeq_dict.npy', "rb"))
Pos_dict = pickle.load(open(githubpath + 'Features/Pos_seqandInterfaceF_dict.npy', "rb"))

# Merge dictionaries
Complex_AllFeatures = dict(list(Pos_dict.items()) + list(Ubench5_dict.items()))

# Create simplified mapping
ComplexInterfaceFeatures = {}
for key in Complex_AllFeatures:
    if '_' in str(key):
        # Handle naming format: "2XA0_chains" -> "2XA0"
        compname = key.split('_')[0]
        ComplexInterfaceFeatures[compname] = Complex_AllFeatures[key]
    else:
        ComplexInterfaceFeatures[key] = Complex_AllFeatures[key]

print(f"Loaded {len(ComplexInterfaceFeatures)} protein complexes")
```

**Pre-extracted Feature Files**:
1. **Pos_seqandInterfaceF_dict.npy**: Interface features for positive examples
2. **NewUbench5InterfaceandSeq_dict.npy**: Interface features for Ubench5 benchmark

---

# 3. COMPREHENSIVE PREPROCESSING PIPELINE

## 3.1 Data Preprocessing Stages

### Stage 1: Data Collection
```python
# Load raw data from file
with open('WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt') as f:
    lines = f.readlines()

# Parse examples
for line in lines:
    parts = line.split()
    pdb_id = parts[0]
    complex_name = parts[1]
    smiles = ' '.join(parts[2:-1])  # Handle spaces in SMILES
    label = float(parts[-1])
```

### Stage 2: Ligand Preprocessing
```python
# For each compound:
# 1. Extract SMILES
# 2. Generate Morgan fingerprint
# 3. Validate (filter None values)
# 4. Convert to numpy array

valid_ligands = []
for smiles in all_smiles:
    fp = getFP(smiles)  # Returns None if invalid
    if fp is not None:
        valid_ligands.append(fp)
```

### Stage 3: Protein Complex Preprocessing
```python
# Option A: Process from PDB files
External_Protein_GNN_Data_dict = processProtein(UniqueProtein, Pdbloc)

# Option B: Load pre-extracted features
ComplexInterfaceFeatures = pickle.load(...)
```

### Stage 4: Feature Scaling
```python
from sklearn.preprocessing import StandardScaler

# Initialize scalers
Pscaler = StandardScaler()
Cscaler = StandardScaler()

# Fit on training data
Ptr_scaled = Pscaler.fit_transform(Ptr)    # Interface features
Ctr_scaled = Cscaler.fit_transform(Ctr)    # Compound features

# Standardization removes mean and divides by std dev
# Formula: X_scaled = (X - mean) / std
```

**Why StandardScaler?**
- Centers features around 0
- Normalizes variance to 1
- Improves MLP training stability
- Prevents features with large ranges from dominating

### Stage 5: Tensor Conversion & GPU Transfer
```python
# Create dictionaries for fast lookup during training
Ptrdict = dict(zip(Ptrname, torch.FloatTensor(Ptr_scaled).cuda()))
Ctrdict = dict(zip(Ctrname, torch.FloatTensor(Ctr_scaled).cuda()))

# Format:
# Ptrdict[complex_id] -> tensor of shape [1, features]
# Ctrdict[compound_id] -> tensor of shape [1, features]
```

---

# 4. MODEL ARCHITECTURE

## 4.1 GNN Model for Protein Processing

### GNN_First_Layer

```python
class GNN_First_Layer(nn.Module):
    """Combines atom features, residue features, and neighborhood information"""
    
    def __init__(self, filters=512):
        super(GNN_First_Layer, self).__init__()
        self.filters = filters
        
        # Weight matrices for different node/edge types
        self.Wv = nn.Parameter(torch.randn(13, filters))          # Atom features
        self.Wr = nn.Parameter(torch.randn(21, filters))          # Residue features
        self.Wsr = nn.Parameter(torch.randn(13, filters))         # Same-residue neighbors
        self.Wdr = nn.Parameter(torch.randn(13, filters))         # Different-residue neighbors
    
    def forward(self, x):
        atoms, residues, same_neigh, diff_neigh = x
        
        # Project features through weight matrices
        node_signals = atoms @ self.Wv                    # Atom feature projection
        residue_signals = residues @ self.Wr              # Residue feature projection
        neigh_signals_same = atoms @ self.Wsr             # Same-residue neighbor projection
        neigh_signals_diff = atoms @ self.Wdr             # Different-residue neighbor projection
        
        # Aggregate neighbor signals
        unsqueezed_same = (same_neigh > -1).unsqueeze(2)
        unsqueezed_diff = (diff_neigh > -1).unsqueeze(2)
        
        same_neigh_features = neigh_signals_same[same_neigh] * unsqueezed_same
        diff_neigh_features = neigh_signals_diff[diff_neigh] * unsqueezed_diff
        
        # Average neighbor signals
        same_norm = torch.sum(same_neigh > -1).type(torch.float)
        diff_norm = torch.sum(diff_neigh > -1).type(torch.float)
        same_norm = torch.clamp(same_norm, min=1.0)
        diff_norm = torch.clamp(diff_norm, min=1.0)
        
        neigh_same_signal = torch.sum(same_neigh_features, axis=1) / same_norm
        neigh_diff_signal = torch.sum(diff_neigh_features, axis=1) / diff_norm
        
        # Combine all signals with ReLU activation
        final_res = torch.relu(node_signals + residue_signals + 
                               neigh_same_signal + neigh_diff_signal)
        
        return final_res, same_neigh, diff_neigh  # Shape: [num_atoms, 512]
```

## 4.2 Integration MLP: IPPI_MLP_Net

```python
class IPPI_MLP_Net(nn.Module):
    """
    Multi-task Integration Network
    Combines three feature streams into a single prediction
    """
    
    def __init__(self):
        super(IPPI_MLP_Net, self).__init__()
        
        # Feature dimensions:
        # - GNN output: 512 (from GNN model)
        # - Interface: 211 (amino acid pair counts, pre-processed)
        # - Compound: 2048 (Morgan fingerprint)
        # Total input: 512 + 211 + 2048 = 2771 (or variant)
        
        # MLP layers
        self.fc1 = nn.Linear(2840, 1024)    # Input layer
        self.fc2 = nn.Linear(1024, 512)     # Hidden layer 1
        self.fc3 = nn.Linear(512, 100)      # Hidden layer 2
        self.fc4 = nn.Linear(100, 1)        # Output layer
    
    def forward(self, gnn_features, compound_features, interface_features):
        """
        Args:
            gnn_features: [batch_size, 512] - GNN output
            compound_features: [batch_size, 2048] - Morgan fingerprint
            interface_features: [batch_size, 211] - AA pair features
        
        Returns:
            predictions: [batch_size, 1] - Inhibition score
        """
        # Concatenate all feature streams
        x = torch.hstack((gnn_features, interface_features, compound_features))
        # Shape: [batch_size, 2771]
        
        # Hidden layers with tanh activation
        x = torch.tanh(self.fc1(x))         # [batch_size, 1024]
        x = torch.tanh(self.fc2(x))         # [batch_size, 512]
        
        # Output layer with relu then linear
        x = torch.relu(self.fc3(x))         # [batch_size, 100]
        x = self.fc4(x)                      # [batch_size, 1]
        
        return x  # Raw prediction score (typically used with sigmoid loss)
```

**Architecture Rationale**:
- **Tanh in middle layers**: Non-linear transformations, bounded output (-1, 1)
- **ReLU in output preparation**: Allows for better gradient flow
- **Single output neuron**: Binary classification (inhibitor vs non-inhibitor)

---

# 5. COMPLETE DATA FLOW DIAGRAM

```
INPUT: WriteAllexamplesRandomBindersIdsAll_24JAN_Binary.txt
  │
  ├─ Line format: PDB_ID  COMPLEX_NAME  LIGAND_SMILES  LABEL
  │
  ├─────────────────────────────────────────────────────────────────┐
  │                                                                 │
  ├─→ COMPOUND STREAM                                              │
  │   │                                                             │
  │   ├─ getFP(SMILES)                                             │
  │   │   ├─ Parse SMILES string                                   │
  │   │   ├─ Generate Morgan fingerprint (r=3, 2048 bits)          │
  │   │   └─ Return: np.array of shape (2048,)                     │
  │   │                                                             │
  │   ├─ StandardScaler.fit_transform()                            │
  │   │   └─ Normalize: (X - mean) / std                           │
  │   │                                                             │
  │   └─ torch.FloatTensor(Ctr).cuda()                             │
  │       └─ Create GPU tensor: shape (batch, 2048)                │
  │                                                                 │
  ├─→ PROTEIN STREAM                                               │
  │   │                                                             │
  │   ├─ Load PDB files                                            │
  │   │   ├─ Parse structure with BioPython                        │
  │   │   ├─ Extract one-hot atom encodings                        │
  │   │   ├─ Extract one-hot residue encodings                     │
  │   │   └─ Extract connectivity graphs                           │
  │   │                                                             │
  │   ├─ Convert to PyTorch tensors                                │
  │   │   └─ [one_hot_atom, one_hot_res, neigh_same, neigh_diff]  │
  │   │                                                             │
  │   └─ GNN_model() forward pass                                  │
  │       ├─ GNN_First_Layer: Combines features with neighbors     │
  │       ├─ GNN_Layers: Propagates information                    │
  │       └─ Returns: Protein embedding (batch, 512)               │
  │                                                                 │
  ├─→ INTERFACE STREAM                                             │
  │   │                                                             │
  │   ├─ Load pre-extracted features from pickle files:            │
  │   │   ├─ Pos_seqandInterfaceF_dict.npy                         │
  │   │   └─ NewUbench5InterfaceandSeq_dict.npy                    │
  │   │                                                             │
  │   ├─ chainLabel() extraction (if computed dynamically):        │
  │   │   ├─ getDist(): Calculate inter-chain distances            │
  │   │   │  └─ Threshold: 8.0 Angstroms                           │
  │   │   │                                                         │
  │   │   ├─ generate_pair_features(): Count AA pairs              │
  │   │   │  └─ 20 AA × 20 AA = 210 pairs + 1 special = 211       │
  │   │   │                                                         │
  │   │   ├─ extract_feats(): Order features consistently          │
  │   │   │  └─ Load from prote_letter_pair_keys.npy               │
  │   │   │                                                         │
  │   │   └─ normalize(InterfaceF, norm='l2')                      │
  │   │      └─ L2 normalization                                   │
  │   │                                                             │
  │   ├─ StandardScaler.fit_transform()                            │
  │   │   └─ Normalize: (X - mean) / std                           │
  │   │                                                             │
  │   └─ torch.FloatTensor(Ptr).cuda()                             │
  │       └─ Create GPU tensor: shape (batch, 211)                 │
  │                                                                 │
  └──────────────────────────────────────────────────────────────┘
           │              │                 │
           ↓              ↓                 ↓
      [2048-dim]    [512-dim]          [211-dim]
      Compound      Protein (GNN)      Interface
      Features      Features           Features
           │              │                 │
           └──────────┬───┴────────────────┘
                      ↓
              torch.hstack(gnn, interface, compound)
                      ↓
                Shape: (batch, 2771)
                      ↓
              ┌───────────────────────────┐
              │   IPPI_MLP_Net            │
              │                           │
              │  FC1: 2771 → 1024 (tanh)  │
              │  FC2: 1024 → 512 (tanh)   │
              │  FC3: 512 → 100 (relu)    │
              │  FC4: 100 → 1 (linear)    │
              │                           │
              └───────────────────────────┘
                      ↓
            Prediction Score (Real value)
                      ↓
            Binary Cross-Entropy Loss
                      ↓
            Backpropagation & Optimization
```

---

# 6. KEY TECHNICAL DETAILS & NOTES

## 6.1 Data Processing Characteristics

| Aspect | Details |
|---|---|
| **Total Complexes** | 23 protein complexes |
| **Total Inhibitors** | 714 positive examples |
| **Negative Examples** | Generated randomly or from external databases |
| **Train/Test Split** | Leave-One-Complex-Out (LOCO) cross-validation |
| **Feature Validation** | Invalid SMILES automatically filtered |
| **Preprocessing** | StandardScaler with fit on training data |

## 6.2 Feature Engineering Decisions

1. **Morgan Fingerprints**
   - Chosen for chemical similarity representation
   - Radius 3 captures local neighborhood information
   - 2048 bits provide good resolution without over-fitting

2. **GNN for Protein Structures**
   - Captures 3D geometric information
   - Propagates information through atomic connectivity
   - Produces fixed 512-dimensional embeddings

3. **Interface Features**
   - Pre-computed for efficiency
   - 8Å threshold balances sensitivity/specificity
   - L2 normalization emphasizes pattern over magnitude
   - Consistent ordering prevents non-determinism

4. **StandardScaler Application**
   - Applied after all feature extractions
   - Fitted on training set only
   - Prevents data leakage between train/test splits

## 6.3 Important Implementation Notes

**File Handling**:
```python
# Handle whitespace in SMILES
parts = line.split()
smiles = ' '.join(parts[2:-1])  # Join middle tokens
label = float(parts[-1])
```

**Validation**:
```python
# Filter invalid SMILES
if getFP(smiles) is not None:
    # Process this example
```

**Feature Ordering**:
```python
# Load pre-defined key ordering to ensure consistency
key_list = np.load('prote_letter_pair_keys.npy')
# This ensures same feature order across all runs
```

**GPU Management**:
```python
# Move tensors to GPU
use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")

# All models moved to GPU
GNN_model.cuda()
IPPI_Net.cuda()
```

## 6.4 Potential Preprocessing Issues

1. **SMILES with spaces**: Handled by joining middle tokens
2. **Non-standard amino acids**: Mapped to special entry ('_', '_')
3. **Missing PDB files**: Automatically skipped with warning
4. **Invalid structures**: BioPython handles gracefully

---

# SUMMARY

The complete PPI Inhibitors pipeline implements a three-stream deep learning architecture:

1. **Compound Stream**: SMILES → Morgan Fingerprints (2048-dim)
2. **Protein Stream**: PDB → GNN Processing → Embeddings (512-dim)
3. **Interface Stream**: Inter-chain distances → AA pair patterns (211-dim)

All streams are standardized, concatenated, and fed through a 4-layer MLP for binary classification of protein-inhibitor interactions.

Key preprocessing aspects:
- Data loaded from tab-separated text file
- Features extracted and validated independently
- StandardScaler applied for normalization
- PyTorch tensors created for GPU-accelerated training
- Pre-computed interface features loaded for efficiency

