# Data

One of the goals of this competition is to explore and compare many different ways of representing molecules. Small molecules have been represented with SMILES, graphs, 3D structures, and more, including more esoteric methods such as spherical convolutional neural nets. We encourage competitors to explore not only different methods of making predictions but also to try different ways of representing the molecules.

??? note "Citation"

    This is directly from the [Kaggle competition](https://www.kaggle.com/competitions/leash-BELKA) [^kaggle-comp]

## Format

`[train/test].[csv/parquet]` - The train or test data, available in both the csv and parquet formats.

-   `id`: A unique example_id that we use to identify the molecule-binding target pair.
-   `buildingblock1_smiles`: The structure, in SMILES, of the first building block
-   `buildingblock2_smiles`: The structure, in SMILES, of the second building block
-   `buildingblock3_smiles`: The structure, in SMILES, of the third building block
-   `molecule_smiles`: The structure of the fully assembled molecule, in SMILES. This includes the three building blocks and the triazine core. Note we use a [Dy] as the stand-in for the DNA linker.
-   `protein_name`: The protein target name
-   `binds`: The target column. A binary class label of whether the molecule binds to the protein. Not available for the test set.

## SMILES

We provide the molecules in SMILES format.

SMILES is a concise string notation used to represent the structure of chemical molecules. It encodes the molecular graph, including atoms, bonds, connectivity, and stereochemistry as a linear sequence of characters, by traversing the molecule graph. SMILES is widely used in machine learning applications for chemistry, such as molecular property prediction, drug discovery, and materials design, as it provides a standardized and machine-readable format for representing and manipulating chemical structures.

The SMILES in this dataset should be sufficient to be translated into any other chemical representation format that you want to try. A simple way to perform some of these translations is with RDKit.

## Targets

Proteins are encoded in the genome, and names of the genes encoding those proteins are typically bestowed by their discoverers and regulated by the Hugo Gene Nomenclature Committee. The protein products of these genes can sometimes have different names, often due to the history of their discovery.

We screened three protein targets for this competition.

### EPHX2 (sEH)

The first target, epoxide hydrolase 2, is encoded by the EPHX2 genetic locus, and its protein product is commonly named “soluble epoxide hydrolase”, or abbreviated to sEH.
Hydrolases are enzymes that catalyze certain chemical reactions, and EPHX2/sEH also hydrolyzes certain phosphate groups. EPHX2/sEH is a potential drug target for high blood pressure and diabetes progression, and small molecules inhibiting EPHX2/sEH from earlier DEL efforts made it to clinical trials.

<div id="seh-view" class="mol-container"></div>
<script>
var uri = 'https://files.rcsb.org/view/3I28.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#seh-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        // viewer.center();
        viewer.setView([ -54.14339497666153, -0.520440320071149, -74.35662880640125, -87.04864632953522, -0.7375537098595212, 0.06842126065013904, 0.39277791256246963, 0.5450307950626022 ]);
        viewer.setClickable({}, true, function(atom,viewer,event,container) {
            console.log(viewer.getView());
        });
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

EPHX2/sEH was also screened with DELs, and hits predicted with ML approaches, in a recent study but the screening data were not published. We included EPHX2/sEH to allow contestants an external gut check for model performance by comparing to these previously-published results.

We screened EPHX2/sEH purchased from Cayman Chemical, a life sciences commercial vendor. For those contestants wishing to incorporate protein structural information in their submissions, the amino sequence is positions 2-555 from UniProt entry P34913, the crystal structure can be found in PDB entry 3I28, and predicted structure can be found in AlphaFold2 entry 34913. Additional EPHX2/sEH crystal structures with ligands bound can be found in PDB.

### BRD4

The second target, bromodomain 4, is encoded by the BRD4 locus and its protein product is also named BRD4. Bromodomains bind to protein spools in the nucleus that DNA wraps around (called histones) and affect the likelihood that the DNA nearby is going to be transcribed, producing new gene products. Bromodomains play roles in cancer progression and a number of drugs have been discovered to inhibit their activities.

<div id="BRD4-view" class="mol-container"></div>
<script>
var uri = 'https://alphafold.ebi.ac.uk/files/AF-O60885-F1-model_v4.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#BRD4-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        // viewer.center();
        viewer.setView([ 1.8962023531609569, -5.7543611915211486, -5.216623587636537, -442.4397691230421, 0.08929251685511222, 0.09354281047578777, -0.04145215887468375, 0.9907362452068654 ]);
        viewer.setClickable({}, true, function(atom,viewer,event,container) {
            console.log(viewer.getView());
        });
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

BRD4 has been screened with DEL approaches previously but the screening data were not published. We included BRD4 to allow contestants to evaluate candidate molecules for oncology indications.

We screened BRD4 purchased from Active Motif, a life sciences commercial vendor. For those contestants wishing to incorporate protein structural information in their submissions, the amino acid sequence is positions 44-460 from [UniProt entry O60885-1](https://www.uniprot.org/uniprotkb/O60885/entry#O60885-1), the crystal structure (for a single domain) can be found in [PDB entry 7USK](https://www.rcsb.org/structure/7USK) and predicted structure can be found in [AlphaFold2 entry O60885](https://alphafold.ebi.ac.uk/entry/O60885).
Additional BRD4 crystal structures with ligands bound can be found in PDB.

### ALB (HSA)

The third target, serum albumin, is encoded by the ALB locus and its protein product is also named ALB. The protein product is sometimes abbreviated as HSA, for “human serum albumin”. ALB, the most common protein in the blood, is used to drive osmotic pressure (to bring fluid back from tissues into blood vessels) and to transport many ligands, hormones, fatty acids, and more.

<div id="ALB-view" class="mol-container"></div>
<script>
var uri = 'https://files.rcsb.org/view/1AO6.pdb';
jQuery.ajax( uri, {
    success: function(data) {
        // https://3dmol.org/doc/GLViewer.html
        let viewer = $3Dmol.createViewer(
            document.querySelector('#ALB-view'),
            { backgroundAlpha: '0.0' }
        );
        viewer.addModel( data, 'pdb' );
        viewer.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
        viewer.setStyle({chain: 'B'}, {});
        // viewer.center({chain: 'A'});
        viewer.setView([ -29.55014298131246, -31.851400043459382, -23.55865710560621, -143.61579985812273, -0.6529007355218823, 0.24661957208370378, -0.4024504947276198, 0.5923959955247251 ]);
        viewer.setClickable({}, true, function(atom,viewer,event,container) {
            console.log(viewer.getView());
        });
        viewer.render();
    },
    error: function(hdr, status, err) {
        console.error( "Failed to load " + uri + ": " + err );
    },
});
</script>

Albumin, being the most abundant protein in the blood, often plays a role in absorbing candidate drugs in the body and sequestering them from their target tissues. Adjusting candidate drugs to bind less to albumin and other blood proteins is a strategy to help these candidate drugs be more effective.

ALB has been screened with DEL approaches previously but the screening data were not published. We included ALB to allow contestants to build models that might have a larger impact on drug discovery across many disease types. The ability to predict ALB binding well would allow drug developers to improve their candidate small molecule therapies much more quickly than physically manufacturing many variants and testing them against ALB empirically in an iterative process.

We screened ALB purchased from Active Motif. For those contestants wishing to incorporate protein structural information in their submissions, the amino acid sequence is positions 25 to 609 from UniProt entry P02768, the crystal structure can be found in PDB entry 1AO6, and predicted structure can be found in AlphaFold2 entry P02768. Additional ALB crystal structures with ligands bound can be found in PDB.

<!-- REFERENCES -->

[^kaggle-comp]: Andrew Blevins, Ian K Quigley, Brayden J Halverson, Nate Wilkinson, Rebecca S Levin, Agastya Pulapaka, Walter Reade, Addison Howard. (2024). Leash Bio - Predict New Medicines with BELKA. Kaggle. https://kaggle.com/competitions/leash-BELKA
