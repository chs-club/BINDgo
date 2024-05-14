# Experiments

In this competition, you’ll develop machine learning (ML) models to predict the binding affinity of small molecules to specific protein targets&mdash;a critical step in drug development for the pharmaceutical industry that would pave the way for more accurate drug discovery.
You’ll help predict which drug-like small molecules (chemicals) will bind to three possible protein targets.

??? note "Citation"

    This is directly from the [Kaggle competition](https://www.kaggle.com/competitions/leash-BELKA) [^kaggle-comp]

## DELs are libraries of small molecules with unique DNA barcodes covalently attached

Traditional high-throughput screening requires keeping individual small molecules in separate, identifiable tubes and demands a lot of liquid handling to test each one of those against the protein target of interest in a separate reaction. The logistical overhead of these efforts tends to restrict screening collections, called libraries, to 50K-5M small molecules. A scalable solution to this problem, DNA-encoded chemical libraries, was described in 2009. As DNA sequencing got cheaper and cheaper, it became clear that DNA itself could be used as a label to identify, and deconvolute, collections of molecules in a complex mixture. DELs leverage this DNA sequencing technology.

These barcoded small molecules are in a pool (many in a single tube, rather than one tube per small molecule) and are exposed to the protein target of interest in solution. The protein target of interest is then rinsed to remove small molecules in the DEL that don’t bind the target, and the remaining binders are collected and their DNA sequenced.

## DELs are manufactured by combining different building blocks

An intuitive way to think about DELs is to imagine a Mickey Mouse head as an example of a small molecule in the DEL. We attach the DNA barcode to Mickey’s chin. Mickey’s left ear is connected by a zipper; Mickey’s right ear is connected by velcro. These attachment points of zippers and velcro are analogies to different chemical reactions one might use to construct the DEL.

We could purchase ten different Mickey Mouse faces, ten different zipper ears, and ten different velcro ears, and use them to construct our small molecule library. By creating every combination of these three, we’ll have 1,000 small molecules, but we only needed thirty building blocks (faces and ears) to make them. This combinatorial approach is what allows DELs to have so many members: the library in this competition is composed of 133M small molecules. The 133M small molecule library used here, AMA014, was provided by AlphaMa. It has a triazine core and superficially resembles the DELs described here.

![](https://www.googleapis.com/download/storage/v1/b/kaggle-user-content/o/inbox%2F1095143%2F1901c6caa0c6c011617f4dec525d7bbe%2FKaggle%20v2%20(1).png?generation=1712179256934503&alt=media)

<!-- REFERENCES -->

[^kaggle-comp]: Andrew Blevins, Ian K Quigley, Brayden J Halverson, Nate Wilkinson, Rebecca S Levin, Agastya Pulapaka, Walter Reade, Addison Howard. (2024). Leash Bio - Predict New Medicines with BELKA. Kaggle. https://kaggle.com/competitions/leash-BELKA
