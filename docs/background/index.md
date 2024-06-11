# Background

In this competition, you’ll develop machine learning (ML) models to predict the binding affinity of small molecules to specific protein targets&mdash;a critical step in drug development for the pharmaceutical industry that would pave the way for more accurate drug discovery.
You’ll help predict which drug-like small molecules (chemicals) will bind to three possible protein targets.

??? note "Citation"

    This is directly from the [Kaggle competition](https://www.kaggle.com/competitions/leash-BELKA) [^kaggle-comp]

## Description

Small molecule drugs are chemicals that interact with cellular protein machinery and affect the functions of this machinery in some way.
Often, drugs are meant to inhibit the activity of single protein targets, and those targets are thought to be involved in a disease process.
A classic approach to identify such candidate molecules is to physically make them, one by one, and then expose them to the protein target of interest and test if the two interact.
This can be a fairly laborious and time-intensive process.

The US Food and Drug Administration (FDA) has approved roughly 2,000 novel molecular entities in its entire history.
However, the number of chemicals in druglike space has been estimated to be 10^60, a space far too big to physically search.
There are likely effective treatments for human ailments hiding in that chemical space, and better methods to find such treatments are desirable to us all.

To evaluate potential search methods in small molecule chemistry, competition host Leash Biosciences physically tested some 133M small molecules for their ability to interact with one of three protein targets using DNA-encoded chemical library (DEL) technology.
This dataset, the Big Encoded Library for Chemical Assessment (BELKA), provides an excellent opportunity to develop predictive models that may advance drug discovery.

Datasets of this size are rare and restricted to large pharmaceutical companies.
The current best-curated public dataset of this kind is perhaps bindingdb, which, at 2.8M binding measurements, is much smaller than BELKA.

This competition aims to revolutionize small molecule binding prediction by harnessing ML techniques.
Recent advances in ML approaches suggest it might be possible to search chemical space by inference using well-trained computational models rather than running laboratory experiments.
Similar progress in other fields suggest using ML to search across vast spaces could be a generalizable approach applicable to many domains.
We hope that by providing BELKA we will democratize aspects of computational drug discovery and assist the community in finding new lifesaving medicines.

Here, you’ll build predictive models to estimate the binding affinity of unknown chemical compounds to specified protein targets.
You may use the training data provided; alternatively, there are a number of methods to make small molecule binding predictions without relying on empirical binding data (e.g., DiffDock, and this contest was designed to allow for such submissions).

Your work will contribute to advances in small molecule chemistry used to accelerate drug discovery.

<!-- REFERENCES -->

[^kaggle-comp]: Andrew Blevins, Ian K Quigley, Brayden J Halverson, Nate Wilkinson, Rebecca S Levin, Agastya Pulapaka, Walter Reade, Addison Howard. (2024). Leash Bio - Predict New Medicines with BELKA. Kaggle. https://kaggle.com/competitions/leash-BELKA
