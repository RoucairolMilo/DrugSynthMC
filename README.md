DrugSynthMC is a molecule generation software.

As is, it generates molecules by building the SMILES character by character to produce Drug-like molecules. It stops once the molecule is drug-like.
You can modify the score function to suit your needs.

--- How to use ---

- Clone the repository.
- Install the libraries using the toml.
- Run main (don't forget the --release option).

--- How to modify ---

The structure of that program is as follows:
- It generates the molecules using search model named SMILESgen.rs in src\models. This model defines the search space.
- That model is called by various algorithms in the src\methods folder.
- ngrams are in the ngrams folder, neural networks are in the Neural folder

To use different algorithms, change the one called in the main.
To use different priors, replace them in the ngrams and Neural folders.
