# CROWN-seq

**C**onversion **R**esistance detection **O**n **W**hole-transcriptomic **N6**,2′-O-dimethyladenosine transcription-start sites by **seq**uencing (**CROWN-seq**), is a single base resolution m6Am detection and quantification method developed by Jianheng Liu (Fox) @ **Jaffrey Lab** (https://www.jaffreylab.org/). By this method, one can profile the m6Am stoichiometry at single transcription-start nucleotide level. This Github Repo is for the computational pipeline for CROWN-Seq.

## Note

This repo is currently for illustration only. The detailed experimental protocol, and some other documents will be uploaded as soon as manuscript accepted.

## How to

The CROWN-Seq pipeline consists of three steps:

1. Obtain accurate transcription-start sites (TSS) by ReCappable-Seq or CROWN-Seq.
2. Merge all TSSs found.
3. Measure the m6Am levels. 

Read the notebooks under `notebooks` folder for details. 

Scripts for metadata generation are under `scripts` folder (and also under the `notebooks` folder).

## Contact

Please email Jianheng Liu (Fox) if you have any question:

jil4026@med.cornell.edu or jhfoxliu@gmail.com

## License

MIT.

