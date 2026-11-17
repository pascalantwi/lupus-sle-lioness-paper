# References

<!-- The project's reference list. This project keeps no Zotero library: this file is the only source
for anything cited in a deliverable from the proposal onward (workflow/rules.md, section 2).

- The agent writes each entry from a record it opened (the DOI's Crossref record or the publisher's
  page, PubMed, arXiv, bioRxiv, CRAN or Bioconductor for software) and names that record, with the
  date, in the `record` field. Nothing is written from memory.
- `references.bib` is generated from this file by
  `python3 -I manuscript/build-bib.py references/references.md manuscript/references.bib --citations manuscript/citations.csv`
  and is never edited by hand. Faults go to `reference-fixes.md`.
- The quoted line under each heading is rewritten by build-bib.py on every run: edit the fields, not
  the quoted line.
- The passage each citation rests on is recorded in `citations.csv`, not here.
- `full_text` says where the full text was read: a URL (an open copy), `references/pdfs/<key>.pdf`
  (a PDF the user added; the folder is git-ignored), `abstract only` (allowed for a supporting
  citation), or `pending` (a core citation still waiting for its full text).
- The user reviews this file once, at the end (workflow/03-manuscript.md, step 6, the reference review).

Entry form (fields: type, author, editor, title, journal, booktitle, publisher, address, series,
edition, institution, school, year, volume, number, pages, doi, url, note, record, full_text):

### key2024

- type: article
- author: Family, Given and Family, Given
- title: Title in sentence case, with {R} or {LIONESS} in braces to keep the capitals
- journal: Journal Name
- year: 2024
- volume: 1
- number: 2
- pages: 10-20
- doi: 10.xxxx/xxxxx
- record: Crossref record for the DOI, opened 2026-11-03
- full_text: https://... (open copy), or references/pdfs/key2024.pdf, or abstract only, or pending
-->

## Seed entries

Cited by every manuscript of this workflow (the ADEMP sentence and the reproducibility paragraph).

### morris2019

> Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine*, 38(11), 2074-2102. https://doi.org/10.1002/sim.8086

- type: article
- author: Morris, Tim P. and White, Ian R. and Crowther, Michael J.
- title: Using simulation studies to evaluate statistical methods
- journal: Statistics in Medicine
- year: 2019
- volume: 38
- number: 11
- pages: 2074-2102
- doi: 10.1002/sim.8086
- record: seed entry carried from the research-workflow template's references.bib
- full_text: abstract only

### marwick2018

> Marwick, B., Boettiger, C., & Mullen, L. (2018). Packaging data analytical work reproducibly using R (and friends). *The American Statistician*, 72(1), 80-88. https://doi.org/10.1080/00031305.2017.1375986

- type: article
- author: Marwick, Ben and Boettiger, Carl and Mullen, Lincoln
- title: Packaging data analytical work reproducibly using {R} (and friends)
- journal: The American Statistician
- year: 2018
- volume: 72
- number: 1
- pages: 80-88
- doi: 10.1080/00031305.2017.1375986
- record: seed entry carried from the research-workflow template's references.bib
- full_text: abstract only

## Project entries

<!-- Entries added from stage 1 onward, in the order they were added. -->

## Project entries

Carried on 2026-11-13 from the first LIONESS project (`work-flow-for-lioness-paper`), where each record was opened on the date its entry gives, and new entries for this project.

### kuijjer2019lioness

> Kuijjer, M. L., Tung, M. G., Yuan, G., Quackenbush, J., & Glass, K. (2019). Estimating sample-specific regulatory networks. *iScience*, 14, 226-240. https://doi.org/10.1016/j.isci.2019.03.021

- type: article
- author: Kuijjer, M. L. and Tung, M. G. and Yuan, G. and Quackenbush, J. and Glass, K.
- title: Estimating sample-specific regulatory networks
- journal: iScience
- year: 2019
- volume: 14
- pages: 226-240
- doi: 10.1016/j.isci.2019.03.021
- record: PubMed record PMID 30981959, opened 2026-11-03
- full_text: https://pmc.ncbi.nlm.nih.gov/articles/PMC6463816/ (read 2026-11-03)

### kuijjer2019lionessr

> Kuijjer, M. L., Hsieh, P.-H., Quackenbush, J., & Glass, K. (2019). lionessR: single sample network inference in R. *BMC Cancer*, 19(1), 1003. https://doi.org/10.1186/s12885-019-6235-7

- type: article
- author: Kuijjer, Marieke L. and Hsieh, Ping-Han and Quackenbush, John and Glass, Kimberly
- title: {lionessR}: single sample network inference in {R}
- journal: BMC Cancer
- year: 2019
- volume: 19
- number: 1
- pages: 1003
- doi: 10.1186/s12885-019-6235-7
- record: Crossref record for the DOI, opened 2026-11-03
- full_text: https://pmc.ncbi.nlm.nih.gov/articles/PMC6815019/ (read 2026-11-03)

### devlin1975robust

> Devlin, S. J., Gnanadesikan, R., & Kettenring, J. R. (1975). Robust estimation and outlier detection with correlation coefficients. *Biometrika*, 62(3), 531-545. https://doi.org/10.1093/biomet/62.3.531

- type: article
- author: Devlin, Susan J. and Gnanadesikan, R. and Kettenring, J. R.
- title: Robust estimation and outlier detection with correlation coefficients
- journal: Biometrika
- year: 1975
- volume: 62
- number: 3
- pages: 531-545
- doi: 10.1093/biomet/62.3.531
- record: Crossref record for the DOI, opened 2026-11-03
- full_text: pending

### yu2017iena

> Yu, X., Zhang, J., Sun, S., Zhou, X., Zeng, T., & Chen, L. (2017). Individual-specific edge-network analysis for disease prediction. *Nucleic Acids Research*, 45(20), e170. https://doi.org/10.1093/nar/gkx787

- type: article
- author: Yu, Xiangtian and Zhang, Jingsong and Sun, Shaoyan and Zhou, Xin and Zeng, Tao and Chen, Luonan
- title: Individual-specific edge-network analysis for disease prediction
- journal: Nucleic Acids Research
- year: 2017
- volume: 45
- number: 20
- pages: e170
- doi: 10.1093/nar/gkx787
- record: Crossref record for the DOI, opened 2026-11-03
- full_text: https://pmc.ncbi.nlm.nih.gov/articles/PMC5714249/ (read 2026-11-03)

### deschildre2024evaluation

> Deschildre, J., Vandemoortele, B., Loers, J. U., De Preter, K., & Vermeirssen, V. (2024). Evaluation of single-sample network inference methods for precision oncology. *npj Systems Biology and Applications*, 10(1), 18. https://doi.org/10.1038/s41540-024-00340-w

- type: article
- author: Deschildre, Joke and Vandemoortele, Boris and Loers, Jens Uwe and De Preter, Katleen and Vermeirssen, Vanessa
- title: Evaluation of single-sample network inference methods for precision oncology
- journal: npj Systems Biology and Applications
- year: 2024
- volume: 10
- number: 1
- pages: 18
- doi: 10.1038/s41540-024-00340-w
- record: Crossref record for the DOI, opened 2026-11-03
- full_text: https://pmc.ncbi.nlm.nih.gov/articles/PMC10869342/ (read 2026-11-03)

### kuijjer2026challenges

> Kuijjer, M., De Marzio, M., & Glass, K. (2026). Challenges and opportunities in single-sample network modeling. bioRxiv preprint. https://doi.org/10.64898/2026.02.27.708608

- type: misc
- author: Kuijjer, Marieke and De Marzio, Margherita and Glass, Kimberly
- title: Challenges and opportunities in single-sample network modeling
- year: 2026
- doi: 10.64898/2026.02.27.708608
- note: bioRxiv preprint
- record: Crossref record for the DOI, opened 2026-11-03; PubMed PMID 41867852, opened 2026-11-03
- full_text: https://pmc.ncbi.nlm.nih.gov/articles/PMC13001488/ (read 2026-11-03); bioRxiv PDF of the main paper supplied by the user, read 2026-11-06 (references/pdfs/kuijjer2026challenges.pdf; Figure 1, the recast equations, and the Methods, which place the derivations in Supplemental File 1); Supplement 1 (Supplemental File 1, 7 pages) supplied by the user and read 2026-11-06 (references/pdfs/kuijjer2026challenges-supp1.pdf)

### quenouille1956bias

> Quenouille, M. H. (1956). Notes on bias in estimation. *Biometrika*, 43(3-4), 353-360. https://doi.org/10.1093/biomet/43.3-4.353

- type: article
- author: Quenouille, M. H.
- title: Notes on bias in estimation
- journal: Biometrika
- year: 1956
- volume: 43
- number: 3-4
- pages: 353-360
- doi: 10.1093/biomet/43.3-4.353
- record: Crossref record for the DOI, opened 2026-11-06
- full_text: not read; pre-1980 source cited as definitional (the jackknife)

### miller1974jackknife

> Miller, R. G. (1974). The jackknife: a review. *Biometrika*, 61(1), 1-15. https://doi.org/10.1093/biomet/61.1.1

- type: article
- author: Miller, Rupert G.
- title: The jackknife: a review
- journal: Biometrika
- year: 1974
- volume: 61
- number: 1
- pages: 1-15
- doi: 10.1093/biomet/61.1.1
- record: Crossref record for the DOI, opened 2026-11-06
- full_text: not read; pre-1980 source cited as definitional (jackknife pseudo-values)

### friedman2010glmnet

> Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software*, 33(1). https://doi.org/10.18637/jss.v033.i01

- type: article
- author: Friedman, Jerome and Hastie, Trevor and Tibshirani, Robert
- title: Regularization paths for generalized linear models via coordinate descent
- journal: Journal of Statistical Software
- year: 2010
- volume: 33
- number: 1
- doi: 10.18637/jss.v033.i01
- record: Crossref record for the DOI, opened 2026-11-04
- full_text: not needed (definitional)

### nadeau2003inference

> Nadeau, C., & Bengio, Y. (2003). Inference for the generalization error. *Machine Learning*, 52(3), 239-281. https://doi.org/10.1023/A:1024068626366

- type: article
- author: Nadeau, Claude and Bengio, Yoshua
- title: Inference for the generalization error
- journal: Machine Learning
- year: 2003
- volume: 52
- number: 3
- pages: 239-281
- doi: 10.1023/A:1024068626366
- record: Crossref record for the DOI, opened 2026-11-03
- full_text: not needed (definitional)

### schuirmann1987tost

> Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics*, 15(6), 657-680. https://doi.org/10.1007/BF01068419

- type: article
- author: Schuirmann, Donald J.
- title: A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability
- journal: Journal of Pharmacokinetics and Biopharmaceutics
- year: 1987
- volume: 15
- number: 6
- pages: 657-680
- doi: 10.1007/BF01068419
- record: Crossref record for the DOI, opened 2026-11-04
- full_text: not needed (definitional)

### mai2021ifn5

> Mai, L., Asaduzzaman, A., Noamani, B., Fortin, P. R., Gladman, D. D., Touma, Z., Urowitz, M. B., & Wither, J. (2021). The baseline interferon signature predicts disease severity over the subsequent 5 years in systemic lupus erythematosus. *Arthritis Research & Therapy*, 23(1), 29. https://doi.org/10.1186/s13075-021-02414-0

- type: article
- author: Mai, Lloyd and Asaduzzaman, Arundip and Noamani, Babak and Fortin, Paul R. and Gladman, Dafna D. and Touma, Zahi and Urowitz, Murray B. and Wither, Joan
- title: The baseline interferon signature predicts disease severity over the subsequent 5 years in systemic lupus erythematosus
- journal: Arthritis Research & Therapy
- year: 2021
- volume: 23
- number: 1
- pages: 29
- doi: 10.1186/s13075-021-02414-0
- record: Crossref record for the DOI, opened 2026-11-03
- full_text: https://pmc.ncbi.nlm.nih.gov/articles/PMC7811214/ (read 2026-11-03)

### lien2024pregnancy

> Lien, H. J. T., Pedersen, T. T., Jakobsen, B., Flatberg, A., Chawla, K., S\aetrom, P., & Fenstad, M. H. (2024). Single-cell resolution of longitudinal blood transcriptome profiles in rheumatoid arthritis, systemic lupus erythematosus and healthy control pregnancies. *Annals of the Rheumatic Diseases*, 83(3), 300-311. https://doi.org/10.1136/ard-2023-224644

- type: article
- author: Lien, Hilde Johanne Tokle and Pedersen, Tina T. and Jakobsen, Bente and Flatberg, Arnar and Chawla, Konika and S{\ae}trom, P{\aa}l and Fenstad, Mona H.
- title: Single-cell resolution of longitudinal blood transcriptome profiles in rheumatoid arthritis, systemic lupus erythematosus and healthy control pregnancies
- journal: Annals of the Rheumatic Diseases
- year: 2024
- volume: 83
- number: 3
- pages: 300-311
- doi: 10.1136/ard-2023-224644
- record: PubMed record PMID 38049980, opened 2026-11-13; GEO series GSE235508, opened 2026-11-13
- full_text: abstract only (PubMed; PMC10894842 is open); cited as supporting for the dataset and its findings

### banchereau2016personalized

> Banchereau, R., Hong, S., Cantarel, B., Baldwin, N., Baisch, J., Edens, M., Cepika, A. M., Acs, P., Turner, J., Anguiano, E., Vinod, P., Kahn, S., Obermoser, G., Blankenship, D., Wakeland, E., Nassi, L., Gotte, A., Punaro, M., Liu, Y. J., ..., & Pascual, V. (2016). Personalized immunomonitoring uncovers molecular networks that stratify lupus patients. *Cell*, 165(3), 551-565. https://doi.org/10.1016/j.cell.2016.03.008

- type: article
- author: Banchereau, R. and Hong, S. and Cantarel, B. and Baldwin, N. and Baisch, J. and Edens, M. and Cepika, A. M. and Acs, P. and Turner, J. and Anguiano, E. and Vinod, P. and Kahn, S. and Obermoser, G. and Blankenship, D. and Wakeland, E. and Nassi, L. and Gotte, A. and Punaro, M. and Liu, Y. J. and Banchereau, J. and Rossello-Urgell, J. and Wright, T. and Pascual, V.
- title: Personalized immunomonitoring uncovers molecular networks that stratify lupus patients
- journal: Cell
- year: 2016
- volume: 165
- number: 3
- pages: 551-565
- doi: 10.1016/j.cell.2016.03.008
- record: Europe PMC record PMID 27040498, opened 2026-11-16; GEO series GSE65391, opened 2026-11-15
- full_text: abstract only (PMC5426482 is open); cited as supporting for the dataset

### petri2019stability

> Petri, M., Fu, W., Ranger, A., Allaire, N., Cullen, P., Magder, L. S., & Zhang, Y. (2019). Association between changes in gene signatures expression and disease activity among patients with systemic lupus erythematosus. *BMC Medical Genomics*, 12(1), 4. https://doi.org/10.1186/s12920-018-0468-1

- type: article
- author: Petri, M. and Fu, W. and Ranger, A. and Allaire, N. and Cullen, P. and Magder, L. S. and Zhang, Y.
- title: Association between changes in gene signatures expression and disease activity among patients with systemic lupus erythematosus
- journal: BMC Medical Genomics
- year: 2019
- volume: 12
- number: 1
- pages: 4
- doi: 10.1186/s12920-018-0468-1
- record: Europe PMC record PMID 30626389, opened 2026-11-16
- full_text: abstract only (PMC6327466 is open); cited as supporting

### farahbod2020composition

> Farahbod, M., & Pavlidis, P. (2020). Untangling the effects of cellular composition on coexpression analysis. *Genome Research*, 30(6), 849-859. https://doi.org/10.1101/gr.256735.119

- type: article
- author: Farahbod, M. and Pavlidis, P.
- title: Untangling the effects of cellular composition on coexpression analysis
- journal: Genome Research
- year: 2020
- volume: 30
- number: 6
- pages: 849-859
- doi: 10.1101/gr.256735.119
- record: Europe PMC record PMID 32580998, opened 2026-11-16
- full_text: abstract only (PMC7370889 is open); cited as supporting

### whitney2003individuality

> Whitney, A. R., Diehn, M., Popper, S. J., Alizadeh, A. A., Boldrick, J. C., Relman, D. A., & Brown, P. O. (2003). Individuality and variation in gene expression patterns in human blood. *Proceedings of the National Academy of Sciences of the United States of America*, 100(4), 1896-1901. https://doi.org/10.1073/pnas.252784499

- type: article
- author: Whitney, A. R. and Diehn, M. and Popper, S. J. and Alizadeh, A. A. and Boldrick, J. C. and Relman, D. A. and Brown, P. O.
- title: Individuality and variation in gene expression patterns in human blood
- journal: Proceedings of the National Academy of Sciences of the United States of America
- year: 2003
- volume: 100
- number: 4
- pages: 1896-1901
- doi: 10.1073/pnas.252784499
- record: Europe PMC record PMID 12578971, opened 2026-11-16
- full_text: abstract only; cited as supporting

### finn2015fingerprinting

> Finn, E. S., Shen, X., Scheinost, D., Rosenberg, M. D., Huang, J., Chun, M. M., Papademetris, X., & Constable, R. T. (2015). Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. *Nature Neuroscience*, 18(11), 1664-1671. https://doi.org/10.1038/nn.4135

- type: article
- author: Finn, E. S. and Shen, X. and Scheinost, D. and Rosenberg, M. D. and Huang, J. and Chun, M. M. and Papademetris, X. and Constable, R. T.
- title: Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity
- journal: Nature Neuroscience
- year: 2015
- volume: 18
- number: 11
- pages: 1664-1671
- doi: 10.1038/nn.4135
- record: Europe PMC record PMID 26457551, opened 2026-11-16
- full_text: abstract only (PMC5008686 is open); cited as supporting

### mishra2026temporal

> Mishra, N., Kimmig, F., Vandeputte, D., Talevi, V., De Commer, L., Verspecht, C., Vich Vila, A., El-Sayed Moustafa, J. S., Kreft, L., Botzki, A., El Darzi, Y., Proost, S., Devolder, L., Wang, D., Bernardes, J. P., Aziz, N. A., SYSCID consortium, Franke, A., Schreiber, S., ..., & Rosenstiel, P. (2026). Large-scale analysis of temporal gene expression variation in peripheral blood. *Nature Communications*, 17(1), 6992. https://doi.org/10.1038/s41467-026-73218-6

- type: article
- author: Mishra, N. and Kimmig, F. and Vandeputte, D. and Talevi, V. and De Commer, L. and Verspecht, C. and Vich Vila, A. and El-Sayed Moustafa, J. S. and Kreft, L. and Botzki, A. and El Darzi, Y. and Proost, S. and Devolder, L. and Wang, D. and Bernardes, J. P. and Aziz, N. A. and {SYSCID consortium} and Franke, A. and Schreiber, S. and Dermitzakis, E. T. and Vieira-Silva, S. and Falony, G. and Small, K. S. and Breteler, M. M. B. and Schultze, J. L. and Raes, J. and Rosenstiel, P.
- title: Large-scale analysis of temporal gene expression variation in peripheral blood
- journal: Nature Communications
- year: 2026
- volume: 17
- number: 1
- pages: 6992
- doi: 10.1038/s41467-026-73218-6
- record: Europe PMC record PMID 42215465, opened 2026-11-16
- full_text: abstract only (PMC13392353 is open); cited as supporting

### altman2021bloodgen3

> Altman, M. C., Rinchai, D., Baldwin, N., Toufiq, M., Whalen, E., Garand, M., Syed Ahamed Kabeer, B., Alfaki, M., Presnell, S. R., Khaenam, P., Ayllón-Benítez, A., Mougin, F., Thébault, P., Chiche, L., Jourde-Chiche, N., Phillips, J. T., Klintmalm, G., O'Garra, A., Berry, M., ..., & Chaussabel, D. (2021). Development of a fixed module repertoire for the analysis and interpretation of blood transcriptome data. *Nature Communications*, 12(1), 4385. https://doi.org/10.1038/s41467-021-24584-w

- type: article
- author: Altman, M. C. and Rinchai, D. and Baldwin, N. and Toufiq, M. and Whalen, E. and Garand, M. and Syed Ahamed Kabeer, B. and Alfaki, M. and Presnell, S. R. and Khaenam, P. and Ayllón-Benítez, A. and Mougin, F. and Thébault, P. and Chiche, L. and Jourde-Chiche, N. and Phillips, J. T. and Klintmalm, G. and O'Garra, A. and Berry, M. and Bloom, C. and Wilkinson, R. J. and Graham, C. M. and Lipman, M. and Lertmemongkolchai, G. and Bedognetti, D. and Thiebaut, R. and Kheradmand, F. and Mejias, A. and Ramilo, O. and Palucka, K. and Pascual, V. and Banchereau, J. and Chaussabel, D.
- title: Development of a fixed module repertoire for the analysis and interpretation of blood transcriptome data
- journal: Nature Communications
- year: 2021
- volume: 12
- number: 1
- pages: 4385
- doi: 10.1038/s41467-021-24584-w
- record: Crossref record for the DOI, opened 2026-11-17; Europe PMC record PMID 34282143, opened 2026-11-17
- full_text: https://europepmc.org/article/PMC/PMC8289976 (open, CC BY 4.0), read 2026-11-17

### rinchai2021bloodgen3module

> Rinchai, D., Roelands, J., Toufiq, M., Hendrickx, W., Altman, M. C., Bedognetti, D., & Chaussabel, D. (2021). BloodGen3Module: blood transcriptional module repertoire analysis and visualization using R. *Bioinformatics*, 37(16), 2382-2389. https://doi.org/10.1093/bioinformatics/btab121

- type: article
- author: Rinchai, D. and Roelands, J. and Toufiq, M. and Hendrickx, W. and Altman, M. C. and Bedognetti, D. and Chaussabel, D.
- title: BloodGen3Module: blood transcriptional module repertoire analysis and visualization using R
- journal: Bioinformatics
- year: 2021
- volume: 37
- number: 16
- pages: 2382-2389
- doi: 10.1093/bioinformatics/btab121
- record: Crossref record for the DOI, opened 2026-11-17; Europe PMC record PMID 33624743, opened 2026-11-17; Bioconductor package page (BloodGen3Module 1.20.0, GPL-2), opened 2026-11-17
- full_text: abstract only (PMC8388021 is open); cited as definitional for the module lists

### liberzon2015hallmark

> Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database Hallmark gene set collection. *Cell Systems*, 1(6), 417-425. https://doi.org/10.1016/j.cels.2015.12.004

- type: article
- author: Liberzon, A. and Birger, C. and Thorvaldsdóttir, H. and Ghandi, M. and Mesirov, J. P. and Tamayo, P.
- title: The Molecular Signatures Database Hallmark gene set collection
- journal: Cell Systems
- year: 2015
- volume: 1
- number: 6
- pages: 417-425
- doi: 10.1016/j.cels.2015.12.004
- record: Crossref record for the DOI, opened 2026-11-17; Europe PMC record PMID 26771021, opened 2026-11-17
- full_text: abstract only (PMC4707969 is open); cited as definitional for the gene sets

### horvath2008geometric

> Horvath, S., & Dong, J. (2008). Geometric interpretation of gene coexpression network analysis. *PLoS Computational Biology*, 4(8), e1000117. https://doi.org/10.1371/journal.pcbi.1000117

- type: article
- author: Horvath, S. and Dong, J.
- title: Geometric interpretation of gene coexpression network analysis
- journal: PLoS Computational Biology
- year: 2008
- volume: 4
- number: 8
- pages: e1000117
- doi: 10.1371/journal.pcbi.1000117
- record: Crossref record for the DOI, opened 2026-11-17; Europe PMC record PMID 18704157, opened 2026-11-17
- full_text: https://europepmc.org/article/PMC/PMC2446438 (open, CC BY 4.0), read 2026-11-17

### langfelder2008wgcna

> Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9, 559. https://doi.org/10.1186/1471-2105-9-559

- type: article
- author: Langfelder, P. and Horvath, S.
- title: WGCNA: an R package for weighted correlation network analysis
- journal: BMC Bioinformatics
- year: 2008
- volume: 9
- pages: 559
- doi: 10.1186/1471-2105-9-559
- record: Crossref record for the DOI, opened 2026-11-17; Europe PMC record PMID 19114008, opened 2026-11-17
- full_text: abstract only (PMC2631488 is open); cited as supporting
