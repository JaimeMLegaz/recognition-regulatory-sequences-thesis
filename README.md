### Deep recognition of regulatory sequences

This is the project corresponding to the final thesis of the student Jaime Martínez Legaz at the University of Murcia. Here you can find the Markdown scripts that were used in order to work in this thesis and generate the results that appear in it.

By downloading all these Markdowns, anyone can replicate what we did and get the results we got. You only need to instal Rstudio and all the packages required by the libraries we work with, in addition to the package developed by said student that will be used when treating the sequences. To use that R package, it is only needed to write in the R console the following commands:

```
devtools::install_github("JaimeMLegaz/setup-sequences")
library(SetupSequences)
```
## How to work with these Markdowns

We consider everyone reading this knows the basics of R Markdowns. If not, feel free to check these tutorials: https://rmarkdown.rstudio.com/lesson-1.html

You should have the exact same distributions of folders that appear in this project, since the markdowns will check for some files with a specific route. For example, the file "Promoters_ObtainSequences.Rmd" will need a file called "Data/allGenes", so the folder "Data" will have to be in its directory, containing said file "allGenes". As said, not modifying the directory structure of this project should be enough.

The Markdowns are named with this convention: "experiment"\_"step". For example, "Promoters_ObtainSequences.Rmd" refers to the step of obtaining the sequences of the Promoters experiment. There are three different experiments (Promoters, UTR and PromvsUTR) and three different steps (ObtainSequences, TreatingSequences, and TrainNetwork).

The Markdowns should be executed in this order: ObtainSequences -> TreatingSequences -> TrainNetwork, since the first two generate important data required for the next one. That data can be found partially in the folders of this project that refer to these three experiments. Since GitHub doesn't let us to upload files bigger than 25MB, some files couldn't be uploaded. This, however, should not be ever a problem if you just execute the Markdowns in the recomended order, since you will be generating those files.

## Already generated HTML Markdowns

We took the liberty of generating all the HTML documents that show the result of the markdowns. This way, these scripts would only need to be executed if you want to check if the work is legit or are interested on modifying the parameters used.

These HTML documents can be found in the HTML Markdowns folder.


