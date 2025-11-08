# restructure repo for cleaner github installation

In CL

```bash
cd /Users/rhalenathomas/GITHUB/CelltypeR
mv CelltypeR/* ./
rm -r CelltypeR # remove empty folder
```

The R package files are now out in the main repo

# populate package updates to NAMESPACE
This must be done each time the library is updated for the new changes to be in the library


In CL
```bash
cd /Users/rhalenathomas/GITHUB/CelltypeR
# open R
R
```

Inside R in the repo library folder

```r
devtools::document() #  R will see the documentation
# verify directory
getwd()
# check the files
list.files()
# run  the update
devtools::check()
# apply 
devtools::install()

```
- takes a while

# add and push changes

CL

```bash
cd /Users/rhalenathomas/GITHUB/CelltypeR
git add -A
git commit -m "library import fixed"
```


Other notes:
I have functions in progress in other files but these are in the R folder and still accessed. 
The last file to be loaded alphabetically is the file from which a function will be selected
I change the make_seu function that was updated for seurat 5 to make_seu5 to avoid issues.


Check examples

```bash

grep -R "@examples" R/
```

need to put don't run for those that load an object or run something intensive

#' @examples
#' \dontrun{
#' fsc_to_fs(input_folder_fsc, downsample = "none")
#' }
