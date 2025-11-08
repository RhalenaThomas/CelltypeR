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

# branch protection is on need to make pull request
git checkout -b fix-package-imports # switches to a new branch
git push -u origin fix-package-imports


```


After I accept the merge request

```
# Switch to main
git checkout main

# Update it with the merged changes
git pull origin main

# (Optional) delete your local feature branch
git branch -d fix-package-imports

```



Other notes:
I have functions in progress in other files but these are in the R folder and still accessed. 
The last file to be loaded alphabetically is the file from which a function will be selected
I change the make_seu function that was updated for seurat 5 to make_seu5 to avoid issues.


Check examples

```bash

grep -R "@examples" R/
grep -R "\dontrun" R/
```

need to put don't run for those that load an object or run something intensive

