target=master
default_src='develop'
source=${1:-$default_src}

echo "merging $source into $target using branchTEMP"
# in case branchA is not our current branch
git checkout $target

# make merge commit but without conflicts!!
# the contents of 'ours' will be discarded later
git merge -s ours $source

# make temporary branch to merged commit
git branch branchTEMP

# get contents of working tree and index to the one of branchB
git reset --hard $source

# reset to our merged commit but 
# keep contents of working tree and index
git reset --soft branchTEMP

# change the contents of the merged commit
# with the contents of branchB
git commit --amend 

# get rid off our temporary branch
git branch -D branchTEMP

# verify that the merge commit contains only contents of branchB
git diff HEAD ${source}
