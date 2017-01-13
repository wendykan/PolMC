# Very simple script to automake the documentation

# Make all the documentation from header into an html file
echo Processing headerdoc info

find . -name "*.h" -exec headerdoc2html -C -o doc/ {} \;
echo Making table of contents
cd doc

# Annoying problem with gatherheaderdoc and subversion, which I have not solved yet:
# gatherHeaderdoc keeps insisting on going inside the .svn directories.

find . -name ".svn" -exec chmod u-rwx {} \;
gatherheaderdoc .
find . -name ".svn" -exec chmod u+rwx {} \;

