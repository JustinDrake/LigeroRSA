# this script returns version to the build process
#GIT_HASH=`git rev-parse --short HEAD`
GIT_HASH=123
echo $GIT_HASH
echo $GIT_HASH>version.out
