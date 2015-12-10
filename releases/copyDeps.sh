#!/bin/bash



function usage()
{
    cat << EOU
Usage: bash $0 <path to binary> <path to library directory>
EOU
exit 1
}


#Validate the inputs
[[ $# < 2 ]] && usage

#Check if paths are valid
[[ ! -e $1 ]] && echo "Not a valid input $1" && exit 1
[[ -d $2 ]] || echo "Directory $2 doesn't exist. Creating..."&& mkdir -p "$2"

#Get the library dependencies
echo "Collecting shared library dependencies for $1..."
deps=$(ldd $1 | awk 'BEGIN{ORS=" "}$1\
~/^\//{print $1}$3~/^\//{print $3}'\
 | sed 's/,$/\n/')
echo "Copying dependencies to $2"

#Copy the deps
for dep in $deps
do
    echo "Copying $dep to $2"
    cp "$dep" "$2"
done

echo "Done!"
